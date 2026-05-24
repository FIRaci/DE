function [F_comm, F_sens, best_fit] = solve_DE(H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, F_warm_start)
% SOLVE_DE: Royal Bloodline + Expansion Strategy
% Cải tiến: Mở rộng vùng tìm kiếm quanh Vua để không bị kẹt

    generations = 200;   
    pop_size = 50;      
    F_weight = 0.5;     % Tăng nhẹ F để nhảy xa hơn
    CR = 0.9;           
    
    [U, M, N] = size(H_comm);
    dim = 2 * U * M * N; 
    
    %% 2. Initialization (Expanded Royal Mode)
    x_warm = [real(F_warm_start(:)); imag(F_warm_start(:))]';
    population = zeros(pop_size, dim);
    
    % --- 1. THE KING (Vua - Giữ nguyên) ---
    population(1, :) = x_warm;
    
    % --- 2. ROYAL GUARDS (Cận vệ - Bám sát) ---
    % 30% dân số chỉ biến động nhẹ để tinh chỉnh
    num_guards = round(pop_size * 0.3);
    for i = 2:num_guards
        population(i, :) = x_warm + 0.005 * randn(1, dim); 
    end
    
    % --- 3. EXPEDITIONARY FORCE (Lính viễn chinh - Mở rộng) ---
    % 70% dân số đi xa hơn để tìm kiếm cơ hội (như bọn Original)
    % Noise tăng lên 0.1 - 0.2 để thoát khỏi vùng cục bộ của JSC
    for i = (num_guards + 1):pop_size
        population(i, :) = x_warm + 0.15 * randn(1, dim); 
    end
    
    pop_obj = zeros(pop_size, 1);
    pop_vio = zeros(pop_size, 1);
    
    for i = 1:pop_size
        [pop_obj(i), pop_vio(i)] = evaluate_candidate(population(i, :), H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
    end
    
    best_idx = 1;
    for i = 2:pop_size
        if check_deb_wins(pop_obj(i), pop_vio(i), pop_obj(best_idx), pop_vio(best_idx))
            best_idx = i;
        end
    end
    best_sol = population(best_idx, :);
    best_fit = pop_obj(best_idx); 
    
    %% 3. Main Loop
    for g = 1:generations
        for i = 1:pop_size
            idxs = [1:i-1, i+1:pop_size];
            r = idxs(randperm(length(idxs), 3));
            mutant = population(r(1), :) + F_weight * (population(r(2), :) - population(r(3), :));
            cross_points = rand(1, dim) < CR;
            if ~any(cross_points), cross_points(randi(dim)) = true; end
            trial = population(i, :); trial(cross_points) = mutant(cross_points);
            
            [trial_obj, trial_vio] = evaluate_candidate(trial, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
            
            if check_deb_wins(trial_obj, trial_vio, pop_obj(i), pop_vio(i))
                population(i, :) = trial;
                pop_obj(i) = trial_obj;
                pop_vio(i) = trial_vio;
                
                [b_obj, b_vio] = evaluate_candidate(best_sol, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
                if check_deb_wins(trial_obj, trial_vio, b_obj, b_vio)
                    best_sol = trial;
                    best_fit = trial_obj;
                end
            end
        end
    end
    [F_comm, F_sens] = reconstruct_matrices(best_sol, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle);
end

%% --- HELPER FUNCTIONS ---
function [obj, vio] = evaluate_candidate(x, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N)
    [F_comm, F_sens] = reconstruct_matrices(x, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle);
    ssnr = compute_sensing_SNR_handle(sigmasq_radar_rcs, sensing_beamsteering, F_comm, F_sens);
    SINR = compute_SINR_handle(H_comm, F_comm, F_sens, sigmasq_comm);
    min_sinr = min(SINR);
    
    % Comm là Bố (100x)
    obj = 100 * min_sinr + ssnr; 
    
    target_threshold = gamma_target - 1e-4;
    if min_sinr >= target_threshold, vio = 0; else, vio = target_threshold - min_sinr; end
end

function wins = check_deb_wins(obj_A, vio_A, obj_B, vio_B)
    if (vio_A == 0) && (vio_B > 0), wins = true;
    elseif (vio_A > 0) && (vio_B > 0), wins = (vio_A < vio_B);
    elseif (vio_A == 0) && (vio_B == 0), wins = (obj_A > obj_B);
    else, wins = false; end
end

function [F_comm, F_sens] = reconstruct_matrices(x, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle)
    num_elements = U * M * N; F_real = x(1:num_elements); F_imag = x(num_elements+1:end);
    F_comm = reshape(F_real + 1j * F_imag, [U, M, N]);
    pow_per_ap_comm = sum(abs(F_comm).^2, [1, 3]); 
    for m = 1:M, if pow_per_ap_comm(1, m, 1) > P_comm + 1e-4, F_comm(:, m, :) = F_comm(:, m, :) * sqrt(P_comm / pow_per_ap_comm(1, m, 1)); end; end
    F_sens_base = beam_mapping_handle(H_comm, sensing_beamsteering); 
    curr_p = max(sum(abs(F_sens_base).^2, [1, 3]), [], 'all'); if curr_p > 0, F_sens = F_sens_base * sqrt(P_sensing / curr_p); else, F_sens = F_sens_base; end
end