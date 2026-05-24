function [F_comm, F_sens, best_fit] = solve_SHADE(H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, F_warm_start)
% SOLVE_SHADE: Royal Bloodline + Expansion

    generations = 200; pop_size = 50; H_mem_size = 50;
    [U, M, N] = size(H_comm); dim = 2 * U * M * N;    
    
    x_warm = [real(F_warm_start(:)); imag(F_warm_start(:))]';
    population = zeros(pop_size, dim);
    
    % --- 1. THE KING ---
    population(1, :) = x_warm; 
    
    % --- 2. ROYAL GUARDS (30%) ---
    num_guards = round(pop_size * 0.3);
    for i = 2:num_guards
        population(i, :) = x_warm + 0.005 * randn(1, dim);
    end
    
    % --- 3. EXPEDITIONARY FORCE (70% - Noise lớn) ---
    for i = (num_guards + 1):pop_size
        population(i, :) = x_warm + 0.15 * randn(1, dim); 
    end
    
    pop_obj = zeros(pop_size, 1); pop_vio = zeros(pop_size, 1);
    for i = 1:pop_size
        [pop_obj(i), pop_vio(i)] = evaluate_candidate(population(i, :), H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
    end
    best_idx = 1; best_sol = population(1, :); best_fit = pop_obj(1);
    
    M_cr = 0.5 * ones(H_mem_size, 1); M_f = 0.5 * ones(H_mem_size, 1); archive = []; k_mem = 1;
    
    for g = 1:generations
        S_cr = []; S_f = []; S_imp = []; new_pop = population;
        sorted_idx = sort_population_deb(pop_obj, pop_vio);
        
        for i = 1:pop_size
            ri = randi(H_mem_size); CRi = max(0, min(1, normrnd(M_cr(ri), 0.1)));
            while true, Fi = M_f(ri) + 0.1 * tan(pi * (rand - 0.5)); if Fi > 0, break; end, end; Fi = min(1, Fi);
            p = rand() * 0.15 + 0.05; top_p = max(1, round(p * pop_size)); pbest_idx = sorted_idx(randi(top_p)); x_pbest = population(pbest_idx, :);
            idxs = [1:i-1, i+1:pop_size]; r1 = population(idxs(randi(length(idxs))), :);
            pop_archive = [population; archive]; r2 = pop_archive(randi(size(pop_archive, 1)), :);
            
            mutant = population(i, :) + Fi * (x_pbest - population(i, :)) + Fi * (r1 - r2);
            cross_points = rand(1, dim) < CRi; if ~any(cross_points), cross_points(randi(dim)) = true; end
            trial = population(i, :); trial(cross_points) = mutant(cross_points);
            
            [trial_obj, trial_vio] = evaluate_candidate(trial, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
            
            if check_deb_wins(trial_obj, trial_vio, pop_obj(i), pop_vio(i))
                new_pop(i, :) = trial;
                if pop_vio(i) > 0, improvement = max(0, pop_vio(i) - trial_vio); else, improvement = max(0, trial_obj - pop_obj(i)); end
                pop_obj(i) = trial_obj; pop_vio(i) = trial_vio;
                S_cr = [S_cr; CRi]; S_f = [S_f; Fi]; S_imp = [S_imp; improvement];
                archive = [archive; population(i, :)]; if size(archive, 1) > pop_size, archive(randi(size(archive, 1)), :) = []; end
                
                [b_obj, b_vio] = evaluate_candidate(best_sol, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
                if check_deb_wins(trial_obj, trial_vio, b_obj, b_vio)
                    best_sol = trial; best_fit = trial_obj;
                end
            end
        end
        population = new_pop;
        if ~isempty(S_cr)
            w = S_imp / (sum(S_imp) + 1e-10); M_cr(k_mem) = sum(w .* S_cr); sum_w_f = sum(w .* S_f); if sum_w_f == 0, sum_w_f = 1e-10; end; M_f(k_mem) = sum(w .* (S_f.^2)) / sum_w_f; k_mem = mod(k_mem, H_mem_size) + 1;
        end
    end
    [F_comm, F_sens] = reconstruct_matrices(best_sol, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle);
end

function idx = sort_population_deb(obj, vio)
    scores = vio * 1e9 - obj; [~, idx] = sort(scores, 'ascend');
end

function [obj, vio] = evaluate_candidate(x, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N)
    [F_comm, F_sens] = reconstruct_matrices(x, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle);
    ssnr = compute_sensing_SNR_handle(sigmasq_radar_rcs, sensing_beamsteering, F_comm, F_sens);
    SINR = compute_SINR_handle(H_comm, F_comm, F_sens, sigmasq_comm);
    min_sinr = min(SINR);
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