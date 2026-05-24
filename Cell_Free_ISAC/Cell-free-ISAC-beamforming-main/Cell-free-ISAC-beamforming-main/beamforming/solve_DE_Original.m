function [F_comm, F_sens, best_fit] = solve_DE_Original(H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle)
% SOLVE_DE_ORIGINAL: Chạy ngẫu nhiên + Phạt nặng

    generations = 200; pop_size = 50; F_weight = 0.5; CR = 0.8; % Gen ít thôi cho nhanh
    [U, M, N] = size(H_comm); dim = 2 * U * M * N; 
    population = -1 + 2 * rand(pop_size, dim); 
    
    fitness = zeros(pop_size, 1);
    for i = 1:pop_size
        fitness(i) = calculate_fitness_soft(population(i, :), H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
    end
    [best_fit, best_idx] = max(fitness); best_sol = population(best_idx, :);
    
    for g = 1:generations
        for i = 1:pop_size
            idxs = [1:i-1, i+1:pop_size]; r = idxs(randperm(length(idxs), 3));
            mutant = population(r(1), :) + F_weight * (population(r(2), :) - population(r(3), :));
            cross_points = rand(1, dim) < CR; if ~any(cross_points), cross_points(randi(dim)) = true; end
            trial = population(i, :); trial(cross_points) = mutant(cross_points);
            
            f_trial = calculate_fitness_soft(trial, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
            
            if f_trial > fitness(i)
                population(i, :) = trial; fitness(i) = f_trial;
                if f_trial > best_fit, best_fit = f_trial; best_sol = trial; end
            end
        end
    end
    [F_comm, F_sens] = reconstruct_matrices(best_sol, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle);
end

function val = calculate_fitness_soft(x, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N)
    [F_comm, F_sens] = reconstruct_matrices(x, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle);
    ssnr = compute_sensing_SNR_handle(sigmasq_radar_rcs, sensing_beamsteering, F_comm, F_sens);
    SINR = compute_SINR_handle(H_comm, F_comm, F_sens, sigmasq_comm);
    min_sinr = min(SINR);
    
    % PHẠT CỰC NẶNG NẾU KHÔNG ĐẠT (ÂM VÔ CỰC)
    if min_sinr < gamma_target
        val = -1e9; 
    else
        val = ssnr;
    end
end

function [F_comm, F_sens] = reconstruct_matrices(x, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle)
    num_elements = U * M * N; F_real = x(1:num_elements); F_imag = x(num_elements+1:end);
    F_comm = reshape(F_real + 1j * F_imag, [U, M, N]);
    pow_per_ap_comm = sum(abs(F_comm).^2, [1, 3]); 
    for m = 1:M, if pow_per_ap_comm(1, m, 1) > 0, F_comm(:, m, :) = F_comm(:, m, :) * sqrt(P_comm / pow_per_ap_comm(1, m, 1)); end; end
    F_sens_base = beam_mapping_handle(H_comm, sensing_beamsteering); 
    curr_p = max(sum(abs(F_sens_base).^2, [1, 3]), [], 'all'); if curr_p > 0, F_sens = F_sens_base * sqrt(P_sensing / curr_p); else, F_sens = F_sens_base; end
end