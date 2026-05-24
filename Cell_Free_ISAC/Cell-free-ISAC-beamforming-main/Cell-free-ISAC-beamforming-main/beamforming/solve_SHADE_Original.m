function [F_comm, F_sens, best_fit] = solve_SHADE_Original(H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle)
% SOLVE_SHADE_ORIGINAL: Thuần chủng + Phạt nặng

    generations = 200; pop_size = 50; H_mem_size = 50;
    [U, M, N] = size(H_comm); dim = 2 * U * M * N;    
    
    population = -1 + 2 * rand(pop_size, dim); 
    fitness = zeros(pop_size, 1);
    for i = 1:pop_size
        fitness(i) = calculate_fitness_soft(population(i, :), H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
    end
    [best_fit, best_idx] = max(fitness); best_sol = population(best_idx, :);
    
    M_cr = 0.5 * ones(H_mem_size, 1); M_f = 0.5 * ones(H_mem_size, 1); archive = []; k_mem = 1;
    
    for g = 1:generations
        S_cr = []; S_f = []; S_imp = []; new_pop = population;
        [~, sorted_idx] = sort(fitness, 'descend'); 
        for i = 1:pop_size
            ri = randi(H_mem_size); CRi = max(0, min(1, normrnd(M_cr(ri), 0.1)));
            while true, Fi = M_f(ri) + 0.1 * tan(pi * (rand - 0.5)); if Fi > 0, break; end, end; Fi = min(1, Fi);
            p = rand() * 0.15 + 0.05; top_p = max(1, round(p * pop_size)); pbest_idx = sorted_idx(randi(top_p)); x_pbest = population(pbest_idx, :);
            idxs = [1:i-1, i+1:pop_size]; r1 = population(idxs(randi(length(idxs))), :);
            pop_archive = [population; archive]; r2 = pop_archive(randi(size(pop_archive, 1)), :);
            mutant = population(i, :) + Fi * (x_pbest - population(i, :)) + Fi * (r1 - r2);
            cross_points = rand(1, dim) < CRi; if ~any(cross_points), cross_points(randi(dim)) = true; end
            trial = population(i, :); trial(cross_points) = mutant(cross_points);
            f_trial = calculate_fitness_soft(trial, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N);
            if f_trial > fitness(i)
                new_pop(i, :) = trial; S_cr = [S_cr; CRi]; S_f = [S_f; Fi]; S_imp = [S_imp; f_trial - fitness(i)];
                archive = [archive; population(i, :)]; if size(archive, 1) > pop_size, archive(randi(size(archive, 1)), :) = []; end
                fitness(i) = f_trial; if f_trial > best_fit, best_fit = f_trial; best_sol = trial; end
            end
        end
        population = new_pop;
        if ~isempty(S_cr)
            w = S_imp / (sum(S_imp) + 1e-10); M_cr(k_mem) = sum(w .* S_cr); sum_w_f = sum(w .* S_f); if sum_w_f == 0, sum_w_f = 1e-10; end; M_f(k_mem) = sum(w .* (S_f.^2)) / sum_w_f; k_mem = mod(k_mem, H_mem_size) + 1;
        end
    end
    [F_comm, F_sens] = reconstruct_matrices(best_sol, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle);
end

function val = calculate_fitness_soft(x, H_comm, sensing_beamsteering, P_comm, P_sensing, sigmasq_comm, gamma_target, sigmasq_radar_rcs, beam_mapping_handle, compute_SINR_handle, compute_sensing_SNR_handle, U, M, N)
    [F_comm, F_sens] = reconstruct_matrices(x, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle);
    ssnr = compute_sensing_SNR_handle(sigmasq_radar_rcs, sensing_beamsteering, F_comm, F_sens);
    SINR = compute_SINR_handle(H_comm, F_comm, F_sens, sigmasq_comm);
    min_sinr = min(SINR);
    if min_sinr < gamma_target, val = -1e9; else, val = ssnr; end
end

function [F_comm, F_sens] = reconstruct_matrices(x, H_comm, sensing_beamsteering, P_comm, P_sensing, U, M, N, beam_mapping_handle)
    num_elements = U * M * N; F_real = x(1:num_elements); F_imag = x(num_elements+1:end);
    F_comm = reshape(F_real + 1j * F_imag, [U, M, N]);
    pow_per_ap_comm = sum(abs(F_comm).^2, [1, 3]); 
    for m = 1:M, if pow_per_ap_comm(1, m, 1) > 0, F_comm(:, m, :) = F_comm(:, m, :) * sqrt(P_comm / pow_per_ap_comm(1, m, 1)); end; end
    F_sens_base = beam_mapping_handle(H_comm, sensing_beamsteering); 
    curr_p = max(sum(abs(F_sens_base).^2, [1, 3]), [], 'all'); if curr_p > 0, F_sens = F_sens_base * sqrt(P_sensing / curr_p); else, F_sens = F_sens_base; end
end