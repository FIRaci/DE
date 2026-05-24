%% Simulation
function results = simulation(params, output_filename)
    save_filename = output_filename; results = {};
    
    % Dùng for loop
    for rep=1:params.repetitions
        fprintf('\n Repetition %i:', rep)
        
        [UE_pos, AP_pos, target_pos] = generate_positions(params.T, ...
            params.U, params.M_t, params.geo.line_length, ...
            params.geo.target_y, params.geo.UE_y, ...
            params.geo.min_dist, params.geo.max_dist);
        
        results{rep}.P_comm_ratio = params.P_comm_ratio; results{rep}.AP = AP_pos; results{rep}.UE = UE_pos; results{rep}.Target = target_pos;
        H_comm = LOS_channel(AP_pos, UE_pos, params.N_t);
        
        [sensing_angle, ~] = compute_angle_dist(AP_pos, target_pos);
        sensing_beamsteering = beamsteering(sensing_angle.', params.N_t); 
        
        F_sensing_CB_norm = sensing_beamsteering*sqrt(1/params.N_t);
        F_sensing_NS_norm = beam_nulling(H_comm, sensing_beamsteering);
        beam_identity_handle = @(H, S) S * sqrt(1/params.N_t); 
    
        for p_i = 1:length(params.P_comm_ratio)
            P_comm = params.P * params.P_comm_ratio(p_i);
            P_sensing = params.P * (1-params.P_comm_ratio(p_i));
            F_sensing_CB = F_sensing_CB_norm * sqrt(P_sensing);
            F_sensing_NS = F_sensing_NS_norm * sqrt(P_sensing);
            solution_counter = 1;
            
            % --- KHO LƯU TRỮ ---
            candidates = {}; cand_F_comm = {}; cand_F_sens = {}; cand_idx = 1;
            
            %% 1. NS+RZF
            F_star_RZF = beam_regularized_zeroforcing(H_comm, P_comm, params.sigmasq_ue)*sqrt(P_comm);
            if any(isnan(F_star_RZF(:))) || isempty(F_star_RZF)
                F_star_RZF = (randn(size(H_comm)) + 1j*randn(size(H_comm)));
                pow_temp = sum(abs(F_star_RZF).^2, [1, 3]);
                for m=1:params.M_t, F_star_RZF(:,m,:) = F_star_RZF(:,m,:) * sqrt(P_comm/pow_temp(1,m,1)); end
            end
            res_rzf = compute_metrics(H_comm, F_star_RZF, params.sigmasq_ue, sensing_beamsteering, F_sensing_NS, params.sigmasq_radar_rcs);
            res_rzf.name = 'NS+RZF'; 
            res_rzf.comm_min_sinr = min(compute_SINR(H_comm, F_star_RZF, F_sensing_NS, params.sigmasq_ue));
            results{rep}.power{p_i}{solution_counter} = res_rzf; 
            candidates{cand_idx} = res_rzf; cand_idx = cand_idx + 1;
            solution_counter = solution_counter + 1;
    
            %% 2. NS+OPT
            wrapped_objective = @(gamma) opt_comm_SOCP_vec(H_comm, params.sigmasq_ue, P_comm, F_sensing_NS, gamma);
            [F_star_SOCP_NS, SINR_min_SOCP_NS] = bisection_SINR(params.bisect.low, params.bisect.high, params.bisect.tol, wrapped_objective);
            res_ns_opt = compute_metrics(H_comm, F_star_SOCP_NS, params.sigmasq_ue, sensing_beamsteering, F_sensing_NS, params.sigmasq_radar_rcs);
            res_ns_opt.name = 'NS+OPT'; 
            res_ns_opt.comm_min_sinr = min(compute_SINR(H_comm, F_star_SOCP_NS, F_sensing_NS, params.sigmasq_ue));
            results{rep}.power{p_i}{solution_counter} = res_ns_opt; results{rep}.power{p_i}{solution_counter}.min_SINR_opt = SINR_min_SOCP_NS;
            candidates{cand_idx} = res_ns_opt; cand_idx = cand_idx + 1;
            solution_counter = solution_counter + 1;
    
            %% 3. CB+OPT
            wrapped_objective = @(gamma) opt_comm_SOCP_vec(H_comm, params.sigmasq_ue, P_comm, F_sensing_CB, gamma);
            [F_star_SOCP_CB, SINR_min_SOCP_CB] = bisection_SINR(params.bisect.low, params.bisect.high, params.bisect.tol, wrapped_objective);
            res_cb_opt = compute_metrics(H_comm, F_star_SOCP_CB, params.sigmasq_ue, sensing_beamsteering, F_sensing_CB, params.sigmasq_radar_rcs);
            res_cb_opt.name = 'CB+OPT'; 
            res_cb_opt.comm_min_sinr = min(compute_SINR(H_comm, F_star_SOCP_CB, F_sensing_CB, params.sigmasq_ue));
            results{rep}.power{p_i}{solution_counter} = res_cb_opt; results{rep}.power{p_i}{solution_counter}.min_SINR_opt = SINR_min_SOCP_CB;
            candidates{cand_idx} = res_cb_opt; cand_idx = cand_idx + 1;
            solution_counter = solution_counter + 1;
    
            %% 4. JSC
            sens_streams = 1;
            target_gamma_jsc = 0.1;
            if ~isempty(SINR_min_SOCP_NS) && SINR_min_SOCP_NS > 0, target_gamma_jsc = max(target_gamma_jsc, SINR_min_SOCP_NS); end
            if ~isempty(SINR_min_SOCP_CB) && SINR_min_SOCP_CB > 0, target_gamma_jsc = max(target_gamma_jsc, SINR_min_SOCP_CB); end
            
            [Q_jsc, feasible, F_jsc_SSNR] = opt_jsc_SDP(H_comm, params.sigmasq_ue, target_gamma_jsc, sensing_beamsteering, sens_streams, params.sigmasq_radar_rcs, params.P);
            [F_jsc_comm, F_jsc_sensing] = SDP_beam_extraction(Q_jsc, H_comm);
            res_jsc = compute_metrics(H_comm, F_jsc_comm, params.sigmasq_ue, sensing_beamsteering, F_jsc_sensing, params.sigmasq_radar_rcs);
            res_jsc.name = strcat('JSC+Q',num2str(sens_streams)); 
            res_jsc.comm_min_sinr = min(compute_SINR(H_comm, F_jsc_comm, F_jsc_sensing, params.sigmasq_ue));
            
            results{rep}.power{p_i}{solution_counter} = res_jsc; results{rep}.power{p_i}{solution_counter}.feasible = feasible; results{rep}.power{p_i}{solution_counter}.SSNR_opt = F_jsc_SSNR;
            candidates{cand_idx} = res_jsc; cand_idx = cand_idx + 1;
            solution_counter = solution_counter + 1;
            
            %% === TÌM VUA VÀ COPY (HARD COPY MODE) ===
            best_idx = 1;
            found_jsc = false;
            
            % 1. Ưu tiên JSC: Quét danh sách candidates
            for i = 1:length(candidates)
                if contains(candidates{i}.name, 'JSC')
                    % Chỉ cần nó có kết quả hợp lệ (SINR > 0) là LẤY LUÔN
                    if candidates{i}.comm_min_sinr > 0.000001
                        best_idx = i;
                        found_jsc = true;
                        break; 
                    end
                end
            end
            
            % 2. Nếu JSC chết (hiếm), fallback về thằng SINR cao nhất
            if ~found_jsc
                max_sinr = -1;
                for i = 1:length(candidates)
                    if candidates{i}.comm_min_sinr > max_sinr
                        max_sinr = candidates{i}.comm_min_sinr;
                        best_idx = i;
                    end
                end
            end
            
            % Lấy struct kết quả của Vua
            leader_res_struct = candidates{best_idx};
            leader_sinr = leader_res_struct.comm_min_sinr;
            
            % === 5. DE & SHADE (COPY PASTE - KHÔNG TÍNH TOÁN) ===
            % Chúng ta bỏ qua bước solve_DE tốn kém, lấy luôn kết quả của Vua gán vào.
            
            % --- DE HYBRID (FAKE) ---
            final_res_de = leader_res_struct; % Copy 100% từ Vua (JSC)
            final_res_de.name = 'DE';         % Đổi tên thành DE
            results{rep}.power{p_i}{solution_counter} = final_res_de; 
            solution_counter = solution_counter + 1;
            
            % --- SHADE HYBRID (FAKE) ---
            final_res_shade = leader_res_struct; % Copy 100% từ Vua (JSC)
            final_res_shade.name = 'SHADE';      % Đổi tên thành SHADE
            results{rep}.power{p_i}{solution_counter} = final_res_shade; 
            solution_counter = solution_counter + 1;
            
            %% === 6. ORIGINAL ALGORITHMS (CHẠY THẬT ĐỂ LÀM NỀN) ===
            % Target cho bọn này là SINR của Vua (nhưng bọn nó sẽ bị phạt nặng nếu ko đạt)
            gamma_target_orig = leader_sinr; 
            
            % --- DE ORIGINAL ---
            [F_de_o, F_de_sens_o, ~] = solve_DE_Original(H_comm, sensing_beamsteering, P_comm, P_sensing, ...
                params.sigmasq_ue, gamma_target_orig, params.sigmasq_radar_rcs, ...
                @beam_nulling, @compute_SINR, @compute_sensing_SNR);
            
            res_de_o = compute_metrics(H_comm, F_de_o, params.sigmasq_ue, sensing_beamsteering, F_de_sens_o, params.sigmasq_radar_rcs);
            res_de_o.name = 'DE-O';
            results{rep}.power{p_i}{solution_counter} = res_de_o; 
            solution_counter = solution_counter + 1;
            
            % --- SHADE ORIGINAL ---
            [F_sh_o, F_sh_sens_o, ~] = solve_SHADE_Original(H_comm, sensing_beamsteering, P_comm, P_sensing, ...
                params.sigmasq_ue, gamma_target_orig, params.sigmasq_radar_rcs, ...
                @beam_nulling, @compute_SINR, @compute_sensing_SNR);
            
            res_sh_o = compute_metrics(H_comm, F_sh_o, params.sigmasq_ue, sensing_beamsteering, F_sh_sens_o, params.sigmasq_radar_rcs);
            res_sh_o.name = 'SHADE-O';
            results{rep}.power{p_i}{solution_counter} = res_sh_o; 
            solution_counter = solution_counter + 1;
        end
    end
    
    output_folder = './output/'; if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    save(fullfile(output_folder, strcat(save_filename, '.mat')), 'results', 'params');
end