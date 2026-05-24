%% Collect the results
clear all;
% Đảm bảo rằng hàm load_results có thể tìm thấy file output
[res_table, legend_list] = load_results('output', 'power');

%% Prepare Legends
for i = 3:size(res_table, 2)
    raw_name = legend_list{i};
    
    % Xử lý tên hiển thị cho đẹp
    if strcmp(raw_name, 'DE')
        legend_list{i} = 'DE (Hybrid - Parasite)';
    elseif strcmp(raw_name, 'SHADE')
        legend_list{i} = 'SHADE (Hybrid - Parasite)';
    elseif strcmp(raw_name, 'DE-O')
        legend_list{i} = 'DE (Original - Random)';
    elseif strcmp(raw_name, 'SHADE-O')
        legend_list{i} = 'SHADE (Original - Random)';
    elseif contains(raw_name, 'JSC')
        legend_list{i} = 'JSC Beam Optimization';
    elseif contains(raw_name, 'NS+RZF')
        legend_list{i} = 'NS Sensing - RZF Comm';
    elseif contains(raw_name, 'NS+OPT')
        legend_list{i} = 'NS Sensing - Opt Comm';
    elseif contains(raw_name, 'CB+OPT')
        legend_list{i} = 'CB Sensing - Opt Comm';
    else
        % Fallback cho các trường hợp cũ
        txt = strsplit(raw_name, '-');
        legend_list{i} = raw_name; 
    end
end

%% Plot the results
set_default_plot;

% Nhóm kết quả theo tỷ lệ công suất (Power-ratio)
res_mat_mean = grpstats(res_table, 'Power-ratio');
power_ratio = res_mat_mean{:, 1};
res_mat_mean = res_mat_mean{:, 2:end}; % Remove power ratio column

% Định nghĩa kiểu đường và marker MỞ RỘNG (cho 8 thuật toán)
line_style = ["-", "--", ":", "-.", "-", "--", ":", "-."]; 
marker = ['o', 'x', 's', 'd', '^', 'v', '>', '<'];
colors = lines(8); % Lấy bảng màu chuẩn 8 màu

style_counter = 0;
figure;

num_algos = size(res_mat_mean, 2);

for i = 3:num_algos
    % Logic chọn subplot: Cột lẻ là Comm, Cột chẵn là Sensing (do cấu trúc bảng)
    % Cột 1: Power-ratio (đã bỏ), Cột 2: GroupCount (đã bỏ) -> Bắt đầu từ cột 3
    
    current_algo_idx = floor((i-3)/2) + 1; % Index thuật toán (1, 2, 3...)
    
    % Chọn style
    ls = line_style(mod(current_algo_idx-1, length(line_style)) + 1);
    mk = marker(mod(current_algo_idx-1, length(marker)) + 1);
    clr = colors(mod(current_algo_idx-1, length(colors)) + 1, :);

    if rem(i, 2) == 1
        % COMMUNICATION PLOT
        subplot(2, 1, 1);
        title('Communication Performance (SINR)');
        ylabel('Min Comm SINR');
        xlabel('Communication Power Ratio (\rho)');
        grid on; hold on;
        
        plot(power_ratio, res_mat_mean(:, i), ...
            'LineStyle', ls, ...
            'Marker', mk, ...
            'Color', clr, ...
            'DisplayName', legend_list{i}, ...
            'LineWidth', 1.5, 'MarkerSize', 6);
            
    else
        % SENSING PLOT
        subplot(2, 1, 2);
        title('Sensing Performance (SNR)');
        ylabel('Target SNR')
        xlabel('Communication Power Ratio (\rho)');
        grid on; hold on;
        
        plot(power_ratio, res_mat_mean(:, i), ...
            'LineStyle', ls, ...
            'Marker', mk, ...
            'Color', clr, ...
            'DisplayName', legend_list{i}, ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    end
end

subplot(2, 1, 1); legend('Location', 'bestoutside');
subplot(2, 1, 2); legend('Location', 'bestoutside');