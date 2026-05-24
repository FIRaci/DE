%% Collect the results
clear all;
[res_table, legend_list] = load_results('output', 'dist');

%% Prepare Legends
for i = 3:size(res_table, 2)
    raw_name = legend_list{i};
    
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
        legend_list{i} = raw_name;
    end
end

%% Plot the results
set_default_plot;

% Nhóm kết quả theo khoảng cách
% Cột 1 là Distance
dist_groups = discretize(res_table.('Min UE-Target Distance'), [0:5:50]);
mean_dist = splitapply(@mean, res_table(:, 1), dist_groups);

% Định nghĩa kiểu đường và marker MỞ RỘNG
line_style = ["-", "--", ":", "-.", "-", "--", ":", "-."]; 
marker = ['o', 'x', 's', 'd', '^', 'v', '>', '<'];
colors = lines(8);

figure;
num_cols = size(res_table, 2);

for i = 3:num_cols
    current_algo_idx = floor((i-3)/2) + 1;
    
    ls = line_style(mod(current_algo_idx-1, length(line_style)) + 1);
    mk = marker(mod(current_algo_idx-1, length(marker)) + 1);
    clr = colors(mod(current_algo_idx-1, length(colors)) + 1, :);

    data = res_table(:, i);
    mean_val = splitapply(@mean, data, dist_groups);

    if rem(i, 2) == 1
        % COMM PLOT
        subplot(2, 1, 1);
        title('Communication Performance vs Distance');
        ylabel('Min Comm SINR');
        xlabel('Target-Closest UE Distance (m)');
        grid on; hold on;
        
        plot(mean_dist, mean_val, ...
            'LineStyle', ls, ...
            'Marker', mk, ...
            'Color', clr, ...
            'DisplayName', legend_list{i}, ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    else
        % SENSING PLOT
        subplot(2, 1, 2);
        title('Sensing Performance vs Distance');
        ylabel('Target SNR')
        xlabel('Target-Closest UE Distance (m)');
        grid on; hold on;
        
        plot(mean_dist, mean_val, ...
            'LineStyle', ls, ...
            'Marker', mk, ...
            'Color', clr, ...
            'DisplayName', legend_list{i}, ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    end
end

subplot(2, 1, 1); legend('Location', 'bestoutside');
subplot(2, 1, 2); legend('Location', 'bestoutside');