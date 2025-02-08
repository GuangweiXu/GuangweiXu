function analyze_network_results_inh(base_path)
    % 获取所有i_stren值对应的文件夹
    i_stren_values = 0:0.1:1;
    all_results = cell(length(i_stren_values), 1);
    
    % 定义神经元数量
    Ne1=5; Ni1=2; Ne2=8; Ni2=3; Ne3=9; Ni3=3;exc3=Ne3;inh3=Ni3;
    % 计算第三层兴奋性和抑制性神经元的索引
    exc3_start = Ne1 + Ni1 + Ne2 + Ni2 + 1;
    exc3_end = exc3_start + Ne3 - 1;
    exc3_indices = exc3_start:exc3_end;
    inh3_start = Ne1 + Ni1 + Ne2 + Ni2 + Ne3 + 1;
    inh3_end = inh3_start + Ni3 - 1;
    inh3_indices = inh3_start:inh3_end;    
    % 添加同步性分析的存储变量
    coherence_values = zeros(length(i_stren_values), 1);
    sync_index = zeros(length(i_stren_values), 1);
    mean_xcorr = zeros(length(i_stren_values), 1);
    % 加载所有结果
    for i = 1:length(i_stren_values)
        i_stren = i_stren_values(i);
        folder_path = fullfile(base_path, ['neuron_' num2str(i_stren)]);
        file_path = fullfile(folder_path, ['avg_firing_rates_' num2str(i_stren) '.mat']);
        
        if exist(file_path, 'file')
            loaded_data = load(file_path);
            all_results{i} = loaded_data;
        else
            warning(['File not found: ' file_path]);
        end
    end
    
    % 初始化存储数组
    num_neurons = 500; 
    num_sweeps = 1; 
    num_exc3_pairs = nchoosek(length(exc3_indices), 2);
    
    % 存储平均发放率和噪音相关性
    mean_inh3_rates = zeros(length(i_stren_values), length(inh3_indices));
    std_inh3_rates = zeros(length(i_stren_values), length(inh3_indices));
    mean_exc3_rates = zeros(length(i_stren_values), exc3);
    noise_correlations = zeros(length(i_stren_values), num_exc3_pairs);
    mean_noise_corr = zeros(length(i_stren_values), 1);
    std_noise_corr = zeros(length(i_stren_values), 1);
    % 在初始化存储数组部分添加
    mean_synchrony = zeros(length(i_stren_values), 1);
    std_synchrony = zeros(length(i_stren_values), 1);
     % 添加共同输入分析的存储变量
    common_input_strength = zeros(length(i_stren_values), 1);
    common_input_correlation = zeros(length(i_stren_values), 1);
    neuron1_influence = zeros(length(i_stren_values), exc3);
    % 处理数据
    for i = 1:length(i_stren_values)
        current_results = all_results{i};
        
        % 收集第三层抑制性神经元的发放率
        for j = 1:length(inh3_indices)
            neuron_idx = inh3_indices(j);
            rate_field = ['avg_firing_rate' num2str(neuron_idx)];
            if isfield(current_results, rate_field)
                rates = current_results.(rate_field);
                mean_inh3_rates(i,j) = mean(rates(:));
                std_inh3_rates(i,j) = std(rates(:));
            end
        end
        
        % 收集inh3神经元的发放率
        for j = 1:inh3
            neuron_idx = inh3_indices(j);
            rate_field = ['avg_firing_rate' num2str(neuron_idx)];
            if isfield(current_results, rate_field)
                rates = current_results.(rate_field);
                mean_inh3_rates(i,j) = mean(rates(:));
            end
        end
        % 收集exc3神经元的发放率
        for j = 1:exc3
            neuron_idx = exc3_indices(j);
            rate_field = ['avg_firing_rate' num2str(neuron_idx)];
            if isfield(current_results, rate_field)
                rates = current_results.(rate_field);
                mean_exc3_rates(i,j) = mean(rates(:));
            end
        end
        
        pair_idx = 1;
        temp_noise_corr = zeros(num_neurons, num_sweeps, num_exc3_pairs);
        
        for n1_idx = 1:length(exc3_indices)
            for n2_idx = (n1_idx+1):length(exc3_indices)
                n1 = exc3_indices(n1_idx);
                n2 = exc3_indices(n2_idx);              
                
                for neuron = 1:num_neurons
                    for sweep = 1:num_sweeps                        
                          
                        if ~isempty(current_results.fine_grained_rates{neuron, sweep})
                            % 获取两个神经元的发放率时间序列
                            rates1 = current_results.fine_grained_rates{neuron, sweep}(n1, :);
                            rates2 = current_results.fine_grained_rates{neuron, sweep}(n2, :);                            
                            % 检查数据有效性
                            if ~isempty(rates1) && ~isempty(rates2) && ...
                               ~any(isnan(rates1)) && ~any(isnan(rates2)) && ...
                               ~all(rates1 == 0) && ~all(rates2 == 0)  % 确保不是全零序列
                                 % 计算标准化的噪声（z-score）
                                std1 = std(rates1, 'omitnan');
                                std2 = std(rates2, 'omitnan');
                                % 移除平均响应以获得噪音
                                noise1 = (rates1 - mean(rates1, 'omitnan'))/std1;
                                noise2 = (rates2 - mean(rates2, 'omitnan'))/std2;
                                
                                % 计算相关系数
                                if std(noise1) > 0 && std(noise2) > 0  % 确保有变化
                                    temp_corr = corrcoef(noise1, noise2, 'rows', 'complete');
                                    if ~isempty(temp_corr) && ~isnan(temp_corr(1,2))
                                        temp_noise_corr(neuron, sweep, pair_idx) = temp_corr(1,2);
                                    end
                                end
                            end                        
                        end
                    end
                end
                pair_idx = pair_idx + 1;
            end
        end
        
        % 计算平均噪音相关性，忽略NaN值
        for pair = 1:num_exc3_pairs
            valid_corrs = temp_noise_corr(:,:,pair);
            valid_corrs = valid_corrs(~isnan(valid_corrs));
            if ~isempty(valid_corrs)
                noise_correlations(i, pair) = mean(valid_corrs, 'omitnan');
            else
                noise_correlations(i, pair) = NaN;
            end
        end
        
        % 计算总体平均和标准差，忽略NaN值
        valid_pairs = ~isnan(noise_correlations(i, :));
        if any(valid_pairs)
            mean_noise_corr(i) = mean(noise_correlations(i, valid_pairs), 'omitnan');
            std_noise_corr(i) = std(noise_correlations(i, valid_pairs), 'omitnan');
        else
            mean_noise_corr(i) = NaN;
            std_noise_corr(i) = NaN;
        end
        % 处理同步性度量
    all_sync_values = [];
    if isfield(current_results, 'synchrony_measure')
        for neuron = 1:num_neurons
            for sweep = 1:num_sweeps
                if ~isempty(current_results.synchrony_measure{neuron, sweep})
                    % 确保我们获取的是一维数组
                    current_sync = current_results.synchrony_measure{neuron, sweep};
                    if ~isempty(current_sync)
                        % 如果current_sync是多维的，取平均值
                        if numel(current_sync) > 1
                            current_sync = mean(current_sync(:), 'omitnan');
                        end
                        all_sync_values = [all_sync_values; current_sync];
                    end
                end
            end
        end
    end
     % 计算平均同步性和标准差
        if ~isempty(all_sync_values)
            mean_synchrony(i) = mean(all_sync_values, 'omitnan');
            std_synchrony(i) = std(all_sync_values, 'omitnan');
        else
            mean_synchrony(i) = NaN;
            std_synchrony(i) = NaN;
        end

        % 分析神经元1的活动
                neuron1_rate_field = 'avg_firing_rate1';
                if isfield(current_results, neuron1_rate_field)
                    neuron1_rates = current_results.(neuron1_rate_field);
                    
                    % 计算神经元1与每个exc3神经元的相关性
                    for exc3_idx = 1:length(exc3_indices)
                        target_neuron = exc3_indices(exc3_idx);
                        
                        for neuron = 1:num_neurons
                            for sweep = 1:num_sweeps
                                if ~isempty(current_results.fine_grained_rates{neuron, sweep})
                                    % 获取神经元1和目标exc3神经元的发放率时间序列
                                    rates1 = current_results.fine_grained_rates{neuron, sweep}(1, :);  % 神经元1
                                    rates2 = current_results.fine_grained_rates{neuron, sweep}(target_neuron, :);
                                    
                                    % 检查数据有效性
                                    if ~isempty(rates1) && ~isempty(rates2) && ...
                                       ~any(isnan(rates1)) && ~any(isnan(rates2)) && ...
                                       ~all(rates1 == 0) && ~all(rates2 == 0)  % 确保不是全零序列
                                        
                                        % 计算标准化的噪声（z-score）
                                        std1 = std(rates1, 'omitnan');
                                        std2 = std(rates2, 'omitnan');
                                        
                                        % 移除平均响应以获得噪音
                                        noise1 = (rates1 - mean(rates1, 'omitnan'))/std1;
                                        noise2 = (rates2 - mean(rates2, 'omitnan'))/std2;
                                        
                                        % 计算相关系数
                                        if std(noise1) > 0 && std(noise2) > 0  % 确保有变化
                                            temp_corr = corrcoef(noise1, noise2, 'rows', 'complete');
                                            if ~isempty(temp_corr) && ~isnan(temp_corr(1,2))
                                                temp_neuron1_corr(neuron, sweep, exc3_idx) = temp_corr(1,2);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    % 计算每个exc3神经元与神经元1的平均相关性
                    neuron1_influence = zeros(1, length(exc3_indices));
                    for exc3_idx = 1:length(exc3_indices)
                        valid_corrs = temp_neuron1_corr(:,:,exc3_idx);
                        valid_corrs = valid_corrs(~isnan(valid_corrs));
                        if ~isempty(valid_corrs)
                            neuron1_influence(exc3_idx) = mean(valid_corrs, 'omitnan');
                        else
                            neuron1_influence(exc3_idx) = NaN;
                        end
                    end
                    
                    % 存储结果
                    common_input_correlation(i) = mean(neuron1_influence, 'omitnan');
                    neuron1_influence_matrix(i,:) = neuron1_influence;  % 存储详细的影响模式
                    
                    % 计算共同输入强度（使用神经元1的平均发放率）
                    if isfield(current_results, 'avg_firing_rate1')
                        common_input_strength(i) = mean(current_results.avg_firing_rate1(:), 'omitnan');
                    end
                end


    end
    % 确保所有数据都是列向量
mean_synchrony = mean_synchrony(:);  % 转换为列向量
mean_noise_corr = mean_noise_corr(:);  % 转换为列向量
    % % 调用绘图函数
    % plot_analysis(i_stren_values, mean_inh3_rates, std_inh3_rates, ...
    %     mean_exc3_rates, inh3_indices, mean_noise_corr, std_noise_corr, ...
    %     noise_correlations);


% 计算同步性指标
        all_exc3_rates = zeros(length(exc3_indices), size(current_results.fine_grained_rates{1,1}, 2));
        for j = 1:length(exc3_indices)
            all_exc3_rates(j,:) = current_results.fine_grained_rates{1,1}(exc3_indices(j),:);
        end
        
       

        % 调用修改后的绘图函数
    plot_analysis(i_stren_values, mean_inh3_rates, std_inh3_rates, ...
    mean_exc3_rates, inh3_indices, mean_noise_corr, std_noise_corr, ...
    noise_correlations, mean_synchrony, std_synchrony, ...
    common_input_strength, common_input_correlation, neuron1_influence);

end



function sync_idx = calculate_sync_index(spike_rates)
    % 计算瞬时群体同步性指数
    T = size(spike_rates, 2);
    N = size(spike_rates, 1);
    
    % 计算归一化的群体活动方差
    pop_rate = mean(spike_rates, 1);
    pop_var = var(pop_rate);
    individual_var = mean(var(spike_rates, [], 2));
    
    sync_idx = (pop_var - individual_var/N)/(individual_var*(N-1)/N);
end

function mean_cc = calculate_mean_xcorr(spike_rates)
    % 计算平均互相关
    N = size(spike_rates, 1);
    cc_sum = 0;
    count = 0;
    
    for i = 1:N
        for j = (i+1):N
            [c,~] = xcorr(spike_rates(i,:), spike_rates(j,:), 'coeff');
            cc_sum = cc_sum + max(c);
            count = count + 1;
        end
    end
    
    mean_cc = cc_sum/count;
end

function plot_analysis(i_stren_values, mean_inh3_rates, std_inh3_rates, ...
    mean_exc3_rates, inh3_indices, mean_noise_corr, std_noise_corr, ...
    noise_correlations, mean_synchrony, std_synchrony, ...
    common_input_strength, common_input_correlation, neuron1_influence)
    
    % 图1：INH3抑制强度与抑制性神经元活动的关系
    figure('Position', [100, 100, 1200, 400]);
    
    subplot(1,3,1);
    mean_rates = mean(mean_exc3_rates, 2, 'omitnan');
    std_rates = std(mean_exc3_rates, 0, 2, 'omitnan');
    
    errorbar(i_stren_values, mean_rates, std_rates, 'b-', 'LineWidth', 2);
    hold on;
    plot(i_stren_values, mean_rates, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    
    xlabel('PV Inhibition Strength');
    ylabel('Mean Firing Rate (Hz)');
    title('Layer 3 Excitatory Neurons Response');
    grid on;
    
    subplot(1,3,2);
    imagesc(i_stren_values, 1:size(mean_exc3_rates,2), mean_exc3_rates');
    colormap('jet');
    colorbar;
    xlabel('PV Inhibition Strength');
    ylabel('Exc Neuron Index');
    title('Individual Exc Neuron Responses');
    
    subplot(1,3,3);
    plot(i_stren_values, mean_inh3_rates, 'LineWidth', 2);
    xlabel('PV Inhibition Strength');
    ylabel('Firing Rate (Hz)');
    title('inhibitory Neuron Activity');
    grid on;
    
    saveas(gcf, 'inhibitory_neuron_inhibition_analysis.png');
    
    % 图2：噪音相关性分析
    figure('Position', [100, 500, 1200, 400]);
    
    % 子图1：平均噪音相关性随INH3抑制强度的变化
    subplot(1,3,1);
    
    % 处理NaN值
    valid_idx = ~isnan(mean_noise_corr) & ~isnan(std_noise_corr);
    valid_x = i_stren_values(valid_idx);
    valid_mean = mean_noise_corr(valid_idx);
    valid_std = std_noise_corr(valid_idx);
    
    errorbar(valid_x, valid_mean, valid_std, 'b-', 'LineWidth', 2);
    hold on;
    plot(valid_x, valid_mean, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    
    % 只对有效数据进行拟合
    if length(valid_x) > 2  % 确保有足够的点进行拟合
        try
            [fit_curve, gof] = fit(valid_x', valid_mean', 'smoothingspline');
            plot(valid_x, fit_curve(valid_x), 'r--', 'LineWidth', 1.5);
            % 添加R²值
            text(0.7*max(valid_x), 0.9*max(valid_mean), ...
                ['R^2 = ' num2str(gof.rsquare, '%.3f')], ...
                'FontSize', 12);
        catch
            warning('无法完成曲线拟合，可能是数据点不足或存在其他问题');
        end
    end
    
    xlabel('PV Inhibition Strength');
    ylabel('Mean Noise Correlation');
    title('Layer 3 Excitatory Neurons Noise Correlation');
    grid on;
    
    % 子图2：噪音相关性矩阵热图
    subplot(1,3,2);
    % 处理热图数据中的NaN值
    noise_corr_plot = noise_correlations;
    noise_corr_plot(isnan(noise_corr_plot)) = 0;  % 或者使用其他合适的值
    
    imagesc(i_stren_values, 1:size(noise_corr_plot,2), noise_corr_plot');

    colormap('jet');
    colorbar;
    xlabel('PV Inhibition Strength');
    ylabel('Neuron Pair Index');
    title('Pairwise Noise Correlations');
    
    % 子图3：噪音相关性分布
    subplot(1,3,3);
    valid_corr = noise_correlations(~isnan(noise_correlations));
    if ~isempty(valid_corr)
        histogram(valid_corr, 50, 'Normalization', 'probability');
    end
    xlabel('Noise Correlation Value');
    ylabel('Probability');
    title('Distribution of Noise Correlations');
    grid on;
    
    saveas(gcf, 'noise_correlation_analysis.png');
    
    % 输出统计信息
    fprintf('\nAnalysis Results:\n');
    fprintf('Mean noise correlation: %.3f ± %.3f\n', ...
        mean(mean_noise_corr, 'omitnan'), std(mean_noise_corr, 'omitnan'));
    fprintf('Correlation range: [%.3f, %.3f]\n', ...
        min(valid_corr), max(valid_corr));
    
    % INH3强度与噪音相关性的相关系数
    if length(valid_x) > 1
        inh3_corr = corrcoef(valid_x, valid_mean, 'rows', 'complete');
        fprintf('Correlation between inhibitory strength and noise correlation: %.3f\n', inh3_corr(1,2));
        if exist('gof', 'var')
            fprintf('R-squared of noise correlation fit: %.3f\n', gof.rsquare);
        end
    end

figure('Position', [100, 500, 1200, 400]);
i_stren_values=i_stren_values(:);
% 子图1：细粒度同步性随时间的变化
subplot(1,3,1);
errorbar(i_stren_values, mean_synchrony, std_synchrony, 'k-', 'LineWidth', 2);
hold on;
plot(i_stren_values, mean_synchrony, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
xlabel('PV Inhibition Strength');
ylabel('Fine-grained Synchrony');
title('Time-resolved Synchrony Measure');
grid on;

% 子图2：同步性与噪音相关性的关系
subplot(1,3,2);
valid_idx = ~isnan(mean_synchrony(:)) & ~isnan(mean_noise_corr(:));
scatter(mean_synchrony(valid_idx), mean_noise_corr(valid_idx), 50, 'filled');
xlabel('Synchrony Measure');
ylabel('Noise Correlation');
title('Synchrony vs. Noise Correlation');
grid on;

% 添加相关系数
if sum(valid_idx) > 1
    [r, p] = corrcoef(mean_synchrony(valid_idx), mean_noise_corr(valid_idx));
    text(0.1, 0.9, sprintf('r = %.3f\np = %.3f', r(1,2), p(1,2)), ...
        'Units', 'normalized');
end

% 子图3：同步性分布
subplot(1,3,3);
histogram(mean_synchrony(~isnan(mean_synchrony)), 20, 'Normalization', 'probability');
xlabel('Synchrony Value');
ylabel('Probability');
title('Distribution of Synchrony Values');
grid on;

saveas(gcf, 'detailed_synchrony_analysis.png');

% 添加统计信息输出
fprintf('\nDetailed Synchrony Analysis Results:\n');
fprintf('Mean fine-grained synchrony: %.3f ± %.3f\n', ...
    mean(mean_synchrony, 'omitnan'), std(mean_synchrony, 'omitnan'));
if sum(valid_idx) > 1
    fprintf('Correlation between synchrony and noise correlation: r = %.3f, p = %.3f\n', ...
        r(1,2), p(1,2));
end

% 添加共同输入分析图
    figure('Position', [100, 500, 1200, 400]);
    
    % 子图1：共同输入强度与抑制强度的关系
    subplot(1,3,1);
    % plot(i_stren_values, common_input_strength, 'b-o', 'LineWidth', 2);
    % hold on;
    plot(i_stren_values, common_input_correlation, 'r-o', 'LineWidth', 2);
    xlabel('PV Inhibition Strength');
    ylabel('Correlation');
    legend('Neuron 1 Rate', 'Common Input Correlation');
    title('Common Input Analysis');
    grid on;
    
    % 子图2：神经元1对各个exc3神经元的影响
    subplot(1,3,2);
    imagesc(i_stren_values, 1:size(neuron1_influence,2), neuron1_influence');
    colormap('jet');
    colorbar;
    xlabel('PV Inhibition Strength');
    ylabel('Exc3 Neuron Index');
    title('Neuron 1 Influence on Exc3');
    
    % 子图3：相关性分解
    subplot(1,3,3);
    plot(i_stren_values, mean_noise_corr, 'k-', 'LineWidth', 2);
    hold on;
    plot(i_stren_values, common_input_correlation, 'r--', 'LineWidth', 2);
    plot(i_stren_values, mean_noise_corr - common_input_correlation, 'b--', 'LineWidth', 2);
    xlabel('PV Inhibition Strength');
    ylabel('Correlation');
    legend('Total Correlation', 'Common Input', 'Residual');
    title('Correlation Decomposition');
    grid on;
    
    saveas(gcf, 'common_input_analysis.png');
    
    % 添加统计信息输出
    fprintf('\nCommon Input Analysis:\n');
    fprintf('Mean common input correlation: %.3f ± %.3f\n', ...
        mean(common_input_correlation, 'omitnan'), std(common_input_correlation, 'omitnan'));
    
    % 计算共同输入与总相关性的关系
    valid_idx = ~isnan(common_input_correlation) & ~isnan(mean_noise_corr);
    if sum(valid_idx) > 1
        [r, p] = corrcoef(common_input_correlation(valid_idx), ...
            mean_noise_corr(valid_idx));
        fprintf('Correlation between common input and total correlation: r = %.3f, p = %.3f\n', ...
            r(1,2), p(1,2));
    end
    
    % 分析抑制对共同输入效应的影响
    common_input_reduction = (common_input_correlation(1) - ...
        common_input_correlation(end)) / common_input_correlation(1) * 100;
    fprintf('Common input correlation reduction: %.1f%%\n', common_input_reduction);

end