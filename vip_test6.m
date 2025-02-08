% vip_test_7.m
clear,clc,close all
global avg_firing_rate8 avg_firing_rate9 avg_firing_rate10 
global_time = 0;

% Network architecture parameters
Ne1=5; Ni1=2; Ne2=8; Ni2=3; Ne3=9; Ni3=3; Nvip=4; 
N = Ne1+Ni1+Ne2+Ni2+Ne3+Ni3+Nvip;
vip = (Ne1+Ni1+Ne2+Ni2+Ne3+Ni3+1):N;

% Define neuron groups
group1 = 1:(Ne1+Ni1);
group2 = (Ne1+Ni1+1):(Ne1+Ni1+Ne2+Ni2);
group3 = (Ne1+Ni1+Ne2+Ni2+1):N;

% Define excitatory and inhibitory indices
exc1 = 1:Ne1;
inh1 = (Ne1+1):(Ne1+Ni1);
exc2 = (Ne1+Ni1+1):(Ne1+Ni1+Ne2);
inh2 = (Ne1+Ni1+Ne2+1):(Ne1+Ni1+Ne2+Ni2);
exc3 = (Ne1+Ni1+Ne2+Ni2+1):(Ne1+Ni1+Ne2+Ni2+Ne3);
inh3 = (Ne1+Ni1+Ne2+Ni2+Ne3+1):(Ne1+Ni1+Ne2+Ni2+Ne3+Ni3);

% Simulation parameters
conn_prob_within = 0.7;
conn_prob_between = 0.5;
conn_prob_vip_inh=1;
num_neuron = 500;
i_stren_values = 0:0.1:1;
results = cell(length(i_stren_values), 1);

% Start parallel pool
if isempty(gcp('nocreate'))
    parpool;
end
nsweep=1;
parfor i_stren_idx = 1:length(i_stren_values)
    i_stren = i_stren_values(i_stren_idx);
    current_results = struct();
    current_results.avg_firing_rate1 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate19 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate20 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate21 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate22 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate23 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate24 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate25 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate26 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate27 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate28 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate29 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate30 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate31 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate32 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate33 = zeros(num_neuron, nsweep);
    current_results.avg_firing_rate34 = zeros(num_neuron, nsweep);
    current_results.fine_grained_rates = cell(num_neuron, nsweep);
    % 在主循环中添加spike时间序列的存储
current_results.spike_times = cell(N, num_neuron, nsweep); % 存储每个神经元的spike时间
    for ineuron = 1:num_neuron
        % Generate parameters for current neuron
        re = rand(Ne1 + Ne2 + Ne3, 1);
        ri = rand(Ni1 + Ni2 + Ni3, 1);
        
        % Initialize neuron parameters
        a = zeros(N, 1);
        b = zeros(N, 1);
        c = zeros(N, 1);
        d = zeros(N, 1);
        
        % Set excitatory neuron parameters
        exc_indices = [1:Ne1, (Ne1+Ni1+1):(Ne1+Ni1+Ne2), (Ne1+Ni1+Ne2+Ni2+1):(Ne1+Ni1+Ne2+Ni2+Ne3)];
        a(exc_indices) = 0.02;
        b(exc_indices) = 0.2;
        c(exc_indices) = -65 + 15 * re.^2;
        d(exc_indices) = 8 - 6 * re.^2;
        
        % Set inhibitory neuron parameters
        inh_indices = [(Ne1+1):(Ne1+Ni1), (Ne1+Ni1+Ne2+1):(Ne1+Ni1+Ne2+Ni2), (Ne1+Ni1+Ne2+Ni2+Ne3+1):N-4];
        a(inh_indices) = 0.08 + 0.03 * ri;
        b(inh_indices) = 0.25 - 0.05 * ri;
        c(inh_indices) = -65;
        d(inh_indices) = 2;
        
        % Set VIP neuron parameters
        a(vip) = 0.1 + 0.05 * rand(Nvip, 1);
        b(vip) = 0.25 - 0.05 * rand(Nvip, 1);
        c(vip) = -65 + 10 * randn(Nvip, 1);
        d(vip) = 2 + rand(Nvip, 1);
        
        % Initialize membrane potentials
        v = -65 * ones(N,1);
        u = b.*v;
        v(group3) = -65 + 15*rand(length(group3), 1);
        u(group3) = b(group3) .* v(group3);
        v(vip) = -65 + 15*rand(Nvip, 1);
        u(vip) = b(vip) .* v(vip);
        
        % Initialize connectivity matrix
        S_raw = zeros(N, N);
        S_raw = setup_connectivity(S_raw, N, exc1, exc2, exc3, inh1, inh2, inh3, group1, group2, group3, vip, conn_prob_within, conn_prob_between,conn_prob_vip_inh, i_stren);        
        % Run sweeps
        firings_ulity = cell(1,nsweep);
       
        for isweep = 1:nsweep
             [firings, rates, spike_times] = run_sweep(v, u, a, b, c, d, S_raw, N, group3, vip, inh3);
        
            % 存储spike时间序列
            for n = 1:N
                current_results.spike_times{n, ineuron, isweep} = spike_times{n};
            end
            firings_ulity{isweep} = firings;
            
            % Store results
            current_results.avg_firing_rate1(ineuron, isweep) = rates(1);
            current_results.avg_firing_rate19(ineuron, isweep) = rates(19);
            current_results.avg_firing_rate20(ineuron, isweep) = rates(20);
            current_results.avg_firing_rate21(ineuron, isweep) = rates(21);
            current_results.avg_firing_rate22(ineuron, isweep) = rates(22);
            current_results.avg_firing_rate23(ineuron, isweep) = rates(23);
            current_results.avg_firing_rate24(ineuron, isweep) = rates(24);
            current_results.avg_firing_rate25(ineuron, isweep) = rates(25);
            current_results.avg_firing_rate26(ineuron, isweep) = rates(26);
            current_results.avg_firing_rate27(ineuron, isweep) = rates(27);
            current_results.avg_firing_rate28(ineuron, isweep) = rates(28);
            current_results.avg_firing_rate29(ineuron, isweep) = rates(29);
            current_results.avg_firing_rate30(ineuron, isweep) = rates(30);
            current_results.avg_firing_rate31(ineuron, isweep) = rates(31);
            current_results.avg_firing_rate32(ineuron, isweep) = rates(32);
            current_results.avg_firing_rate33(ineuron, isweep) = rates(33);
            current_results.avg_firing_rate34(ineuron, isweep) = rates(34);
            % current_results.fine_grained_rates{ineuron, isweep} = calculate_fine_rates(firings, N);
             % 计算发放率和同步性
            [fine_rates, sync_measure] = calculate_fine_rates(firings, spike_times, N);
            current_results.fine_grained_rates{ineuron, isweep} = fine_rates;
            current_results.synchrony_measure{ineuron, isweep} = sync_measure;
        end
        current_results.firings_ulity{ineuron} = firings_ulity;
    end
    results{i_stren_idx} = current_results;
end

% 处理和保存结果
for i = 1:length(i_stren_values)    
    i_stren = i_stren_values(i);
    save_path = ['F:\neural model\neuron_VIP\neuron_' num2str(i_stren)];
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    save_filename = ['avg_firing_rates_' num2str(i_stren) '.mat'];
    parsave(fullfile(save_path, save_filename), results{i});
    disp(['Results for i_stren = ' num2str(i_stren) ' saved to: ' fullfile(save_path, save_filename)]);
end
base_path = 'F:\neural model\neuron_VIP';
analyze_network_results(base_path);
% 辅助函数
function parsave(fname, data)
    save(fname, '-struct', 'data');
end

% Helper functions (defined in the same file)
function S_raw = setup_connectivity(S_raw, N, exc1, exc2, exc3, inh1, inh2, inh3, group1, group2, group3, vip, conn_prob_within, conn_prob_between,conn_prob_vip_inh, i_stren)
    % Initialize connection matrix
    S_raw = zeros(N, N);
    
    % Setup connections for each neuron pair
    for i = 1:N
        for j = 1:N
            if i == j
                S_raw(i,j) = 0;
                continue;
            end

            % Handle excitatory connections
            if ismember(j, [exc1, exc2, exc3])
                if j == 1
                    if ismember(i, [group1, group2])
                        S_raw(i,j) = (rand() <= conn_prob_between) * rand_weight(50, 10);
                    else
                        S_raw(i,j) = 0;
                    end
                elseif (any(i == group1) && any(j == exc1)) || ...
                       (any(i == group2) && any(j == exc2)) || ...
                       (any(i == group3) && any(j == exc3))
                    S_raw(i,j) = (rand() <= conn_prob_within) * rand_weight(8, 2) * (0.8 + 0.4*rand());
                elseif any(i == group2) && any(j == group1)
                    S_raw(i,j) = (rand() <= conn_prob_between) * rand_weight(50, 10);
                elseif any(i == group3) && any(j == group2)
                    S_raw(i,j) = (rand() <= conn_prob_between) * rand_weight(8, 2) * (0.8 + 0.4*rand());
                elseif any(i == group3) && any(j == group1)
                    S_raw(i,j) = (rand() <= 0.3) * rand_weight(2, 0.5);  % 弱连接
                else
                    S_raw(i,j) = 0;
                end
            
            % Handle inhibitory connections
            elseif ismember(j, [inh1, inh2, inh3])
                if j == 1
                    S_raw(i,j) = 0;
                elseif ismember(i, [exc1, exc2, exc3]) && ...
                       ((any(i == exc1) && any(j == group1)) || ...
                        (any(i == exc2) && any(j == group2)) || ...
                        (any(i == exc3) && any(j == group3)))
                    S_raw(i,j) = (rand() <= conn_prob_within) * -rand_weight(8, 2);
                elseif i == 1 
                    S_raw(i,j) = 0;
                else
                    S_raw(i,j) = 0;
                end
            end
        end
    end
    
    % Setup VIP neuron connections
    % VIP neurons receive input from exc2 and exc3
    for i = vip
        for j = [exc2, exc3]
            S_raw(i,j) = (rand() <= conn_prob_between) * rand_weight(8, 2);
        end
    end

    for i = inh3
        for j = [exc3,exc2]
            S_raw(i,j) = (rand() <= conn_prob_between) * rand_weight(2, 0.1);
        end
        
    end

    % VIP neurons inhibit inh3 neurons with non-linear saturation
    for i = inh3
        for j = vip
            % base_weight = rand_weight(50, 10);
            base_weight=30;
            S_raw(i,j) = (rand() <= conn_prob_vip_inh) * -base_weight * (i_stren);
            % Add non-linear saturation to VIP inhibition
            % S_raw(i,j) = (rand() <= conn_prob_within) * -base_weight * (1 / (1 + exp(-(i_stren-1))));
        end
    end
    
    % % Scale VIP output strength
    % for j = vip
    %     S_raw(:,j) = S_raw(:,j) * i_stren;
    % end
end

function w = rand_weight(mean_val, std_dev)
    w = mean_val + std_dev * randn();
    w = max(0, w);  % Ensure non-negative weights
end

function [firings, rates, spike_times] = run_sweep(v, u, a, b, c, d, S_raw, N, group3, vip, inh3)
    T = 20000;
    dt = 0.5;
    firings = [];
    spike_times = cell(N, 1); % 为每个神经元创建一个cell来存储spike时间
    % Initialize variables for this sweep
    S_dynamic = S_raw;
    adaptation = zeros(N, 1);
    tau_adapt = 150;
    noise_scale = 10;
    adaptation_strength=0.1;
    % Generate noise for the entire simulation
    noise_streams = randn(N, T);
    
    for t = 1:T
        % Calculate input current
        I = calculate_input(t, N, noise_streams(:,t), noise_scale, vip, group3, inh3);
        
        % Find fired neurons
        fired = find(v >= 30);
        
        if ~isempty(fired)
            firings = [firings; t+0*fired, fired];
            for i = 1:length(fired)
                neuron_idx = fired(i);
                spike_times{neuron_idx} = [spike_times{neuron_idx}; t*dt]; % 转换为实际时间（ms）
            end
            v(fired) = c(fired);
            u(fired) = u(fired) + d(fired);
            
            % Update adaptation
            adaptation(fired) = adaptation(fired) + adaptation_strength;
        end
        
        % Update membrane potential and recovery variable
        I = I - adaptation;
        % Synaptic transmission with probabilistic release
        synaptic_release_prob = 0.8;
        % effective_S = S_dynamic .* (rand(size(S_dynamic)) < synaptic_release_prob);
        effective_S = S_dynamic ;
        I = I + sum(effective_S(:,fired), 2);
        v = v + 0.5 * (0.04 * v.^2 + 5 * v + 140 - u + I);
        v = v + 0.5 * (0.04 * v.^2 + 5 * v + 140 - u + I);
        u = u + a.*(b.*v - u);
    

        adaptation = adaptation + dt * (-adaptation / tau_adapt);
    end
    
    % Calculate rates
    rates = zeros(N, 1);
    for i = 1:N
        rates(i) = sum(firings(:,2) == i) / (T/1000);
    end
end

function I = calculate_input(t, N, noise, noise_scale, vip, group3, inh3)
    % Calculate base input
    I1 = 3 * sin(50*pi*t/10000) + 3;
    I1 = abs(I1);
    I_weak = poissrnd(2, [N-1, 1]);    
    I= [I1; I_weak];
    
    
    % Add noise
    I = I + noise * noise_scale;
    % I(group3) = I(group3) + randn() * noise_scale;    
    

     % 增强VIP神经元的驱动
    % I(vip) = I(vip) + 5 + 2 * randn(length(vip), 1); % 增加基础电流
    
    % 为inh3添加稳定的背景输入
    background_current = 0;
    I(inh3) = I(inh3) + background_current;
    
    % 使用较小的噪声尺度
    noise_scale_inh = 0;
    I(inh3) = I(inh3) + noise(inh3) * noise_scale_inh;
end

function [rates, sync_measure] = calculate_fine_rates(firings, spike_times, N)
    window_size = 100;
    num_windows = 20;
    rates = zeros(N, num_windows);
    sync_measure = zeros(1, num_windows);
    for window = 1:num_windows
        start_time = (window-1)*window_size + 1;
        end_time = window*window_size;
        for neuron = 1:N
            spike_times_neuron = spike_times{neuron};
            spikes_in_window = sum(spike_times_neuron >= start_time & ...
                                 spike_times_neuron < end_time);
            spikes_in_window_1 = sum(firings(:,2) == neuron & firings(:,1) >= start_time & firings(:,1) < end_time);
            rates(neuron, window) = spikes_in_window_1 / (window_size/1000);
        end
        % 计算同步性度量
        all_spikes_in_window = [];
        for neuron = 19:N-6
            spikes = spike_times{neuron};
            window_spikes = spikes(spikes >= start_time & spikes < end_time);
            all_spikes_in_window = [all_spikes_in_window; window_spikes];
        end
        
        % 使用ISI CV作为同步性度量
        
        if length(all_spikes_in_window) > 1
            ISIs = diff(sort(all_spikes_in_window));
            sync_measure(window) = std(ISIs) / mean(ISIs);
        end
    end
end

function save_results(results, i_stren_values)
    for i = 1:length(i_stren_values)
        i_stren = i_stren_values(i);
        save_path = ['F:\neural model\neuron_VIP\neuron_' num2str(i_stren)];
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        save_filename = ['avg_firing_rates_' num2str(i_stren) '.mat'];
        save(fullfile(save_path, save_filename), '-struct', 'results{i}');
    end
end