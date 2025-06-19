%@brief 利用粒子群算法整定PID参数并进行结果可视化
%@author 23010341 杳泽
clc; clear; close all;

% 加载系统模型
load("./gs_result.mat");
sys_gs = P1D_fit_by_matlab;

% 验证系统模型存在
if ~exist('sys_gs', 'var')
    error('系统传递函数sys_gs必须在工作区中');
end

% 系统参数设置
target_temp = 35;     % 目标温度(°C)
Ts = 0.5;             % 采样时间(秒)
t_sim = 5000;         % 延长仿真时间至5000秒
err_band = 0.02 * target_temp; % 2%误差带宽度

% 显示系统模型特性
fprintf('===== 系统模型信息 =====\n');
disp(sys_gs);

% 提取模型参数
model_Kp = sys_gs.Kp;
model_Tp1 = sys_gs.Tp1;
model_Td = sys_gs.Td;

fprintf('\n关键模型参数:\n');
fprintf('增益 Kp = %.4f\n', model_Kp);
fprintf('时间常数 Tp1 = %.4f 秒\n', model_Tp1);
fprintf('延迟时间 Td = %.4f 秒\n', model_Td);

% ===== PID目标函数 =====
function [cost, t, y] = pid_objective(K, sys, target_temp, Ts, t_sim)
    % 创建PID控制器
    Kp = K(1); Ki = K(2); Kd = K(3);
    C = pid(Kp, Ki, Kd);
    
    % 从idproc构建传递函数（考虑延迟）
    sys_tf = tf(sys.Kp, [sys.Tp1, 1], 'InputDelay', sys.Td);
    
    % 构建闭环系统
    sys_cl = feedback(C * sys_tf, 1);
    
    % 仿真时间设置（考虑延迟）
    t = (0:Ts:t_sim)';
    r = zeros(size(t));
    r(t>=1) = target_temp; % 1秒后施加阶跃输入
    
    % 初始化变量
    cost = 1e6; % 默认大代价
    y = zeros(size(t));
    
    try
        % 仿真系统响应
        y = lsim(sys_cl, r, t);
        
        % 验证仿真结果
        if any(isnan(y)) || any(isinf(y)) || max(y) > 500
            error('仿真结果异常');
        end
        
        % 计算误差
        err = r - y;
        
        % 计算超调量
        [peak, idx] = max(y);
        if peak > target_temp
            overshoot = (peak - target_temp) / target_temp * 100;
        else
            overshoot = 0;
        end
        
        % 计算上升时间(10%-90%)
        idx_10 = find(y >= 0.1 * target_temp, 1, 'first');
        idx_90 = find(y >= 0.9 * target_temp, 1, 'first');
        if isempty(idx_10) || isempty(idx_90)
            rise_time = t_sim; % 惩罚值
        else
            rise_time = t(idx_90) - t(idx_10);
        end
        
        % 计算调节时间（改进方法：从后向前查找）
        abs_error = abs(y - target_temp);
        in_band = abs_error <= (0.02 * target_temp);
        
        % 找到最后一个超出误差带的点
        last_out_index = find(~in_band, 1, 'last');
        
        if isempty(last_out_index)
            % 系统从未超出：取首次进入时间
            settle_index = find(in_band, 1, 'first');
            if isempty(settle_index)
                settle_time = t_sim; % 从未进入误差带
            else
                settle_time = t(settle_index);
            end
        elseif last_out_index < length(t)
            % 检查最后一个超出点后的稳定性
            post_band = in_band(last_out_index+1:end);
            
            if all(post_band)
                % 调节时间取最后一个超出点后的第一个采样点
                settle_time = t(last_out_index+1);
            else
                settle_time = t_sim; % 后续仍有超出
            end
        else
            settle_time = t_sim; % 结束时仍超出
        end
        
        % ITAE指标（时间加权绝对误差）
        itae = trapz(t, t.*abs(err));
        
        % 综合成本函数（考虑系统特性）
        cost = itae + ...              % 主要误差指标
               500 * abs(overshoot) + ...  % 超调惩罚
               10 * settle_time + ...     % 调节时间惩罚
               rise_time;                 % 上升时间惩罚
        
        % 额外惩罚：系统不稳定或响应过慢
        if settle_time > 0.8 * t_sim
            cost = cost * 2;
        end
        
    catch ME
        % 显示错误信息
        fprintf('仿真错误: %s\n', ME.message);
        % 保持大代价
    end
end

% ===== 粒子群优化算法实现 =====
function [best_params, best_cost] = pso_pid_optimization(objective_func, lb, ub, options)
    % 初始化粒子群
    n_particles = options.SwarmSize;
    n_dims = length(lb);
    
    % 初始化粒子位置和速度
    particles = zeros(n_particles, n_dims);
    velocity = zeros(n_particles, n_dims);
    
    % 基于系统特性的初始猜测
    base_Kp = 0.5 / options.model_Kp; % 基于系统增益的初始估计
    base_Ki = 0.1 / (options.model_Tp1 + options.model_Td);
    base_Kd = 0.1 * (options.model_Tp1 + options.model_Td);
    
    for i = 1:n_particles
        % 在基础值附近随机初始化
        particles(i, :) = [base_Kp, base_Ki, base_Kd] .* (0.5 + rand(1, 3));
        velocity(i, :) = -0.1*(ub - lb) + 0.2*(ub - lb).*rand(1, n_dims);
    end
    
    % 初始化个体最优和全局最优
    pbest = particles;
    pbest_cost = inf(n_particles, 1);
    gbest = particles(1, :);
    gbest_cost = inf;
    
    % 记录迭代过程
    history = zeros(options.MaxIterations, 1);
    
    % PSO主循环
    for iter = 1:options.MaxIterations
        % 评估所有粒子
        valid_particles = 0;
        for i = 1:n_particles
            [cost, ~, ~] = objective_func(particles(i, :));
            
            % 只考虑有效解
            if cost < 1e6
                valid_particles = valid_particles + 1;
                
                % 更新个体最优
                if cost < pbest_cost(i)
                    pbest(i, :) = particles(i, :);
                    pbest_cost(i) = cost;
                end
                
                % 更新全局最优
                if cost < gbest_cost
                    gbest = particles(i, :);
                    gbest_cost = cost;
                end
            end
        end
        
        % 如果没有有效粒子，重新初始化部分粒子
        if valid_particles == 0
            fprintf('迭代 %d: 无有效粒子，重新初始化...\n', iter);
            for i = 1:ceil(n_particles/3)
                idx = randi(n_particles);
                particles(idx, :) = lb + rand(1, n_dims) .* (ub - lb);
                velocity(idx, :) = -0.1*(ub - lb) + 0.2*(ub - lb).*rand(1, n_dims);
            end
            continue;
        end
        
        % 记录全局最优
        history(iter) = gbest_cost;
        
        % 更新粒子速度和位置
        w = options.InertiaRange(1) + (options.InertiaRange(2) - options.InertiaRange(1)) * ...
            (1 - iter/options.MaxIterations); % 线性递减惯性权重
        
        for i = 1:n_particles
            % 更新速度
            r1 = rand(1, n_dims);
            r2 = rand(1, n_dims);
            
            velocity(i, :) = w * velocity(i, :) + ...
                options.SelfAdjustment * r1 .* (pbest(i, :) - particles(i, :)) + ...
                options.SocialAdjustment * r2 .* (gbest - particles(i, :));
            
            % 速度限制
            velocity(i, :) = min(max(velocity(i, :), -options.VelocityLimit), options.VelocityLimit);
            
            % 更新位置
            particles(i, :) = particles(i, :) + velocity(i, :);
            
            % 位置边界处理
            particles(i, :) = max(lb, min(ub, particles(i, :)));
        end
        
        % 显示迭代信息
        if mod(iter, 5) == 0 || iter == 1
            fprintf('迭代 %d/%d: 最佳代价 = %.4f, 有效粒子 = %d/%d\n', ...
                iter, options.MaxIterations, gbest_cost, valid_particles, n_particles);
        end
    end
    
    % 返回结果
    best_params = gbest;
    best_cost = gbest_cost;
    
    % 绘制收敛曲线
    figure;
    plot(1:options.MaxIterations, history, 'b-', 'LineWidth', 1.5);
    xlabel('迭代次数');
    ylabel('最佳适应度值');
    title('PSO收敛曲线');
    grid on;
end

% ===== 计算控制信号 =====
function [t, u] = compute_control_signal(sys, params, target_temp, Ts)
    % 创建PID控制器
    Kp = params(1); Ki = params(2); Kd = params(3);
    C = pid(Kp, Ki, Kd);
    
    % 从idproc构建传递函数（考虑延迟）
    sys_tf = tf(sys.Kp, [sys.Tp1, 1], 'InputDelay', sys.Td);
    
    % 构建闭环系统
    sys_cl = feedback(C * sys_tf, 1);
    
    % 仿真时间
    t = (0:Ts:5000)'; % 延长仿真时间
    
    % 创建参考信号
    r = zeros(size(t));
    r(t>=1) = target_temp;
    
    % 仿真系统响应
    y = lsim(sys_cl, r, t);
    
    % 计算误差
    err = r - y;
    
    % 积分误差项
    integ_err = cumtrapz(t, err);
    
    % PID输出 = Kp*e + Ki*∫e + Kd*de/dt
    u = zeros(size(err));
    for i = 2:length(t)
        de = (err(i) - err(i-1)) / Ts;
        u(i) = Kp * err(i) + Ki * integ_err(i) + Kd * de;
        % 限幅0-10V
        u(i) = max(0, min(10, u(i)));
    end
end

% ===== 绘制响应曲线 =====
function plot_pid_response(params, performance, target_temp)
    Kp = params(1); Ki = params(2); Kd = params(3);
    overshoot = performance.overshoot;
    settle_time = performance.settle_time;
    rise_time = performance.rise_time;
    peak_time = performance.peak_time;
    t = performance.t;
    y = performance.y;
    u = performance.u;
    err = performance.err;
    
    % 创建图形
    figure('Position', [100, 100, 1000, 800], 'Color', 'w');
    
    % 温度响应
    subplot(3, 1, 1);
    plot(t, y, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(t, target_temp*ones(size(t)), 'r--', 'LineWidth', 1.5);
    
    % 2%误差带
    err_band = 0.02 * target_temp;
    plot([t(1), t(end)], [target_temp+err_band, target_temp+err_band], 'k:');
    plot([t(1), t(end)], [target_temp-err_band, target_temp-err_band], 'k:');
    
    % 标记调节时间
    if isfinite(settle_time) && settle_time < t(end)
        settle_idx = find(t >= settle_time, 1);
        plot([settle_time, settle_time], [min(y), max(y)], 'g--', 'LineWidth', 1);
        plot(settle_time, y(settle_idx), 'go', 'MarkerSize', 8, 'LineWidth', 1.5);
        text(settle_time, max(y)*1.4, sprintf('调节时间: %.1f s', settle_time), ...
            'HorizontalAlignment', 'right', 'BackgroundColor', 'white');
    end
    
    % 标记超调
    if overshoot > 0
        [peak, idx] = max(y);
        plot(peak_time, peak, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
        text(peak_time, peak*1.2, sprintf('超调: %.1f%%', overshoot), ...
            'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
    end
    
    % 标记上升时间
    if isfinite(rise_time)
        text(t(end)*0.7, target_temp*0.5, sprintf('上升时间: %.1f s', rise_time), ...
            'BackgroundColor', 'white');
    end
    
    title(sprintf('温度响应曲线 (K_p=%.3f, K_i=%.3f, K_d=%.3f)', Kp, Ki, Kd));
    xlabel('时间 (秒)');
    ylabel('温度 (°C)');
    legend('系统响应', '目标温度', '2%误差带', 'Location', 'southeast');
    grid on;
    ylim([0, max(50, 1.2*max(y))]);
    
    % 控制信号
    subplot(3, 1, 2);
    plot(t, u, 'm-', 'LineWidth', 1.5);
    title('控制信号 (电压)');
    xlabel('时间 (秒)');
    ylabel('电压 (V)');
    grid on;
    ylim([0, 10.5]);
    
    % 误差曲线
    subplot(3, 1, 3);
    plot(t, err, 'c-', 'LineWidth', 1.5);
    hold on;
    plot([t(1), t(end)], [err_band, err_band], 'k:');
    plot([t(1), t(end)], [-err_band, -err_band], 'k:');
    title('温度误差');
    xlabel('时间 (秒)');
    ylabel('误差 (°C)');
    grid on;
end

% ================= 主程序 =================
fprintf('\n===== 开始优化PSO算法PID整定 - 目标温度 %.1f°C =====\n', target_temp);

% 设置PSO参数（考虑系统特性）
pso_options = struct(...
    'SwarmSize', 40, ...           % 增加粒子数量
    'MaxIterations', 100, ...      % 最大迭代次数
    'InertiaRange', [0.4, 0.9], ...% 惯性权重范围
    'SelfAdjustment', 1.2, ...     % 个体学习因子
    'SocialAdjustment', 1.5, ...   % 社会学习因子
    'VelocityLimit', 0.1, ...      % 减小速度限制
    'model_Kp', sys_gs.Kp, ...     % 添加模型参数
    'model_Tp1', sys_gs.Tp1, ...   % 添加模型参数
    'model_Td', sys_gs.Td ...      % 添加模型参数
);

% 参数范围 [Kp, Ki, Kd]（基于模型特性调整）
lb = [0.01, 0.0001, 0.001];   % 参数下限（更保守）
ub = [5, 0.05, 5];            % 参数上限（考虑大时间常数）

fprintf('\n参数范围:\n');
fprintf('Kp: [%.4f, %.4f]\n', lb(1), ub(1));
fprintf('Ki: [%.4f, %.4f]\n', lb(2), ub(2));
fprintf('Kd: [%.4f, %.4f]\n', lb(3), ub(3));

% 创建目标函数
objective_func = @(K) pid_objective(K, sys_gs, target_temp, Ts, t_sim);

% 运行PSO优化
fprintf('\n===== 运行改进的粒子群优化算法 =====\n');
[params_opt, best_cost] = pso_pid_optimization(objective_func, lb, ub, pso_options);

% 显示优化结果
fprintf('\n===== PSO优化完成 =====\n');
fprintf('最优成本值: %.4f\n', best_cost);
fprintf('最优PID参数: Kp=%.4f, Ki=%.4f, Kd=%.4f\n', ...
        params_opt(1), params_opt(2), params_opt(3));

% 使用优化参数进行最终仿真
fprintf('\n===== 进行最终仿真 =====\n');
[~, t, y] = objective_func(params_opt);
[t_u, u] = compute_control_signal(sys_gs, params_opt, target_temp, Ts);

% 确保时间向量一致
if length(t) ~= length(t_u)
    t = t_u;
end

% 计算性能指标
err = target_temp - y;
[peak, idx] = max(y);
if peak > target_temp
    overshoot = (peak - target_temp) / target_temp * 100;
else
    overshoot = 0;
end
peak_time = t(idx);

% 上升时间(10%-90%)
idx_10 = find(y >= 0.1 * target_temp, 1, 'first');
idx_90 = find(y >= 0.9 * target_temp, 1, 'first');
if isempty(idx_10) || isempty(idx_90)
    rise_time = NaN;
else
    rise_time = t(idx_90) - t(idx_10);
end

% 调节时间计算（最终验证）
abs_error = abs(y - target_temp);
in_band = abs_error <= (0.02 * target_temp);
last_out_index = find(~in_band, 1, 'last');
if isempty(last_out_index)
    settle_index = find(in_band, 1, 'first');
    if isempty(settle_index)
        settle_time = Inf;
    else
        settle_time = t(settle_index);
    end
elseif last_out_index < length(t)
    post_band = in_band(last_out_index+1:end);
    if all(post_band)
        settle_time = t(last_out_index+1);
    else
        settle_time = Inf;
    end
else
    settle_time = Inf;
end

% 显示性能指标
fprintf('\n===== 系统性能指标 =====\n');
fprintf('超调量: %.2f%%\n', overshoot);
fprintf('峰值时间: %.2f 秒\n', peak_time);
fprintf('上升时间: %.2f 秒\n', rise_time);
fprintf('调节时间(2%%误差带): %.2f 秒\n', settle_time);
fprintf('稳态误差: %.4f°C\n', mean(abs_error(end-50:end)));

% 封装性能数据
performance = struct(...
    'overshoot', overshoot, ...
    'settle_time', settle_time, ...
    'rise_time', rise_time, ...
    'peak_time', peak_time, ...
    't', t, ...
    'y', y, ...
    'u', u(1:length(t)), ... % 确保长度一致
    'err', err ...
);

% 绘制响应曲线
plot_pid_response(params_opt, performance, target_temp);

% 保存结果
save('pso_optimized_pid.mat', 'params_opt', 'overshoot', 'rise_time', 'settle_time',"peak_time");
fprintf('\n结果已保存到 pso_optimized_pid.mat\n');