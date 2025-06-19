% @brief  通过智能优化算法，针对得到的传递函数进行PID整定并可视化显示
clc; clear; 

load("./gs_result.mat")
sys_gs = P1D_fit_by_matlab;

% 验证系统模型存在
if ~exist('sys_gs', 'var')
    error('系统传递函数sys_gs必须在工作区中');
end

% 设定目标温度
target_temp = 35;  % 目标温度35°C
Ts = 0.5; % 采样时间（在函数外部定义）



% 显示系统模型特性
fprintf('===== 系统模型信息 =====\n');
disp(sys_gs);

% 定义PID整定目标函数（修复输出参数问题）
function [cost, t, y] = pid_objective(K, sys, target_temp, Ts)
    % 创建PID控制器
    Kp = K(1); Ki = K(2); Kd = K(3);
    C = pid(Kp, Ki, Kd);
    
    % 构建闭环系统
    sys_cl = feedback(C * sys, 1);
    
    % 使用连续时间仿真避免离散化问题
    t = (0:Ts:1000)'; % 仿真时间800秒
    r = zeros(size(t));
    r(t>=1) = target_temp; % 1秒后施加阶跃输入到目标温度
    
    % 初始化y，避免未定义
    y = zeros(size(t));
    cost = 1e8; % 默认大代价
    
    try
        % 使用lsim进行连续时间仿真
        y = lsim(sys_cl, r, t);
        
        % 检查仿真结果是否有效
        if any(isnan(y)) || any(isinf(y))
            error('仿真结果包含NaN或Inf值');
        end
        
        % 计算超调量
        [peak, idx] = max(y);
        overshoot = max(0, (peak - target_temp) / target_temp * 100);
        
        % 计算调节时间 (最后一次进入2%误差带)
        idx_settle = find(abs(y - target_temp) <= 0.02 * target_temp, 1, 'last');
        if isempty(idx_settle)
            settle_time = Inf;
        else
            settle_time = t(idx_settle);
        end
        
        % 计算上升时间
        idx_10 = find(y >= 0.1 * target_temp, 1, 'first');
        idx_90 = find(y >= 0.9 * target_temp, 1, 'first');
        if isempty(idx_10) || isempty(idx_90)
            rise_time = Inf;
        else
            rise_time = t(idx_90) - t(idx_10);
        end
        
        % 计算峰值时间
        peak_time = t(idx);
        
        % ITAE指标（时间加权绝对误差）
        err = r - y;
        itae = trapz(t, t.*abs(err));
        
        % 综合成本函数（包含多个性能指标）
        cost = 0.4 * itae + 0.3 * abs(overshoot) + 0.2 * settle_time + 0.1 * rise_time;
        
        % 添加约束惩罚
        if overshoot > 15 % 超调量约束
            cost = cost + 100 * (overshoot - 15);
        end
        if settle_time > 400 % 调节时间约束
            cost = cost + 10 * (settle_time - 400);
        end
        
    catch ME
        % 如果出错，保持大代价
        fprintf('仿真错误: %s\n', ME.message);
    end
end

% 优化PID函数（修复输出处理问题）
function [best_params, performance] = optimize_pid(sys, target_temp, Ts, algorithm)
    % 定义PID参数范围 [Kp, Ki, Kd]
    lb = [0.1, 0.001, 0.001]; % 参数下限
    ub = [10, 2, 2];          % 参数上限
    
    % 初始猜测
    initial_guess = mean([lb; ub], 1);
    
    % 优化选项
    options = optimoptions('ga', ...
        'Display', 'iter', ...
        'MaxGenerations', 50, ...
        'PopulationSize', 30, ...
        'FunctionTolerance', 1e-3, ...
        'PlotFcn', {@gaplotbestf});
    
    % 优化函数（确保输出为一个标量）
    objective_func = @(K) pid_objective(K, sys, target_temp, Ts);
    
    % 执行优化
    [params_opt, fval, ~, output] = ga(objective_func, 3, [], [], [], [], lb, ub, [], options);
    
    % 输出优化信息
    fprintf('\n优化算法: %s\n', algorithm);
    fprintf('优化结果值: %.4f\n', fval);
    fprintf('迭代次数: %d\n', output.generations);
    fprintf('函数评估次数: %d\n', output.funccount);
    
    % 重新计算性能指标
    [~, t, y] = pid_objective(params_opt, sys, target_temp, Ts);
    
    % 计算性能指标
    [peak, idx] = max(y);
    overshoot = max(0, (peak - target_temp) / target_temp * 100);
    idx_settle = find(abs(y - target_temp) <= 0.02 * target_temp, 1, 'first');
    if isempty(idx_settle)
        settle_time = Inf;
    else
        settle_time = t(idx_settle);
    end
    peak_time = t(idx);
    
    % 封装性能指标
    performance = struct(...
        'Overshoot', overshoot, ...
        'SettleTime', settle_time, ...
        'PeakTime', peak_time, ...
        'Response', [t, y] ...
    );
    
    best_params = struct('Kp', params_opt(1), 'Ki', params_opt(2), 'Kd', params_opt(3));
end

% 绘制结果
function plot_results(params, performance, title_str, sys, target_temp)
    % 创建图形
    figure('Position', [100, 100, 800, 600], 'Color', 'w');
    
    % 提取响应数据
    t = performance.Response(:, 1);
    y = performance.Response(:, 2);
    
    % 绘制温度响应
    plot(t, y, 'b-', 'LineWidth', 2);
    hold on;
    
    % 绘制目标温度线
    plot([t(1), t(end)], [target_temp, target_temp], 'r-', 'LineWidth', 1.5);
    
    % 绘制2%误差带
    err_band = 0.02 * target_temp;
    plot([t(1), t(end)], [target_temp-err_band, target_temp-err_band], 'k:');
    plot([t(1), t(end)], [target_temp+err_band, target_temp+err_band], 'k:');
    
    % 标记调节时间
    if isfinite(performance.SettleTime)
        plot([performance.SettleTime, performance.SettleTime], [0, max(y)], 'g--');
    end
    
    % 标记峰值点
    [~, idx] = max(y);
    plot(t(idx), y(idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    
    % 设置图形属性
    grid on;
    title(sprintf('%s PID控制 (Kp=%.3f, Ki=%.3f, Kd=%.3f)', title_str, ...
        params.Kp, params.Ki, params.Kd));
    xlabel('时间 (s)');
    ylabel('温度 (°C)');
    ylim([0, max(45, target_temp*1.3)]);
    
    % 添加图例
    legend_items = {'系统响应', '目标温度', '2%误差带', '2%误差带', '调节时间', '峰值点'};
    legend(legend_items, 'Location', 'best');
    
    % 添加性能指标文本
    perf_text = {sprintf('超调量: %.2f%%', performance.Overshoot), ...
                sprintf('峰值时间: %.2f s', performance.PeakTime), ...
                sprintf('调节时间: %.2f s', performance.SettleTime)};
    
    text(0.05, 0.95, perf_text, 'Units', 'normalized', ...
         'FontSize', 10, 'VerticalAlignment', 'top', 'BackgroundColor', [1, 1, 1, 0.7]);
end

% 计算控制信号
function [t, u] = compute_control_signal(sys, params, target_temp, Ts)
    % 创建PID控制器
    Kp = params.Kp; Ki = params.Ki; Kd = params.Kd;
    C = pid(Kp, Ki, Kd);
    
    % 构建闭环系统
    sys_cl = feedback(C * sys, 1);
    
    % 仿真时间
    t = (0:Ts:1000)';
    
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
        de = (err(i) - err(i-1)) / (t(i) - t(i-1));
        u(i) = Kp * err(i) + Ki * integ_err(i) + Kd * de;
        % 限幅0-10V
        u(i) = max(0, min(10, u(i)));
    end
end

% ================= 主程序 =================
fprintf('\n===== 开始PID整定 - 目标温度 %.1f°C =====\n', target_temp);

% 使用遗传算法优化PID
fprintf('\n===== 使用遗传算法优化PID =====\n');
[params_opt, perf] = optimize_pid(sys_gs, target_temp, Ts, '遗传算法');

% 绘制结果
plot_results(params_opt, perf, '遗传算法优化', sys_gs, target_temp);

% 计算并绘制控制信号
[t, u] = compute_control_signal(sys_gs, params_opt, target_temp, Ts);

figure('Position', [100, 100, 800, 400], 'Color', 'w');
plot(t, u, 'LineWidth', 2);
grid on;
title('PID控制器输出');
xlabel('时间 (s)');
ylabel('控制电压 (V)');
yline(0, 'k');
yline(10, 'r');
legend('控制信号', '电压下限', '电压上限', 'Location', 'best');
ylim([-1, 11]);

% 输出性能指标
fprintf('\n===== PID控制器性能指标 =====\n');
fprintf('Kp = %.4f\n', params_opt.Kp);
fprintf('Ki = %.4f\n', params_opt.Ki);
fprintf('Kd = %.4f\n', params_opt.Kd);
fprintf('超调量: %.2f%%\n', perf.Overshoot);
fprintf('峰值时间: %.2f s\n', perf.PeakTime);
fprintf('调节时间: %.2f s\n', perf.SettleTime);