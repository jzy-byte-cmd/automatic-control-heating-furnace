%@brief 利用遗传算法整定PID参数并进行结果可视化
%@author 23010341 杳泽
% PID优化整定 - 温度控制系统
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
t_sim = 1000;         % 仿真时间(秒)
err_band = 0.02 * target_temp; % 2%误差带宽度

% 显示系统模型特性
fprintf('===== 系统模型信息 =====\n');
disp(sys_gs);

% PID优化目标函数
function cost = pid_objective(K, sys, target_temp, Ts, t_sim)
    % 创建PID控制器
    Kp = K(1); Ki = K(2); Kd = K(3);
    C = pid(Kp, Ki, Kd);
    
    % 构建闭环系统
    sys_cl = feedback(C * sys, 1);
    
    % 仿真时间设置
    t = (0:Ts:t_sim)';
    r = zeros(size(t));
    r(t>=1) = target_temp; % 1秒后施加阶跃输入
    
    % 初始化变量
    cost = 1e6; % 默认大代价
    
    try
        % 仿真系统响应
        y = lsim(sys_cl, r, t);
        
        % 验证仿真结果
        if any(isnan(y)) || any(isinf(y))
            error('仿真结果包含NaN或Inf值');
        end
        
        % 计算误差
        err = r - y;
        
        % ITAE指标（时间加权绝对误差）
        itae = trapz(t, t.*abs(err));
        
        % 计算超调量
        [peak, idx] = max(y);
        if peak > target_temp
            overshoot = (peak - target_temp) / target_temp * 100;
        else
            overshoot = 0;
        end
        
        % 计算峰值时间
        peak_time = t(idx);
        
        % 计算上升时间(10%-90%)
        y_target = target_temp;
        idx_10 = find(y >= 0.1 * y_target, 1, 'first');
        idx_90 = find(y >= 0.9 * y_target, 1, 'first');
        if isempty(idx_10) || isempty(idx_90)
            rise_time = t_sim; % 惩罚值
        else
            rise_time = t(idx_90) - t(idx_10);
        end
        
        % 计算调节时间（改进方法：从后向前查找）
        err_band = 0.02 * target_temp;
        abs_error = abs(y - target_temp);
        in_band = abs_error <= err_band;
        
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
        
        % 综合成本函数（权重可调整）
        cost = itae + 100 * abs(overshoot) + 10 * settle_time + rise_time;
        
    catch
        % 保持大代价
    end
end

% 遗传算法优化PID参数
fprintf('\n===== 开始遗传算法优化PID参数 =====\n');

% 参数范围 [Kp, Ki, Kd]
lb = [0.1, 0.001, 0.001]; % 参数下限
ub = [10, 4, 80];          % 参数上限

% 遗传算法选项
options = optimoptions('ga', ...
    'Display', 'iter', ...
    'MaxGenerations', 100, ...
    'PopulationSize', 100, ...
    'EliteCount', 5, ...
    'CrossoverFraction', 0.8, ...
    'FunctionTolerance', 1e-3, ...
    'PlotFcn', {@gaplotbestf});

% 目标函数
objective_func = @(K) pid_objective(K, sys_gs, target_temp, Ts, t_sim);

% 运行遗传算法
[params_opt, fval, ~, output] = ga(objective_func, 3, [], [], [], [], lb, ub, [], options);

% 显示优化结果
fprintf('\n===== 优化完成 =====\n');
fprintf('最优成本值: %.4f\n', fval);
fprintf('迭代次数: %d\n', output.generations);
fprintf('函数评估次数: %d\n', output.funccount);
fprintf('最优PID参数: Kp=%.4f, Ki=%.4f, Kd=%.4f\n', params_opt(1), params_opt(2), params_opt(3));

% 使用优化参数进行最终仿真
Kp = params_opt(1); Ki = params_opt(2); Kd = params_opt(3);
C = pid(Kp, Ki, Kd);
sys_cl = feedback(C * sys_gs, 1);

t = (0:Ts:t_sim)';
r = zeros(size(t));
r(t>=1) = target_temp;

y = lsim(sys_cl, r, t);

% 计算控制信号
u = zeros(size(t));
err = r - y;
integ_err = cumtrapz(t, err); % 积分误差
for i = 2:length(t)
    de = (err(i) - err(i-1)) / Ts; % 微分误差
    u(i) = Kp * err(i) + Ki * integ_err(i) + Kd * de;
    u(i) = max(0, min(10, u(i))); % 限幅0-10V
end

% 计算性能指标
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
in_band = abs_error <= err_band;
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

% 绘制响应曲线
figure('Position', [100, 100, 1000, 800], 'Color', 'w');

% 温度响应
subplot(3, 1, 1);
plot(t, y, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, r, 'r--', 'LineWidth', 1.5);
plot([t(1), t(end)], [target_temp*(1+0.02), target_temp*(1+0.02)], 'k:');
plot([t(1), t(end)], [target_temp*(1-0.02), target_temp*(1-0.02)], 'k:');

% 标记调节时间
if isfinite(settle_time)
    plot([settle_time, settle_time], [min(y), max(y)], 'g--', 'LineWidth', 1);
    plot(settle_time, y(t == settle_time), 'go', 'MarkerSize', 8, 'LineWidth', 1.5);
    text(settle_time, max(y)-5, sprintf('调节时间: %.2f s', settle_time), ...
        'HorizontalAlignment', 'right');
end

% 标记超调和峰值时间
if overshoot > 0
    plot(peak_time, peak, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot([peak_time, peak_time], [min(y), max(y)], 'g--', 'LineWidth', 1);
    plot(peak_time, y(t == peak_time), 'go', 'MarkerSize', 8, 'LineWidth', 1.5);
    text(peak_time, peak+2, sprintf('超调: %.2f%%', overshoot), ...
        'HorizontalAlignment', 'right');
    text(peak_time, peak+7, sprintf('峰值时间: %.2fs', peak_time), ...
        'HorizontalAlignment', 'right');
end

title(sprintf('温度响应曲线 (K_p=%.3f, K_i=%.3f, K_d=%.3f)', Kp, Ki, Kd));
xlabel('时间 (秒)');
ylabel('温度 (°C)');
legend('系统响应', '目标温度', '2%误差带', 'Location', 'southeast');
grid on;
ylim([0, max(50, peak*1.1)]);

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

% 保存结果
save('optimized_pid.mat', 'params_opt', 'overshoot', 'rise_time', 'settle_time');