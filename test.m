clc; clear; 
%% 数据导入与预处理
data = readtable("./data.csv");
time = data.time;
tem = data.temperature;
volte = data.volte;

% 检测阶跃起始点（假设阶跃变化始于温度开始上升时）
[~, startIdx] = max(diff(tem) > mean(diff(tem))*2); % 找到温度显著上升的位置
step_time = time(startIdx); % 阶跃发生时刻
step_value = mean(volte(startIdx:end)); % 阶跃输入值

% 重构阶跃输入信号
input_signal = zeros(size(time));
input_signal(time >= step_time) = step_value;

% 输出信号预处理：去趋势
output_signal = detrend(tem, 'linear');

% 截取阶跃后的有效数据
validIdx = time >= (step_time - 10); % 保留阶跃前10秒数据用于初始条件估计
valid_time = time(validIdx);
valid_input = input_signal(validIdx);
valid_output = output_signal(validIdx);

% 数据重采样（确保等间隔）
Ts = mean(diff(valid_time)); % 平均采样周期
time_vec = (valid_time(1):Ts:valid_time(end))';
interp_input = interp1(valid_time, valid_input, time_vec, 'linear', 'extrap');
interp_output = interp1(valid_time, valid_output, time_vec, 'linear', 'extrap');

% 创建系统辨识数据对象
sys_data = iddata(interp_output, interp_input, Ts);

%% 初始参数估计
% 估算初始时滞（交叉相关法）
[c, lags] = xcorr(interp_output, interp_input - min(interp_input));
[~, maxIdx] = max(abs(c));
initial_delay = lags(maxIdx) * Ts; % 初始时滞估计

% 估算稳态增益
steady_output = mean(interp_output(end-round(10/Ts):end)); % 最后10秒平均值
transient_output = mean(interp_output(1:round(1/Ts))); % 初始1秒平均值
initial_gain = (steady_output - transient_output) / step_value; % 比例增益估计

% 估算时间常数（63%响应法）
target_value = transient_output + 0.63*(steady_output - transient_output);
[~, tcIdx] = min(abs(interp_output - target_value));
initial_Tp = time_vec(tcIdx) - step_time; % 时间常数估计

%% 传递函数辨识
% 配置辨识选项
opt = procestOptions;
opt.InitialCondition = 'estimate'; % 估计初始条件
opt.SearchOptions.MaxIterations = 50; % 增加迭代次数
opt.Display = 'on'; % 显示优化过程
opt.Focus = 'simulation'; % 优化仿真性能

% 设置合理参数范围
opt.K.Minimum = 0.8 * initial_gain; % 增益下限
opt.K.Maximum = 1.2 * initial_gain; % 增益上限
opt.Tp1.Minimum = 0.5 * initial_Tp; % 时间常数下限
opt.Tp1.Maximum = 2 * initial_Tp; % 时间常数上限
opt.Td.Minimum = max(0, 0.5 * initial_delay); % 时滞下限
opt.Td.Maximum = 2 * initial_delay; % 时滞上限

% 使用初始估计值
initialParams = struct(...
    'K', initial_gain, ...      % 初始增益
    'Tp1', initial_Tp, ...      % 初始时间常数
    'Td', initial_delay);       % 初始时滞

% 执行传递函数辨识
fprintf('\n=== 开始传递函数辨识 ===\n');
fprintf('初始估计参数:\n');
fprintf('  比例增益 K = %.4f\n', initial_gain);
fprintf('  时间常数 Tp = %.2f s\n', initial_Tp);
fprintf('  时滞 Td = %.2f s\n', initial_delay);

% 使用PROCEST进行模型辨识
model = procest(sys_data, 'P1D', initialParams, opt);

% 显示最终模型参数
fprintf('\n=== 辨识结果 ===\n');
fprintf('最终比例增益 K = %.4f\n', model.K);
fprintf('最终时间常数 Tp = %.2f s\n', model.Tp1);
fprintf('最终时滞 Td = %.2f s\n', model.Td);
fprintf('拟合优度 = %.1f%%\n', model.Report.Fit.FitPercent);

%% 模型验证与分析
% 模型与实际输出对比
figure('Position', [100, 100, 800, 500]);
compare(sys_data, model);
title('模型与实际输出对比');
grid on;

% 残差分析
figure('Position', [100, 100, 800, 400]);
resid(sys_data, model);
title('残差分析');
grid on;

% 模型阶跃响应
figure('Position', [100, 100, 800, 400]);
step(model);
title('模型阶跃响应');
grid on;

% 传递函数表达式
[num, den] = tfdata(model);
s = tf('s');
tf_model = model.K * exp(-model.Td * s) / (model.Tp1 * s + 1);
disp(' ');
disp('传递函数表达式:');
tf_model

% 零点极点图
figure('Position', [100, 100, 800, 400]);
pzmap(tf_model);
title('极点零点图');
grid on;

% 波特图分析
figure('Position', [100, 100, 800, 500]);
bode(tf_model);
title('系统波特图');
grid on;