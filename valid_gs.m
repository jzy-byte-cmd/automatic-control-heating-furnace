% @brief  对比工作区中的传递函数的拟合度 MSE等数据 并作可视化显示  valid_gs.m
% @author 23010341 杳泽

clc;clear;

load("./gs_result.mat");    % 加载传递函数和数据到工作空间
sys1 = P1D_ref;             % 给予的参考传递函数
sys2 = P1D_2dots;           % 实际两点法计算所得
sys3 = P1D_fit_by_matlab;   % MATLAB 拟合出的函数

% 创建iddata对象
Ts = 0.5;  % 采样周期
data = iddata(tem, vlote, Ts, ...
              'Name', 'Heating_Furnace_Experiment', ...
              'OutputName', 'Temperature', ...
              'OutputUnit', '°C', ...
              'InputName', 'Voltage', ...
              'InputUnit', 'V', ...
              'TimeUnit', 'seconds');

% 验证数据对象
fprintf('===== 实验数据信息 =====\n');
fprintf('数据点数: %d\n', length(tem));
fprintf('总时长: %.2f秒\n', max(time));
fprintf('输入范围: %.1fV to %.1fV\n', min(vlote), max(vlote));
fprintf('温度范围: %.1f°C to %.1f°C\n\n', min(tem), max(tem));

% 确认传递函数对象存在并验证
if ~exist('sys1', 'var') || ~exist('sys2', 'var') || ~exist('sys3', 'var')
    error('传递函数sys1, sys2和sys3必须在工作区中');
end

% 设置一致的输入输出名称
sys1.InputName = {'Voltage'};
sys1.OutputName = {'Temperature'};
sys2.InputName = {'Voltage'};
sys2.OutputName = {'Temperature'};
sys3.InputName = {'Voltage'};
sys3.OutputName = {'Temperature'};

% ================= 计算三个模型的指标 =================
% 1. 拟合度 (Fit to estimation data)
[~, fit1] = compare(data, sys1);
[~, fit2] = compare(data, sys2);
[~, fit3] = compare(data, sys3);

% 2. 均方误差 (MSE)
[y1, ~] = compare(data, sys1);
mse1 = mean((data.OutputData - y1.OutputData).^2);

[y2, ~] = compare(data, sys2);
mse2 = mean((data.OutputData - y2.OutputData).^2);

[y3, ~] = compare(data, sys3);
mse3 = mean((data.OutputData - y3.OutputData).^2);

% 3. 最终预测误差 (FPE)
fpe1 = fpe(sys1);
fpe2 = fpe(sys2);
fpe3 = fpe(sys3);

% 4. 模型复杂性参数
params1 = numel(sys1.Report.Parameters.ParVector);  % 参数数量
params2 = numel(sys2.Report.Parameters.ParVector);
params3 = numel(sys3.Report.Parameters.ParVector);

% ================= 结果展示 =================
fprintf('\n===== 模型性能分析结果 =====\n');
fprintf('模型1 (参考模型):\n');
fprintf('  拟合度 (Fit): %.2f%%\n', fit1);
fprintf('  MSE: %.6f (°C)²\n', mse1);
fprintf('  FPE: %.6f\n', fpe1);
fprintf('  参数数量: %d\n', params1);

fprintf('\n模型2 (两点法):\n');
fprintf('  拟合度 (Fit): %.2f%%\n', fit2);
fprintf('  MSE: %.6f (°C)²\n', mse2);
fprintf('  FPE: %.6f\n', fpe2);
fprintf('  参数数量: %d\n', params2);

fprintf('\n模型3 (MATLAB拟合):\n');
fprintf('  拟合度 (Fit): %.2f%%\n', fit3);
fprintf('  MSE: %.6f (°C)²\n', mse3);
fprintf('  FPE: %.6f\n', fpe3);
fprintf('  参数数量: %d\n', params3);

% ================= 综合性能比较 =================
fprintf('\n===== 综合性能排名 =====\n');
% 按拟合度排序
[~, fit_rank] = sort([fit1, fit2, fit3], 'descend');
fprintf('拟合度排名: 模型%d > 模型%d > 模型%d\n', fit_rank(1), fit_rank(2), fit_rank(3));

% 按MSE排序
[~, mse_rank] = sort([mse1, mse2, mse3]);
fprintf('MSE排名: 模型%d < 模型%d < 模型%d\n', mse_rank(1), mse_rank(2), mse_rank(3));

% 按FPE排序
[~, fpe_rank] = sort([fpe1, fpe2, fpe3]);
fprintf('FPE排名: 模型%d < 模型%d < 模型%d\n', fpe_rank(1), fpe_rank(2), fpe_rank(3));

% ================= 可视化代码 =================

figure('Position', [100, 100, 1200, 800], 'Color', 'w')

% 1. 实际数据与模型输出对比
subplot(2, 2, [1, 2])
compare(data, sys1, sys2,sys3);
title('模型与实际数据对比')
grid on


% 2. 误差分布
subplot(2, 2, 3)
hold on

% 计算每个模型的误差
err1 = data.OutputData - y1.OutputData;
err2 = data.OutputData - y2.OutputData;
err3 = data.OutputData - y3.OutputData;

% 动态设置Bin范围（确保所有模型可见）
all_errors = [err1; err2; err3];
max_err = max(abs(all_errors)) * 1.1; % 增加10%边距
bin_edges = linspace(-max_err, max_err, 50); % 增加bins数量

histogram(err1, bin_edges, 'FaceColor', 'b', 'FaceAlpha', 0.6)
histogram(err2, bin_edges, 'FaceColor', 'r', 'FaceAlpha', 0.6)
histogram(err3, bin_edges, 'FaceColor', 'g', 'FaceAlpha', 0.6)
title('模型误差分布')
xlabel('误差 (°C)')
ylabel('频率')
legend({'参考模型', '两点法模型', 'MATLAB拟合模型'}, 'Location', 'best')
grid on
hold off

% 3. 残差序列
subplot(2, 2, 4)
hold on
plot(err1, 'b')
plot(err2, 'r')
plot(err3, 'g')
title('残差序列')
xlabel('样本点')
ylabel('残差值')
legend({'参考模型', '两点法模型', 'MATLAB拟合模型'}, 'Location', 'best') % 完整的图例标签
yline(0, 'k--', 'LineWidth', 1.5)
grid on
hold off