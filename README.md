# ECUST 自动控制原理 实验报告
## 任务说明（任务B）
1.	利用提供的temperature.csv，使用经典辨识方法（两点法等）辨识出加热炉对象的输入输出模型，并验证辨识方法的准确性。
2.	根据辨识出的模型，设计 PID 控制器对系统进行闭环控制，并使用智能优化算法调整PID控制器参数，使加热炉温度稳定控制在35℃。同时，计算并分析动态指标和稳态指标。
3.	提交的材料需包含完整的问题分析、算法设计过程以及最终的控制效果。

## 使用工具
matlab system identification toolbox, Global Optimization tool box

python math库,matplotlib库

## 各个文件功能即操作步骤说明

1. `read_data.m` 读取数据到工作区内为time,tem,vlote
2. 利用matlab system identification 导入工作区的time tem vlote，并进行传递函数的拟合，完成后导出回工作区并save("./gs_result.mat");
3. `valid_gs.m`利用工作区内的传递函数和数据（运行时会load gs_result.mat），计算MSE，FPE，和拟合度，并可视化对比不同传递函数的准确性
4. `inheritance_pid.m` 利用遗传算法整定PID参数，并计算超调量、调节时间、峰值时间等参数，并可视化；结果保存为`inheritance_pid*.mat`
5. `PSO.m`利用粒子群算法整定PID参数，并计算超调量、调节时间、峰值时间等参数，并可视化;结果保存为`pso_optimized-*.mat`
