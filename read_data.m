% @brief  读取阶跃响应数据到工作区  read_data.m    
% @author 23010341 杳泽

clc;clear;

data=readtable("./data.csv") ;   %读取csv文件
time = data.time;
tem=data.temperature;
vlote=data.volte;   %时间 温度和电压


