%% Combine and Standardize TWT
clear; close all; clc;
%% Load all .csv data
T1 = readtable('E:\MCS\MCS122123\GPR\processed\MCS122123-TWT.csv');
T2 = readtable('E:\MCS\MCS021324\GPR\processed\GPR-TWT.csv');
T3 = readtable('E:\MCS\MCS031524\GPR\processed\MCS031524-TWT.csv');
T4 = readtable('E:\MCS\MCS031824\GPR\processed\MCS031824-TWT.csv');
T5 = readtable('E:\MCS\MCS032024\GPR\processed\MCS032024-TWT.csv');
T6 = readtable('E:\MCS\MCS041824\GPR\processed\MCS041824-TWT.csv');
%% Standardize Each DataSet
% T1
T1.zTWT = (T1.TWT-mean(T1.TWT))./std(T1.TWT);
T1.meanTWT = ones(height(T1),1).*mean(T1.TWT);
T1.stdTWT =  ones(height(T1),1).*std(T1.TWT);
% T2
T2.zTWT = (T2.TWT-mean(T2.TWT))./std(T2.TWT);
T2.meanTWT = ones(height(T2),1).*mean(T2.TWT);
T2.stdTWT =  ones(height(T2),1).*std(T2.TWT);
% T3
T3.zTWT = (T3.TWT-mean(T3.TWT))./std(T3.TWT);
T3.meanTWT = ones(height(T3),1).*mean(T3.TWT);
T3.stdTWT =  ones(height(T3),1).*std(T3.TWT); 
% T4
T4.zTWT = (T4.TWT-mean(T4.TWT))./std(T4.TWT);
T4.meanTWT = ones(height(T4),1).*mean(T4.TWT);
T4.stdTWT =  ones(height(T4),1).*std(T4.TWT);
% T5
T5.zTWT = (T5.TWT-mean(T5.TWT))./std(T5.TWT);
T5.meanTWT = ones(height(T5),1).*mean(T5.TWT);
T5.stdTWT =  ones(height(T5),1).*std(T5.TWT);
% T6
T6.zTWT = (T6.TWT-mean(T6.TWT))./std(T6.TWT);
T6.meanTWT = ones(height(T6),1).*mean(T6.TWT);
T6.stdTWT =  ones(height(T6),1).*std(T6.TWT);

%% Combine all Tables
T = [T1;T2;T3;T4;T5;T6];
%% Write Output
writetable(T,'E:\MCS\MCS2024-TWT.csv');
%% Test
synthTWT = (T.zTWT.*T2.stdTWT(1))+T2.meanTWT(1);