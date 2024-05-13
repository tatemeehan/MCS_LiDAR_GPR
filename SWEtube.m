%% Process SWE Tube
clear; close all; clc
addpath('C:\Users\RDCRLTGM\Desktop\git-repository\preprocessing')
% Load tubular data
T = readtable('E:\SNOWWI\SWEtube\Gm0324_swetube.xlsx');
% Convert to UTM
[x,y,zone] = deg2utm(T.Latitude,T.Longitude);
% Calculate Density
area = 30; % Cross Sectional Area 30 cm^2
% volume = area.*T.Depth;
volume = area.*T.Depth_cm_;
% density = T.snowtube_g_./volume; %g/cc
density = T.sweTube_g_./volume; %g/cc
density = density.*1000; % kg/m^2
% Cacluate SWE
% SWE = (T.Depth./100).*density; % mm
SWE = (T.Depth_cm_./100).*density; % mm

T.Date_Time = datetime(T.Date_Time,'Format','yyyy-MM-dd HH:mm:ss');
isTableCol = @(t, colname) ismember(colname, t.Properties.VariableNames);
isFt = isTableCol(T,'Altitude_ft_');
if isFt
    % Convert to meters
    T.Altitude_m_ = T.Altitude_ft_.*0.3048;
end
% Write output .csv
% Tout = table(T.Date_Time,T.Longitude,T.Latitude,T.Altitude_m_,x,y,zone,T.Depth,density,SWE,'VariableNames',{'Date_Time','Longitude_DD','Latitude_DD','Elevation_masl','UTM_X','UTM_Y','UTM_Zone','Depth_cm','Density_kgm3','SWE_mm'});
Tout = table(T.Date_Time,T.Longitude,T.Latitude,T.Altitude_m_,x,y,zone,T.Depth_cm_,density,SWE,'VariableNames',{'Date_Time','Longitude_DD','Latitude_DD','Elevation_masl','UTM_X','UTM_Y','UTM_Zone','Depth_cm','Density_kgm3','SWE_mm'});

Tout = sortrows(Tout,'Date_Time');
writetable(Tout,'E:\SNOWWI\SWEtube\GM_2024-03-25_2024-03-28_SWE.csv')