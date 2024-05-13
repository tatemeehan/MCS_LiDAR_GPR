%% Process Depth Probe
clear; close all; clc
addpath C:\Users\RDCRLTGM\Desktop\git-repository\preprocessing\
% Load tubular data
% T = readtable('E:\MCS\MCS021324\Probe\MCS021324_5points.csv');
T = readtable('E:\SNOWWI\SWEtube\Gm0324_depth.xlsx');
T.Depth = T.Depth_cm_;
% Convert to UTM
[x,y,zone] = deg2utm(T.Latitude,T.Longitude);
% Convert to meters
Elevation = T.Altitude_ft_.*0.3048;
% Check if Time Zone is incorrect!
isChangeTZ = 0;
if isChangeTZ
    T.Date_Time = datetime(T.Date_Time,'Format','yyyy-MM-dd HH:mm:ss','TimeZone','America/Anchorage');
    T.Date_Time.TimeZone = 'America/Denver';
else
        T.Date_Time = datetime(T.Date_Time,'Format','yyyy-MM-dd HH:mm:ss','TimeZone','America/Denver');
end

% Write Table
Tout = table(T.Date_Time,T.Longitude,T.Latitude,Elevation,x,y,zone,T.Depth,T.Description,'VariableNames',{'Date_Time','Longitude_DD','Latitude_DD','Elevation_masl','UTM_X','UTM_Y','UTM_Zone','Depth_cm','Comments'});
% Ensure Comments Section is Type Cell
Tout.Comments = cell(Tout.Comments);
% Sort Data by Time of Day
Tout = sortrows(Tout,'Date_Time');

writetable(Tout,'E:\SNOWWI\SWEtube\GM0324_depth.csv')

%% Combine multiple files
isCombine = 0;
f1 = 'E:\MCS\MCS021324\Probe\MCS20240213_DEPTH_1.csv';
f2 = 'E:\MCS\MCS021324\Probe\MCS20240213_DEPTH_2.csv';
% f3 = 'E:\MCS\MCS021324\Probe\MCS20240213_DEPTH_3.csv';
if isCombine
    T1 = readtable(f1);
    % Ensure Comments Section is String
    tmp = cellstr(string(T1.Comments));
    T1.Comments = tmp;
    T2 = readtable(f2);
    % For no good reason matlab reads empty text cells as NaN (numerical)
    tmp = cellstr(string(T2.Comments));
    T2.Comments = tmp;
%     T3 = readtable(f3);
%     % Quite annoying matlab
%     tmp = cellstr(string(T3.Comments));
%     T3.Comments = tmp;
    % Write Combined Table
    Tout = [T1;T2];
    Tout = sortrows(Tout,'Date_Time');
    writetable(Tout,'E:\MCS\MCS021324\Probe\MCS20240213_DEPTH.csv')
end
