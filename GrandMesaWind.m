%% Grand Mesa Weather Data
%https://nsidc.org/data/snex_met/versions/1#anchor-2
cmap = csvread('D:\git-repository\GreenTrACS_MxRadar\colorMaps\RdYlBu.csv');
cmap = flipud(cmap);
%% Wind Data
addpath('D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\snotel\Windrose');
% winD = 'CAGMSWx_2019-10-01-2020-09-30.csv';
winDir = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\snotel\Wx';
f1 = 'SNEX_Met_MW_final_output.csv';
f2 = 'SNEX_Met_MM_final_output.csv';

W = readtable([winDir,'\',f1]);
WM = readtable([winDir,'\',f2]);
% Filter by Date
d1 = datetime(2019,10,01); d2 = datetime(2020,02,01);
rmvIx = isbetween(W.TIMESTAMP,d1,d2);
W = W(rmvIx,:);
WM = WM(rmvIx,:);

WMwindSpeed = WM.WSms_10ft_Avg(1:2250);
WwindSpeed = W.WSms_10ft_Avg(1:2250);
% Scale Factor 2.18
rmvIx = find(WMwindSpeed<0);
WMwindSpeed(rmvIx) = [];
mean(WwindSpeed)./mean(WMwindSpeed)
% Filter Bad Data
% W(3178:3181,:) = [];
% Correct Wind Direction
RefN = 0; RefE = 90;
dir              = mod((RefN-W.WindDir_10ft_D1_WVT)/(RefN-RefE)*90,360);  
%% Li & Pomeroy Wind Transport Threshold
Ut = 9.43 + 0.18.*W.AirTC_10ft_Avg + 0.0033.*W.AirTC_10ft_Avg.^2;
wix = find(W.WSms_10ft_Avg >= Ut - 1.97);
% wix = find(W.WSms_10ft_Avg >= Ut);


median(W.WSms_10ft_Avg(wix))
median(dir(wix))
mean(dir(wix))
%% Wind Rose
% properties = {'cmap',cmap,'titlestring',{'Mesa West (MW) Snow Transport Wind Rose';'Oct. 01, 2019 - Feb. 12, 2020'},'lablegend','Wind Speed (m/s)','anglenorth',0,'angleeast',90};

properties = {'cmap',cmap,'titlestring',{'Mesa West (MW) Snow Transport Wind Rose';'Oct. 01, 2019 - Feb. 12, 2020'},'lablegend','Wind Speed (m/s)','anglenorth',0,'angleeast',90,'vwinds',[0:10]};

wrh = WindRose(W.WindDir_10ft_D1_WVT(wix),W.WSms_10ft_Avg(wix),properties);


    % 'titlestring'         Cell/String.    {'Wind Rose';' '}       Figure title. It is recommended to include an empty line below the main string.
    % 'lablegend'           String.         'Wind speeds in m/s'    String that will appear at the top of the legend. Can be empty.
    % 'legendvariable'      String.         'W_S'  
    
set(gca,'fontsize',18,'fontweight','bold','fontname','serif');

% All Wind
properties = {'cmap',cmap,'titlestring',{'Mesa West (MW) Wind Rose';'Oct. 01, 2019 - Feb. 12, 2020'},'lablegend','Wind Speed (m/s)','anglenorth',0,'angleeast',90};


wrh = WindRose(W.WindDir_10ft_D1_WVT,W.WSms_10ft_Avg,properties);

    % 'titlestring'         Cell/String.    {'Wind Rose';' '}       Figure title. It is recommended to include an empty line below the main string.
    % 'lablegend'           String.         'Wind speeds in m/s'    String that will appear at the top of the legend. Can be empty.
    % 'legendvariable'      String.         'W_S'  
    
set(gca,'fontsize',18,'fontweight','bold','fontname','serif');