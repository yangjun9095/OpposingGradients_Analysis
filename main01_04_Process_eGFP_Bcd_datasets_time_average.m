function main01_04_Process_eGFP_Bcd_datasets_time_average
%% DESCRIPTION
% This script is for processing the Bcd datasets.
% 1) Synchronize the datasets with the beginning of NC13
% 2) Background subtraction
% 3) Average across multiple embryos
% 3-1) time-average for certain time windows.
% 4) Save these for future usage.

%% Load the datasets
% I have two choices, 
% (1) My own datasets that I took with Paul, 40% AP bin
% (2) Liz & Jonathan's datasets

% For now, let's use my datasets since they were taken with pretty much the
% same condition as 
BcdData = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Bcd-Averaged.mat')

%% Extract useful fields
BcdTime = BcdData.ElapsedTime;
BcdNC13 = BcdData.nc13;
BcdNC14 = BcdData.nc14;

BcdFluo = BcdData.MeanVectorAP;
BcdFluoSD = BcdData.SDVectorAP;

%% Time-average for different time-windows
% First, let's define different time windows for averaging
% 0-5min, 0-10min, 0-20 min into NC14 as the first round.

tWindows{1} = 0:5; % min
tWindows{2} = 0:10; % min
tWindows{3} = 0:20; % min


end