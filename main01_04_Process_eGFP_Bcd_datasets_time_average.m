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

end