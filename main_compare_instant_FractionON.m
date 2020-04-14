function main_compare_instant_FractionON(DataType)
%% Description
% This script is to analyze the instantaneous fraction on for each DataType,
% then, save into .mat structure that could be used to generate different
% types of plots, heatmaps, etc.
% The goal is to grab all datasets from each DataType(construct), then
% synchronize the Fraction of active nuclei, then average it over time.
% Use the datasets processed by AverageDatasets.m

%% Average datasets using AverageDatasets.m
%% Load datasets
% This assumes that the DataType is the name of the constructs in the
% DataStatus.xlsx tab, for example, r3-new

DropboxPath = 'S:/YangJoon/Dropbox';
filePath = [DropboxPath,filesep,'OpposingGradient/OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'];

r0Data = load([filePath, filesep, 'r0-new.mat']);
r1Data = load([filePath, filesep, 'r1-new.mat']);
r2Data = load([filePath, filesep, 'r2-new.mat']);
r3Data = load([filePath, filesep, 'r3-new.mat']);

r1closeData = load([filePath, filesep, 'r1-close.mat']);
r1midData = load([filePath, filesep, 'r1-mid.mat']);

r2closeData = load([filePath, filesep, 'r2_1+2.mat']);
r2farData = load([filePath, filesep, 'r2_1+3.mat']);

% Runt null datasets ([0,0,0] and [1,1,1])
r0NullData = load([filePath, filesep, 'r0_RuntNull.mat']);
r3NullData = load([filePath, filesep, 'r3_RuntNull.mat']);
%% Put them in a structure for convenience
AveragedData{1} = r0Data;
AveragedData{2} = r1Data;
AveragedData{3} = r2Data;
AveragedData{4} = r3Data;
AveragedData{5} = r1closeData;
AveragedData{6} = r1midData;
AveragedData{7} = r2closeData;
AveragedData{8} = r2farData;
AveragedData{9} = r0NullData;
AveragedData{10} = r3NullData;

%% Info
% The filed "FractionON" is already averaged
% FractionON_individual
%% Generate plots

end