function main12_Compare_Averaged_SpotFluo_AllConstructs(varargin)
%% DESCRIPTION
% This script is generate a structure with averaged datasets per construct, using AverageDatasets.m
% We will run AverageDatasets for all DataTypes, then load into a master
% structure, then save into a variabile, so that we can load at any
% plotting script.


%% Average datasets using AverageDatasets.m

% Define the file path
DropboxPath = 'S:\YangJoon\Dropbox\OpposingGradient';

AverageDatasets('r0-new','NC',13,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r1-new','NC',13,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r2-new','NC',13,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r3-new','NC',13,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);

% Other variants
%% Average datasets using AverageDatasets.m (Runt nulls)
% Note that we mostly have NC14 only,
AverageDatasets('r0_RuntNull','NC',14,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r1_RuntNull','NC',14,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r1_mid_RuntNull','NC',14,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r1_close_RuntNull','NC',14,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r2_RuntNull','NC',14,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r2_close_RuntNull','NC',14,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r2_far_RuntNull','NC',14,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
AverageDatasets('r3_RuntNull','NC',14,'savePath',[DropboxPath,filesep,'OpposingGradients_ProcessedData']);
%% Load datasets into a structure (a master one to save)

% This assumes that the DataType is the name of the constructs in the
% DataStatus.xlsx tab, for example, r0-new

DropboxPath = 'S:/YangJoon/Dropbox/OpposingGradient';
filePath = [DropboxPath,filesep,'OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'];

AveragedData{1,1} = 'DataType';
AveragedData{1,2} = 'averagedData';

DataTypesForFit = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};

% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};
                
                
for i=1:length(DataTypesForFit)
    AveragedData{i+1,1} = DataTypesForFit{i}; % DataType
    AveragedData{i+1,2} = load([filePath, filesep, DataTypesForFit{i} ,'.mat']);
end


%% Save the structure, compiledData for future usage (in plotting scripts)
save('S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedData.mat',...
        'AveragedData')


end