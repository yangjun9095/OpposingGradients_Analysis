function main12_02_Compare_Averaged_Spotfluo_overALLnuclei_Allconstructs
%% DESCRIPTION
% This script is to compare 
% 1) Mean spot fluo over ALL nuclei

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

%% Color map
% Use cbrewer (https://colorbrewer2.org/#)
% Create 10 distintive colors for now.
% I might need to revisit these color schemes so that we can recognize
% these better.
[colormap] = cbrewer('qual', 'Set1', 12);
%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
purple = [171,133,172]/255;
colorDict.purple = (4*purple - [1,1,1])/3;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
brown = [179,155,142]/255;
colorDict.brown = (2*brown - [1,1,1])/1;

colorDict.darkgreen = [126,157,144]/255;
colorDict.lightgreen = [205,214,209]/255;
thickpink = [132,27,69]/255;
colorDict.thickpink = (3*thickpink + [1,1,1]) / 4; % adding white

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.magenta; colorDict.thickpink;...
                colormap(10,:); colormap(11,:)]; 
            
%% Calculate the averaged spot fluorescence over ALL nuclei


end