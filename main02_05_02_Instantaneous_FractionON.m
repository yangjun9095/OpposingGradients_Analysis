function main02_05_02_Instantaneous_FractionON
% This is a script for investigating how the instantaneous fraction on
% changes over time and AP bins (for different constructs).

% Caveats : The spot segmentaiton is not "perfect". So, we'd need to be
% careful in terms of interpreting this.

%% Step1. Process the datasets
% I'm using AverageDatasets.m script to synchronize and average multiple
% embryos from one construct. This will average instantaneous fraction on
% as FractionON
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
dataTypes = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3'};

% Run the AverageDatasets for all Data Types.
for i=4:length(dataTypes)
    AverageDatasets(dataTypes{i},'savePath',filePath)
    pause(1);
end

%% Step2. Load the processed datasets
filePath = '/Users/yangjoonkim/Dropbox/OpposingGradient/OpposingGradients_ProcessedData/AveragedDatasets_Feb2020';
r0Data = load([filePath, filesep,dataTypes{1},'.mat']);

r3Data = load([filePath, filesep,dataTypes{4},'.mat']);

r1mid_Data = load([filePath, filesep,dataTypes{6},'.mat']);

Data = r1mid_Data;

% Extract the instantaneous Fraction ON (for individual embryos)
FractionON_individual = Data.FractionON_individual;

tempVar = Data.MeanVectorAP_individual;
numEmbryos = length(tempVar(1,1,:));
ElapsedTime = Data.ElapsedTime;

% Average using nanmean
FractionON_average = nanmean(FractionON_individual,3);
FractionON_SEM = nanstd(FractionON_individual,0,3)./sqrt(numEmbryos);

% Heatmap 
% First, generate strings of APbins and time points.
% X labels
for i=1:41
    xLabels{i} = num2str((i-1)*2.5);
end
% Y labels
for j=1:length(ElapsedTime)
    yLabels{j} = num2str(ElapsedTime(j));
end
figure(2)
h = heatmap(FractionON_average,'Colormap',viridis)
xlabel('AP axis (EL)')
ylabel('time(min)')

%% Sanity check for the r1-mid
MaxFractionON_individual = squeeze(nanmax(FractionON_individual(Data.nc14+3:end,:,:),[],1));

hold on
for i=1:numEmbryos
    plot(0:0.025:1, MaxFractionON_individual(:,i))
    pause
end

% embryo #3,4,5,6 seem to have reasonable fraction on
end