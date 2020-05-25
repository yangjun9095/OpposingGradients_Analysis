function main_compare_instant_MS2_fluo(DataType)
%% Description
% This script is to analyze the instantaneous MS2 fluorescence for each DataType,
% then, save into .mat structure that could be used to generate different
% types of plots, heatmaps, etc.
% The goal is to grab all datasets from each DataType(construct), then
% synchronize the Fraction of active nuclei, then average it over time.
% Use the datasets processed by AverageDatasets.m

% Define the file path
DropboxPath = 'S:/YangJoon/Dropbox';
filePath = [DropboxPath,filesep,'OpposingGradient/OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'];
%% Average datasets using AverageDatasets.m
% Use the AverageDatasets.m script
%[000],[111] with Runt has nc12,13, and 14, thus we will average from nc13.
% AverageDatasets('r0-new','NC',13,'savePath',filePath)
% AverageDatasets('r3-new','NC',13,'savePath',filePath)
% % For the [000],[111] Runt nulls, we mostly have only nc14, so let's just
% % average nc14 only.
AverageDatasets('r0_RuntNull','NC',14,'savePath',filePath)
% AverageDatasets('r3_RuntNull','NC',14,'savePath',filePath)
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

% Runt null datasets ([0,0,0] and [1,1,1]) at nc13
% I'm loading these datsets separately just because we have only one
% dataset per construct that has intact nc13 at this point...(before the
% quarantine period).
r0NullData_nc13 = load([filePath, filesep, 'r0_RuntNull_NC13.mat']);
r3NullData_nc13 = load([filePath, filesep, 'r3_RuntNull_NC13.mat']);
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
% FractionON_individual has dimension of (Frames x APbins x embryos)
% We can use this FractionON_individual to calculate the error bar

%% Part1. compare different constructs
% First, we want to compare the [000] with or without Runt protein
% Compare the instantaneous fraction on in NC14 for [000] and [000], Runt
% null. 
% Caveats : [000], Runt null doesn't have nuclear marker, thus its fraction
% on is basically spot density (number of spots / APbin area)
%% Make sure that the synchronization and the run of AverageDatasets is done for the correct DataStatus.xlsx


% check for [000], Runt null
AverageDatasets('r0_RuntNull','NC',14,'savePath',filePath)
r0NullData = load([filePath, filesep, 'r0_RuntNull.mat']);
AveragedData{9} = r0NullData;

%%
data = AveragedData{9};
APbin = 10;
[~,~,numEmbryos] = size(data.MeanVectorAP_individual);

for APbin=9:16
    clf
    hold on
    for i=1:numEmbryos
    errorbar(data.ElapsedTime, data.MeanVectorAP_individual(:,APbin,i),data.SDVectorAP_individual(:,APbin,i))
    %pause
    end
    legend('embryo1','embryo2')%,'embryo3')
	APbin
    pause
end
%% Extract useful fields and calculate the mean and STD of fraction on over individual embryos

% From the AveragedData structure, note that I'd need to pull specific
% elements, then the calculation is repetitive. 

% flag the null datasets (as they might need somewhat different treatment)
nullFlag = [0,0,0,0,0,0,0,0,1,1];
nameFlag = {'[000]','[100]','[011]','[111]','[001]','[010]','[110]','[101]',...
            '[000], Runt null','[111], Runt null'};
        
elongateFlag = [0,0,0,0,0,0,0,0,1,0]; % elongate the time vector for visualization

% Define the AP bins to plot
APaxis = 0:0.025:1;
APAxis = [0.15:0.025:0.6]*100;
APrange = 7:25;
nAPbins = length(APaxis);

k=1; % counter for subplot

for i=[1,9,4,10] % [000],[111],[000]-Runt null, [111]-Runt null
    % extract useful field, do calculation, then plot the heatmaps using the subplot
    % APbin x Time)
    data = AveragedData{i};
    % time, nc14, MeanVectorAP, SDVectorAP, etc.
    Time = data.ElapsedTime;
    nc14 = data.nc14;
    fluo_mean = nanmean(data.MeanVectorAP_individual,3);
    fluo_STD = nanstd(data.MeanVectorAP_individual,0,3);
    [~,~,numEmbryos] = size(data.MeanVectorAP_individual);
    fluo_SEM = fluo_STD./sqrt(numEmbryos);
    
    % some processing for visualization
    % 1) filter out the nans to be zeros
    fluo_mean(isnan(fluo_mean)) = 0;
    % 2) [000], Runt nulls only up to 30 minutes.
    % add data points of nans up to 45 minutes
    % We only have up until ~30 minutes for the [000], Runt nullfor the plotting purpose, let's add some zero values after the time point
    % in here, let's hardcode for the [000] datasets, ~40 min into nc14
    % we can either use the elongateFlag or tMax

    if elongateFlag(i)
        tMax = 48; % min
        tRes = median(diff(Time));
        nLength = tMax/tRes;
        nLength_add = ceil(nLength - length(Time));
        [~,nAPbins] = size(data.MeanVectorAP);
        vec_add = nan(nLength_add,nAPbins);
        tvec = tRes*[1:nLength_add] + Time(end);

        Time = [Time, tvec];
        fluo_mean = [fluo_mean; vec_add];
    end

    
    % Plot the heatmap (using the subplot)
%     subplot(2,2,k)
%     pcolor(APAxis,Time(nc14:end) - Time(nc14),fluo_mean(nc14:end,APrange))
%     colormap(viridis)
%     colorbar
%     caxis([0 600])
%     xlabel('position (% embryo length)')
%     set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
%     ylabel('time into nc14 (min)')
%     set(gca,'ytick',0:10:Time(end)-Time(nc14),'yticklabels',[0:10:Time(end)-Time(nc14)])
%     title(nameFlag{i})
%     StandardFigure(gcf,gca)
    
    % Save the calculated fields to a cell for future usage.
    Fluo_mean{k} = fluo_mean;
    time{k} = Time;
    NC14{k} = nc14;
    k=k+1;
end

%% Calculate the ratio (fold-change)
% Rationale : I'd like to see the fold-change of the instantaneous MS2 fluo for the
% [111] construct with/without Runt. Then, I want to see if there's any
% clear trend that we can connect with the Runt protein concentration
% dynamics.

% Use the fields calculated above
% Fluo_mean{k}, time{k}, k=3,4
% step1 : trim the matrices such that they have the same frame length.

tLength1 = length(time{3}) - NC14{3};
tLength2 = length(time{4}) - NC14{4};

frameLength = min(tLength1, tLength2);

r3_fluo_mean = Fluo_mean{3}(NC14{3}:NC14{3}+frameLength,:);
r3Null_fluo_mean = Fluo_mean{4}(NC14{4}:NC14{4}+frameLength,:);

% How to put an error bar in here, using some confidence interval?
fc_fluo_111 = r3_fluo_mean./r3Null_fluo_mean;

pcolor(APAxis,time{4}(1:frameLength),...
        fc_fractionON_111(1:frameLength,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc14 (min)')
%set(gca,'ytick',0:10:r3Null_Time(end),'yticklabels',[0:10:r3Null_Time(end)])
set(gca,'ytick',0:10:45,'yticklabels',[0:10:45])
caxis([0 1])
title('[111]-fold-change')
StandardFigure(gcf,gca)

% save the plot
FigPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\Figures_OpposingGradients\Data\AveragedMS2_TimeTracesV2\Heatmaps_RuntNulls';
saveas(gcf, [FigPath, filesep, 'fold_change_meanfluo_ONnuclei_[111]_nc14.tif'])
saveas(gcf, [FigPath, filesep, 'fold_change_meanfluo_ONnuclei_[111]_nc14.pdf'])

%% heatmap for the fold-change of averaged MS2 fluo over ALL nuclei
% As this quantity is a multiplication of averaged fluo(over ON nuclei) and
% fraction on
fc_fluo_AllNuclei_111 = fc_fluo_111.*fc_fractionON_111;

pcolor(APAxis,time{4}(1:frameLength),...
        fc_fluo_AllNuclei_111(1:frameLength,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc14 (min)')
%set(gca,'ytick',0:10:r3Null_Time(end),'yticklabels',[0:10:r3Null_Time(end)])
set(gca,'ytick',0:10:45,'yticklabels',[0:10:45])
caxis([0 1])
title('[111]-fold-change')
StandardFigure(gcf,gca)

% save the plot
FigPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\Figures_OpposingGradients\Data\AveragedMS2_TimeTracesV2\Heatmaps_RuntNulls';
saveas(gcf, [FigPath, filesep, 'fold_change_meanfluo_AllNuclei_[111]_nc14.tif'])
saveas(gcf, [FigPath, filesep, 'fold_change_meanfluo_AllNuclei_[111]_nc14.pdf'])
%% plot the Runt concentration gradient dynamics
% load the dataset (post-BG subtraction)
RuntData = load([filePath, filesep, 'Runt-1min-200Hz-mixed_BGsubtracted-Averaged.mat']);

% plot a heatmap of Runt concentration (AU) over AP and time in nc14
Runt_Time = RuntData.ElapsedTime;
Runt_nc14 = RuntData.nc14;
Runt_fluo_mean = RuntData.MeanVectorAP_BGsubtracted;
Runt_fluo_SEM = RuntData.SEVectorAP;


% plot the heatmap
% for plotting, let's plot until 45 min, to match with the fraction on
% heatmap.
frameEnd = Runt_nc14 + 45;
pcolor(APAxis,Runt_Time(Runt_nc14:frameEnd) - Runt_Time(Runt_nc14),...
            Runt_fluo_mean(Runt_nc14:frameEnd,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc14 (min)')
set(gca,'ytick',0:10:Runt_Time(frameEnd)-Runt_Time(Runt_nc14),'yticklabels',[0:10:Runt_Time(frameEnd)-Runt_Time(Runt_nc14)])
caxis([100 400])

title('Runt protein')
StandardFigure(gcf,gca)

% save the plot
saveas(gcf, [FigPath, filesep, 'Runt_protein_BGsubtracted_nc14.tif'])
saveas(gcf, [FigPath, filesep, 'Runt_protein_BGsubtracted_nc14.pdf'])

end