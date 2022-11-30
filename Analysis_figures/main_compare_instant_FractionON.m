function main_compare_instant_FractionON(DataType)
%% Description
% This script is to analyze the instantaneous fraction on for each DataType,
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
AverageDatasets('r0-new','NC',13,'savePath',filePath)
AverageDatasets('r3-new','NC',13,'savePath',filePath)
% For the [000],[111] Runt nulls, we mostly have only nc14, so let's just
% average nc14 only.
AverageDatasets('r0_RuntNull','NC',14,'savePath',filePath)
AverageDatasets('r3_RuntNull','NC',14,'savePath',filePath)
%% Load datasets
% This assumes that the DataType is the name of the constructs in the
% DataStatus.xlsx tab, for example, r3-new

DropboxPath = 'S:/YangJoon/Dropbox/OpposingGradient';
filePath = [DropboxPath,filesep,'OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'];

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
% AveragedData{1} = r0Data;
% AveragedData{2} = r1Data;
% AveragedData{3} = r2Data;
% AveragedData{4} = r3Data;
% AveragedData{5} = r1closeData;
% AveragedData{6} = r1midData;
% AveragedData{7} = r2closeData;
% AveragedData{8} = r2farData;
% AveragedData{9} = r0NullData;
% AveragedData{10} = r3NullData;

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
%% Extract useful fields and calculate the mean and STD of fraction on over individual embryos
% [000]
r0_Time = r0Data.ElapsedTime;
r0_nc14 = r0Data.nc14;
%r0_FractionON_mean = r0Data.FractionON;
r0_FractionON_mean = nanmean(r0Data.FractionON_individual,3);
r0_FractionON_STD = nanstd(r0Data.FractionON_individual,0,3);
[~,~,r0_numEmbryos] = size(r0Data.FractionON_individual);
r0_FractionON_SEM = r0_FractionON_STD./sqrt(r0_numEmbryos);

% [111]
r3_Time = r3Data.ElapsedTime;
r3_nc14 = r3Data.nc14;
r3_FractionON_mean = nanmean(r3Data.FractionON_individual,3);
r3_FractionON_STD = nanstd(r3Data.FractionON_individual,0,3);
[~,~,r3_numEmbryos] = size(r3Data.FractionON_individual);
r3_FractionON_SEM = r3_FractionON_STD./sqrt(r3_numEmbryos);

% [000], Runt null
r0Null_Time = r0NullData.ElapsedTime;
r0Null_nc14 = r0NullData.nc14;
r0Null_FractionON_mean = nanmean(r0NullData.FractionON_individual,3);
r0Null_FractionON_STD = nanstd(r0NullData.FractionON_individual,0,3);
[~,~,r0Null_numEmbryos] = size(r0NullData.FractionON_individual);
r0Null_FractionON_SEM = r0Null_FractionON_STD./sqrt(r0Null_numEmbryos);

% Since the fraction on is #spots/APbinArea, we want to rescale this using
% the 95% percentile value of FractionON
[row,col] = size(r0Null_FractionON_mean);
r0Null_scale = prctile(reshape(r0Null_FractionON_mean,[1,row*col]),99);
r0Null_FractionON_mean_rescaled = r0Null_FractionON_mean./r0Null_scale;
r0Null_FractionON_SEM_rescaled = r0Null_FractionON_SEM./r0Null_scale;

% We only have up until ~30 minutes for the [000], Runt nullfor the plotting purpose, let's add some zero values after the time point
% in here, let's hardcode for the [000] datasets, ~40 min into nc14
r0Null_tMax = 45; % min
tRes = median(diff(r0Null_Time));
nLength = r0Null_tMax/tRes;
nLength_add = ceil(nLength - length(r0Null_Time));
[~,nAPbins] = size(r0NullData.FractionON);
vec_add = nan(nLength_add,nAPbins);
tvec = tRes*[1:nLength_add] + r0Null_Time(end);

r0Null_Time_elongated = [r0Null_Time, tvec];
r0Null_FractionON_mean_rescaled_elongated = [r0Null_FractionON_mean_rescaled; vec_add];


% [111], Runt null
r3Null_Time = r3NullData.ElapsedTime;
r3Null_nc14 = r3NullData.nc14;
r3Null_FractionON_mean = nanmean(r3NullData.FractionON_individual,3);
r3Null_FractionON_STD = nanstd(r3NullData.FractionON_individual,0,3);
[~,~,r3Null_numEmbryos] = size(r3NullData.FractionON_individual);
r3Null_FractionON_SEM = r3Null_FractionON_STD./sqrt(r3Null_numEmbryos);

% Since the fraction on is #spots/APbinArea, we want to rescale this using
% the 95% percentile value of FractionON
[row,col] = size(r3Null_FractionON_mean);
r3Null_scale = prctile(reshape(r3Null_FractionON_mean,[1,row*col]),99);
r3Null_FractionON_mean_rescaled = r3Null_FractionON_mean./r3Null_scale;
r3Null_FractionON_SEM_rescaled = r3Null_FractionON_SEM./r3Null_scale;
%% Plot the different Fraction ONs(heatmap of Time x APbins)

% Set the x-axis limits (AP bins)
APbin_vec = 15:2.5:60;
APbin_lb = 15;   % % of EL
xlim_lb = find(APbin_vec == APbin_lb);
APbin_ub = 60;
xlim_ub = find(APbin_vec == APbin_ub);   % of EL

% We will compare the nc14 for now.
APaxis = 0:0.025:1;
APAxis = [0.15:0.025:0.6]*100;
APrange = 7:25;
nAPbins = length(APaxis);

subplot(2,2,1)
pcolor(APAxis,r0_Time(r0_nc14:end) - r0_Time(r0_nc14),r0_FractionON_mean(r0_nc14:end,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc14 (min)')
set(gca,'ytick',0:10:r0_Time(end)-r0_Time(r0_nc14),'yticklabels',[0:10:r0_Time(end)-r0_Time(r0_nc14)])
title('[000]')
StandardFigure(gcf,gca)


subplot(2,2,2)
pcolor(APAxis,r0Null_Time_elongated ,r0Null_FractionON_mean_rescaled_elongated(:,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc14 (min)')
%set(gca,'ytick',0:10:r3Null_Time(end),'yticklabels',[0:10:r3Null_Time(end)])
set(gca,'ytick',0:10:45,'yticklabels',[0:10:45])
title('[000], Runt null')
StandardFigure(gcf,gca)


subplot(2,2,3)
pcolor(APAxis,r3_Time(r3_nc14:end) - r3_Time(r3_nc14),r3_FractionON_mean(r3_nc14:end,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc14 (min)')
set(gca,'ytick',0:10:r3_Time(end)-r3_Time(r3_nc14),'yticklabels',[0:10:r3_Time(end)-r3_Time(r0_nc14)])
title('[111]')
StandardFigure(gcf,gca)


subplot(2,2,4)
pcolor(APAxis,r3Null_Time ,r3Null_FractionON_mean_rescaled(:,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc14 (min)')
set(gca,'ytick',0:10:r3Null_Time(end),'yticklabels',[0:10:r3Null_Time(end)])
title('[111], Runt null')
StandardFigure(gcf,gca)

% Save the figure

%% To-do
%% 1) Take specific AP bins, then compare directly for the Runt WT vs. nulls.

%% First, [000] for Runt WT and Runt null
% Let's take an AP bin for now.
APbin1 = 17; %  

% nc14
fractionON_timetrace_r0 = figure(1)
hold on
% [000], Runt WT
tWindow_r0 = r0_nc14:length(r0_Time);
errorbar(r0_Time(tWindow_r0) - r0_Time(r0_nc14),...
            r0_FractionON_mean(tWindow_r0,APbin1),...
            r0_FractionON_SEM(tWindow_r0,APbin1))
        
% [000], Runt nulls
tWindow_r0Null = r0Null_nc14:length(r0Null_Time);
errorbar(r0Null_Time(tWindow_r0Null) - r0Null_Time(r0Null_nc14),...
            r0Null_FractionON_mean_rescaled(tWindow_r0Null,APbin1),...
            r0Null_FractionON_SEM_rescaled(tWindow_r0Null,APbin1))
        
% plot formatting
ylim([0 1.4])
legend('[000]','[000],Runt null')
xlabel(' time into nc14 (min)')
ylabel('fraction of active nuclei')
title('[000]')
StandardFigure(fractionON_timetrace_r0,fractionON_timetrace_r0.CurrentAxes)

FigPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\Figures_OpposingGradients\Data\FractionON_instant';
saveas(fractionON_timetrace_r0, [FigPath, filesep, 'fractionON_time_NC14_[000]_',num2str((APbin1-1)*2.5),'%_RuntNulls.tif'])
saveas(fractionON_timetrace_r0, [FigPath, filesep, 'fractionON_time_NC14_[000]_',num2str((APbin1-1)*2.5),'%_RuntNulls.pdf'])

%% Second, [111] for Runt WT and Runt null
% Let's take an AP bin for now.
APbin1 = 17; %  

% nc14
fractionON_timetrace_r3 = figure(1)
hold on
% [111], Runt WT
tWindow_r3 = r3_nc14:length(r3_Time);
errorbar(r3_Time(tWindow_r3) - r3_Time(r3_nc14),...
            r3_FractionON_mean(tWindow_r3,APbin1),...
            r3_FractionON_SEM(tWindow_r3,APbin1))
        
% [111], Runt nulls
tWindow_r3Null = r3Null_nc14:length(r3Null_Time);
errorbar(r3Null_Time(tWindow_r3Null) - r3Null_Time(r3Null_nc14),...
            r3Null_FractionON_mean_rescaled(tWindow_r3Null,APbin1),...
            r3Null_FractionON_SEM_rescaled(tWindow_r3Null,APbin1))
        
% plot formatting
ylim([0 1.4])
legend('[111]','[111],Runt null')
xlabel(' time into nc14 (min)')
ylabel('fraction of active nuclei')
title('[111]')
StandardFigure(fractionON_timetrace_r3,fractionON_timetrace_r3.CurrentAxes)

FigPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\Figures_OpposingGradients\Data\FractionON_instant';
saveas(fractionON_timetrace_r3, [FigPath, filesep, 'fractionON_time_NC14_[111]_',num2str((APbin1-1)*2.5),'%_RuntNulls.tif'])
saveas(fractionON_timetrace_r3, [FigPath, filesep, 'fractionON_time_NC14_[111]_',num2str((APbin1-1)*2.5),'%_RuntNulls.pdf'])

%% 2) Fold-change of fraction on?
% Rationale : I'd like to see the fold-change of the fraction on for the
% [111] construct with/without Runt. Then, I want to see if there's any
% clear trend that we can connect with the Runt protein concentration
% dynamics.

% Use the fields calculated above
% r3_FractionON_mean, r3Null_FractionON_mean_rescaled
% step1 : trim the matrices such that they have the same frame length.

frameLength = min(length(tWindow_r3), length(tWindow_r3Null));

r3_FractionON_mean_trimmed = r3_FractionON_mean(1:frameLength,:);
r3Null_FractionON_mean_trimmed = r3Null_FractionON_mean_rescaled(1:frameLength,:);

% How to put an error bar in here, using some confidence interval?
fc_fractionON_111 = r3_FractionON_mean_trimmed./r3Null_FractionON_mean_trimmed;

pcolor(APAxis,r3_Time(1:frameLength),...
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
saveas(gcf, [FigPath, filesep, 'fold_change_[111]_nc14.tif'])
saveas(gcf, [FigPath, filesep, 'fold_change_[111]_nc14.pdf'])
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

%% exploratory plot (1-fraction ON)
pcolor(APAxis,r3_Time(1:frameLength),...
        1-fc_fractionON_111(1:frameLength,APrange))
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
%% 3) Confirm the validity of spot density metric as a proxy for the fraction on
% by comparing the [000]/[111] Runt WTs for two quantificaitons.

%% 4) Fraction ON in nc13 only
%% Extract useful fields and calculate the mean and STD of fraction on over individual embryos
% We have to calculate the mean and std again, as the null datasets are
% different .mat files that are from single embryo.
% To-Do
% Define and use structure to deal with this repetitive steps of
% calculation simpler!


% [000]
r0_Time = r0Data.ElapsedTime;
r0_nc13 = r0Data.nc13;
r0_nc14 = r0Data.nc14;
%r0_FractionON_mean = r0Data.FractionON;
r0_FractionON_mean = nanmean(r0Data.FractionON_individual,3);
r0_FractionON_STD = nanstd(r0Data.FractionON_individual,0,3);

% [111]
r3_Time = r3Data.ElapsedTime;
r3_nc13 = r3Data.nc13;
r0_nc14 = r0Data.nc14;
r3_FractionON_mean = nanmean(r3Data.FractionON_individual,3);
r3_FractionON_STD = nanstd(r3Data.FractionON_individual,0,3);

% [000], Runt null
r0Null_Time = r0NullData_nc13.ElapsedTime;
r0Null_nc13 = r0NullData_nc13.nc13;
r0Null_nc14 = r0NullData_nc13.nc14;
r0Null_FractionON_mean = nanmean(r0NullData_nc13.FractionON_individual,3);
r0Null_FractionON_STD = nanstd(r0NullData_nc13.FractionON_individual,0,3);
% Since the fraction on is #spots/APbinArea, we want to rescale this using
% the 99% percentile value of FractionON
tWindow = r0Null_nc13:r0Null_nc14;
[row,col] = size(r0Null_FractionON_mean(tWindow,:));
r0Null_scale = prctile(reshape(r0Null_FractionON_mean(tWindow,:),[1,row*col]),99);
r0Null_FractionON_mean_rescaled = r0Null_FractionON_mean./r0Null_scale;

% % We only have up until ~30 minutes for the [000], Runt null for the plotting purpose, let's add some zero values after the time point
% % in here, let's hardcode for the [000] datasets, ~40 min into nc14
% r0Null_tMax = 45; % min
% tRes = median(diff(r0Null_Time));
% nLength = r0Null_tMax/tRes;
% nLength_add = ceil(nLength - length(r0Null_Time));
% [~,nAPbins] = size(r0NullData.FractionON);
% vec_add = nan(nLength_add,nAPbins);
% tvec = tRes*[1:nLength_add] + r0Null_Time(end);
% 
% r0Null_Time_elongated = [r0Null_Time, tvec];
% r0Null_FractionON_mean_rescaled_elongated = [r0Null_FractionON_mean_rescaled; vec_add];



% [111], Runt null
r3Null_Time = r3NullData_nc13.ElapsedTime;
r3Null_nc13 = r3NullData_nc13.nc13;
r3Null_nc14 = r3NullData_nc13.nc14;
r3Null_FractionON_mean = nanmean(r3NullData_nc13.FractionON_individual,3);
r3Null_FractionON_STD = nanstd(r3NullData_nc13.FractionON_individual,0,3);
% Since the fraction on is #spots/APbinArea, we want to rescale this using
% the 99% percentile value of FractionON
tWindow = r3Null_nc13:r3Null_nc14;
[row,col] = size(r3Null_FractionON_mean(tWindow,:));
r3Null_scale = prctile(reshape(r3Null_FractionON_mean(tWindow,:),[1,row*col]),99);
r3Null_FractionON_mean_rescaled = r3Null_FractionON_mean./r3Null_scale;

% 15% APbin is actually NaNs.
r3Null_FractionON_mean_rescaled(:,7) = nan;
%% Plot the different Fraction ONs(heatmap of Time x APbins)


% We will compare the nc14 for now.
APaxis = 0:0.025:1;
APAxis = [0.15:0.025:0.6]*100;
APrange = 7:25;
nAPbins = length(APaxis);

subplot(2,2,1)
tWindow = r0_nc13:r0_nc14;
pcolor(APAxis,r0_Time(tWindow) - r0_Time(r0_nc13),r0_FractionON_mean(tWindow,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc13 (min)')
set(gca,'ytick',0:5:15,'yticklabels',[0:5:15])
title('[000]')
StandardFigure(gcf,gca)


subplot(2,2,2)
tWindow = r0Null_nc13:r0Null_nc14;
pcolor(APAxis,r0Null_Time(tWindow)-r0Null_Time(r0Null_nc13) ,...
            r0Null_FractionON_mean_rescaled(tWindow,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc13 (min)')
set(gca,'ytick',0:5:15,'yticklabels',[0:5:15])
title('[000], Runt null')
StandardFigure(gcf,gca)


subplot(2,2,3)
tWindow = r3_nc13:r3_nc14;
pcolor(APAxis,r3_Time(tWindow) - r3_Time(r3_nc13),...
            r3_FractionON_mean(tWindow,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc13 (min)')
set(gca,'ytick',0:5:15,'yticklabels',[0:5:15])
title('[111]')
StandardFigure(gcf,gca)


subplot(2,2,4)
tWindow = r3Null_nc13:r3Null_nc14;
pcolor(APAxis,r3Null_Time(tWindow) - r3Null_Time(r3Null_nc13) ,...
            r3Null_FractionON_mean_rescaled(tWindow,APrange))
colormap(viridis)
colorbar
xlabel('position (% embryo length)')
set(gca,'xtick',[15,20,30,40,50,60],'xticklabels',[15,20,30,40,50,60])
ylabel('time into nc13 (min)')
set(gca,'ytick',0:5:15,'yticklabels',[0:5:15])
title('[111], Runt null')
StandardFigure(gcf,gca)

% Save the figure
end