function main01_04_Process_eGFP_Bcd_datasets_time_average
%% DESCRIPTION
% This script is for processing the Bcd datasets.
% 1) Synchronize the datasets with the beginning of NC13
% 2) Background subtraction
% 3) Average across multiple embryos
% 3-1) time-average for certain time windows.
% 4) Save these for future usage.

%% update : Check the consistency with Liz & Jonathan's eGFP-Bcd datasets
% Last updated : 10/4/2020
DataPath = 'S:\YangJoon\Dropbox\OpposingGradient\eGFP-Bcd-From-Liz-Jonathan';
BcdAnt = load([DataPath, filesep, 'BcdGFPAnt.mat']);
BcdAnt = BcdAnt.DataBcd;

BcdMid = load([DataPath, filesep, 'BcdGFPMid.mat']);
BcdMid = BcdMid.DataBcd;

BcdPos = load([DataPath, filesep, 'BcdGFPPos.mat']);
BcdPos = BcdPos.DataBcd;

%% Synchronize the vectors starting from nc14

% Extract the Mean, SD, and Number of nuclei info, as well as
% nc14-end of measurement

% First, find the length of the NC14 for each dataset
nc14_ant = length(BcdAnt.ElapsedTime) - BcdAnt.nc14 + 1;
nc14_mid = length(BcdMid.ElapsedTime) - BcdMid.nc14 + 1;
nc14_pos = length(BcdPos.ElapsedTime) - BcdPos.nc14 + 1;

nc14_max = max([nc14_ant, nc14_mid, nc14_pos]);

% Second, extract the useful fields from NC14 to the length of nc14_min,
% for compilation. We will combine the Ant/Mid/Pos into one
nAPbins = 41; % this is the default, we can change this if we want finer bins.

ElapsedTime = (1:nc14_max)*0.5; % (min)
MeanVectorAP = nan(nc14_max,41);
SDVectorAP = nan(nc14_max,41);
NParticlesAP  = nan(nc14_max,41);

for i=1:nAPbins
    if  i<15
        MeanVectorAP(1:nc14_ant,i) = BcdAnt.MeanVectorAP(BcdAnt.nc14:end,i);
        SDVectorAP(1:nc14_ant,i)	= BcdAnt.SDVectorAP(BcdAnt.nc14:end,i);
        NParticlesAP(1:nc14_ant,i) = BcdAnt.NParticlesAP(BcdAnt.nc14:end,i);
    elseif i>=15 && i<=24
        MeanVectorAP(1:nc14_mid,i) = BcdMid.MeanVectorAP(BcdMid.nc14:end,i);
        SDVectorAP(1:nc14_mid,i)   = BcdMid.SDVectorAP(BcdMid.nc14:end,i);
        NParticlesAP(1:nc14_mid,i) = BcdMid.NParticlesAP(BcdMid.nc14:end,i);
    else 
        MeanVectorAP(1:nc14_pos,i) = BcdPos.MeanVectorAP(BcdPos.nc14:end,i);
        SDVectorAP(1:nc14_pos,i)   = BcdPos.SDVectorAP(BcdPos.nc14:end,i);
        NParticlesAP(1:nc14_pos,i) = BcdPos.NParticlesAP(BcdPos.nc14:end,i);
    end
        
end

%% Save the compiled structure into Bcd_NC14_compiled.mat
%% Save the structure, compiledData for future usage (in plotting scripts)
DataPath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';
save([DataPath, filesep, 'Bcd_NC14_compiled.mat'],...
        'ElapsedTime', 'MeanVectorAP', 'SDVectorAP', 'NParticlesAP')
    
%% Time-average for the beginning 5-10 minutes
% Note that our time-resolution is 30 sec, so 21st time point is 10 min into
% nc14
tRes = median(diff(ElapsedTime));
tWindow = 5/tRes+1:10/tRes+1;

Bcd_timeAveraged_5_10min_nc14 = nanmean(MeanVectorAP(tWindow,:));
Bcd_timeAveraged_5_10min_nc14_SEM = nanstd(MeanVectorAP(tWindow,:),0,1)./sqrt(length(tWindow));
Bcd_timeAveraged_5_10min_nc14_SD = sqrt(nanmean(SDVectorAP(tWindow,:).^2));
    
errorbar(APaxis, Bcd_timeAveraged_5_10min_nc14, Bcd_timeAveraged_5_10min_nc14_SD)




% Bcd_timeAveraged_10min_nc14 = nanmean(MeanVectorAP(tWindow,:));
% Bcd_timeAveraged_10min_nc14_SEM = nanstd(MeanVectorAP(tWindow,:),0,1)./sqrt(length(tWindow));
% Bcd_timeAveraged_10min_nc14_SD = sqrt(nanmean(SDVectorAP(tWindow,:).^2));
%     
% errorbar(APaxis, Bcd_timeAveraged_10min_nc14, Bcd_timeAveraged_10min_nc14_SD)

%% Save the time-averaged values to a .mat file so that we can load anytime for a modeling input.
% DataPath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';
% save([DataPath, filesep, 'Bcd_NC14_TimeAveraged.mat'],...
%         'Bcd_timeAveraged_10min_nc14', 'Bcd_timeAveraged_10min_nc14_SD')
DataPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\mat files';
save([DataPath, filesep, 'Bcd_NC14_TimeAveraged_5_10min.mat'],...
        'Bcd_timeAveraged_5_10min_nc14', 'Bcd_timeAveraged_5_10min_nc14_SD')
    
%% generate plots of Bicoid dynamics

%% nc13-nc14 (end)
% Pick the Bicoid Anterior dataset, then just plot one AP bin
APpos = 30; % % of embryo length
APbin = APpos/2.5 + 1;

tWindow = BcdAnt.nc13:length(BcdAnt.ElapsedTime);

errorbar(BcdAnt.ElapsedTime(tWindow)- BcdAnt.ElapsedTime(BcdAnt.nc13),...
            BcdAnt.MeanVectorAP(tWindow,APbin),...
            BcdAnt.SDVectorAP(tWindow,APbin))

% Mark the beginning of nc14
%xline(BcdAnt.ElapsedTime(BcdAnt.nc14) - BcdAnt.ElapsedTime(BcdAnt.nc13), '--')
xline(20, '--')
% format
xlabel('time into nc13 (min)')
ylabel('Bicoid concentration (AU)')

ylim([0 max(BcdAnt.MeanVectorAP(tWindow,APbin))+20])
box on

StandardFigure(gcf,gca)
% save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\InputTF';
saveas(gcf,[FigPath,filesep,'Bcd_',num2str(APpos),'%_nc13-14.tif'])
saveas(gcf,[FigPath,filesep,'Bcd_',num2str(APpos),'%_nc13-14.pdf'])

%% nc14 only
% Pick the Bicoid Anterior dataset, then just plot one AP bin
APpos = 30; % % of embryo length
APbin = APpos/2.5 + 1;

tWindow = BcdAnt.nc14:length(BcdAnt.ElapsedTime);

errorbar(BcdAnt.ElapsedTime(tWindow)- BcdAnt.ElapsedTime(BcdAnt.nc14)-2,...
            BcdAnt.MeanVectorAP(tWindow,APbin),...
            BcdAnt.SDVectorAP(tWindow,APbin))

% format
xlabel('time into nc 14 (min)')
ylabel('Bicoid concentration (AU)')

ylim([0 max(BcdAnt.MeanVectorAP(tWindow,APbin))+20])
xlim([0 40])
box on

StandardFigure(gcf,gca)
% save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\InputTF';
saveas(gcf,[FigPath,filesep,'Bcd_',num2str(APpos),'%_nc14.tif'])
saveas(gcf,[FigPath,filesep,'Bcd_',num2str(APpos),'%_nc14.pdf'])
%% %%%%%%%%%%%%%%%%%%%%%%%OLD script with my own datasets %%%%%%%%%%%%%%%%%%%%%%%%%5
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

%% generate plots of time-traces

%% Time-average for different time-windows
% First, let's define different time windows for averaging
% 0-5min, 0-10min, 0-20 min into NC14 as the first round.

tWindows{1} = 0:5; % min
tWindows{2} = 0:10; % min
tWindows{3} = 0:20; % min


%% Part 2. Time-average first, then average over multiple embryos.
% Load the "BcdData" structure (only in nc14) - Liz & Jonathan
% this structure consists of 6 embryos spanning from anterior-posterior
% For each embryo, they have different frame lengths
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\mat files\BcdLevelsNC14AllDatasets.mat')

%% check the synchronization
APbin = 13;

for i=1:length(BcdData)
    hold on
    errorbar(BcdData(i).TimeNC14, BcdData(i).BcdNC14(:,APbin), BcdData(i).SDBcdNC14(:,APbin))
end

%% time-average for the beginning 5-10 minutes into nc14
% frameRate = 0.5; % [min]
% tWindow = 5/frameRate + 1:10/frameRate+1;

% Note that the frame rate actually differs for different embryos, thus we
% need to get different time windows vectors for each embryo.

% calculate the time-average Bicoid concentration for each embryo
for i=1:length(BcdData)
    clear Time
    Time = BcdData(i).TimeNC14;
    tWindow = find(Time>=5 & Time<=10);
    time_averaged_fluo(i,:) = nanmean(BcdData(i).BcdNC14(tWindow,:));
end

%% calculate the average and SE over multiple embryos
time_averaged_fluo_mean = nanmean(time_averaged_fluo);
time_averaged_fluo_sem = nanstd(time_averaged_fluo)/sqrt(length(BcdData));

%% compile the data into a structure
BcdData_processed.timewindow = '5-10min_nc14'
BcdData_processed.time_averaged_fluo_mean = time_averaged_fluo_mean;
BcdData_processed.time_averaged_fluo_mean_sem = time_averaged_fluo_sem;
BcdData_processed.time_averaged_fluo = time_averaged_fluo;


%% save the results
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\InputTF';
save([FilePath, filesep, 'BcdData_processed_nc14.mat'],'BcdData_processed')

%% Part3. Test this idea with different time window of time-averaging
% Just to check whether it's a scaling issue or not.

% 1) 0-10 minutes into nc14
% calculate the time-average Bicoid concentration for each embryo
for i=1:length(BcdData)
    clear Time
    Time = BcdData(i).TimeNC14;
    tWindow = find(Time>=0 & Time<=10);
    time_averaged_fluo_10min(i,:) = nanmean(BcdData(i).BcdNC14(tWindow,:));
end

%% calculate the average and SE over multiple embryos
time_averaged_fluo_10min_mean = nanmean(time_averaged_fluo_10min);
time_averaged_fluo_10min_sem = nanstd(time_averaged_fluo_10min)/sqrt(length(BcdData));
end