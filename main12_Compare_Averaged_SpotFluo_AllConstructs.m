function main12_compare_r0123(varargin)
%% DESCRIPTION
% This script is to compare 
% 1) Mean spot fluo
% 2) integrated mRNA profiles (over AP, at different NC ranges, etc.) 
% between males and females of the same synthetic enhancers.

% The plan is to use this script to compare r0, r1, r2, and r3 

%% Average datasets using AverageDatasets.m
% 7/17/2019
% Old construct, r0, with 1 bp mutation in the eve promoter. I'm using this
% since I have more anterior data point, as well as NC14.
% AverageDatasets('r0','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');

% New datasets, mixed sex
%AverageDatasets('r0-new-mixed','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');
% New datasets, with all females
AverageDatasets('r0-new','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');
AverageDatasets('r1-new','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');
AverageDatasets('r2-new','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');
AverageDatasets('r3-new','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');
%% Load datasets
% This assumes that the DataType is the name of the constructs in the
% DataStatus.xlsx tab, for example, r3-new

% FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
 FilePath = '/Users/yangjoonkim/Dropbox/OpposingGradient/OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'
% This is because our new datasets(r0-new-male, r0-new-female lack NC14)

%r0Data = load([FilePath, filesep, 'r0-new-female']);
r0Data = load([FilePath, filesep, 'r0-new.mat']);
r1Data = load([FilePath, filesep, 'r1-new.mat']);
r2Data = load([FilePath, filesep, 'r2-new']);
r3Data = load([FilePath, filesep, 'r3-new']);

r1closeData = load([FilePath, filesep, 'r1-close.mat']);
r1midData = load([FilePath, filesep, 'r1-mid.mat']);

r2closeData = load([FilePath, filesep, 'r2_1+2.mat']);
r2farData = load([FilePath, filesep, 'r2_1+3.mat']);

%% Put them in a structure for convenience
AveragedData{1} = r0Data;
AveragedData{2} = r1Data;
AveragedData{3} = r2Data;
AveragedData{4} = r3Data;
AveragedData{5} = r1closeData;
AveragedData{6} = r1midData;
AveragedData{7} = r2closeData;
AveragedData{8} = r2farData;
%% 
%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;
colorDict.lightgreen = [205,214,209]/255;
colorDict.thickpink = [132,27,69]/255;
% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.magenta; colorDict.thickpink]; 
%% Plot the time traces of Averaged spot fluo (for NC14)
% for each datatype (since it's averaged already by AverageDatasets.m)
% I'll take NC14 and MeanVectorAP, SEVectorAP then plot them over time
% (ElpasedTime)

for APbin = 9:25 % (APbin-1)*2.5 % of Embryo length
    clf
    hold on
    for i=1:4%length(AveragedData)
        % initialize the variables
        % Extract useful info
        nc13 = AveragedData{i}.nc13;
        nc14 = AveragedData{i}.nc14;
        Time = AveragedData{i}.ElapsedTime;
        tWindow = nc13:nc14;
        %tWindow = nc14:length(ElapsedTime); % NC14
        
        SpotFluo = AveragedData{i}.MeanVectorAP;
        SEMSpotFluo = AveragedData{i}.SEVectorAP;

        errorbar(Time(tWindow)-Time(nc13),...
                 SpotFluo(tWindow,APbin), SEMSpotFluo(tWindow,APbin),'Color',ColorChoice(i,:))
    end
    
    xlim([0 18])
    xlabel('time into NC13(min)')
    ylabel('mean spot fluorescence (AU)')
    legend('0','1','2','3')
    StandardFigure(gcf,gca)
    %pause
    FigPath = '/Users/yangjoonkim/Dropbox/Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces'
    saveas(gcf,[FigPath,filesep,'averaged_MS2_TimeTraces_r0123' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.tif']); 
    saveas(gcf,[FigPath,filesep,'averaged_MS2_TimeTraces_r0123' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.pdf']); 
end
%% individual embryo spot fluo- time trace check
% for r0, since it drops too quickly in NC14, I want to double check if
% there's some spot/particle segmentation issue (systematically)
Time_r0 = r0Data.ElapsedTime;
NC13_r0 = r0Data.nc13;
NC14_r0 = r0Data.nc14;
MeanVectorAP_r0 = r0Data.MeanVectorAP_individual;
SDVectorAP_r0 = r0Data.SDVectorAP_individual;
NParticlesAP_r0 = r0Data.NParticlesAP_individual;

%% plot individual embryos of r0
hold on
APbin = 15;
for embryo=1:5
    errorbar(Time_r0(NC14_r0:end) - Time_r0(NC14),...
                MeanVectorAP_r0(NC14_r0:end, APbin, embryo),...
                SDVectorAP_r0(NC14_r0:end, APbin, embryo)./...
                sqrt(NParticlesAP_r0(NC14_r0:end, APbin, embryo)))
end

%% Plot the time traces of Averaged spot fluo (over ALL nuclei)
% for each datatype (since it's averaged already by AverageDatasets.m)
% I'll take NC13/14 and MeanVectorAP, SEVectorAP then plot them over time
% (ElpasedTime)

for APbin = 9:25 % (APbin-1)*2.5 % of Embryo length
    clf
    hold on
    for i=1:4%length(AveragedData)
        % initialize the variables
        % Extract useful info
        nc13 = AveragedData{i}.nc13;
        nc14 = AveragedData{i}.nc14;
        Time = AveragedData{i}.ElapsedTime;
        %tWindow = nc13:nc14;
        %tInitial = nc13;
        tWindow = nc14:length(Time); % NC14
        tInitial =nc14;
        
        SpotFluo_ind = AveragedData{i}.MeanVectorAP_individual;
        instFractionON_ind = AveragedData{i}.FractionON_individual;
        SpotFluo_overALLnuclei = SpotFluo_ind.*instFractionON_ind;
        SpotFluo_mean_ALLnuclei = nanmean(SpotFluo_overALLnuclei,3); % over embryos
        [~,~,numEmbryos] = size(SpotFluo_ind);
        SEMSpotFluo = nanstd(SpotFluo_overALLnuclei,0,3)./sqrt(numEmbryos); % over embryos
        %SEMSpotFluo = AveragedData{i}.SEVectorAP;

        errorbar(Time(tWindow)-Time(tInitial),...
                 SpotFluo_mean_ALLnuclei(tWindow,APbin), SEMSpotFluo(tWindow,APbin),'Color',ColorChoice(i,:))
    end
    
%     xlim([0 15])%nc13
    xlim([0 40])
%     ylim([0 1200])
    xlabel('time into NC14(min)')
    ylabel('mean spot fluorescence (AU)')
    legend('0','1','2','3')
    StandardFigure(gcf,gca)

    FigPath = '/Users/yangjoonkim/Dropbox/Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_overAllNuclei_TimeTraces'
    saveas(gcf,[FigPath,filesep,'averaged_MS2_TimeTraces_r0123' ,'_AP=',num2str((APbin-1)*2.5),'%_NC14' , '.tif']); 
    saveas(gcf,[FigPath,filesep,'averaged_MS2_TimeTraces_r0123' ,'_AP=',num2str((APbin-1)*2.5),'%_NC14' , '.pdf']); 
end



%% OLD
%% Assign fields to variables

% r0
Time_r0 = r0Data.ElapsedTime;
NC13_r0 = r0Data.nc13;
NC14_r0 = r0Data.nc14;
MeanVectorAP_r0 = r0Data.MeanVectorAP;
SDVectorAP_r0 = r0Data.SDVectorAP;
SEVectorAP_r0 = r0Data.SEVectorAP;
NParticlesAP_r0 = r0Data.NParticlesAP;

% r1
Time_r1 = r1Data.ElapsedTime;
NC13_r1 = r1Data.nc13;
NC14_r1 = r1Data.nc14;
MeanVectorAP_r1 = r1Data.MeanVectorAP;
SDVectorAP_r1 = r1Data.SDVectorAP;
SEVectorAP_r1 = r1Data.SEVectorAP;
NParticlesAP_r1 = r1Data.NParticlesAP;

% r2
Time_r2 = r2Data.ElapsedTime;
NC13_r2 = r2Data.nc13;
NC14_r2 = r2Data.nc14;
MeanVectorAP_r2 = r2Data.MeanVectorAP;
SDVectorAP_r2 = r2Data.SDVectorAP;
SEVectorAP_r2 = r2Data.SEVectorAP;
NParticlesAP_r2 = r2Data.NParticlesAP;

% r3
Time_r3 = r3Data.ElapsedTime;
NC13_r3 = r3Data.nc13;
NC14_r3 = r3Data.nc14;
MeanVectorAP_r3 = r3Data.MeanVectorAP;
SDVectorAP_r3 = r3Data.SDVectorAP;
SEVectorAP_r3 = r3Data.SEVectorAP;
NParticlesAP_r3 = r3Data.NParticlesAP;

%% Optional (Do synchronization more carefully, then average the mean spot fluo)
% %% r0
% Time_r0 = r0Data.ElapsedTime;
% NC13_r0 = r0Data.nc13;
% NC14_r0 = r0Data.nc14;
% MeanVectorAP_individual_r0 = r0Data.MeanVectorAP_individual;
% SDVectorAP_individual_r0 = r0Data.SDVectorAP_individual;
% NParticlesAP_individual_r0 = r0Data.NParticlesAP_individual;
% numEmbryos_r0 = length(r0Data.MeanVectorAP_individual(1,1,:));
% 
% APbin = 12;
% for APbin = 10:21
%     clf
%     hold on
%     for i=1:numEmbryos_r0
%         errorbar(Time_r0, MeanVectorAP_individual_r0(:,APbin,i), SDVectorAP_individual_r0(:,APbin,i))
%     end
%     pause
% end
% 
% %% r1
% Time_r1 = r1Data.ElapsedTime;
% NC13_r1 = r1Data.nc13;
% NC14_r1 = r1Data.nc14;
% MeanVectorAP_individual_r1 = r1Data.MeanVectorAP_individual;
% SDVectorAP_individual_r1 = r1Data.SDVectorAP_individual;
% NParticlesAP_individual_r1 = r1Data.NParticlesAP_individual;
% numEmbryos_r1 = length(r1Data.MeanVectorAP_individual(1,1,:));
% 
% APbin = 12;
% for APbin = 10:21
%     clf
%     hold on
%     for i=1:numEmbryos_r1
%         errorbar(Time_r1, MeanVectorAP_individual_r1(:,APbin,i), SDVectorAP_individual_r1(:,APbin,i))
%     end
%     pause
% end
% %% r2
% Time_r2 = r2Data.ElapsedTime;
% NC13_r2 = r2Data.nc13;
% NC14_r2 = r2Data.nc14;
% MeanVectorAP_individual_r2 = r2Data.MeanVectorAP_individual;
% SDVectorAP_individual_r2 = r2Data.SDVectorAP_individual;
% NParticlesAP_individual_r2 = r2Data.NParticlesAP_individual;
% numEmbryos_r2 = length(r2Data.MeanVectorAP_individual(1,1,:));
% 
% APbin = 12;
% for APbin = 10:21
%     clf
%     hold on
%     for i=1:numEmbryos_r2
%         errorbar(Time_r2, MeanVectorAP_individual_r2(:,APbin,i), SDVectorAP_individual_r2(:,APbin,i))
%     end
%     pause
% end
% 
% %% r3
% Time_r3 = r3Data.ElapsedTime;
% NC13_r3 = r3Data.nc13;
% NC14_r3 = r3Data.nc14;
% MeanVectorAP_individual_r3 = r3Data.MeanVectorAP_individual;
% SDVectorAP_individual_r3 = r3Data.SDVectorAP_individual;
% NParticlesAP_individual_r3 = r3Data.NParticlesAP_individual;
% numEmbryos_r3 = length(r3Data.MeanVectorAP_individual(1,1,:));
% 
% APbin = 12;
% for APbin = 10:21
%     clf
%     hold on
%     for i=1:numEmbryos_r3
%         errorbar(Time_r3, MeanVectorAP_individual_r3(:,APbin,i), SDVectorAP_individual_r3(:,APbin,i))
%     end
%     pause
% end
% 
% %% Re-sync the MeanVectorAP if needed (fine-tuning)
% 
% %% (Optional) Averaging over multiple embryos - average of mean spot fluo, with SEM = SD/sqrt(# of embryos)
% % Caveat : How should I calculate the SEM of mean spot fluo?
% % Should I use the "number of embryos", or "the total number of spots"?
% 
% % r0
% MeanVectorAP_r0 = nanmean(MeanVectorAP_individual_r0,3);
% SDVectorAP_r0 = nanstd(MeanVectorAP_individual_r0,0,3);
% SEVEctorAP_r0 = SDVectorAP_r0./sqrt(numEmbryos_r0);
% 
% % r1
% MeanVectorAP_r1 = nanmean(MeanVectorAP_individual_r1,3);
% SDVectorAP_r1 = nanstd(MeanVectorAP_individual_r1,0,3);
% SEVEctorAP_r1 = SDVectorAP_r1./sqrt(numEmbryos_r1);
% 
% % r2
% MeanVectorAP_r2 = nanmean(MeanVectorAP_individual_r2,3);
% SDVectorAP_r2 = nanstd(MeanVectorAP_individual_r2,0,3);
% SEVEctorAP_r2 = SDVectorAP_r2./sqrt(numEmbryos_r2);
% 
% % r3
% MeanVectorAP_r3 = nanmean(MeanVectorAP_individual_r3,3);
% SDVectorAP_r3 = nanstd(MeanVectorAP_individual_r3,0,3);
% SEVEctorAP_r3 = SDVectorAP_r3./sqrt(numEmbryos_r3);

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow; colorDict.red; colorDict.brown]; % 4 embryos max. it could be extended easily
%% Plot the MeanVectorAP
APbin = 9; % APbin position

% % NC13
% Range_r0 = NC13_r0:NC14_r0;
% Range_r1 = NC13_r1:NC14_r1;
% Range_r2 = NC13_r2:NC14_r2;
% Range_r3 = NC13_r3:NC14_r3;
% tDelay = zeros(1,4);

% NC14
Range_r0 = NC14_r0:length(Time_r0);
Range_r1 = NC14_r1:length(Time_r1);
Range_r2 = NC14_r2:length(Time_r2);
Range_r3 = NC14_r3:length(Time_r3);
tDelay = [Time_r0(NC14_r0), Time_r1(NC14_r1), Time_r2(NC14_r2), Time_r3(NC14_r3)];

MS2TraceFig = figure;
hold on
errorbar(Time_r0(Range_r0)-tDelay(1), MeanVectorAP_r0(Range_r0,APbin),...
            SEVectorAP_r0(Range_r0,APbin),'Color',ColorChoice(1,:))
        
errorbar(Time_r1(Range_r1)-tDelay(2), MeanVectorAP_r1(Range_r1,APbin),...
            SEVectorAP_r1(Range_r1,APbin),'Color',ColorChoice(2,:))
        
errorbar(Time_r2(Range_r2)-tDelay(3), MeanVectorAP_r2(Range_r2,APbin),...
            SEVectorAP_r2(Range_r2,APbin),'Color',ColorChoice(3,:))
        
errorbar(Time_r3(Range_r3)-tDelay(4), MeanVectorAP_r3(Range_r3,APbin),...
            SEVectorAP_r3(Range_r3,APbin),'Color',ColorChoice(4,:))
% xlim, ylim
%xlim([0 20])
ylim([0 max(MeanVectorAP_r0(Range_r0,APbin)) + 100])
title(['Mean spot fluorescence over time @ AP = ',num2str((APbin-1)*2.5),'%'])
xlabel('Time into NC14 (min)')
ylabel('Mean spot fluorescence (AU)')
legend('r0','r1','r2','r3')
StandardFigure(MS2TraceFig,MS2TraceFig.CurrentAxes)
% standardizeFigure_YJK(gca,legend,[])

% Save Figure
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\hbP2-r0123-Averaged_MS2_traces_multipleEmbryos';
saveas(MS2TraceFig,[FigPath,filesep, 'Mean_MS2_traces at AP=', num2str((APbin-1)*2.5),'%' , '_NC14' , '.tif']); 
saveas(MS2TraceFig,[FigPath,filesep, 'Mean_MS2_traces at AP=', num2str((APbin-1)*2.5),'%' , '_NC14' , '.pdf']); 

%% 
%% Section 1. Mean Spot fluo
%% Mean spot fluorescence between constructs

%% Section 2. Accumulated mRNA
%  %%  Plot the Accumulated mRNA
% % APaxis = 0:0.025:1;
% % 
% % 
% % AccumulatedmRNA_All_r0(isnan(AccumulatedmRNA_All_SD_r0)) = nan;
% % AccumulatedmRNA_All_r1(isnan(AccumulatedmRNA_All_SD_r1)) = nan;
% % AccumulatedmRNA_All_r2(isnan(AccumulatedmRNA_All_SD_r2)) = nan;
% % AccumulatedmRNA_All_r3(isnan(AccumulatedmRNA_All_SD_r3)) = nan;
% % 
% % % NC13
% % AccumulatedmRNA_NC13_figure = figure
% % hold on
% % errorbar(APaxis, AccumulatedmRNA_All_r0(NC14_r0,:), AccumulatedmRNA_All_SD_r0(NC14_r0,:))
% % errorbar(APaxis, AccumulatedmRNA_All_r1(NC14_r1,:), AccumulatedmRNA_All_SD_r1(NC14_r1,:))
% % errorbar(APaxis, AccumulatedmRNA_All_r2(NC14_r2,:), AccumulatedmRNA_All_SD_r2(NC14_r2,:))
% % errorbar(APaxis, AccumulatedmRNA_All_r3(NC14_r3,:), AccumulatedmRNA_All_SD_r3(NC14_r3,:))
% % 
% % title('Accumulated mRNA over AP @ NC13')
% % xlabel('AP axis (EL)')
% % ylabel('Accumulated mRNA (AU)')
% % legend('r0','r1','r2','r3')
% % 
% % AccumulatedmRNA_NC14_figure = figure
% % hold on
% % errorbar(APaxis, AccumulatedmRNA_All_r0(end,:), AccumulatedmRNA_All_SD_r0(end,:))
% % errorbar(APaxis, AccumulatedmRNA_All_r1(end,:), AccumulatedmRNA_All_SD_r1(end,:))
% % errorbar(APaxis, AccumulatedmRNA_All_r2(end,:), AccumulatedmRNA_All_SD_r2(end,:))
% % errorbar(APaxis, AccumulatedmRNA_All_r3(end,:), AccumulatedmRNA_All_SD_r3(end,:))
% % 
% % title('Accumulated mRNA over AP @ NC14')
% % xlabel('AP axis (EL)')
% % ylabel('Accumulated mRNA (AU)')
% % legend('r0','r1','r2','r3')
% 
 %% Alternative calculation for the Accumulated mRNA
% % Since my script, AverageDatasets didn't take into account of the
% % APbinArea for calculating the total mRNA
% % I'll try IntegratemRNA.m script for this.
% 
% % 1) Load the datasets
% r0Data = LoadMS2Sets('r0')
% r1Data = LoadMS2Sets('r1-new-female')
% r2Data = LoadMS2Sets('r2-new-female')
% r3Data = LoadMS2Sets('r3-new-female')
% 
% % 2) IntegratemRNA
% [TotalProd_r0,TotalProdError_r0,TotalProdN_r0,...
%     MeanTotalProd_r0,SDTotalProd_r0,SETotalProd_r0]=IntegratemRNA(r0Data,1,2)
% 
% [TotalProd_r1,TotalProdError_r1,TotalProdN_r1,...
%     MeanTotalProd_r1,SDTotalProd_r1,SETotalProd_r1]=IntegratemRNA(r1Data,1,2)
% 
% [TotalProd_r2,TotalProdError_r2,TotalProdN_r2,...
%     MeanTotalProd_r2,SDTotalProd_r2,SETotalProd_r2]=IntegratemRNA(r2Data,1,2)
% 
% [TotalProd_r3,TotalProdError_r3,TotalProdN_r3,...
%     MeanTotalProd_r3,SDTotalProd_r3,SETotalProd_r3]=IntegratemRNA(r3Data,1,2)
% 
% %% Plot the Accumulated mRNA
% APaxis = 0:0.025:1;
% % Figure Path
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Transcription-Output\hbP2-r0123-AccumulatedmRNA-females';
% % 
% % AccumulatedmRNA_All_r0(isnan(AccumulatedmRNA_All_SD_r0)) = nan;
% % AccumulatedmRNA_All_r1(isnan(AccumulatedmRNA_All_SD_r1)) = nan;
% % AccumulatedmRNA_All_r2(isnan(AccumulatedmRNA_All_SD_r2)) = nan;
% % AccumulatedmRNA_All_r3(isnan(AccumulatedmRNA_All_SD_r3)) = nan;
% 
% % NC13
% NC= 13;
% AccumulatedmRNA_NC13_figure = figure
% hold on
% errorbar(APaxis, MeanTotalProd_r0(:,NC), SETotalProd_r0(:,NC))
% errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
% errorbar(APaxis, MeanTotalProd_r2(:,NC), SETotalProd_r2(:,NC))
% errorbar(APaxis, MeanTotalProd_r3(:,NC), SETotalProd_r3(:,NC))
% xlim([0.2 0.5])
% 
% title('Accumulated mRNA over AP @ NC13')
% xlabel('AP axis (EL)')
% ylabel('Accumulated mRNA (AU)')
% legend('r0','r1','r2','r3')
% 
% StandardFigure(AccumulatedmRNA_NC13_figure,AccumulatedmRNA_NC13_figure.CurrentAxes)
% %saveas(AccumulatedmRNA_NC13_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC13','_SE','.tif'])
% %saveas(AccumulatedmRNA_NC13_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC13','_SE','.pdf'])
% 
% % NC14
% NC = 14;
% AccumulatedmRNA_NC14_figure = figure
% hold on
% errorbar(APaxis, MeanTotalProd_r0(:,NC), SETotalProd_r0(:,NC))
% errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
% errorbar(APaxis, MeanTotalProd_r2(:,NC), SETotalProd_r2(:,NC))
% errorbar(APaxis, MeanTotalProd_r3(:,NC), SETotalProd_r3(:,NC))
% 
% xlim([0.2 0.5])
% 
% title('Accumulated mRNA over AP @ NC14')
% xlabel('AP axis (EL)')
% ylabel('Accumulated mRNA (AU)')
% legend('r0','r1','r2','r3')
% 
% StandardFigure(AccumulatedmRNA_NC14_figure,AccumulatedmRNA_NC14_figure.CurrentAxes)
% %saveas(AccumulatedmRNA_NC14_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC14','_SE','.tif'])
% %saveas(AccumulatedmRNA_NC14_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC14','_SE','.pdf'])
% 
 %% Side thing (males)
% % Compare the total mRNA profile for diffrent sex
% % 1) Load the datasets
% 
% r1Data_male = LoadMS2Sets('r1-new-male')
% r2Data_male = LoadMS2Sets('r2-new-male')
% r3Data_male = LoadMS2Sets('r3-new-male')
% 
% % 2) IntegratemRNA
% 
% [TotalProd_r1_male,TotalProdError_r1_male,TotalProdN_r1_male,...
%     MeanTotalProd_r1_male,SDTotalProd_r1_male,SETotalProd_r1_male]=IntegratemRNA(r1Data_male,1,2)
% 
% [TotalProd_r2_male,TotalProdError_r2_male,TotalProdN_r2_male,...
%     MeanTotalProd_r2_male,SDTotalProd_r2_male,SETotalProd_r2_male]=IntegratemRNA(r2Data_male,1,2)
% 
% [TotalProd_r3_male,TotalProdError_r3_male,TotalProdN_r3_male,...
%     MeanTotalProd_r3_male,SDTotalProd_r3_male,SETotalProd_r3_male]=IntegratemRNA(r3Data_male,1,2)
% 
% %% Plot for different sexes
% NC= 13;
% AccumulatedmRNA_NC13_figure = figure
% hold on
% % females (r0 is mixed)
% errorbar(APaxis, MeanTotalProd_r0(:,NC), SETotalProd_r0(:,NC))
% errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
% errorbar(APaxis, MeanTotalProd_r2(:,NC), SETotalProd_r2(:,NC))
% errorbar(APaxis, MeanTotalProd_r3(:,NC), SETotalProd_r3(:,NC))
% % males
% errorbar(APaxis, MeanTotalProd_r1_male(:,NC), SETotalProd_r1_male(:,NC))
% errorbar(APaxis, MeanTotalProd_r2_male(:,NC), SETotalProd_r2_male(:,NC))
% errorbar(APaxis, MeanTotalProd_r3_male(:,NC), SETotalProd_r3_male(:,NC))
% xlim([0.2 0.5])
% 
% title('Accumulated mRNA over AP @ NC13')
% xlabel('AP axis (EL)')
% ylabel('Accumulated mRNA (AU)')
% legend('r0','r1','r2','r3','r1-male','r2-male','r3-male')

%% Section 2. Accumulated mRNA (using AccumulatedmRNA.m script)
%% Accumulate mRNA for individual embryos, then average with proper SEM
% (SD/sqrt(number of embryos))
% MinParticles = 2;
AccumulatedmRNA('r0',2);
AccumulatedmRNA('r1-new-female',2);
AccumulatedmRNA('r2-new-female',2);
AccumulatedmRNA('r3-new-female',2);

%% load the Accumulated mRNA
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

%% Calculate the 
%% Section 3. Fraction ON
% Use the 
%% 
end