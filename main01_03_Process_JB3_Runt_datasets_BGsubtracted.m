function main01_02_Process_JB3_Runt_datasets
%% DESCRIPTION
% This script is for processing the Runt datasets.
% 1) Synchronize the datasets with the beginning of NC13
% 2) Background subtraction
% 3) Average for each sex
% 4) Save these for future usage.
% 

%% Step1. Sync the datasets using AverageDatasets_LlamaTaggedProtein.m
% Female
AverageDatasets_LlamaTaggedProtein('Runt-1min-200Hz-Female','savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData')
% Male
AverageDatasets_LlamaTaggedProtein ('Runt-1min-200Hz-Male','savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData')
% NoNB
% AverageDatasets_NuclearProtein('NoNB','savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData')
%% Load the datasets that are synchronized above.
% File path for calling the datasets
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

Runt_female = load([FilePath, filesep, 'Runt-1min-200Hz-Female_BGsubtracted-Averaged.mat']);

Runt_male = load([FilePath, filesep, 'Runt-1min-200Hz-Male_BGsubtracted-Averaged.mat']);

%% Extract useful fields
% Time, MeanVectorAP_indiviual, SDVector_individual,
% NParticlesAP_individual, etc.

% Female
MeanFluo_female = Runt_female.MeanVectorAP_BGsubtracted_individual;
SDFluo_female = Runt_female.SDVectorAP_BGsubtracted_individual;
NNuclei_female = Runt_female.NParticlesAP_BGsubtracted_individual;
SEFluo_female = SDFluo_female./sqrt(NNuclei_female); % for SEM of individual embryos.
Time_female = Runt_female.ElapsedTime;
NC13_female = Runt_female.nc13;
NC14_female = Runt_female.nc14;

% Male
MeanFluo_male = Runt_male.MeanVectorAP_BGsubtracted_individual;
SDFluo_male = Runt_male.SDVectorAP_BGsubtracted_individual;
NNuclei_male = Runt_male.NParticlesAP_BGsubtracted_individual;
SEFluo_male = SDFluo_male./sqrt(NNuclei_male); % for SEM of individual embryos.
Time_male = Runt_male.ElapsedTime;
NC13_male = Runt_male.nc13;
NC14_male = Runt_male.nc14;

%% Plot for checking
% As of 11/15/2019, I don't have any anterior data point.
% So, I'll process more embryos for this, but for now, I'll average the
% data points across AP bins as an approximation of anterior points.

% Female
% APaxis = 0:0.025:1;
% EmbryoIndex = 4;
% 
% hold on
% for i=1:length(MeanFluo_female(:,1,EmbryoIndex))
%     errorbar(0:0.025:1, MeanFluo_female(i,:,EmbryoIndex), SEFluo_female(i,:,EmbryoIndex))
% end

% Male
% APaxis = 0:0.025:1;
% EmbryoIndex = 3;
% 
% hold on
% for i=1:length(MeanFluo_male(:,1,EmbryoIndex))
%     errorbar(0:0.025:1, MeanFluo_male(i,:,EmbryoIndex), SEFluo_male(i,:,EmbryoIndex))
%     pause
% end

%% Plot for checking (over Time)
tLength = min([length(Time_female), length(Time_male)]);
APbin = 15;
hold on
% female
for embryo=1:length(MeanFluo_female(1,1,:))
    errorbar(Time_female(1:tLength), MeanFluo_female(1:tLength,APbin,embryo), SDFluo_female(1:tLength,APbin,embryo),'r')
    pause
end

% male
for embryo=1:length(MeanFluo_male(1,1,:))
    errorbar(Time_male(1:tLength), MeanFluo_male(1:tLength,APbin,embryo), SDFluo_male(1:tLength,APbin,embryo),'b')
    pause
end

%% Some post-processing for Male dataset : Prefix = '2018-12-02-RuntJB3-vasa-eGFP-Pos5'
% since it has one frame that the reflection dominates the nuclear fluo.
% 24th frame would be an average between 23rd, and 25th, same goes for the
% SD.
% MAKE SURE TO CHECK THIS INDEXING AT THE DATASTATUS.XLSX!!!
%  i=2; % 2019-07-09, YJK, 2nd dataset for now. I need a better way for
% % dealing with this.
% MeanFluo_male(24,:,i) = (MeanFluo_male(23,:,i) + MeanFluo_male(25,:,i))/2;
% SDFluo_male(24,:,i) = (SDFluo_male(23,:,i) + SDFluo_male(25,:,i))/2;

%% Averaging over multiple embryos (for each sex)
% Female
numEmbryos_female = length(MeanFluo_female(1,1,:));
Averaged_Fluo_female = nanmean(MeanFluo_female, 3);
SD_Fluo_female = nanstd(MeanFluo_female, 0 , 3);
SE_Fluo_female = SD_Fluo_female./sqrt(numEmbryos_female); % # of embryos
Time_female = Time_female(1:tLength); % New time

% Male
numEmbryos_male = length(MeanFluo_male(1,1,:));
Averaged_Fluo_male = nanmean(MeanFluo_male, 3);
SD_Fluo_male = nanstd(MeanFluo_male, 0 , 3);
SE_Fluo_male = SD_Fluo_male./sqrt(numEmbryos_male); % # of embryos
Time_male = Time_male(1:tLength); % New time
%% Plot to check (Runt over Time - male vs female)
APbin = 21;
APpos = (APbin - 1) * 2.5;

hold on
errorbar(Time_male,...
        Averaged_Fluo_male(1:tLength,APbin), SE_Fluo_male(1:tLength,APbin))
errorbar(Time_female,...
        Averaged_Fluo_female(1:tLength,APbin), SE_Fluo_female(1:tLength,APbin))

% xlim([0 max(Time_female)])
xlim([0 50])
ylim([100 max(Averaged_Fluo_female(:,APbin)) + 50])
title('Runt protein concentration over Time')
xlabel('Time (min)')
ylabel('Runt protein conc. (AU)')
legend('male','female','Location','NorthWest')
StandardFigure(gcf,gca)

% Save the plots
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\BGsubtracted_by_Cyto';
saveas(gcf,[FigPath,filesep,'Runt_Time_NC13-14_MaleFemale @ AP =',num2str(APpos),'%','.tif']); 
saveas(gcf,[FigPath,filesep,'Runt_Time_NC13-14_MaleFemale @ AP =',num2str(APpos),'%','.pdf']); 

%% Runt conc. over Time (male & female) - NC14 only
% Pick one AP bin
APbin = 21;
APpos = (APbin - 1) * 2.5;

hold on
errorbar(Time_male(NC14_male:end) - Time_male(NC14_male),...
        Averaged_Fluo_male(NC14_male:tLength,APbin), SE_Fluo_male(NC14_male:tLength,APbin))
errorbar(Time_female(NC14_female:end)- Time_female(NC14_female),...
        Averaged_Fluo_female(NC14_female:tLength,APbin), SE_Fluo_female(NC14_female:tLength,APbin))

xlim([0 25])
%ylim([0 300])
title('Runt protein concentration over Time')
xlabel('Time (min)')
ylabel('Runt protein conc. (AU)')
legend('male','female','Location','NorthWest')
StandardFigure(gcf,gca)

% Save the figures
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\BGsubtracted_by_Cyto';
saveas(gcf,[FigPath,filesep,'Runt_Time_NC14_MaleFemale @ AP =',num2str(APpos),'%','.tif']); 
saveas(gcf,[FigPath,filesep,'Runt_Time_NC14_MaleFemale @ AP =',num2str(APpos),'%','.pdf']); 

%% Save the processed fields (after the background subtraction & averaging)
BGsubtracted = 'CytoFluo';
Averaged = false;

% Save female data
savedVariables_female_preprocessed = {};
savedVariables_female_preprocessed = [savedVariables_female_preprocessed,...
                                    'MeanFluo_female',...
                                    'SDFluo_female',...
                                    'SEFluo_female',...
                                    'BGsubtracted','Averaged',...
                                    'Time_female','NC13_female','NC14_female'];
                     
% Save male data
savedVariables_male_preprocessed = {};
savedVariables_male_preprocessed = [savedVariables_male_preprocessed,...
                                    'MeanFluo_male',...
                                    'SDFluo_male',...
                                    'SEFluo_male',...
                                    'BGsubtracted','Averaged',...
                                    'Time_male','NC13_male','NC14_male'];
                                
% Save the variables into .mat file.
save([FilePath,filesep,'Runt_preprocessed_female.mat'],...
    savedVariables_female_preprocessed{:},'-v7.3');

save([FilePath,filesep,'Runt_preprocessed_male.mat'],...
    savedVariables_male_preprocessed{:},'-v7.3');

%% Save the processed fields
BGsubtracted = 'CytoFluo';
Averaged = true;

% Save female data
savedVariables_female = {};
savedVariables_female = [savedVariables_female, 'Averaged_Fluo_female', ...
                        'SD_Fluo_female','SE_Fluo_female',...
                        'BGsubtracted', 'numEmbryos_female','Time_female'];
% Save male data
savedVariables_male = {};
savedVariables_male = [savedVariables_male, 'Averaged_Fluo_male',...
                        'SD_Fluo_male','SE_Fluo_male',...
                        'BGsubtracted', 'numEmbryos_male',  'Time_male'];


save([FilePath,filesep,'Runt_processed_female.mat'],...
    savedVariables_female{:},'-v7.3');

save([FilePath,filesep,'Runt_processed_male.mat'],...
    savedVariables_male{:},'-v7.3');

%% For time-averaging, refer to main01_05_Time_Average_RuntProfile.m script.
%% Section 3. (Optional)
% Plot the Runt profile in different ways.
% %% (1) Runt conc. - individual vs. mean (before the background-subtraction)
% % Female
% % Cut the time frames, such that the time points that I'm not sure about
% % the Z-centering would not be included.
% tIndices = 1:31; % 0 to 30 minutes
% 
% % Pick one AP bin.
% APbin = 17;
% 
% RuntFigure_individual_mean = figure;
% hold on
% for i=1:numEmbryos_female
%     errorbar(Time_female(tIndices), MeanFluo_female_BGsubtracted(tIndices,APbin,i), SEFluo_female_BGsubtracted(tIndices,APbin,i))
% end
%     errorbar(Time_female(tIndices), Averaged_Fluo_female(tIndices,APbin), SE_Fluo_female(tIndices,APbin))
%     
% % Find the y Maximum, for plotting    
% yMax1 = max(max(MeanFluo_female_BGsubtracted(tIndices,APbin,:)));
% yMax2 = max(Averaged_Fluo_female(tIndices,APbin));
% yMax = max([yMax1 yMax2]);
% xlim([0 30])
% ylim([0 yMax + 20])
% 
% % Format for figures    
% title([{'Runt nuclear concentration over time';'(Background-subtracted)'}])
% xlabel('Time (min)')
% ylabel('Nuclear concentration (AU)')
% legend('embryo1', 'embryo2', 'embryo3', 'Average','Location','northwest')
% 
% StandardFigure(RuntFigure_individual_mean, RuntFigure_individual_mean.CurrentAxes)
% 
% % Save figure
% % FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\Runt_background_subtracted_07082019';
% % saveas(gcf,[FigPath,filesep,'Runt_BGsubtracted_over_Time','_individual_vs_Mean_SEM',' AP = ',num2str((APbin-1)*2.5),'%','.tif'])
% % saveas(gcf,[FigPath,filesep,'Runt_BGsubtracted_over_Time','_individual_vs_Mean_SEM',' AP = ',num2str((APbin-1)*2.5),'%','.pdf'])
% 
% %% (1) Runt conc. - individual vs. mean (before the background-subtraction)
% % Male
% 
% % Cut the time frames, such that the time points that I'm not sure about
% % the Z-centering would not be included.
% tIndices = 1:31; % 0 to 30 minutes
% 
% APbin = 17;
% RuntFigure_individual_mean = figure;
% hold on
% for i=1:numEmbryos_male
%     errorbar(Time_male(tIndices), MeanFluo_male_BGsubtracted(tIndices,APbin,i), SEFluo_male_BGsubtracted(tIndices,APbin,i))
% end
%     errorbar(Time_male(tIndices), Averaged_Fluo_male(tIndices,APbin), SE_Fluo_male(tIndices,APbin))
%     
% % Find the y Maximum, for plotting    
% yMax1 = max(max(MeanFluo_male_BGsubtracted(tIndices,APbin,:)));
% yMax2 = max(Averaged_Fluo_male(tIndices,APbin));
% yMax = max([yMax1 yMax2]);
% xlim([0 30])
% ylim([0 yMax + 50])
% 
% % Format for figures    
% title([{'Runt nuclear concentration over time';'(Background-subtracted)'}])
% xlabel('Time (min)')
% ylabel('Nuclear concentration (AU)')
% legend('embryo1', 'embryo2', 'embryo3', 'Average','Location','northwest')
% 
% StandardFigure(RuntFigure_individual_mean, RuntFigure_individual_mean.CurrentAxes)
% 
% % Save figure
% % FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\Runt_background_subtracted_07082019';
% % saveas(gcf,[FigPath,filesep,'Runt_BGsubtracted_over_Time','_individual_vs_Mean_SEM',' AP = ',num2str((APbin-1)*2.5),'%','.tif'])
% % saveas(gcf,[FigPath,filesep,'Runt_BGsubtracted_over_Time','_individual_vs_Mean_SEM',' AP = ',num2str((APbin-1)*2.5),'%','.pdf'])
% 
% 
% %% (2) Runt conc. over time (BG-subtracted, then Aveaged)
% 
% % Cut the time frames, such that the time points that I'm not sure about
% % the Z-centering would not be included.
% tIndices = 1:31; % 0 to 30 minutes
% 
% APbin = 14;
% 
% hold on
% errorbar(Time_male(tIndices), Averaged_Fluo_male(tIndices,APbin), SE_Fluo_male(tIndices,APbin))
% errorbar(Time_female(tIndices), Averaged_Fluo_female(tIndices,APbin), SE_Fluo_female(tIndices,APbin))
% 
% title([{'Runt nuclear concentration over time';'(Background-subtracted)'}])
% xlabel('Time (min)')
% ylabel('Nuclear concentration (AU)')
% legend('male','female')
% StandardFigure(gcf,gca)
% 
% % Save figure
% % FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\Runt_background_subtracted_07082019';
% % saveas(gcf,[FigPath,filesep,'Runt_BGsubtracted_over_Time','_Averaged_Male_Female',' AP = ',num2str((APbin-1)*2.5),'%','.tif'])
% % saveas(gcf,[FigPath,filesep,'Runt_BGsubtracted_over_Time','_Averaged_Male_Female',' AP = ',num2str((APbin-1)*2.5),'%','.pdf'])
% 
% %% (3) Runt conc. over AP axis (@ different time points)
% APaxis = 0:0.025:1;
% 
% for tPoint = 1:30%length(Time_male)
%     clf
% %     RuntFigure_APaxis = figure;
%     hold on
%     errorbar(APaxis, Averaged_Fluo_male(tPoint,:), SE_Fluo_male(tPoint,:))
%     errorbar(APaxis, Averaged_Fluo_female(tPoint,:), SE_Fluo_female(tPoint,:))
% 
%     xlim([0.2 0.7])
%     ylim([0 350])
%     title([{'Runt nuclear concentration over AP';'(Background-subtracted, Averaged)'}])
%     xlabel('AP axis (Embryo Length)')
%     ylabel('Nuclear concentration (AU)')
%     legend('male','female','Location','northwest')
%     StandardFigure(gcf,gca)
%     pause
%     % Save figure
% %     FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\Runt_background_subtracted_07082019';
% %     saveas(gcf,[FigPath,filesep,'Runt_BGsubtracted_AP','_Averaged_Male_Female',' Time = ',num2str(tPoint),'min','.tif'])
% %     saveas(gcf,[FigPath,filesep,'Runt_BGsubtracted_AP','_Averaged_Male_Female',' Time = ',num2str(tPoint),'min','.pdf'])
% 
% end
% 
% % Saving for all time points could make a movie using ImageJ.


% %% (4) Runt conc. (time-averaged for different time windows)
% % NC13
% % tWindow1
% tWindow1 = 4:6; % 3-5 minutes, when the MS2 spot shows the initial rise.
% % tWindow2
% tWindow2 = 1:11; % 0-10 minutes, when the Runt keeps increasing.
% % tWindow3
% tWindow3 = 1:16; % 0-15 minutes, whole NC13.
% 
% tWindow{1} = tWindow1;
% tWindow{2} = tWindow2;
% tWindow{3} = tWindow3;
% 
% for i=1:length(tWindow)
%     % Male
%     time_averaged_Fluo_male(i,:) = nanmean(Averaged_Fluo_male(tWindow{i},:));
%     time_averaged_SE_Fluo_male(i,:) = nanmean(SE_Fluo_male(tWindow{i},:));
%     % Female
%     time_averaged_Fluo_female(i,:) = nanmean(Averaged_Fluo_female(tWindow{i},:));
%     time_averaged_SE_Fluo_female(i,:) = nanmean(SE_Fluo_female(tWindow{i},:));
% end
% 
% 
% %% plot to check
% hold on
% 
% for i=1:3%length(tWindow)
%     errorbar(APaxis, time_averaged_Fluo_female(i,:), time_averaged_SE_Fluo_female(i,:))
%     errorbar(APaxis, time_averaged_Fluo_male(i,:), time_averaged_SE_Fluo_male(i,:))
% end
% 
% % Formatting
% xlim([0.2 0.7])
% %ylim([0 250])
% title([{'Runt nuclear concentration over AP';'(Time-averaged)'}])
% xlabel('AP axis (Embryo Length)')
% ylabel('Nuclear concentration (AU)')
% legend('3-5 min','0-10 min','0-15 min','Location','northwest')
% StandardFigure(gcf,gca)
% 
% % Save figure
% % FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\Runt_background_subtracted_07082019';
% % saveas(gcf,[FigPath,filesep,'Runt_time_averaged_AP','_female_','NC13','.tif'])
% % saveas(gcf,[FigPath,filesep,'Runt_time_averaged_AP','_female_','NC13','.pdf'])
% 
% %% Plot for individual time points vs time-averaged
% hold on
% errorbar(APaxis, Averaged_Fluo_female(4,:), SE_Fluo_female(4,:))
% errorbar(APaxis, Averaged_Fluo_female(6,:), SE_Fluo_female(6,:))
% errorbar(APaxis, Averaged_Fluo_female(11,:), SE_Fluo_female(11,:))
% 
% errorbar(APaxis, time_averaged_Fluo_female(2,:), time_averaged_SE_Fluo_female(2,:))
% 
% xlim([0.2 0.7])
% %ylim([0 250])
% title([{'Runt nuclear concentration over AP';'(Time-averaged)'}])
% xlabel('AP axis (Embryo Length)')
% ylabel('Nuclear concentration (AU)')
% legend('3 min','5 min','10 min','0-10 min','Location','northwest')
% StandardFigure(gcf,gca)
% 
% % Save figure
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\Runt_background_subtracted_07082019';
% saveas(gcf,[FigPath,filesep,'Runt_AP_@timepoints','_female_','NC13','.tif'])
% saveas(gcf,[FigPath,filesep,'Runt_AP_@timepoints','_female_','NC13','.pdf'])
% %% time-averaged Runt profile over AP
% % NC14
% % tWindow1
% tWindow4 = 19:24; % 1-6 minutes of NC14.
% % tWindow2
% tWindow5 = 19:30; % 1-10 minutes, when the Runt keeps increasing.
% % tWindow3
% tWindow6 = 19:40; % 1-15 minutes, whole NC13.
% 
% tWindow{4} = tWindow4;
% tWindow{5} = tWindow5;
% tWindow{6} = tWindow6;
% 
% for i=4:length(tWindow)
%     % Male
%     time_averaged_Fluo_male(i,:) = nanmean(Averaged_Fluo_male(tWindow{i},:));
%     time_averaged_SE_Fluo_male(i,:) = nanmean(SE_Fluo_male(tWindow{i},:));
%     % Female
%     time_averaged_Fluo_female(i,:) = nanmean(Averaged_Fluo_female(tWindow{i},:));
%     time_averaged_SE_Fluo_female(i,:) = nanmean(SE_Fluo_female(tWindow{i},:));
% end
% %% plot to check
% hold on
% for i=4:length(tWindow)
%     errorbar(APaxis, time_averaged_Fluo_female(i,:), time_averaged_SE_Fluo_female(i,:))
% end
% 
% % Formatting
% xlim([0.2 0.7])
% %ylim([0 250])
% title([{'Runt nuclear concentration over AP';'(Time-averaged)'}])
% xlabel('AP axis (Embryo Length)')
% ylabel('Nuclear concentration (AU)')
% legend('0-5 min','0-10 min','0-20 min','Location','northwest')
% StandardFigure(gcf,gca)
% 
% % Save figure
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\Runt_background_subtracted_07082019';
% saveas(gcf,[FigPath,filesep,'Runt_time_averaged_AP','_female_','NC14','.tif'])
% saveas(gcf,[FigPath,filesep,'Runt_time_averaged_AP','_female_','NC14','.pdf'])
% 
% %% Save the time-averaged vectors with tWindow
% Runt_female_time_averaged.tWindow = tWindow;
% Runt_female_time_averaged.Averaged_Fluo = time_averaged_Fluo_female;
% Runt_female_time_averaged.SE_Fluo = time_averaged_SE_Fluo_female;
% 
% Runt_male_time_averaged.tWindow = tWindow;
% Runt_male_time_averaged.Averaged_Fluo = time_averaged_Fluo_male;
% Runt_male_time_averaged.SE_Fluo = time_averaged_SE_Fluo_male;
% 
% 
% % Save male data
% save([FilePath,filesep,'Runt_time_averaged_female.mat'],...
%    '-struct', 'Runt_female_time_averaged','-v7.3');
% save([FilePath,filesep,'Runt_time_averaged_male.mat'],...
%    '-struct', 'Runt_male_time_averaged','-v7.3');
end