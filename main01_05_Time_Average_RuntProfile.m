function main01_05_Time_Average_RuntProfile
%% DESCRIPTION
% This script is for post-processing the Runt data that has already been
% (1) Background subtracted, (2) Synchronized (from NC13)
% to get Time-averaged Runt profile over AP (over different time windows)
% This is for generating the input TF spatial profile across the AP axis.

%% Import the pre-processed data 
% Either from two scripts below.
% main01_02_process_JB3_Runt_datasets
% main01_03_process_JB3_Runt_datasets_BGsubtracted


%% Load the pre-processed data (before averaging over multiple embryos)
FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

RuntFemale = load([FilePath, filesep,'Runt_preprocessed_female.mat'])
   
RuntMale = load([FilePath, filesep,'Runt_preprocessed_male.mat'])

%% Extract useful fields
% Female datasets
MeanFluo_female = RuntFemale.MeanFluo_female;
SDFluo_female = RuntFemale.SDFluo_female;
NC13_female = RuntFemale.NC13_female;
NC14_female = RuntFemale.NC14_female;
Time_female = RuntFemale.Time_female;

% Male datasets
MeanFluo_male = RuntMale.MeanFluo_male;
SDFluo_male = RuntMale.SDFluo_male;
NC13_male = RuntMale.NC13_male;
NC14_male = RuntMale.NC14_male;
Time_male = RuntMale.Time_male;

%% Calculate the time-averaged Runt profile (spatial profile)

% First, let's define the time window (for averaging)
% 0-10 minutes of NC14
% The datasets were taken with frame rate of 1 frame / 1 min, so we can
% easily define the time window.
tWindow = NC14_female:NC14_female +9 ; % beginning 10 minutes of NC14
% In future, this wWindow can be fed into this function as an input.

% In case indexing is different for males and females
% if NC14_female == NC14_male
%     tWindow = NC14_female:NC14_female + 9; % beginning 10 minutes of NC14
% else
%     tWindow(1,:) = NC14_female:NC14_female + 9;
%     tWindow(2,:) = NC14_male:NC14_male + 9;
% end

% Second, let's calculate the temporal average (over AP axis) for each
% embryo.
tAveraged_Fluo_female = squeeze(nanmean(MeanFluo_female(tWindow,:,:),1));
tAveraged_SDFluo_female = squeeze(nanmean(SDFluo_female(tWindow,:,:),1));

tAveraged_Fluo_male = squeeze(nanmean(MeanFluo_male(tWindow,:,:),1));
tAveraged_SDFluo_male = squeeze(nanmean(SDFluo_male(tWindow,:,:),1));
%% Plot to check individual time-averaged Runt spatial profile.
APaxis = 0:0.025:1;

nEmbryos_female = length(MeanFluo_female(1,1,:));
nEmbryos_male = length(MeanFluo_male(1,1,:));

hold on
for embryo=1:nEmbryos_female
    errorbar(APaxis,tAveraged_Fluo_female(:,embryo), tAveraged_SDFluo_female(:,embryo),'r')
    pause
end

for embryo=1:nEmbryos_male
    errorbar(APaxis,tAveraged_Fluo_male(:,embryo), tAveraged_SDFluo_male(:,embryo),'b')
    pause
end
hold off

%% Take the average of "time-averaged spatial Runt profile" over embryos.
AveragedFluo_tAveraged_female = nanmean(tAveraged_Fluo_female,2);
SDFluo_tAveraged_female = nanstd(tAveraged_Fluo_female,0,2);
%nEmbryos_female
SEFluo_tAveraged_female = SDFluo_tAveraged_female./sqrt(nEmbryos_female);

% Filter out weird embryos, by weird embryos, I mean the ones that I'm not
% sure whether they capture the center of the nuclei, or whether it has the
% right time...meaning whether it was cut in the middle in which case it's
% not synchronized with other datasets well.
embryos_male = [1,2,3]; % indexing the embryos from the DataStatus
nEmbryos_male = length(embryos_male);

AveragedFluo_tAveraged_male = nanmean(tAveraged_Fluo_male(:,embryos_male),2);
SDFluo_tAveraged_male = nanstd(tAveraged_Fluo_male(:,embryos_male),0,2);
%nEmbryos_male
SEFluo_tAveraged_male = SDFluo_tAveraged_male./sqrt(nEmbryos_male);

%% Plot to check the averaged (both over embryo ,and over time window) Runt spatial profile.
hold on
errorbar(APaxis, AveragedFluo_tAveraged_male, SEFluo_tAveraged_male)
errorbar(APaxis, AveragedFluo_tAveraged_female, SEFluo_tAveraged_female)

title({'Runt protein spatial profile';'(time-averaged) : 0-10 min into NC14'})
xlabel('AP axis (EL)')
ylabel('Runt concentration (AU)')
legend('male','female','Location','NorthWest')
StandardFigure(gcf,gca)

% Save the plots
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\BGsubtracted_by_Cyto';
saveas(gcf,[FigPath,filesep,'Runt_APprofile_NC14_0-10min_sex_comparison','.tif']); 
saveas(gcf,[FigPath,filesep,'Runt_APprofile_NC14_0-10min_sex_comparison','.pdf']); 

%% Section2. 
% Alright, from the plot above, we can conclude that the Runt profile looks
% more or less same between sex.
% Let's compile male/female data together, to get somewhat averaged
% profile of Runt across AP, then we can use this data to change the AP
% axis as Runt concentration!

% Combine the two time-averaged fluo matrices. (concatenate)
% tAveraged_Fluo_female, tAveraged_Fluo_male
tAveraged_Fluo_mixed = cat(2, tAveraged_Fluo_female, tAveraged_Fluo_male(:,embryos_male));

% Get the average over embryos (all embryos regardless of their sex)
AveragedFluo_tAveraged_mixed = nanmean(tAveraged_Fluo_mixed,2);
SDFluo_tAveraged_mixed = nanstd(tAveraged_Fluo_mixed,0,2);
SEFluo_tAveraged_mixed = SDFluo_tAveraged_mixed./sqrt(length(tAveraged_Fluo_mixed(1,:)));

%% Plot the averaged Runt profile over AP (over time (0-10 min, NC14, and over embryos)
APaxis = 0:0.025:1;

errorbar(APaxis, AveragedFluo_tAveraged_mixed, SEFluo_tAveraged_mixed)
title('Runt concentration over AP')
xlabel('AP axis')
ylabel('Runt concentration (AU)')
legend('Runt-time-averaged')

%% Show that the Runt spatial profile looks the same regardless the averaging time-window.
% Define different time windows
% In this case, we have the same index for NC14 for both male and female
% So, I'll just grab female value. This should be dealt with carefully when
% it's not true.
NC14 = NC14_female;

tWindow1 = 0:5; %0-5 minutes
tWindow2 = 0:10;  %0-10 minutes
tWindow3 = 0:20; % 0-20 minutes

tWindows = {tWindow1 + NC14; tWindow2 + NC14 ; tWindow3 + NC14} ;

% Time-average for female, and male (in case the matrix dimension is not
% the same), then concatenate the matrices.
for i=1:length(tWindows)
    tAveraged_Fluo_female_temp(i,:,:) = squeeze(nanmean(MeanFluo_female(tWindows{i},:,:),1));
    tAveraged_Fluo_male_temp(i,:,:) = squeeze(nanmean(MeanFluo_male(tWindows{i},:,:),1));
end

% Concatenate the matrices (dimension : tWindows x APbins x embryos)
tAveraged_Fluo_mixed = cat(3, tAveraged_Fluo_female_temp, tAveraged_Fluo_male_temp);

% Get the average over embryos (all embryos regardless of their sex)
AveragedFluo_tAveraged_mixed = nanmean(tAveraged_Fluo_mixed,3);
SDFluo_tAveraged_mixed = nanstd(tAveraged_Fluo_mixed,0,3);
SEFluo_tAveraged_mixed = SDFluo_tAveraged_mixed./sqrt(length(tAveraged_Fluo_mixed(1,1,:)));

%% Save the useful fields into mat files.
% Define the file path
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

% Sort out the fields to be saved
savedVariables_processed = {};
savedVariables_processed = [savedVariables_processed,...
                                    'AveragedFluo_tAveraged_mixed',...
                                    'SDFluo_tAveraged_mixed',...
                                    'SEFluo_tAveraged_mixed',...
                                    'tWindows', 'NC14',...
                                    'tAveraged_Fluo_mixed'];
% Save the variables into .mat file.
save([FilePath,filesep,'Runt_TimeAveraged_mixedSex_NC14.mat'],...
    savedVariables_processed{:},'-v7.3');

%% Check if the saving worked fine.
A = load([FilePath, filesep, 'Runt_TimeAveraged_mixedSex_NC14.mat'])

%% Plot the spatial gradient for different time windows of temporal averaging
hold on
for i=1:length(tWindows)
    errorbar(APaxis, AveragedFluo_tAveraged_mixed(i,:), SEFluo_tAveraged_mixed(i,:))
end
xlim([0.1 0.7])
xticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7])

title('Runt gradient (time-averaged)')
xlabel('AP axis (EL)')
ylabel('Runt concentration (AU)')
legend('0-5 min','0-10 min','0-20 min','Location','NorthWest')
StandardFigure(gcf,gca)

% Save the figure
% Save the figures
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\BGsubtracted_by_Cyto';
saveas(gcf,[FigPath,filesep,'Runt_time_averaged_NC14_multiwindows_mixedSex','.tif']); 
saveas(gcf,[FigPath,filesep,'Runt_time_averaged_NC14_multiwindows_mixedSex','.pdf']); 
end