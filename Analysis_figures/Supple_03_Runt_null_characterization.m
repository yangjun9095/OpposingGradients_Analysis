function Supple_03_Runt_null_characterization
% DESCRIPTION 
% This script is for comparing the expression of JB3-Runt;r0 with that of
% run3 (Runt null);r3 (with MCP-GFP(4F).

% Last updated : 1/6/2020

%% Using MCP-GFP (4F) (instead of MCP-mCherry)
% I'm going to use MCP-GFP (4F) for imaging rN-MS2 reporter, and also for
% detection of JB3-Runt (we only keep the data where there's no Runt
% pattern (stripe in the mid-late NC14).

% First, let's load datasets with WT Runt, as a baseline for a comparison.
% r0
Data_r0 = LoadMS2Sets('r0-new','dontCompare');

% r0 (old construct) - the only difference is a 1 bp mutation in the eve
% promoter
%Data_r0_old = LoadMS2Sets('r0','dontCompare');

% r3
Data_r3 = LoadMS2Sets('r3-new','dontCompare');

% Runt nulls
% r0, Runt null
Data_r0_RuntNull = LoadMS2Sets('r0_RuntNull');
% r3, Runt null
Data_r3_RuntNull = LoadMS2Sets('r3_RuntNull');

%% Extract the initial slopes (fitted using FitMeanAPAsymmetric.m)
% r0
[fittedRate_r0,fittedRateSD_r0,fittedTon_r0] = Extract_Fields_MeanFits(Data_r0,'Asymmetric');
% r3
[fittedRate_r3,fittedRateSD_r3,fittedTon_r3] = Extract_Fields_MeanFits(Data_r3,'Asymmetric');

% Nulls
% r0, Runt null
[fittedRate_r0_RuntNull,fittedRateSD_r0_RuntNull,fittedTon_r0_RuntNull] = Extract_Fields_MeanFits(Data_r0_RuntNull,'Asymmetric');

% r3, Runt null
[fittedRate_r3_RuntNull,fittedRateSD_r3_RuntNull,fittedTon_r3_RuntNull] = Extract_Fields_MeanFits(Data_r3_RuntNull,'Asymmetric');

% Side thing
%[fittedRate_r0_old,fittedRateSD_r0_old,fittedTon_r0_old] = Extract_Fields_MeanFits(Data_r0_old,'Asymmetric');

%% Average over embryos (per construct)
average_fittedRate_r0 = nanmean(fittedRate_r0,3);
SEM_fittedRate_r0 = nanstd(fittedRate_r0,0,3)/sqrt(length(Data_r0));

average_fittedRate_r3 = nanmean(fittedRate_r3,3);
SEM_fittedRate_r3 = nanstd(fittedRate_r3,0,3)/sqrt(length(Data_r3));

% Nulls
average_fittedRate_r0_RuntNull = nanmean(fittedRate_r0_RuntNull,3);
SEM_fittedRate_r0_RuntNull = nanstd(fittedRate_r0_RuntNull,0,3)/sqrt(length(Data_r0_RuntNull));

average_fittedRate_r3_RuntNull = nanmean(fittedRate_r3_RuntNull,3);
SEM_fittedRate_r3_RuntNull = nanstd(fittedRate_r3_RuntNull,0,3)/sqrt(length(Data_r3_RuntNull));

% side
% average_fittedRate_r0_old = nanmean(fittedRate_r0_old,3);
% SEM_fittedRate_r0_old = nanstd(fittedRate_r0_old,0,3)/sqrt(length(Data_r0_old));
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
colorDict.darkBlue = [16,16,55]/255;
colorDict.darkOrange = [207,134,41]/255;

ColorChoice = [colorDict.blue; colorDict.red; colorDict.darkBlue; colorDict.darkOrange];

%% Compare the r0 expression from old vs new constructs (just 1 bp mutation in the eve promoter)
% NC = 14; % nc14
% NC = NC-11; % For indexing in the FitResults
% 
% hold on
% errorbar(0:0.025:1, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC))
% errorbar(0:0.025:1, average_fittedRate_r0_old(:,NC), SEM_fittedRate_r0_old(:,NC))
%% plot averaged initial slopes (fitted by Asymmetric fit)
NC = 14; % nc14
NC = NC-11; % For indexing in the FitResults

hold on
errorbar(0:0.025:1, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC),'Color',ColorChoice(1,:)) % Blue
errorbar(0:0.025:1, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC),'Color',ColorChoice(2,:)) % Red

errorbar(0:0.025:1, average_fittedRate_r0_RuntNull(:,NC), SEM_fittedRate_r0_RuntNull(:,NC),'Color',ColorChoice(3,:)) % darkBlue
errorbar(0:0.025:1, average_fittedRate_r3_RuntNull(:,NC), SEM_fittedRate_r3_RuntNull(:,NC),'Color',ColorChoice(4,:)) % dark Orange

legend('r0','r3','r0, Runt null','r3, Runt null')

xlim([0.15 0.6])

xlabel('AP Position')
ylabel('Initial rate (AU/min)')
title('Initial rate over AP axis (nc 14)')
StandardFigure(gcf, gca)

% Save the figure and data
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\';

saveas(gcf,[FigPath 'InitialRate_AsymmetricFit_r0_r3_withRuntNulls_MCP-GFP4F' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath 'InitialRate_AsymmetricFit_r0_r3_withRuntNulls_MCP-GFP4F' , '_NC14' , '.pdf']); 

%% Save the fields for future usage.
% Define a variable with everything.
savedVariables = {};
savedVariables = [savedVariables,...
                  'fittedRate_r0','fittedRate_r0_RuntNull',...
                  'fittedRate_r3','fittedRate_r3_RuntNull',...
                  'fittedRateSD_r0','fittedRateSD_r0_RuntNull',...
                  'fittedRateSD_r3','fittedRateSD_r3_RuntNull',...
                  'average_fittedRate_r0','average_fittedRate_r0_RuntNull',...
                  'average_fittedRate_r3','average_fittedRate_r3_RuntNull',...
                  'SEM_fittedRate_r0','SEM_fittedRate_r0_RuntNull',...
                  'SEM_fittedRate_r3','SEM_fittedRate_r3_RuntNull'];   

filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
save([filePath,filesep,'InitialSlopes_ONnuclei_RuntNull_r0_r3.mat'],...
    savedVariables{:},'-v7.3');
%% DEPRECATED (MCP-mCherry)

% Here, I need to check the dosage (in case of MCP-mCherry), 
% which I can quickly check with the offset values.

% %% Load the datasets
% r0_RuntJB3 = LoadMS2Sets('r0-RuntJB3-MCP-mCherry','dontCompare')
% r3_RuntNull = LoadMS2Sets('r3-Run3-MCP-mCherry','dontCompare')
% 
% %r0_RuntJB3 = r0_RuntJB3.Particles;
% % r3_RuntNull = r3_RuntNull.Particles;
% %% Extract the fitted values from the datasets
% % Caveat : 
% % 1) Asymmetric fit : MeanFitAPAsymmetric.m
% [fittedRate_r0,fittedRateSD_r0,fittedTon_r0] = Extract_Fields_MeanFits(r0_RuntJB3,'Asymmetric');
% [fittedRate_r3_RuntNull,fittedRateSD_r3_RuntNull,fittedTon_r3_RuntNull] = Extract_Fields_MeanFits(r3_RuntNull,'Asymmetric');
% 
% 
% %% Plot individuals for a sanity check
% NC=2 % NC14
% hold on
% for i=1:3
%     errorbar(0:0.025:1, fittedRate_r0(:,NC,i),fittedRateSD_r0(:,NC,i))
%     pause
% end
% 
% figure(2)
% hold on
% for i=1:3
%     errorbar(0:0.025:1, fittedRate_r3_RuntNull(:,NC,i),fittedRateSD_r3_RuntNull(:,NC,i))
%     pause
% end
% 
% %% Calculate the average using nanmean, nanstd
% % r0, Runt-JB3
% average_fittedRate_r0 = nanmean(fittedRate_r0,3);
% SEM_fittedRate_r0 = nanstd(fittedRate_r0,0,3)/sqrt(length(r0_RuntJB3));
% 
% % r3, Runt null (run3)
% average_fittedRate_r3_RuntNull = nanmean(fittedRate_r3_RuntNull,3);
% SEM_fittedRate_r3_RuntNull = nanstd(fittedRate_r3_RuntNull,0,3)/sqrt(length(r3_RuntNull));
% 
% %% Color definition
% % This is defining the line color
% colorDict = struct();
% colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
% colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
% colorDict.yellow = [234,194,100]/255;
% colorDict.cyan = [108,188,233]/255;
% colorDict.magenta = [208,109,171]/255;
% colorDict.lightBlue = [115,142,193]/255;
% colorDict.purple = [171,133,172]/255;
% colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
% colorDict.brown = [179,155,142]/255;
% colorDict.darkgreen = [126,157,144]/255;
% 
% ColorChoice = [colorDict.purple; colorDict.green; colorDict.yellow; colorDict.red; colorDict.brown];
% %% Plot the averaged fittedRate (Asymmetric) (initial rate of RNAP loading), and SEM
% 
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\';
% 
% % NC13
% InitialRate_NC13_figure = figure;
% nc = 2; % NC13
% hold on
% errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc),'Color',ColorChoice(1,:))
% 
% errorbar(0:0.025:1,average_fittedRate_r3_RuntNull(:,nc),SEM_fittedRate_r3_RuntNull(:,nc),'Color',ColorChoice(2,:))
% 
% xlim([0.15 0.5])
% %ylim([0 400])
% 
% legend('r0','r3, Runt null')
% xlabel('AP Position')
% ylabel('Initial rate (AU/min)')
% title('Initial rate of RNAP loading along AP axis, at NC 13')
% StandardFigure(InitialRate_NC13_figure, InitialRate_NC13_figure.CurrentAxes)
% 
% saveas(InitialRate_NC13_figure,[FigPath 'InitialRate_AsymmetricFit_r0_vs_r3_run3' , '_NC13' , '.tif']); 
% saveas(InitialRate_NC13_figure,[FigPath 'InitialRate_AsymmetricFit_r0_vs_r3_run3' , '_NC13' , '.pdf']); 
% 
% % NC14
% InitialRate_NC14_figure = figure;
% nc = 3; % NC14
% hold on
% errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc),'Color',ColorChoice(1,:))
% 
% errorbar(0:0.025:1,average_fittedRate_r3_RuntNull(:,nc),SEM_fittedRate_r3_RuntNull(:,nc),'Color',ColorChoice(2,:))
% 
% xlim([0.15 0.5])
% %ylim([0 400])
% 
% legend('r0','r3, Runt null')
% xlabel('AP Position')
% ylabel('Initial rate (AU/min)')
% title('Initial rate of RNAP loading along AP axis, at NC 14')
% StandardFigure(InitialRate_NC14_figure, InitialRate_NC14_figure.CurrentAxes)
% 
% saveas(InitialRate_NC14_figure,[FigPath 'InitialRate_AsymmetricFit_r0_vs_r3_run3' , '_NC14' , '.tif']); 
% saveas(InitialRate_NC14_figure,[FigPath 'InitialRate_AsymmetricFit_r0_vs_r3_run3' , '_NC14' , '.pdf']); 
% 
% %% Extract the fitted values from the datasets
% % Caveat : 
% 
% % 2) Linear fit : MeanFitAPLinearSlope.m
% [fittedRate_r0,fittedRateSD_r0,fittedTon_r0] = Extract_Fields_MeanFits(r0_RuntJB3,'Linear');
% [fittedRate_r3_RuntNull,fittedRateSD_r3_RuntNull,fittedTon_r3_RuntNull] = Extract_Fields_MeanFits(r3_RuntNull,'Linear');
% 
% %% Calculate the average using nanmean, nanstd
% % r0, Runt-JB3
% average_fittedRate_r0 = nanmean(fittedRate_r0,3);
% SEM_fittedRate_r0 = nanstd(fittedRate_r0,0,3)/sqrt(length(r0_RuntJB3));
% 
% % r3, Runt null (run3)
% average_fittedRate_r3_RuntNull = nanmean(fittedRate_r3_RuntNull,3);
% SEM_fittedRate_r3_RuntNull = nanstd(fittedRate_r3_RuntNull,0,3)/sqrt(length(r3_RuntNull));
% %% Plot the averaged fittedRate (Linear) (initial rate of RNAP loading), and SEM
% 
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_linear\';
% 
% % NC13
% InitialRate_NC13_figure = figure;
% nc = 2; % NC13
% hold on
% errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc),'Color',ColorChoice(1,:))
% 
% errorbar(0:0.025:1,average_fittedRate_r3_RuntNull(:,nc),SEM_fittedRate_r3_RuntNull(:,nc),'Color',ColorChoice(2,:))
% 
% xlim([0.15 0.5])
% %ylim([0 400])
% 
% legend('r0','r3, Runt null')
% xlabel('AP Position')
% ylabel('Initial rate (AU/min)')
% title('Initial rate of RNAP loading along AP axis, at NC 13')
% StandardFigure(InitialRate_NC13_figure, InitialRate_NC13_figure.CurrentAxes)
% 
% saveas(InitialRate_NC13_figure,[FigPath 'InitialRate_LinearFit_r0_vs_r3_run3' , '_NC13' , '.tif']); 
% saveas(InitialRate_NC13_figure,[FigPath 'InitialRate_LinearFit_r0_vs_r3_run3' , '_NC13' , '.pdf']); 
% 
% % NC14
% InitialRate_NC14_figure = figure;
% nc = 3; % NC14
% hold on
% errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc),'Color',ColorChoice(1,:))
% 
% errorbar(0:0.025:1,average_fittedRate_r3_RuntNull(:,nc),SEM_fittedRate_r3_RuntNull(:,nc),'Color',ColorChoice(2,:))
% 
% xlim([0.15 0.5])
% %ylim([0 400])
% 
% legend('r0','r3, Runt null')
% xlabel('AP Position')
% ylabel('Initial rate (AU/min)')
% title('Initial rate of RNAP loading along AP axis, at NC 14')
% StandardFigure(InitialRate_NC14_figure, InitialRate_NC14_figure.CurrentAxes)
% 
% saveas(InitialRate_NC14_figure,[FigPath 'InitialRate_LinearFit_r0_vs_r3_run3' , '_NC14' , '.tif']); 
% saveas(InitialRate_NC14_figure,[FigPath 'InitialRate_LinearFit_r0_vs_r3_run3' , '_NC14' , '.pdf']); 
% 
% %% Save the fitted initial rate and SEM 
% 
% % Define the values
% 
% save('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedInitialRate_AsymmetricFit.mat',...
%         'average_fittedRate_r0','SEM_fittedRate_r0',...
%         'average_fittedRate_r1_female','SEM_fittedRate_r1_female',...
%         'average_fittedRate_r2_female','SEM_fittedRate_r2_female',...
%         'average_fittedRate_r3_female','SEM_fittedRate_r3_female')
% %% Check the MCP-mCherry dosage
% 
% % r0
% for i=1:length(r0_RuntJB3)
%     Offset1(i) = nanmean(r0_RuntJB3(i).MeanOffsetVector);
% end
% 
% % r3, Runt null
% for i=1:length(r3_RuntNull)
%     Offset2(i) = nanmean(r3_RuntNull(i).MeanOffsetVector);
% end
% 
% MCP_dosage_Fig = figure;
% hold on
% plot([1 1 1], Offset1,'o')
% plot([2 2 2 ], Offset2,'o')
% title('MCP-mCherry dosage')
% xlabel('Constructs')
% ylabel('MCP-mCherry dosage (Offset)')
% %legend('','')
% xlim([0 3])
% StandardFigure(MCP_dosage_Fig,MCP_dosage_Fig.CurrentAxes)
% 
% % save the plots
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\';
% 
% saveas(MCP_dosage_Fig,[FigPath 'MCP-mCherry-dosage_r0_r3(run3)' , '.tif']); 
% saveas(MCP_dosage_Fig,[FigPath 'MCP-mCherry-dosage_r0_r3(run3)' , '.pdf']); 5
% 
% % Save the data



end