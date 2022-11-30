function main02_04_compare_Initial_Loading_Rates_diff_Fits
% Here, the goal is to compare the inital rate of RNAP loading, acquired
% with different methods of fitting.
% For example, I've used FitMeanAPAsymmetric.m, and FitMeanAPLinearSlope.m
% to extract the initial slope. 
% For me, the linear slope method is actually kinda noisy...I'll just
% compare as it is for now.

%% Load the datasets
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
% FitMeanAPAsymmetric
FittedRates_Asymm = load([filePath, filesep, 'AveragedInitialRate_AsymmetricFit.mat']);
% FitMeanAPLinearSlope
FittedRates_Linear = load([filePath, filesep, 'AveragedInitialRate_LinearFit.mat']);

%% Define fields for future usage
% First, Asymmetric fit
% 1) Mean
average_fittedRate_r0_Asymm = FittedRates_Asymm.average_fittedRate_r0;
average_fittedRate_r1_Asymm = FittedRates_Asymm.average_fittedRate_r1_female;
average_fittedRate_r2_Asymm = FittedRates_Asymm.average_fittedRate_r2_female;
average_fittedRate_r3_Asymm = FittedRates_Asymm.average_fittedRate_r3_female;
% 2) SEM
SEM_fittedRate_r0_Asymm = FittedRates_Asymm.SEM_fittedRate_r0;
SEM_fittedRate_r1_Asymm = FittedRates_Asymm.SEM_fittedRate_r1_female;
SEM_fittedRate_r2_Asymm = FittedRates_Asymm.SEM_fittedRate_r2_female;
SEM_fittedRate_r3_Asymm = FittedRates_Asymm.SEM_fittedRate_r3_female;

% Second, Linear fit
average_fittedRate_r0_Linear = FittedRates_Linear.average_fittedRate_r0;
average_fittedRate_r1_Linear = FittedRates_Linear.average_fittedRate_r1_female;
average_fittedRate_r2_Linear = FittedRates_Linear.average_fittedRate_r2_female;
average_fittedRate_r3_Linear = FittedRates_Linear.average_fittedRate_r3_female;
% 2) SEM
SEM_fittedRate_r0_Linear = FittedRates_Linear.SEM_fittedRate_r0;
SEM_fittedRate_r1_Linear = FittedRates_Linear.SEM_fittedRate_r1_female;
SEM_fittedRate_r2_Linear = FittedRates_Linear.SEM_fittedRate_r2_female;
SEM_fittedRate_r3_Linear = FittedRates_Linear.SEM_fittedRate_r3_female;

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
%% Plot altogether


% Second, Linear fit

% NC13
InitialRate_NC13_figure = figure;
nc = 2; % NC13
hold on
% Asymmetric
errorbar(0:0.025:1,average_fittedRate_r0_Asymm(:,nc),SEM_fittedRate_r0_Asymm(:,nc),'Color',ColorChoice(1,:))
% Females
errorbar(0:0.025:1,average_fittedRate_r1_Asymm(:,nc),SEM_fittedRate_r1_Asymm(:,nc),'Color',ColorChoice(2,:))
errorbar(0:0.025:1,average_fittedRate_r2_Asymm(:,nc),SEM_fittedRate_r2_Asymm(:,nc),'Color',ColorChoice(3,:))
errorbar(0:0.025:1,average_fittedRate_r3_Asymm(:,nc),SEM_fittedRate_r3_Asymm(:,nc),'Color',ColorChoice(4,:))

% Linear
errorbar(0:0.025:1,average_fittedRate_r0_Linear(:,nc),SEM_fittedRate_r0_Linear(:,nc),'--','Color',ColorChoice(1,:))
% Females
errorbar(0:0.025:1,average_fittedRate_r1_Linear(:,nc),SEM_fittedRate_r1_Linear(:,nc),'--','Color',ColorChoice(2,:))
errorbar(0:0.025:1,average_fittedRate_r2_Linear(:,nc),SEM_fittedRate_r2_Linear(:,nc),'--','Color',ColorChoice(3,:))
errorbar(0:0.025:1,average_fittedRate_r3_Linear(:,nc),SEM_fittedRate_r3_Linear(:,nc),'--','Color',ColorChoice(4,:))

xlim([0.15 0.5])
ylim([0 400])

legend('r0-Asymm','r1-female-Asymm','r2-female-Asymm','r3-female-Asymm',...
        'r0-Linear','r1-female-Linear','r2-female-Linear','r3-female-Linear')
xlabel('AP Position')
ylabel('Initial rate (AU/min)')
title('Initial rate of RNAP loading along AP axis, at NC 13')
StandardFigure(InitialRate_NC13_figure, InitialRate_NC13_figure.CurrentAxes)
%standardizeFigure_YJK(gca,legend)
% saveas(InitialRate_NC13_figure,[FigPath 'InitialRate_AsymmetricFit_r0123_female' , '_NC13' , '.tif']); 
% saveas(InitialRate_NC13_figure,[FigPath 'InitialRate_AsymmetricFit_r0123_female' , '_NC13' , '.pdf']); 

% NC14
InitialRate_NC14_figure = figure;
nc = 3; % NC13
hold on
% Asymmetric
errorbar(0:0.025:1,average_fittedRate_r0_Asymm(:,nc),SEM_fittedRate_r0_Asymm(:,nc),'Color',ColorChoice(1,:))
% Females
errorbar(0:0.025:1,average_fittedRate_r1_Asymm(:,nc),SEM_fittedRate_r1_Asymm(:,nc),'Color',ColorChoice(2,:))
errorbar(0:0.025:1,average_fittedRate_r2_Asymm(:,nc),SEM_fittedRate_r2_Asymm(:,nc),'Color',ColorChoice(3,:))
errorbar(0:0.025:1,average_fittedRate_r3_Asymm(:,nc),SEM_fittedRate_r3_Asymm(:,nc),'Color',ColorChoice(4,:))

% Linear
errorbar(0:0.025:1,average_fittedRate_r0_Linear(:,nc),SEM_fittedRate_r0_Linear(:,nc),'--','Color',ColorChoice(1,:))
% Females
errorbar(0:0.025:1,average_fittedRate_r1_Linear(:,nc),SEM_fittedRate_r1_Linear(:,nc),'--','Color',ColorChoice(2,:))
errorbar(0:0.025:1,average_fittedRate_r2_Linear(:,nc),SEM_fittedRate_r2_Linear(:,nc),'--','Color',ColorChoice(3,:))
errorbar(0:0.025:1,average_fittedRate_r3_Linear(:,nc),SEM_fittedRate_r3_Linear(:,nc),'--','Color',ColorChoice(4,:))

xlim([0.15 0.5])
ylim([0 400])

legend('r0-Asymm','r1-female-Asymm','r2-female-Asymm','r3-female-Asymm',...
        'r0-Linear','r1-female-Linear','r2-female-Linear','r3-female-Linear')
xlabel('AP Position')
ylabel('Initial rate (AU/min)')
title('Initial rate of RNAP loading along AP axis, at NC 14')
StandardFigure(InitialRate_NC13_figure, InitialRate_NC13_figure.CurrentAxes)
%standardizeFigure_YJK(gca,legend)
% saveas(InitialRate_NC14_figure,[FigPath 'InitialRate_AsymmetricFit_r0123_female' , '_NC14' , '.tif']); 
% saveas(InitialRate_NC14_figure,[FigPath 'InitialRate_AsymmetricFit_r0123_female' , '_NC14' , '.pdf']); 
end