function main03_01_calculate_Fold_Change_InitialSlopes
% This script is for calculating the fold-change of initial slopes
% Fold-change (FC) = Rate(R>0) / Rate(R=0) ~ Rate(R>0)/Rate(nR=0)
% We will make an assumption that the basal rate of Txn for each construct
% is equivalent to the rate of Txn in the absence of Repressor binding
% sites.

%% Load the datasets
% Define the file paths (Dropbox)
% YJK Notes.
% I need to edit these functions such that they can read from a defined
% directories, as well as using a MovieDatabase in that directory.

% Opposing Gradients Results folder
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

% Load the fitted initial slopes
% InitialSlope = load([filePath, filesep, 'AveragedInitialRate_AsymmetricFit.mat']);

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
ColorChoice = [colorDict.magenta; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.blue; colorDict.purple; colorDict.thickpink]; 
%% Calculate the fold-change
% We need to learn how to calculate the error bars of division (both
% numereator and denominator have error range).
NC = 14;
NC = NC-11; % NC12, 13 and 14 in order

BasalActivity = average_fittedRate_r0(:,NC);

% r1
FC_r1 = average_fittedRate_r1(:,NC)./BasalActivity;
FC_r1_close = average_fittedRate_r1_close(:,NC)./BasalActivity;
FC_r1_mid = average_fittedRate_r1_mid(:,NC)./BasalActivity;

% r2
FC_r2 = average_fittedRate_r2(:,NC)./BasalActivity;
FC_r2_close = average_fittedRate_r2_close(:,NC)./BasalActivity;
FC_r2_far = average_fittedRate_r2_far(:,NC)./BasalActivity;

FC_r3 = average_fittedRate_r3(:,NC)./BasalActivity;

% Error (SEM) divided by the basal activity (r0)
% r1_error
FC_r1_SEM = SEM_fittedRate_r1(:,NC)./BasalActivity;
FC_r1_close_SEM = SEM_fittedRate_r1_close(:,NC)./BasalActivity;
FC_r1_mid_SEM = SEM_fittedRate_r1_mid(:,NC)./BasalActivity;

% r2_error
FC_r2_SEM = SEM_fittedRate_r2(:,NC)./BasalActivity;
FC_r2_close_SEM = SEM_fittedRate_r2_close(:,NC)./BasalActivity;
FC_r2_far_SEM = SEM_fittedRate_r2_far(:,NC)./BasalActivity;

% r3 error
FC_r3_SEM = SEM_fittedRate_r3(:,NC)./BasalActivity;

%% plot for checking
APaxis = 0:0.025:1;

%load('Runt-1min-200Hz-mixed_BGsubtracted-Averaged.mat')
load([filePath, filesep, 'Runt_TimeAveraged_mixedSex_NC14.mat'])
%tAveraged_Fluo_mixed
%SEFluo_tAveraged_mixed
% 2nd row is the 0-10 minutes of the time window
Runt_NC14_0_10min = AveragedFluo_tAveraged_mixed(2,:);
Runt_NC14_0_10min_SE = SEFluo_tAveraged_mixed(2,:);

%% r1 - different positions
hold on
errorbar(APaxis, FC_r1, FC_r1_SEM,'Color',ColorChoice(2,:))
errorbar(APaxis, FC_r1_close, FC_r1_close_SEM,'Color',ColorChoice(5,:))
errorbar(APaxis, FC_r1_mid, FC_r1_mid_SEM,'Color',ColorChoice(6,:))

% Limit the AP range to the region that makes sense.
% The region posterior to 35% (or 40%?) might not be appropriate for FC calculation,
% since the fitting was more sensitive and more error-prone. 
% xlim([0.2 0.35])
% xticks([0.2 0.25 0.3 0.35])

title('Fold-Change over AP axis- 1 Runt binding site')
xlabel('AP axis (EL)')
ylabel('Fold-Change')
legend('original','close','mid')

StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fold_Change\';
saveas(gcf,[FigPath 'FoldChange_r1_variants' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath 'FoldChange_r1_variants' , '_NC14' , '.pdf']); 
%% r1 - different positions (X-axis is Runt conc.)
% Let's try to limit this to 20%-35% of the AP bins (9th-15th AP bins)
AP1 = 20; % (% of EL)
AP2 = 35; % (% of EL)
APrange = (AP1/2.5 + 1):(AP2/2.5 - 1);

hold on
errorbar(Runt_NC14_0_10min(APrange), FC_r1(APrange), FC_r1_SEM(APrange))
errorbar(Runt_NC14_0_10min(APrange), FC_r1_close(APrange), FC_r1_close_SEM(APrange))
errorbar(Runt_NC14_0_10min(APrange), FC_r1_mid(APrange), FC_r1_mid_SEM(APrange))


title('Fold-Change over AP axis- 1 Runt binding site')
xlabel('Runt concentration (AU)')
ylabel('Fold-Change')
legend('original','close','mid')

StandardFigure(gcf,gca)
%% r2 - different construts
hold on
errorbar(APaxis, FC_r2, FC_r2_SEM,'Color',ColorChoice(3,:))
errorbar(APaxis, FC_r2_close, FC_r2_close_SEM,'Color',ColorChoice(7,:))
errorbar(APaxis, FC_r2_far, FC_r2_far_SEM,'Color',ColorChoice(8,:))

% Limit the AP range to the region that makes sense.
% The region posterior to 40% might not be appropriate for FC calculation,
% since the fitting was more sensitive and more error-prone. 
% xlim([0.2 0.35])
% xticks([0.2 0.25 0.3 0.35])

title('Fold-Change over AP axis- 2 Runt binding sites')
xlabel('AP axis (EL)')
ylabel('Fold-Change')
legend('original','close','far')

StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fold_Change\';
saveas(gcf,[FigPath 'FoldChange_r2_variants' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath 'FoldChange_r2_variants' , '_NC14' , '.pdf']); 
%% r2 - different configurations (Runt conc. as the X-axis)
% Let's try to limit this to 20%-35% of the AP bins (9th-15th AP bins)
AP1 = 20; % (% of EL)
AP2 = 35; % (% of EL)
APrange = (AP1/2.5 + 1):(AP2/2.5 - 1);

hold on
errorbar(Runt_NC14_0_10min(APrange), FC_r2(APrange), FC_r2_SEM(APrange))
errorbar(Runt_NC14_0_10min(APrange), FC_r2_close(APrange), FC_r2_close_SEM(APrange))
errorbar(Runt_NC14_0_10min(APrange), FC_r2_far(APrange), FC_r2_far_SEM(APrange))

title('Fold-Change over AP axis- 2 Runt binding sites')
xlabel('Runt concentration (AU)')
ylabel('Fold-Change')
legend('original','close','far')

StandardFigure(gcf,gca)
%% r0,1,2,3 
hold on
errorbar(APaxis, FC_r1, FC_r1_SEM,'Color',ColorChoice(2,:))                                                                                                                                                                            
errorbar(APaxis, FC_r2, FC_r2_SEM,'Color',ColorChoice(3,:))
errorbar(APaxis, FC_r3, FC_r3_SEM,'Color',ColorChoice(4,:))

% Limit the AP range to the region that makes sense.
% The region posterior to 40% might not be appropriate for FC calculation,
% since the fitting was more sensitive and more error-prone. 
% xlim([0.2 0.35])
% xticks([0.2 0.25 0.3 0.35])

title('Fold-Change over AP axis')
xlabel('AP axis (EL)')
ylabel('Fold-Change')
legend('1','2','3')

StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fold_Change\';
saveas(gcf,[FigPath 'FoldChange_rN_originals' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath 'FoldChange_rN_originals' , '_NC14' , '.pdf']); 
end