function main_FoldChange_anlaysis_initialSlopes
%% This script is for different types of analyses for initial slopes in nc14
% First, we will start with calculating the fold-change of initial slopes
% for [0,0,0].

%% Load the datasets
DropboxFolder = 'S:\YangJoon\Dropbox\OpposingGradient'
load([DropboxFolder,filesep,'OpposingGradients_ProcessedData',filesep,'InitialSlopes_ONnuclei_AllConstructs.mat'])

%% Define the fields
% average_fittedRate_
% SEM_fittedRate_

FC_000 = average_fittedRate_r0(:,3)./average_fittedRate_r0(:,3);

FC_100 = average_fittedRate_r1(:,3)./average_fittedRate_r0(:,3);
FC_010 = average_fittedRate_r1_mid(:,3)./average_fittedRate_r0(:,3);
FC_001 = average_fittedRate_r1_close(:,3)./average_fittedRate_r0(:,3);

FC_011 = average_fittedRate_r2(:,3)./average_fittedRate_r0(:,3);
FC_101 = average_fittedRate_r2_far(:,3)./average_fittedRate_r0(:,3);
FC_110 = average_fittedRate_r2_close(:,3)./average_fittedRate_r0(:,3);

FC_111 = average_fittedRate_r3(:,3)./average_fittedRate_r0(:,3);

% Errors
FC_000_error = SEM_fittedRate_r0(:,3)./average_fittedRate_r0(:,3);
FC_100_error = SEM_fittedRate_r1(:,3)./average_fittedRate_r0(:,3);
FC_010_error = SEM_fittedRate_r1_mid(:,3)./average_fittedRate_r0(:,3);
FC_001_error = SEM_fittedRate_r1_close(:,3)./average_fittedRate_r0(:,3);

FC_011_error = SEM_fittedRate_r2(:,3)./average_fittedRate_r0(:,3);
FC_110_error = SEM_fittedRate_r2_close(:,3)./average_fittedRate_r0(:,3);
FC_101_error = SEM_fittedRate_r2_far(:,3)./average_fittedRate_r0(:,3);

FC_111_error = SEM_fittedRate_r3(:,3)./average_fittedRate_r0(:,3);

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
                colorDict.magenta; colorDict.purple; colorDict.thickpink;...
                colorDict.darkgreen]; 
%% Plot to check
% 1) r1 variants
APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
errorbar(APaxis, FC_100, FC_100_error,'LineWidth',2, 'Color',ColorChoice(2,:))
errorbar(APaxis, FC_001, FC_001_error,'LineWidth',2, 'Color',ColorChoice(5,:))
errorbar(APaxis, FC_010, FC_010_error,'LineWidth',2, 'Color',ColorChoice(6,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[1,0,0]','[0,0,1]','[0,1,0]','Location','Southeast')

StandardFigure(gcf,gca)

DropboxFolder = 'S:\YangJoon\Dropbox';
figPath = [DropboxFolder,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\Fold_Change'];
saveas(gcf,[figPath,filesep,'fold_change_r1_variants_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_r1_variants_NC14' , '.tif']);  

%%
% 2) r2 variants
APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
errorbar(APaxis, FC_011, FC_011_error,'LineWidth',2, 'Color',ColorChoice(3,:))
errorbar(APaxis, FC_110, FC_110_error,'LineWidth',2, 'Color',ColorChoice(7,:))
errorbar(APaxis, FC_101, FC_101_error,'LineWidth',2, 'Color',ColorChoice(8,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[0,1,1]','[1,1,0]','[1,0,1]','Location','Northeast')

StandardFigure(gcf,gca)

DropboxFolder = 'S:\YangJoon\Dropbox';
figPath = [DropboxFolder,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\Fold_Change'];
saveas(gcf,[figPath,filesep,'fold_change_r2_variants_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_r2_variants_NC14' , '.tif']);  
%% r0->1->2->3
APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2,'Color',ColorChoice(1,:))
errorbar(APaxis, FC_100, FC_100_error,'LineWidth',2,'Color',ColorChoice(2,:))
errorbar(APaxis, FC_011, FC_011_error,'LineWidth',2,'Color',ColorChoice(3,:))
errorbar(APaxis, FC_111, FC_111_error,'LineWidth',2,'Color',ColorChoice(4,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[1,0,0]','[0,1,1]','[1,1,1]','Location','Southeast')

StandardFigure(gcf,gca)

DropboxFolder = 'S:\YangJoon\Dropbox';
figPath = [DropboxFolder,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\Fold_Change'];
saveas(gcf,[figPath,filesep,'fold_change_r0123_originalPositions_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_r0123_originalPositions_NC14' , '.pdf']); 

%% Let's try to build up the prediction of 2 binding sites one by one
% Here, the phenomenological model is that the fold-change is a consequence
% from individual binding sites in an independent manner. Thus, the effect
% of fold-change from one site and the other (site) is multiplicative.

% Prediction is that FC_1+2 = FC1 * FC2
Prediction_FC_101 = FC_100 .* FC_001;
Prediction_FC_101_error = sqrt(FC_100_error.^2 + FC_001_error.^2);

APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
% errorbar(APaxis, FC_100, FC_100_error, 'Color',ColorChoice(2,:))
% errorbar(APaxis, FC_001, FC_001_error, 'Color',ColorChoice(5,:))
errorbar(APaxis, Prediction_FC_101, Prediction_FC_101_error,'--','LineWidth',2, 'Color',ColorChoice(9,:))
errorbar(APaxis, FC_101, FC_101_error, 'LineWidth',2, 'Color',ColorChoice(8,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[1,0,0] x [0,0,1]','Measurement-[1,0,1]','Location','Southeast')

saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[101]_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[101]_NC14' , '.pdf']); 

%% [010] + [001] = [011]?

% Prediction is that FC_1+2 = FC1 * FC2
Prediction_FC_011 = FC_010 .* FC_001;
Prediction_FC_011_error = sqrt(FC_010_error.^2 + FC_001_error.^2);

APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
% shadedErrorBar(APaxis, Prediction_FC_011, Prediction_FC_011_error)
errorbar(APaxis, Prediction_FC_011, Prediction_FC_011_error,'--','LineWidth',2, 'Color',ColorChoice(9,:))
errorbar(APaxis, FC_011, FC_011_error, 'LineWidth',2, 'Color',ColorChoice(3,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[0,1,0] x [0,0,1]','Measurement-[0,1,1]','Location','Northeast')

saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[011]_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[011]_NC14' , '.pdf']); 

%% [100]+[010] = [110]?
% Prediction is that FC_1+2 = FC1 * FC2
Prediction_FC_110 = FC_100 .* FC_010;
Prediction_FC_110_error = sqrt(FC_100_error.^2 + FC_010_error.^2);

APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
% shadedErrorBar(APaxis, Prediction_FC_011, Prediction_FC_011_error)
errorbar(APaxis, Prediction_FC_110, Prediction_FC_110_error,'--','LineWidth',2, 'Color',ColorChoice(9,:))
errorbar(APaxis, FC_110, FC_110_error, 'LineWidth',2, 'Color',ColorChoice(7,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[1,0,0]x[0,1,0]','Measurement-[1,1,0]','Location','Southeast')

saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[110]_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[110]_NC14' , '.pdf']); 

%% 2 sites + 1 site FC prediction
%% [100]+[011] = [111]?
% Prediction is that FC_1+2 = FC1 * FC2
Prediction_FC_111_100_011 = FC_100 .* FC_011;
Prediction_FC_111_100_011_error = sqrt(FC_100_error.^2 + FC_011_error.^2);

APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
% shadedErrorBar(APaxis, Prediction_FC_011, Prediction_FC_011_error)
errorbar(APaxis, Prediction_FC_111_100_011, Prediction_FC_111_100_011_error,'--','LineWidth',2, 'Color',ColorChoice(9,:))
errorbar(APaxis, FC_111, FC_111_error, 'LineWidth',2, 'Color',ColorChoice(4,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[1,0,0]x[0,1,1]','Measurement-[1,1,1]','Location','Southeast')

saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[100+011]_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[100+011]_NC14' , '.pdf']); 

%% [001]+[110] = [111]?
% Prediction is that FC_1+2 = FC1 * FC2
Prediction_FC_111_001_110 = FC_001 .* FC_110;
Prediction_FC_111_001_110_error = sqrt(FC_001_error.^2 + FC_110_error.^2);

APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
% shadedErrorBar(APaxis, Prediction_FC_011, Prediction_FC_011_error)
errorbar(APaxis, Prediction_FC_111_001_110, Prediction_FC_111_001_110_error,'--','LineWidth',2, 'Color',ColorChoice(9,:))
errorbar(APaxis, FC_111, FC_111_error, 'LineWidth',2, 'Color',ColorChoice(4,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[1,1,0]x[0,0,1]','Measurement-[1,1,1]','Location','Northeast')

saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[001+110]_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[001+110]_NC14' , '.pdf']); 

%% [010]+[101] = [111]?
% Prediction is that FC_1+2 = FC1 * FC2
Prediction_FC_111_010_101 = FC_010 .* FC_101;
Prediction_FC_111_010_101_error = sqrt(FC_010_error.^2 + FC_101_error.^2);

APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
% shadedErrorBar(APaxis, Prediction_FC_011, Prediction_FC_011_error)
errorbar(APaxis, Prediction_FC_111_010_101, Prediction_FC_111_010_101_error,'--','LineWidth',2, 'Color',ColorChoice(9,:))
errorbar(APaxis, FC_111, FC_111_error, 'LineWidth',2, 'Color',ColorChoice(4,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[0,1,0]x[1,0,1]','Measurement-[1,1,1]','Location','Northeast')

saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[010+101]_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_[010+101]_NC14' , '.pdf']); 

%% [1,0,0] + [0,1,0] + [0,0,1]
% Prediction is that FC_1+2 = FC1 * FC2
Prediction_FC_111_ind = FC_100 .* FC_010.*FC_001;
Prediction_FC_111_ind_error = sqrt(FC_100_error.^2 + FC_010_error.^2 + FC_001_error.^2);

APaxis = 0:0.025:1;

hold on
errorbar(APaxis, FC_000, FC_000_error,'LineWidth',2, 'Color',ColorChoice(1,:))
% shadedErrorBar(APaxis, Prediction_FC_011, Prediction_FC_011_error)
errorbar(APaxis, Prediction_FC_111_ind, Prediction_FC_111_ind_error,'--','LineWidth',2, 'Color',ColorChoice(9,:))
errorbar(APaxis, FC_111, FC_111_error, 'LineWidth',2, 'Color',ColorChoice(4,:))

% Draw vertical lines to show the region that we can focus on
xline(0.2,'--')
xline(0.4,'--')

xlim([0.2 0.4])
ylim([0 1.6])
xticks([0.2 0.25 0.3 0.35 0.4])

xlabel('AP axis (EL)')
ylabel('fold-change')
title('fold-change across the AP axis')
legend('[0,0,0]','[1,0,0]x[0,1,0]x[0,0,1]','Measurement-[1,1,1]','Location','Northeast')

saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_all3Sites_NC14' , '.tif']); 
saveas(gcf,[figPath,filesep,'fold_change_Prediction_data_all3Sites_NC14' , '.pdf']);
%% Part2. How do we quantify the deviation of FC from the prediction?
% One idea is just summing up the FC over AP axis (range)

end