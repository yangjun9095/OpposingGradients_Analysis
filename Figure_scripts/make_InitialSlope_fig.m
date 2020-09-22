% generate plots of initial slopes, fold-change, etc.
% Name : make_InitialSlope_fig

clear 
close all
addpath('../utilities')
% set ID variables
DropboxFolder = 'S:/YangJoon/Dropbox/OpposingGradient';
FigureRoot = 'S:/YangJoon/Dropbox/OpposingGradientsFigures/PipelineOutput';

% load data structure
load([DropboxFolder,filesep,'OpposingGradients_ProcessedData',filesep,'compiledData.mat'])

FigPath = [FigureRoot, filesep, 'InitialSlopes'];
mkdir(FigPath)

%% Information for the datasets
DataTypesForFit = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};

% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};
                
%% Color module
% This is defining the line color
% We have 8 distinct datasets, with or without Runt protein.
% I think selecting 8 distinguishable color sets, then changing the
% brightness by either adding/subtracting white would be a better idea than
% selecting 16 different color sets.

colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.purple = [171,133,172]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

%colorDict.magenta = [208,109,171]/255;
%colorDict.lightBlue = [115,142,193]/255;
colorDict.lightgreen = [205,214,209]/255;
colorDict.pink = [232,177,157]/255;
colorDict.thickpink = [132,27,69]/255;

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.darkgreen; colorDict.thickpink]; 

% For now, I'll add white (color+[1 1 1])/2 to make thinner color (for the
% Runt nulls)
%% Plot module
% define the figure handle, fig_name
fig_name = figure; 
hold on
errorbar(X, Y, Z, 'LineWidth',2,'Color',ColorChoice(1,:))

% xTicks, yTicks
xticks([0.2 0.3 0.4 0.5 0.6 0.7 0.8])

set(gca,'yticklabel',[])

% no title, no-caps on the axis labels
xlabel('')
ylabel('')

legend('','','Location','SouthWest')

box on

StandardFigure(fig_name, fig_name.CurrentAxes)

% Save the plot
figPath = '';
saveas(gcf,[figPath,filesep,'name','.tif']); 
saveas(gcf,[figPath,filesep,'name','.pdf']); 

%% initial slope over AP axis for each construct with/without Runt protein
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_slope = figure;

for construct=1:8 % total number of constructs (enhancers)
    clf
    hold on
    errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10},'LineWidth',2,'Color',ColorChoice(construct,:))
    errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10},'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2);

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 400])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 100 200 300 400])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('initial slope (AU)')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_slope, fig_slope.CurrentAxes)

    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

%% Initial slope for all Runt nulls
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_slope = figure;

hold on
for construct=1:8 % total number of constructs (enhancers)
    errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10},'LineWidth',2,'Color',(ColorChoice(construct,:)));
end

% xTicks, yTicks
xlim([0.15 0.6])
ylim([0 400])
xticks([0.2 0.3 0.4 0.5 0.6])
yticks([0 100 200 300 400])

%set(gca,'yticklabel',[])

% no title, no-caps on the axis labels
xlabel('embryo length')
ylabel('initial slope (AU)')

legend(constructNames{1:8},'Location','NorthEast')

box on

StandardFigure(fig_slope, fig_slope.CurrentAxes)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_AllConstructs_nullsOnly','.tif']); 
saveas(gcf,[FigPath,filesep,'InitialSlope_AllConstructs_nullsOnly','.pdf']); 

%% Initial Slope (Runt nulls only) - 1 Runt site
fig_slope = figure;

hold on
for construct=[1,2,5,6] % total number of constructs (enhancers)
    errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10},'LineWidth',2,'Color',(ColorChoice(construct,:)));
end

% xTicks, yTicks
xlim([0.15 0.6])
ylim([0 400])
xticks([0.2 0.3 0.4 0.5 0.6])
yticks([0 100 200 300 400])

%set(gca,'yticklabel',[])

% no title, no-caps on the axis labels
xlabel('embryo length')
ylabel('initial slope (AU)')

legend(constructNames{[1,2,5,6]},'Location','NorthEast')

box on

StandardFigure(fig_slope, fig_slope.CurrentAxes)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_AllConstructs_nulls_1RuntSite','.tif']); 
saveas(gcf,[FigPath,filesep,'InitialSlope_AllConstructs_nulls_1RuntSite','.pdf']); 

%% Initial Slope (Runt nulls only) - 2 Runt sites
fig_slope = figure;

hold on
for construct=[1,3,7,8,4] % total number of constructs (enhancers)
    errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10},'LineWidth',2,'Color',(ColorChoice(construct,:)));
end

% xTicks, yTicks
xlim([0.15 0.6])
ylim([0 400])
xticks([0.2 0.3 0.4 0.5 0.6])
yticks([0 100 200 300 400])

%set(gca,'yticklabel',[])

% no title, no-caps on the axis labels
xlabel('embryo length')
ylabel('initial slope (AU)')

legend(constructNames{[1,3,7,8,4]},'Location','NorthEast')

box on

StandardFigure(fig_slope, fig_slope.CurrentAxes)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_AllConstructs_nulls_2RuntSites','.tif']); 
saveas(gcf,[FigPath,filesep,'InitialSlope_AllConstructs_nulls_2RuntSites','.pdf']); 

%% Calculate the fold-change for each construct
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_FC = figure;

for construct=1:8 % total number of constructs (enhancers)
    % clf
    hold on
    FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
    fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
    fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
    FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;
    errorbar(APaxis, FC, FC_error, 'LineWidth',2,'Color',ColorChoice(construct,:))



    %set(gca,'yticklabel',[])
end
    
    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 1.4])
    xticks([0.2 0.3 0.4 0.5 0.6])
    %yticks([0 0.5 1 1.5 2])
    
    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('fold-change')

    legend(constructNames{1:8},'Location','NorthEast')

    box on

    StandardFigure(fig_FC, fig_FC.CurrentAxes)

    % Save the plot
    saveas(gcf,[FigPath,filesep,'fold-change-all','.tif']); 
    saveas(gcf,[FigPath,filesep,'fold-change-all','.pdf']); 
    
%% Fold-change for 1 Runt binding site
    
fig_FC = figure;

for construct=[1,2,5,6,4] % total number of constructs (enhancers)
    % clf
    hold on
    FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
    fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
    fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
    FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;
    
    errorbar(APaxis, FC, FC_error, 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))



    %set(gca,'yticklabel',[])
end
    yline(1,'--','LineWidth',2)
    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 1.4])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])
    
    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('fold-change')

    legend(constructNames{[1,2,5,6,4]},'Location','NorthEast')

    box on
    StandardFigure(fig_FC, fig_FC.CurrentAxes)

    % Save the plot
%     saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_0,3_ref','.tif']); 
%     saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_0,3_ref','.pdf']); 
    
%% Fold-change for 2 Runt binding sites
    
fig_FC = figure;

for construct=[1,3,7,8,4] % total number of constructs (enhancers)
    % clf
    hold on
    FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
    fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
    fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
    FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;
    
    errorbar(APaxis, FC, FC_error, 'LineWidth',2,'Color',ColorChoice(construct,:))



    %set(gca,'yticklabel',[])
end
    
    yline(1,'--','LineWidth',2)
    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 1.4])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])
    
    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('fold-change')

    legend(constructNames{[1,3,7,8,4]},'Location','NorthEast')

    box on

    StandardFigure(fig_FC, fig_FC.CurrentAxes)

    % Save the plot
%     saveas(fig_FC,[FigPath,filesep,'fold-change-2RuntSites_0,3ref','.tif']); 
%     saveas(fig_FC,[FigPath,filesep,'fold-change-2RuntSites_0,3ref','.pdf']); 

%% %% Let's see whether we can make sense of this fold-change data 
%% Load the Runt dataset to replace the AP axis to the Runt concentration at correpsonding APbin
inputPath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';
Runt_tAveraged = load([inputPath, filesep, 'Runt_TimeAveraged_mixedSex_NC14.mat' ]);

% 2nd row is for the 0-10 min time window
Runt = Runt_tAveraged.AveragedFluo_tAveraged_mixed(2,:);

% 
%% Load Bicoid dataset (mine and Liz&Jonathan's)
    
%% Fold-change (theoretical)
r= logspace(-3,3,100);% [R]/K_R

fc = 1./(1+r);

semilogx(r, fc, 'LineWidth', 2)

yticks([0 0.5 1])
xlabel('r ([R]/K_{R})')
ylabel('fold-change')

StandardFigure(gcf,gca)

% Save the plots
saveas(gcf,[FigPath,filesep,'fold-change-6A1R-prediction','.tif']); 
saveas(gcf,[FigPath,filesep,'fold-change-6A1R-prediction','.pdf']); 
%% Fold-change for 1 Runt binding site : fold-change versus Runt

APaxis = 0:0.025:1;

fig_FC = figure;

for construct=[1,2,5,6,4] % total number of constructs (enhancers)
    % clf
    hold on
    FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
    fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
    fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
    FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;
    
    errorbar(Runt, FC, FC_error, 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))



    %set(gca,'yticklabel',[])
end
    yline(1,'--','LineWidth',2)
    % xTicks, yTicks
    xlim([50 250])
    ylim([0 1.4])
    %xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])
    
    % no title, no-caps on the axis labels
    xlabel('repressor (AU)')
    ylabel('fold-change')

    legend(constructNames{[1,2,5,6,4]},'Location','SouthWest')

    box on
    StandardFigure(fig_FC, fig_FC.CurrentAxes)

    % Save the plot
    saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_RuntAxis','.tif']); 
    saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_RuntAxis','.pdf']); 
    
%% Fold-change (for 1 Runt binding site) versus Runt, along with parameter free fit

APaxis = 0:0.025:1;

fig_FC = figure;

% plot the fold-change data
for construct=[1,2,5,6] % total number of constructs (enhancers)
    % clf
    hold on
    FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
    fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
    fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
    FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;
    
    errorbar(Runt, FC, FC_error, 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))



    %set(gca,'yticklabel',[])
end

    % Make a parameter free prediction 
    KR = 10.^[2 3 4]; % binding affinity (AU)
    Runt_range = linspace(0,250,50);
    
    for i=1:length(KR)
        FC_prediction(:,i) =   1./(1+(Runt_range/KR(i)));
    end
    
    % Plot the fold-change prediction
    hold on
    for i=1:length(KR)
        plot(Runt_range, FC_prediction(:,i),'--','Color',([245 126 43]/255),'LineWidth',2)
    end
    
    
    
    yline(1,'--','LineWidth',2)
    % xTicks, yTicks
    xlim([50 250])
    ylim([0 1.4])
    %xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])
    
    %set(gca,'xscale','log')
    
    % no title, no-caps on the axis labels
    xlabel('repressor (AU)')
    ylabel('fold-change')

    legend(constructNames{[1,2,5,6]},'Location','SouthWest')

    box on
    StandardFigure(fig_FC, fig_FC.CurrentAxes)

    % Save the plot
    saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_RuntAxis_parameterFreePrediction','.tif']); 
    saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_RuntAxis_parameterFreePrediction','.pdf']); 
    
%% Let's pick only one construct

%% 001 (fold-change versus Runt)
APaxis = 0:0.025:1;

fig_FC = figure;

% plot the fold-change data
construct = 5 % total number of constructs (enhancers)
% clf
hold on
FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;

errorbar(Runt, FC, FC_error, 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))


% Make a parameter free prediction 
KR = 10.^[2 3 4]; % binding affinity (AU)
Runt_range = linspace(0,250,50);

for i=1:length(KR)
    FC_prediction(:,i) =   1./(1+(Runt_range/KR(i)));
end

% Plot the fold-change prediction
hold on
for i=1:length(KR)
    plot(Runt_range, FC_prediction(:,i),'--','Color',([245 126 43]/255),'LineWidth',2)
end
    
    
    
yline(1,'--','LineWidth',2)
% xTicks, yTicks
xlim([50 250])
ylim([0 1.4])
%xticks([0.2 0.3 0.4 0.5 0.6])
yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])

%set(gca,'xscale','log')

% no title, no-caps on the axis labels
xlabel('repressor (AU)')
ylabel('fold-change')

legend(constructNames{5},'Location','SouthWest')

box on
StandardFigure(fig_FC, fig_FC.CurrentAxes)

% Save the plot
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_001_parameterFreePrediction','.tif']); 
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_001_parameterFreePrediction','.pdf']); 

%% 001, fold-change versus AP axis
APaxis = 0:0.025:1;

fig_FC = figure;

% plot the fold-change data
construct = 5 % total number of constructs (enhancers)
% clf
hold on
FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;

errorbar(APaxis, FC, FC_error, 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))


% Make a parameter free prediction 
clear FC_prediction
KR = 10.^[2 2.5 3 3.5 4]; % binding affinity (AU)
Runt_range = linspace(0,250,50);

for i=1:length(KR)
    FC_prediction(:,i) =   1./(1+(Runt/KR(i)));
end

% Plot the fold-change prediction
hold on
for i=1:length(KR)
    plot(APaxis, FC_prediction(:,i),'--','Color',([245 126 43]/255),'LineWidth',2)
end
    
    
    
yline(1,'--','LineWidth',2)
% xTicks, yTicks
xlim([0.15 0.6])
ylim([0 1.4])
xticks([0.2 0.3 0.4 0.5 0.6])
yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])

%set(gca,'xscale','log')

% no title, no-caps on the axis labels
xlabel('embryo length')
ylabel('fold-change')

legend(constructNames{5},'Location','SouthWest')

box on
StandardFigure(fig_FC, fig_FC.CurrentAxes)

% Save the plot
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_001_paramFreePredicton_APaxis','.tif']); 
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_001_paramFreePredicton_APaxis','.pdf']); 

%% 010, fold-change versus Runt

APaxis = 0:0.025:1;

fig_FC = figure;

% plot the fold-change data
construct = 6 % total number of constructs (enhancers)
% clf
hold on
FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;

errorbar(Runt, FC, FC_error, 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))


% Make a parameter free prediction 
KR = 10.^[2 3 4]; % binding affinity (AU)
Runt_range = linspace(0,250,50);

for i=1:length(KR)
    FC_prediction(:,i) =   1./(1+(Runt_range/KR(i)));
end

% Plot the fold-change prediction
hold on
for i=1:length(KR)
    plot(Runt_range, FC_prediction(:,i),'--','Color',([245 126 43]/255),'LineWidth',2)
end
    
    
    
yline(1,'--','LineWidth',2)
% xTicks, yTicks
xlim([50 250])
ylim([0 1.4])
%xticks([0.2 0.3 0.4 0.5 0.6])
yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])

%set(gca,'xscale','log')

% no title, no-caps on the axis labels
xlabel('repressor (AU)')
ylabel('fold-change')

legend(constructNames{construct},'Location','SouthWest')

box on
StandardFigure(fig_FC, fig_FC.CurrentAxes)

% Save the plot
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_010_parameterFreePrediction','.tif']); 
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_010_parameterFreePrediction','.pdf']); 

%% 010, fold-change versus AP axis
APaxis = 0:0.025:1;

fig_FC = figure;

% plot the fold-change data
construct = 6 % total number of constructs (enhancers)
% clf
hold on
FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;

errorbar(APaxis, FC, FC_error, 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))


% Make a parameter free prediction 
clear FC_prediction
KR = 10.^[2 2.5 3 3.5 4]; % binding affinity (AU)
Runt_range = linspace(0,250,50);

for i=1:length(KR)
    FC_prediction(:,i) =   1./(1+(Runt/KR(i)));
end

% Plot the fold-change prediction
hold on
for i=1:length(KR)
    plot(APaxis, FC_prediction(:,i),'--','Color',([245 126 43]/255),'LineWidth',2)
end
    
    
    
yline(1,'--','LineWidth',2)
% xTicks, yTicks
xlim([0.15 0.6])
ylim([0 1.4])
xticks([0.2 0.3 0.4 0.5 0.6])
yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])

%set(gca,'xscale','log')

% no title, no-caps on the axis labels
xlabel('embryo length')
ylabel('fold-change')

legend(constructNames{5},'Location','SouthWest')

box on
StandardFigure(fig_FC, fig_FC.CurrentAxes)

% Save the plot
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_010_paramFreePredicton_APaxis','.tif']); 
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_010_paramFreePredicton_APaxis','.pdf']); 


%% 000, fold-change versus AP axis
APaxis = 0:0.025:1;

fig_FC = figure;

% plot the fold-change data
construct = 1 % total number of constructs (enhancers)
% clf
hold on
FC = compiledData{construct+1,9}./compiledData{construct+1+8,9};
fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
FC_error = sqrt(fracError1.^2 + fracError2.^2).*FC;

errorbar(APaxis, FC, FC_error, 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))


% Make a parameter free prediction 
clear FC_prediction
KR = 10.^[2 2.5 3 3.5 4]; % binding affinity (AU)
Runt_range = linspace(0,250,50);

for i=1:length(KR)
    FC_prediction(:,i) =   1./(1+(Runt/KR(i)));
end

% Plot the fold-change prediction
hold on
for i=1:length(KR)
    plot(APaxis, FC_prediction(:,i),'--','Color',([245 126 43]/255),'LineWidth',2)
end
    
    
    
yline(1,'--','LineWidth',2)
% xTicks, yTicks
xlim([0.15 0.6])
ylim([0 1.4])
xticks([0.2 0.3 0.4 0.5 0.6])
yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])

%set(gca,'xscale','log')

% no title, no-caps on the axis labels
xlabel('embryo length')
ylabel('fold-change')

legend(constructNames{construct},'Location','SouthWest')

box on
StandardFigure(fig_FC, fig_FC.CurrentAxes)

% Save the plot
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_000_paramFreePredicton_APaxis','.tif']); 
saveas(fig_FC,[FigPath,filesep,'fold-change-1RuntSite_000_paramFreePredicton_APaxis','.pdf']); 
