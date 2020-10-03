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

%% Step1. Checking the credibility of individual rate fitting (MeanFitAPAsymmetric)

%% individual rate fit with their own errors from the fitting
% Start with 000
for construct = 1:16

    % Define the AP axis
    APaxis = 0:0.025:1;

    % number of embryos
    [~,~,numEmbryos] = size(compiledData{construct+1,3});
    fig_slope = figure;

    % Extract the fitted rate and SD (from confidence intervals)
    fittedRate_ind = compiledData{construct+1, 3};
    fittedRate_ind_SD = compiledData{construct+1, 4};

    NC=3; % nc14

    clf
    clear h
    hold on
    for emb=1:numEmbryos % total number of constructs (enhancers)
        h(emb) = errorbar(APaxis, fittedRate_ind(:,NC,emb),fittedRate_ind_SD(:,NC,emb) ,'LineWidth',1,'CapSize',0)%,'Color',ColorChoice(construct,:))
        %pause
    end
    % plot the averaged one with its STD
    h(numEmbryos+1) = errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10},'-k','LineWidth',2,'CapSize',0),%'Color',ColorChoice(construct,:))

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 400])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 100 200 300 400])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('initial slope (AU)')

    legend([h(end)],constructNames{construct},'Location','NorthEast')

    box on

    StandardFigure(fig_slope, fig_slope.CurrentAxes)

    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_individual.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_individual.pdf']); 
    %pause(1)
end
%% coefficient of variation of individual fits
for construct = 1:16

    % Define the AP axis
    APaxis = 0:0.025:1;

    % number of embryos
    [~,~,numEmbryos] = size(compiledData{construct+1,3});
    fig_slope = figure;

    % Extract the fitted rate and SD (from confidence intervals)
    fittedRate_ind = compiledData{construct+1, 3};
    fittedRate_ind_SD = compiledData{construct+1, 4};

    NC=3; % nc14

    clf
    clear h
    hold on
    for emb=1:numEmbryos % total number of constructs (enhancers)
        h(emb) = plot(APaxis, fittedRate_ind_SD(:,NC,emb)./fittedRate_ind(:,NC,emb) ,'LineWidth',1)%,'Color',ColorChoice(construct,:))
    end
    % plot the averaged one with its STD
    h(numEmbryos+1) = plot(APaxis, compiledData{construct+1,10}./compiledData{construct+1,9},'-k','LineWidth',2),%'Color',ColorChoice(construct,:))

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 1])
    xticks([0.2 0.3 0.4 0.5 0.6])
%     yticks([0 100 200 300 400])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('coefficient of variation')

    legend([h(end)],constructNames{construct},'Location','NorthEast')

    box on

    StandardFigure(fig_slope, fig_slope.CurrentAxes)

    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_CV_individual.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_CV_individual.pdf']); 
    %pause(1)
end

%% STD of individual fits
for construct = 1:16

    % Define the AP axis
    APaxis = 0:0.025:1;

    % number of embryos
    [~,~,numEmbryos] = size(compiledData{construct+1,3});
    fig_slope = figure;

    % Extract the fitted rate and SD (from confidence intervals)
    fittedRate_ind = compiledData{construct+1, 3};
    fittedRate_ind_SD = compiledData{construct+1, 4};

    NC=3; % nc14

    clf
    clear h
    hold on
    for emb=1:numEmbryos % total number of constructs (enhancers)
        h(emb) = plot(APaxis, fittedRate_ind_SD(:,NC,emb) ,'LineWidth',1)%,'Color',ColorChoice(construct,:))
    end
    % plot the averaged one with its STD
    h(numEmbryos+1) = plot(APaxis, compiledData{construct+1,10},'-k','LineWidth',2),%'Color',ColorChoice(construct,:))

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 80])
    xticks([0.2 0.3 0.4 0.5 0.6])
%     yticks([0 100 200 300 400])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('coefficient of variation')

    legend([h(end)],constructNames{construct},'Location','NorthEast')

    box on

    StandardFigure(fig_slope, fig_slope.CurrentAxes)

    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_STD_individual.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_STD_individual.pdf']); 
    %pause(1)
end

%% Filter out the initial rate fits with C.V of higher than 0.2 
% ( This part could be moved to the
% main_compile_InitialSlopes_AllConstructs.m script)

% process the compiledData
compiledData{1,15} = 'fittedRate_filtered_mean';
compiledData{1,16} = 'fittedRate_filtered_SEM';
for construct = 1:16
    vars = {'fittedRate_ind','fittedRate_ind_SD',...
            'CV_ind','CV_filter', 'fittedRate_ind_filtered',...
            'fittedRate_filtered_mean', 'fittedRate_filtered_std', 'fittedRate_filtered_sem'};
    clear vars
    % number of embryos
    [~,~,numEmbryos] = size(compiledData{construct+1,3});

    % Extract the fitted rate and SD (from confidence intervals)
    fittedRate_ind = squeeze(compiledData{construct+1, 3}(:,3,:));
    fittedRate_ind_SD = squeeze(compiledData{construct+1, 4}(:,3,:));
    
    % Coefficient of variation of individual embryos (error from the
    % fitting)
    CV_ind = fittedRate_ind_SD./fittedRate_ind; 
    CV_filter = CV_ind <0.2;
    
    fittedRate_ind_filtered = fittedRate_ind .* CV_filter;
    % convert the zeros to NaNs (as we filtered out the values using 0/1
    % matrix, so 0 actually means NaN.
    fittedRate_ind_filtered(fittedRate_ind_filtered==0) = nan;
    
    fittedRate_filtered_mean = nanmean(fittedRate_ind_filtered,2);
    fittedRate_filtered_std = nanstd(fittedRate_ind_filtered,0,2);
    fittedRate_filtered_sem = fittedRate_filtered_std./sqrt(numEmbryos);
    
    % Save into the compiledData structure
    compiledData{construct+1,15} = fittedRate_filtered_mean;
    compiledData{construct+1,16} = fittedRate_filtered_std;
end


%% Plots of initial slope - post filtering (Runt WT versus Runt nulls) - post filtering step

APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_slope = figure;

FigPath = [FigureRoot, filesep, 'InitialSlopes', filesep, 'filtered'];
mkdir(FigPath)

for construct=1:8 % total number of constructs (enhancers)
    clf
    hold on
    errorbar(APaxis, compiledData{construct+1,15}, compiledData{construct+1,16},'LineWidth',2,'Color',ColorChoice(construct,:))
    errorbar(APaxis, compiledData{construct+1+8,15}, compiledData{construct+1+8,16},'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2);

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



%% 
%% Save the structure, compiledData for future usage (in plotting scripts)
save('S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\compiledData.mat',...
        'compiledData')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% from the V1 %%%%%%%%%%%%%%%%%%%%%%%%



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
