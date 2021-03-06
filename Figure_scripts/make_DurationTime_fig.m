% generate plots of initial slopes, fold-change, etc.
% Name : make_InitialSlope_fig

clear 
close all
addpath('../utilities')
% set ID variables
DropboxFolder = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';
FigureRoot = 'S:/YangJoon/Dropbox/OpposingGradientsFigures/PipelineOutput';

% load data structure
load([DropboxFolder,filesep,'compiledData.mat'])

FigPath = [FigureRoot, filesep, 'DurationTime'];
mkdir(FigPath)


%% ID variables
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

%% Duration time over AP axis for each construct with/without Runt protein
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_slope = figure;

for construct=1:8 % total number of constructs (enhancers)
    clf
    hold on
    errorbar(APaxis, compiledData{construct+1,13}, compiledData{construct+1,14},'LineWidth',2,'Color',ColorChoice(construct,:))
    errorbar(APaxis, compiledData{construct+1+8,13}, compiledData{construct+1+8,14},'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2);

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 60])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 10 20 30 40 50 60])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('duration (min)')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_slope, fig_slope.CurrentAxes)
    pause(1)
    % Save the plot
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

%% fold-change : Duration time over AP axis for each construct with/without Runt protein
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_slope = figure;

for construct=1:8 % total number of constructs (enhancers)
    clf
    FC =  compiledData{construct+1,13}./compiledData{construct+1+8,13};
    fracError1 = compiledData{construct+1,14}./compiledData{construct+1,13};
    fracError2 = compiledData{construct+1+8,14}./compiledData{construct+1+8,13};
    FC_SEM = sqrt(fracError1.^2 + fracError2.^2).*FC;
    
    hold on
    errorbar(APaxis, FC, FC_SEM,'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))
    
    % a guide line for FC=1
    yline(1,'--')
    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 1.4])
    xticks([0.2 0.3 0.4 0.5 0.6])
    %yticks([0 10 20 30 40 50 60])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('fold-change')

    legend(constructNames{construct},'Location','NorthEast')

    box on

    StandardFigure(fig_slope, fig_slope.CurrentAxes)
    pause(1)
    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC.pdf']); 
end

%% Part2. generate plots of T_on and T_off over the AP axis for each construct (Runt Null vs. Runt WT)

