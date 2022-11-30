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

FigPath = [FigureRoot, filesep, 'Fraction_competent'];
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

%% Fraction competent over AP (nc14)
% Start with 000
for construct = 1:8

    % Define the AP axis
    APaxis = 0:0.025:1;

    % number of embryos
    [~,~,numEmbryos] = size(compiledData{construct+1,3});
    fig = figure;

    % Extract the fitted rate and SD (from confidence intervals)
%     fractionON = compiledData{construct+1, 16};
%     fractionON_SEM = compiledData{construct+1, 17};

    NC=3; % nc14

    clf 
    hold on
    
    % plot the fraction competent with SEM
    if construct==1 || construct==4
        h(1) = errorbar(APaxis, compiledData{construct+1,16}(:,NC), compiledData{construct+1,17}(:,NC),'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))
        h(2) = errorbar(APaxis, compiledData{construct+1+8,16}, compiledData{construct+1+8,17},'LineWidth',2,'CapSize',0,'Color',(ColorChoice(construct,:)+[1 1 1])/2);
    else
        h(1) = errorbar(APaxis, compiledData{construct+1,16}(:,NC), compiledData{construct+1,17}(:,NC),'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))
        h(2) = errorbar(APaxis, compiledData{construct+1+8,16}(:,NC), compiledData{construct+1+8,17}(:,NC),'LineWidth',2,'CapSize',0,'Color',(ColorChoice(construct,:)+[1 1 1])/2);
    end
    
    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 1.2])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 0.2 0.4 0.6 0.8 1 1.2])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('fraction competent')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')
    box on

    StandardFigure(fig, fig.CurrentAxes)

    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
    %pause(1)
end

%% Part2. Further analyses
%% fold-change of competent nuclei
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_mRNA = figure;


for construct=1:8 % total number of constructs (enhancers)
    clf
    if construct==1 || construct==4
        FC = compiledData{construct+1,16}(:,NC)./compiledData{construct+1+8,16};
        fracError1 = compiledData{construct+1,17}(:,NC)./ compiledData{construct+1,16}(:,NC);
        fracError2 = compiledData{construct+1+8,17}./ compiledData{construct+1+8,16};
    else
        FC = compiledData{construct+1,16}(:,NC)./compiledData{construct+1+8,16}(:,NC);
        fracError1 = compiledData{construct+1,17}(:,NC)./ compiledData{construct+1,16}(:,NC);
        fracError2 = compiledData{construct+1+8,17}(:,NC)./ compiledData{construct+1+8,16}(:,NC);
    end
    FC_SEM = sqrt(fracError1.^2 + fracError2.^2).*FC;
    errorbar(APaxis, FC, FC_SEM,'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))

    yline(1,'--')

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 1.2])
    xticks([0.2 0.3 0.4 0.5 0.6])
    %yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('fold-change')

    legend(constructNames{construct},'Location','SouthEast')

    box on

    StandardFigure(fig_mRNA, fig_mRNA.CurrentAxes)
    pause(1)
    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC','.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC','.pdf']); 
    
    % Optional
    % Calculate the spatially averaged fold-change(repression)
    FC_spatial_mean_fracComp(construct) = nanmean(FC(9:13));
    FC_spatial_error_fracComp(construct) = sqrt(sum(FC_SEM(9:13).^2))/length(9:13);
end
%% plot the spatially averaged fold-change for each construct

construct = [1,2,5,6,3,7,8,4];
order = 1:8;
errorbar(order, FC_spatial_mean_fracComp(construct), FC_spatial_error_fracComp(construct),'o')
yline(1,'--')
xticklabels({'000','100','001','010','011','110','101','111'})

xlabel('construct')
ylabel('fold-change')

xlim([0 9])
ylim([0 1.2])
StandardFigure(gcf,gca)

saveas(gcf,[FigPath,filesep,'FC_20%-30%_averaged''.tif']); 
saveas(gcf,[FigPath,filesep,'FC_20%-30%_averaged''.pdf']); 
%% fold-change(repression) of fraction of competent nuclei
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_mRNA = figure;


for construct=1:8 % total number of constructs (enhancers)
    clf
    if construct==1 || construct==4
        FC = compiledData{construct+1+8,16}./compiledData{construct+1,16}(:,NC);
        fracError1 = compiledData{construct+1,17}(:,NC)./ compiledData{construct+1,16}(:,NC);
        fracError2 = compiledData{construct+1+8,17}./ compiledData{construct+1+8,16};
    else
        FC = compiledData{construct+1+8,16}(:,NC)./compiledData{construct+1,16}(:,NC);
        fracError1 = compiledData{construct+1,17}(:,NC)./ compiledData{construct+1,16}(:,NC);
        fracError2 = compiledData{construct+1+8,17}(:,NC)./ compiledData{construct+1+8,16}(:,NC);
    end
    FC_SEM = sqrt(fracError1.^2 + fracError2.^2).*FC;
    errorbar(APaxis, FC, FC_SEM,'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))

    yline(1,'--')

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 10])
    xticks([0.2 0.3 0.4 0.5 0.6])
    %yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('repression')

    legend(constructNames{construct},'Location','NorthEast')

    box on

    StandardFigure(fig_mRNA, fig_mRNA.CurrentAxes)
    pause(1)
    % Save the plot
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'_repression','.tif']); 
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'_repression','.pdf']); 
    
    % Optional
    % Calculate the spatially averaged fold-change(repression)
    FC_spatial_mean_fracComp(construct) = nanmean(FC(9:13));
    FC_spatial_error_fracComp(construct) = sqrt(sum(FC_SEM(9:13).^2))/length(9:13);
end
%% plot the spatially averaged fold-change for each construct

construct = [1,2,5,6,3,7,8,4];
order = 1:8;
errorbar(order, FC_spatial_mean_fracComp(construct), FC_spatial_error_fracComp(construct),'o')
yline(1,'--')
xticklabels({'000','100','001','010','011','110','101','111'})

xlabel('construct')
ylabel('repression')

StandardFigure(gcf,gca)

saveas(gcf,[FigPath,filesep,'repression_20%-30%_averaged''.tif']); 
saveas(gcf,[FigPath,filesep,'repression_20%-30%_averaged''.pdf']); 
%% Calculate the fold-change (repression)

% plot the fraction competent with SEM
if construct==1 || construct==4
    h(1) = errorbar(APaxis, compiledData{construct+1,16}(:,NC), compiledData{construct+1,17}(:,NC),'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))
    h(2) = errorbar(APaxis, compiledData{construct+1+8,16}, compiledData{construct+1+8,17},'LineWidth',2,'CapSize',0,'Color',(ColorChoice(construct,:)+[1 1 1])/2);
else
    h(1) = errorbar(APaxis, compiledData{construct+1,16}(:,NC), compiledData{construct+1,17}(:,NC),'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))
    h(2) = errorbar(APaxis, compiledData{construct+1+8,16}(:,NC), compiledData{construct+1+8,17}(:,NC),'LineWidth',2,'CapSize',0,'Color',(ColorChoice(construct,:)+[1 1 1])/2);
end