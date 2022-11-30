% generate plots of initial slopes, fold-change, etc.
% Name : make_InitialSlope_fig

%clear 
%close all
%addpath('../utilities')
% set ID variables
DropboxFolder = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\mat files';
FigureRoot = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput';

% load data structure
load([DropboxFolder,filesep,'AccumulatedData.mat'])

FigPath = [FigureRoot, filesep, 'accumulatedmRNA'];
%mkdir(FigPath)

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

%% Accumulated mRNA over AP axis for each construct with/without Runt protein
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_mRNA = figure;

for construct=1:8 % total number of constructs (enhancers)
    clf
    hold on
    errorbar(APaxis, AccumulatedData{construct+1,8}, AccumulatedData{construct+1,9},'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))
    errorbar(APaxis, AccumulatedData{construct+1+8,8}, AccumulatedData{construct+1+8,9},'LineWidth',2,'CapSize',0,'Color',(ColorChoice(construct,:)+[1 1 1])/2);

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 2.5*10^5])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 1 2 3]*10^5)

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('accumulated mRNA (AU)')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_mRNA, fig_mRNA.CurrentAxes)
    pause
    % Save the plot
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

%% initial slope vs Accumulated mRNA (Runt WT)
% For all 8 constructs, we will plot the initial rate vs accumulated mRNA
% for (+/-) Runt in separate panels.

APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
%fig_mRNA = figure;


%     % Runt WT
%     errorbar(compiledData{construct+1,9}, AccumulatedData{construct+1,8},...
%                 AccumulatedData{construct+1,9},...
%                 AccumulatedData{construct+1,9},...
%                 compiledData{construct+1,10},...
%                 compiledData{construct+1,10},...
%                 'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))
%     % Runt null
%     errorbar(compiledData{construct+1+8,9}, AccumulatedData{construct+1+8,8},...
%                 AccumulatedData{construct+1+8,9},...
%                 AccumulatedData{construct+1+8,9},...
%                 compiledData{construct+1+8,10},...
%                 compiledData{construct+1+8,10},...
%                 'LineWidth',2,'CapSize',0,'Color',(ColorChoice(construct,:)+[1 1 1])/2)
for construct=1:8 % total number of constructs (enhancers)
    clf
    hold on
    
    scatter(compiledData{construct+1,9}, AccumulatedData{construct+1,8},...
            'LineWidth',2,'Color',ColorChoice(construct,:))

    scatter(compiledData{construct+1+8,9}, AccumulatedData{construct+1+8,8},...
            'LineWidth',2,'Color',ColorChoice(construct,:))

    % xTicks, yTicks
    xlim([0 400])
    ylim([0 3*10^5])
    xticks([0 100 200 300 400])
    yticks([0 1 2 3]*10^5)

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('initial rate of RNAP loading (AU/min)')
    ylabel('accumulated mRNA (AU)')

    legend(constructNames{construct},'Location','NorthEast')

    box on

    StandardFigure(gcf, gca)
    pause
    % Save the plot
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

%% initial slope vs accumulated mRNA (Runt WT and Runt nulls)

FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\accumulatedmRNA\accmRNA_vs_initial_slope';

xvec = [];
yvec = [];

hold on
for construct=1:8 % total number of constructs (enhancers)
    %clf
    %hold on
    scatter(compiledData{construct+1,9}, AccumulatedData{construct+1,8},...
            70,ColorChoice(construct,:),"filled")

    scatter(compiledData{construct+1+8,9}, AccumulatedData{construct+1+8,8},...
            70,(ColorChoice(construct,:)+[1 1 1])/2,"filled")
    
    % compute Pearson corr. coeff
    xvec = [xvec; compiledData{construct+1,9};compiledData{construct+1+8,9}];
    yvec = [yvec; AccumulatedData{construct+1,8}'; AccumulatedData{construct+1+8,8}'];
    
  

end
  [RHO,PVAL] = corr(xvec, yvec,'rows','complete');
    str=sprintf('r= %1.2f',RHO);
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');

    % xTicks, yTicks
    xlim([0 400])
    ylim([0 3*10^5])
    xticks([0 100 200 300 400])
    yticks([0 1 2 3]*10^5)

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('initial rate of RNAP loading (AU/min)')
    ylabel('accumulated mRNA (AU)')

    %legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(gcf, gca)
    %pause
    % Save the plot
%    saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
%      saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
    saveas(gcf,[FigPath, filesep, 'all_constructs.pdf'])
%% Accumulated mRNA (from individual embryo) with average
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_mRNA = figure;

for construct=1:16 % total number of constructs (enhancers)
    clf
    hold on
    errorbar(APaxis, AccumulatedData{construct+1,8}, AccumulatedData{construct+1,9},'k','LineWidth',2)%,'Color',ColorChoice(construct,:))
    % individual embryos
    [~,~,numEmbryos] = size(AccumulatedData{construct+1,6});
    NC14 = AccumulatedData{construct+1,5};
    tEnd = min(NC14 + 20, length(AccumulatedData{construct+1,6}(:,1,1)));
    
    accumulatedmRNA_temp = AccumulatedData{construct+1,6};
    accumulatedmRNA_SD_temp = AccumulatedData{construct+1,7};
    
    if construct>9
        for i=1:numEmbryos
            errorbar(APaxis,accumulatedmRNA_temp(tEnd,:,i) , accumulatedmRNA_SD_temp(tEnd,:,i))
        end
    else
        for i=1:numEmbryos
            errorbar(APaxis,accumulatedmRNA_temp(tEnd,:,i) - accumulatedmRNA_temp(NC14,:,i) ,...
                        sqrt(accumulatedmRNA_SD_temp(tEnd,:,i).^2 +accumulatedmRNA_SD_temp(NC14,:,i).^2 ))
        end
    end

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 3*10^5])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 1 2 3]*10^5)

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('accumulated mRNA (AU)')

    legend(constructNames{construct},'Location','NorthEast')

    box on

    StandardFigure(fig_mRNA, fig_mRNA.CurrentAxes)
    pause(1)
    % Save the plot
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end
%% Accumulated mRNA over AP axis for each construct with/without Runt protein
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_mRNA = figure;

for construct=1%:8 % total number of constructs (enhancers)
    %clf
    hold on
    errorbar(APaxis, AccumulatedData{construct+1,10}, AccumulatedData{construct+1,11},'LineWidth',2,'Color',ColorChoice(construct,:))
    errorbar(APaxis, AccumulatedData{construct+1+8,10}, AccumulatedData{construct+1+8,11},'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2);

    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 3*10^5])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 1 2 3]*10^5)

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('accumulated mRNA (AU)')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_mRNA, fig_mRNA.CurrentAxes)
    pause(1)
    % Save the plot
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

%% fold-change of accumulated mRNA
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_mRNA = figure;


for construct=1:8 % total number of constructs (enhancers)
    clf
    FC = AccumulatedData{construct+1,8}./AccumulatedData{construct+1+8,8};
    fracError1 = AccumulatedData{construct+1,9}./ AccumulatedData{construct+1,8};
    fracError2 = AccumulatedData{construct+1+8,9}./ AccumulatedData{construct+1+8,8};
    
    FC_SEM = sqrt(fracError1.^2 + fracError2.^2).*FC;
    errorbar(APaxis, FC, FC_SEM,'LineWidth',2,'Color',ColorChoice(construct,:))


    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 1.4])
    xticks([0.2 0.3 0.4 0.5 0.6])
    yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('embryo length')
    ylabel('fold-change')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_mRNA, fig_mRNA.CurrentAxes)
    pause(1)
    % Save the plot
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC','.tif']); 
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC','.pdf']); 
    % Optional
    % Calculate the spatially averaged fold-change(repression)
    FC_spatial_mean_mRNA(construct) = nanmean(FC(9:13));
    FC_spatial_error_mRNA(construct) = sqrt(sum(FC_SEM(9:13).^2))/length(9:13);
end



%% fold-change(repression) of accumulated mRNA
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_mRNA = figure;


for construct=1:8 % total number of constructs (enhancers)
    clf
    FC = AccumulatedData{construct+1+8,8}./AccumulatedData{construct+1,8};
    fracError1 = AccumulatedData{construct+1,9}./ AccumulatedData{construct+1,8};
    fracError2 = AccumulatedData{construct+1+8,9}./ AccumulatedData{construct+1+8,8};
    
    FC_SEM = sqrt(fracError1.^2 + fracError2.^2).*FC;
    errorbar(APaxis, FC, FC_SEM,'LineWidth',2,'CapSize',0,'Color',ColorChoice(construct,:))


    % xTicks, yTicks
    xlim([0.15 0.6])
    ylim([0 30])
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
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_repression','.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'_repression','.pdf']); 
    
    % Optional
    % Calculate the spatially averaged fold-change(repression)
    FC_spatial_mean_mRNA(construct) = nanmean(FC(9:13));
    FC_spatial_error_mRNA(construct) = sqrt(sum(FC_SEM(9:13).^2))/length(9:13);
end
%% plot the spatially averaged fold-change for each construct

construct = [1,2,5,6,3,7,8,4];
order = 1:8;
errorbar(order, FC_spatial_mean_mRNA(construct), FC_spatial_error_mRNA(construct),'o')
yline(1,'--')
xticklabels({'000','100','001','010','011','110','101','111'})

xlabel('construct')
ylabel('repression')

StandardFigure(gcf,gca)

saveas(gcf,[FigPath,filesep,'repression_20%-30%_averaged''.tif']); 
saveas(gcf,[FigPath,filesep,'repression_20%-30%_averaged''.pdf']); 

%% Added for eLife revision_v1
%% Accumulated mRNA over Runt concentration for each construct with/without Runt protein
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_mRNA = figure;
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\accumulatedmRNA\accumulatedmRNA_vs_RuntConc'

for construct=1:8 % total number of constructs (enhancers)
    clf
    hold on
    errorbar(Runt(9:end), AccumulatedData{construct+1,10}(9:end), AccumulatedData{construct+1,11}(9:end),'LineWidth',2,'Color',ColorChoice(construct,:))
    errorbar(RuntNull(9:end), AccumulatedData{construct+1+8,10}(9:end), AccumulatedData{construct+1+8,11}(9:end),'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2);

    % xTicks, yTicks
    xlim([0 500])
    ylim([0 3*10^5])
    xticks([0 100 200 300 400 500])
    yticks([0 1 2 3]*10^5)

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('Runt concentration (AU)')
    ylabel('accumulated mRNA (AU)')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_mRNA, fig_mRNA.CurrentAxes)
    % pause
    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end