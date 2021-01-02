%% generate plots from the global, nonlinear least-squared fitting


%% Define the file paths
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo';
mkdir(FigPath)

FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';

%% Load the fitting results
% FitResult as well as which fitting method was used

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
%% Load the fitting results

%% generate plots
%% (1) Raw initial rates (WT and null) and fits with CI
APaxis = 0:0.025:1;

% counter
k=1;
for construct = [2,5,6]
    clf
    
    % initialize the variables
    Rate = [];
    Rate_SEM = [];
    Rate_null = [];
    
    % Pull the relevant data from the compiledData
    Rate = compiledData{construct+1,9};
    Rate_SEM = compiledData{construct+1,10};

    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};
    
    % pull the fitted result and error bar (Ypred and delta) from the
    % lsqcurvefit fitting.
    APrange = data_model_input(k).APbins;
    Ypred = FitResult(k).Ypred;
    delta = FitResult(k).delta;
    k=k+1;
    
    hold on
    % Runt nulls
    errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
    shadedErrorBar(APrange, Ypred(1:length(APrange)), delta(1:length(APrange)),'lineProps',{'markerfacecolor',ColorChoice(4,:)})
    % Runt WT
    errorbar(APaxis, Rate, Rate_SEM,'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
    shadedErrorBar(APrange, Ypred(length(APrange)+1:end), delta(length(APrange)+1:end),'lineProps',{'markerfacecolor',ColorChoice(1,:)})


    xlim([0.2 0.6])
    xticks([0.2 0.3 0.4 0.5 0.6])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
    ylabel({'initial RNAP', 'loading rate (AU/min)'})

    StandardFigure(gcf,gca)
    pause(0.2)
%   %Save the plot
%     saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.pdf']); 
end
%% (2) Plot the parameters
% plot the inferred parameters with error bar (which we call as STD)
% calculated from the nlparci function, for each construct.

k=1; % counter
hold on
for construct = [2,5,6]
    errorbar(FitResult(k).params_fit, FitResult(k).STD,'o','Color', ColorChoice(construct,:),'MarkerFaceColor', ColorChoice(construct,:))
    k=k+1; % counter
end

% errorbar(params_fit_100, STD_100,'o','Color', ColorChoice(2,:),'MarkerFaceColor', ColorChoice(2,:))
% errorbar(params_fit_010, STD_010,'o','Color', ColorChoice(6,:),'MarkerFaceColor', ColorChoice(6,:))
% errorbar(params_fit_001, STD_001,'o','Color', ColorChoice(5,:),'MarkerFaceColor', ColorChoice(5,:))

set(gca, 'YScale','log')
legend('100','001','010')

xlim([0 8])
xticklabels({'','K_{b}','K_{r}','\omega_{b}','\omega_{bp}','\omega_{rp}','p','R_{max}',''})
xlabel('parameters')
ylabel('fitted values')

box on
StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,'fittedParams_direct','.tif']); 
% saveas(gcf,[FigPath,filesep,'fittedParams_direct','.pdf']); 
