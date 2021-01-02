function modeling_6A1R_DirectRepression_V3(varargin)
% hbP2 + N Runt binding sites constructs modeling with actual Bcd and Runt
% gradient.
% Last updated : 10/4/2020, by YJK

% Assumptions : 
% 1) Bicoid binds cooperatively, and could be desribed using Hill function.
% 2) Runt is not known for its cooperativity, thus we will assume a direct
% repression (where Runt only interacts with RNAP machinery), but not with
% activator binding/activation.

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

%% Step0. Set up the directories
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo';
mkdir(FigPath)
%% Step1. Get the actual input TF and output data
% Here, I'll use the time-averaged Bcd and Runt profiles processed by
% main01_04/05 scripts.
% Actually, the time-averaging "time window" doesn't matter that much,
% since the gradient scales nicely over time. But, we use the 0-10 min into
% nc14.

FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';

Data = load([FilePath, filesep, 'compiledData.mat']);

compiledData = Data.compiledData;

DataTypes = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};

% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};

%% If the above steps are done, then just import the output
TFData = load([FilePath, filesep, 'TFinput.mat']);
TFdata = TFData.TFinput;

Bicoid = TFdata(:,1);
Runt = TFdata(:,2);
RuntNull = TFdata(:,3);

% Make a matrix whose each column is each TF
TF(:,1) = Bicoid;
TF(:,2) = Runt;
TF(:,3) = RuntNull;
%% Generate Predictions using different models

%% Model-type1 : direct repression

%% Fit with the real data 
% Import the real data from the compiledData
% Then, pick one dataset as our starting point.(let's start with one Run site, for example, [001])                

% initialize the variables
data = [];
input_combined =[];

construct = 5; % N-th element in the DataTypes nomenclature structure.

Rate = compiledData{construct+1,9};
Rate_SEM = compiledData{construct+1,10};

Rate_null = compiledData{construct+1+8,9};
Rate_null_SEM = compiledData{construct+1+8,10};

% Define the range of fitting
APpos1 = 20;% [% of embryo length]
APpos2 = 45;% [% of embryo length]

APbin1 = APpos1/2.5 + 1;
APbin2 = APpos2/2.5 + 1;

fitRange = APbin1:APbin2;

% Trim the data of WT and Null into one vector
data = [Rate_null(fitRange); Rate(fitRange)];

input_combined(:,1) = [Bicoid(fitRange); Bicoid(fitRange)];
input_combined(:,2) = [RuntNull(fitRange); Runt(fitRange)];




% Pick a model from different modes
mdl0 = @model_6A1R_direct_repression_V2;
mdl = @(params, TF) model_6A1R_direct_repression_V2(params, TF);

% Set the parameter bounds and initial value for the query
% params = [Kb, Kr, w_a, w_ap, w_rp, p, R_max];
lb = [0.1 0.1 1 1 0 0 50];
ub = [1000 100 100 100 1 1 1000];
% options.Algorithm = 'levenberg-marquardt';

% options = optimoptions('lsqnonlin','Display','iter', 'Algorithm', 'trust-region-reflective')

optimoptions = optimset('TolFun',1E-6, 'MaxIter', 1E6, 'MaxFunEvals', 1E5);
                   

params0 = [1, 5, 2, 2, 0.2, 0.001, 200]; % example

% Fit for the Runt null
% fun = @(params)model_6A1R_direct_repression_V2(params,TF(fitRange,:)) - Rate_null(fitRange);

% fit using the lsqcurvefit
[params_fit,~,Res,~,~,~,Jacobian] =...
                    lsqcurvefit(mdl0, params0,input_combined, data, lb, ub, optimoptions);

%% Calculate the CI of the model fit and parameters
% First, get the CI for the fitted parameters
CI = nlparci(params_fit, Res, 'jacobian', Jacobian);

% Second, calculate the CI for the predicted fit using nlpredci
[Ypred,delta] = nlpredci(mdl,input_combined, params_fit, Res, 'Jacobian', full(Jacobian));

STD = (CI(:,2) - CI(:,1)) /2;
%% Plot for checking
%% First, Runt nulls
APaxis = 0:0.025:1;

hold on
% Runt nulls
errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
shadedErrorBar(APaxis(fitRange), Ypred(1:length(fitRange)), delta(1:length(fitRange)),'lineProps',{'markerfacecolor',ColorChoice(4,:)})
% Runt WT
errorbar(APaxis, Rate, Rate_SEM,'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
shadedErrorBar(APaxis(fitRange), Ypred(length(fitRange)+1:end), delta(length(fitRange)+1:end),'lineProps',{'markerfacecolor',ColorChoice(1,:)})


xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})

StandardFigure(gcf,gca)

% Save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo\lsqcurvefit_V3\Direct';
% saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.tif']); 
% saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.pdf']); 

%% Plot the parameters

% Save the variables into construct specific variables
% params_fit_100 = params_fit;

hold on
plot(params_fit_100, 'o')
plot(params_fit_010, 'o')
plot(params_fit_001, 'o')

set(gca, 'YScale','log')

%% Covariance of parameters
% CovB = inv(J'*J).*MSE;
% covfig = figure;
% cv = @(x, y) sqrt(abs(x)) ./ sqrt((y'*y));
% imagesc(cv(CovB, b));
% colorbar;
% ylabel('parameter 1')
% xlabel('parameter 2')
% title('Covariance matrix of the parameters');

%% First, Runt nulls
hold on

errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
shadedErrorBar(APaxis, Ypred_null, delta_null,'lineProps',{'markerfacecolor',ColorChoice(4,:)})

xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})

StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_null_fit','.tif']); 
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_null_fit','.pdf']); 


%% Second, Runt WT with its model fit         
hold on
% Runt WT
errorbar(APaxis, Rate, Rate_SEM,'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
shadedErrorBar(APaxis, Ypred, delta,'lineProps',{'color',ColorChoice(1,:),'markerfacecolor',ColorChoice(1,:)})


xlim([0.2 0.6])
xticks([ 0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})

StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_WT_fit','.tif']); 
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_WT_fit','.pdf']); 


%% Third, FC from individual fittings
% Rate = compiledData{construct+1,9};
% Rate_SEM = compiledData{construct+1,10};
% 
% Rate_null = compiledData{construct+1+8,9};
% Rate_null_SEM = compiledData{construct+1+8,10};

% calculate the fold-change, FC
FC = Rate./Rate_null;
% calculate the joint error first, using the fractional error
fracError1 = Rate_SEM./Rate;
fracError2 = Rate_null_SEM./Rate_null;
FC_SEM = sqrt(fracError1.^2 + fracError2.^2).*FC;

% FC from the model-fit
FC_model_fit = Ypred./Ypred_null;

fracError1_model = delta./Ypred;
fracError2_model = delta_null./Ypred_null;

FC_model_fit_error = sqrt(fracError1_model.^2 + fracError2_model.^2).*FC_model_fit;

hold on
errorbar(APaxis, FC, FC_SEM,'o','Color',ColorChoice(5,:),'CapSize',0,'Color',ColorChoice(5,:),'MarkerFaceColor',ColorChoice(5,:))
shadedErrorBar(APaxis, FC_model_fit, FC_model_fit_error, 'lineProps',{'color',ColorChoice(5,:),'markerfacecolor',ColorChoice(5,:)})
% plot(APaxis, FC_competition, 'LineWidth',1.5,'Color',ColorChoice(2,:))
% plot(APaxis, FC_quenching, 'LineWidth',1.5,'Color',ColorChoice(1,:))
% plot(APaxis, FC_direct, 'LineWidth',1.5,'Color',ColorChoice(4,:))




xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1 1.2])

xlabel('embryo length')
ylabel('fold-change')

StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC_fit','.tif']); 
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC_fit','.pdf']); 



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end