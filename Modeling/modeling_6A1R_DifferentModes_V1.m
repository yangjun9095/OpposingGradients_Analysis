function modeling_6A1R_DifferentModes_V1
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

BcdData = load([FilePath, filesep, 'Bcd_NC14_TimeAveraged.mat']);

RuntData = load([FilePath, filesep, 'Runt_TimeAveraged_mixedSex_NC14.mat']);

compiledData = load([FilePath, filesep, 'compiledData.mat']);

%% Extract the useful fields from BcdData and RuntData
% Bcd
BcdFluo_tAveraged = BcdData.Bcd_timeAveraged_10min_nc14;
SDBcdFluo_tAveraged = BcdData.Bcd_timeAveraged_10min_nc14_SD;
% Runt
RuntFluo_tAveraged = RuntData.AveragedFluo_tAveraged_mixed(2,:); % 0-10 mins into NC14
SERuntFluo_tAveraged = RuntData.SEFluo_tAveraged_mixed(2,:);% 0-10 mins into NC14

%% Define Bcd and Runt using these time-averaged vectors
Bcd = BcdFluo_tAveraged;
Runt = RuntFluo_tAveraged;

% Transpose to plug in as inputs
Bcd = Bcd';
Runt = Runt';

% Define the Runt null as zeros for the Runt protein input.
RuntNull = zeros(41,1);
%% Generate Predictions using different models

%% Model-type : Make this as an input function such that we can pick different model as we want

%% Round2 : Recycle the parameters from the Runt null fit
%% Fit with the real data 
% Import the real data from the compiledData
% Then, pick one dataset as our starting point.(let's start with one Run site, for example, [001])                

construct = 5; % N-th element in the DataTypes nomenclature structure.

Rate = compiledData{construct+1,9};
Rate_SEM = compiledData{construct+1,10};

Rate_null = compiledData{construct+1+8,9};
Rate_null_SEM = compiledData{construct+1+8,10};

APpos1 = 20;% [% of embryo length]
APpos2 = 45;% [% of embryo length]

APbin1 = APpos1/2.5 + 1;
APbin2 = APpos2/2.5 + 1;

fitRange = APbin1:APbin2;

% set the parameter bounds and initial value for the query
% params =[Kb, Kr, w_b, w_bp, w_br, w_rp, w_brp, p, R_max];
lb = [0.1 0.1 1 1 0 0 0 0 0];
ub = [100 100 10 10 1 1 10 1 1000];
% options.Algorithm = 'levenberg-marquardt';
options = optimoptions('lsqnonlin','Display','iter', 'Algorithm', 'trust-region-reflective');

params0 = [10   1    1.0000    10    0.8    1   0.51    0.1  181.1421]; % example

% Fit for the Runt null
fun = @(params)model_6A1R_combination_all_V1(Bcd(fitRange), RuntNull(fitRange),...
                        params)  - Rate_null(fitRange);
params_fit_null = lsqnonlin(fun, params0, lb, ub, options)

% set the parameter bounds and initial value for the query
% params =[Kb, Kr, w_b, w_bp, w_br, w_rp, w_brp, p, R_max];
lb = [params_fit_null(1) 0.1 params_fit_null(3) params_fit_null(4) 0 0 0 params_fit_null(8) params_fit_null(9)];
ub = [params_fit_null(1) 100 params_fit_null(3) params_fit_null(4) 10 1 1 params_fit_null(8) params_fit_null(9)];

% recycle the parameters for the activation, then start with some range of
% paramter for an initial conditon.
params1 = [params_fit_null(1) 10 params_fit_null(3) params_fit_null(4) 1 0.5 1 params_fit_null(8) params_fit_null(9)];

% Fit for the Runt WT (with more parameters)
fun = @(params)model_6A1R_combination_all_V1(Bcd(fitRange), Runt(fitRange),...
                        params1)  - Rate(fitRange);
params_fit = lsqnonlin(fun, params1, lb, ub, options)

%% Check the fitting result
% generate the predicted rate profile (over AP) based on the fitted
% parameters, for both Runt WT and Runt nulls.

% Generate the model prediction with a set of parameters fitted above.

% params_fit = [100   0.1    1.0000    10    3.5    1.0000    0.51    0.0002  181.1421];
params_fit = [98.1263   1.0000    1.0000    9.9770    0.9    0.5000    0.95    0.0002  366.1176];

Rate_FittedParams = model_6A1R_combination_all_V1(Bcd, Runt,...
                        params_fit);    
% Runt null
Rate_null_FittedParams = model_6A1R_combination_all_V1(Bcd, RuntNull,...
                        params_fit);

                    
% First, Runt nulls
hold on

errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
plot(APaxis, Rate_null_FittedParams, 'LineWidth',1.5,'Color',ColorChoice(4,:))


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

% Second, Runt WT with its model fit         
% figure
hold on
% Runt WT
errorbar(APaxis, Rate, Rate_SEM,'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
plot(APaxis, Rate_FittedParams, 'LineWidth',1.5,'Color',ColorChoice(1,:))


xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
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
FC_model_fit = Rate_FittedParams./Rate_null_FittedParams;

figure
hold on
errorbar(APaxis, FC, FC_SEM,'o','Color',ColorChoice(5,:),'CapSize',0,'MarkerFaceColor',ColorChoice(5,:))
plot(APaxis, FC_model_fit, 'LineWidth',1.5,'Color',ColorChoice(5,:))


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

%% Fit the Fold-change

construct = 5; % 5th element in the DataTypes nomenclature structure.

Rate = compiledData{construct+1,9};
Rate_null = compiledData{construct+1+8,9};

FC_data = Rate./Rate_null;

APpos1 = 20;% [% of embryo length]
APpos2 = 30;% [% of embryo length]

APbin1 = APpos1/2.5 + 1;
APbin2 = APpos2/2.5 + 1;

fitRange = APbin1:APbin2;

% set the parameter bounds and initial value for the query
% params =[Kb, Kr, w_a, w_ap, w_ar, w_rp, p];
lb = [0 0 1 1 0 0 0];
ub = [100 1000 100 100 1 1 100];
%options.Algorithm = 'levenberg-marquardt';

params0 = [20, 20, 2, 2, 0.5, 0.2, 0.001]; % example

% Fit for the Runt null
fun = @(params)model_FC_6A1R_direct_repression_V1(Bcd(fitRange), Runt(fitRange),...
                        params)  - FC_data(fitRange);
params_fit_FC = lsqnonlin(fun, params0, lb, ub)%, options)




%% Check the fitting result
% generate the predicted rate profile (over AP) based on the fitted
% parameters, for both Runt WT and Runt nulls.

% Runt WT
FC_FittedParams = model_FC_6A1R_direct_repression_V1(Bcd, Runt,...
                        params_fit_FC);    

hold on
% Runt WT
plot(APaxis, FC_data)
plot(APaxis, FC_FittedParams)







end