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
compiledData = compiledData.compiledData;
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
RuntNull = zeros(41,1);
%% If the above steps are done, then just import the output
TFData = load([FilePath, filesep, 'TFinput.mat']);
TFdata = TFData.TFinput;

Bcd = TFdata(:,1);
Run = TFdata(:,2);
RunNull = TFdata(:,3);

%% Generate Predictions using different models

%% Model-type1 : direct repression

%% For a set of parameters, calculate the fold-change

%% Fit with the real data 
% Import the real data from the compiledData
% Then, pick one dataset as our starting point.(let's start with one Run site, for example, [001])                

construct = 5; % 5th element in the DataTypes nomenclature structure.

Rate = compiledData{construct+1,9};
Rate_null = compiledData{construct+1+8,9};

APpos1 = 20;% [% of embryo length]
APpos2 = 45;% [% of embryo length]

APbin1 = APpos1/2.5 + 1;
APbin2 = APpos2/2.5 + 1;

fitRange = APbin1:APbin2;

% set the parameter bounds and initial value for the query
% params = [Kb, Kr, w_a, w_ap, w_rp, p, R_max];
lb = [0.1 0.1 1 1 0 0 50];
ub = [100 100 10 10 1 1 1000];

options = optimoptions('lsqnonlin','Display','iter', 'Algorithm', 'trust-region-reflective');

params0 = [1, 5, 2, 2, 0.2, 0.001, 200]; % example

% Fit for the Runt null
fun = @(params)model_6A1R_direct_repression_V1(Bcd(fitRange), RuntNull(fitRange),...
                        params)  - Rate_null(fitRange);
params_fit_null = lsqnonlin(fun,params0, lb, ub, options)

% Fit for the Runt WT (with more parameters)
fun = @(params)model_6A1R_direct_repression_V1(Bcd(fitRange), Runt(fitRange),...
                        params)  - Rate(fitRange);
params_fit = lsqnonlin(fun,params0, lb, ub, options)


%% Check the fitting
% Runt null
fittedRate_RuntNull = model_6A1R_direct_repression_V1(Bcd, RuntNull, params_fit_null);

hold on
plot(APaxis, Rate_null)
plot(APaxis, fittedRate_RuntNull)

% Runt WT
fittedRate_RuntWT = model_6A1R_direct_repression_V1(Bcd, Runt, params_fit);

hold on
plot(APaxis, Rate)
plot(APaxis, fittedRate_RuntWT)


%% Round2 : Recycle the parameters from the Runt null fit
%% Fit with the real data 
% Import the real data from the compiledData
% Then, pick one dataset as our starting point.(let's start with one Run site, for example, [001])                

construct = 6; % 5th element in the DataTypes nomenclature structure.

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
% params = [Kb, Kr, w_a, w_ap, w_rp, p, R_max];
lb = [0.1 0.1 1 1 0 0 50];
ub = [100 1000 10 10 1 1 1000];
% options.Algorithm = 'levenberg-marquardt';
options = optimoptions('lsqnonlin','Display','iter', 'Algorithm', 'trust-region-reflective')

params0 = [1, 1, 2, 2, 0.2, 0.001, 200]; % example

% Fit for the Runt null
fun = @(params)model_6A1R_direct_repression_V1(Bcd(fitRange), RuntNull(fitRange), params) - ...
                                Rate_null(fitRange);

% fit using the lsqnonlin
[params_fit_null,~,residual,~,~,~,jacobian] =...
                    lsqnonlin(fun,params0, lb, ub, options);

CI_null = nlparci(params_fit_null,residual,'jacobian',jacobian);

% set the parameter bounds and initial value for the query
% params = [Kb, Kr, w_a, w_ap, w_rp, p, R_max];
lb = [params_fit_null(1) 0.1 params_fit_null(3) params_fit_null(4) 0.0001 params_fit_null(6) params_fit_null(7)];
ub = [params_fit_null(1) 1000 params_fit_null(3) params_fit_null(4) 1 params_fit_null(6) params_fit_null(7)];
% options.Algorithm = 'levenberg-marquardt';

% recycle the parameters for the activation, then start with some range of
% paramter for an initial conditon.
params1 = [params_fit_null(1) 1 params_fit_null(3) params_fit_null(4) 0.2 params_fit_null(6) params_fit_null(7)];

% Fit for the Runt WT (with more parameters)
fun = @(params)model_6A1R_direct_repression_V1(Bcd(fitRange), Runt(fitRange), params1)-...
                                    Rate(fitRange);
                    
[params_fit,~,residual,~,~,~,jacobian] = ...
                    lsqnonlin(fun, params1, lb, ub, options);

CI = nlparci(params_fit,residual,'jacobian',jacobian);

%% Parameter sensitivity test (optional) - for the [001] coonstruct
% generate the predicted rate profile (over AP) based on the fitted
% parameters, for both Runt WT and Runt nulls.

% Generate the model prediction with a set of parameters fitted above.

% params_fit_direct_001 = [54.2112    5.0000    1.0000   10.2631    0.2000    0.0000  365.2318];

% Test the parameter sensitivity of Bicoid-dependent parameters
% params_fit_null = [99.9645   5.0000  1.0926   8.3891  0.2000    0.0005  373.88];
% params

% % Runt null
% Rate_null_FittedParams = model_6A1R_direct_repression_V1(Bcd, RuntNull,...
%                         params_fit_null);
%                     
% Test the parameter sensitivity of Runt-dependent parameters
% params_fit_direct_001 = params_fit;
% params_fit = [params_fit_direct_001(1), params_fit_direct_001(2),...
%                 params_fit_direct_001(3), params_fit_direct_001(4),...
%                 1, params_fit_direct_001(6),...
%                 params_fit_direct_001(7)];

% % Runt WT
% Rate_FittedParams = model_6A1R_direct_repression_V1(Bcd, Runt,...
%                         params_fit);    


% FC_fromFitting = Rate_FittedParams./Rate_null_FittedParams; 
% 
% hold on
% plot(APaxis, Rate_FittedParams)
% plot(APaxis, FC_fromFitting)
%% Calculate the fit using the fitted parameters
% Runt null
Rate_null_FittedParams = model_6A1R_direct_repression_V1(Bcd, RuntNull,...
                        params_fit_null);

% Runt WT
Rate_FittedParams = model_6A1R_direct_repression_V1(Bcd, Runt,...
                        params_fit);    

%% First, Runt nulls
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


%% Second, Runt WT with its model fit    

% Tune the parameters if needed
params_fit =[85.2969    100.0000    1.0000    7.5320    0.0500    0.0005  193.3277];

% Runt WT
Rate_FittedParams = model_6A1R_direct_repression_V1(Bcd, Runt,...
                        params_fit);    
hold on
% Runt WT
errorbar(APaxis, Rate, Rate_SEM,'o', 'Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
plot(APaxis, Rate_FittedParams, 'LineWidth',1.5,'Color',ColorChoice(1,:))


xlim([0.2 0.6])
xticks([ 0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})

box on
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

figure;
hold on
errorbar(APaxis, FC, FC_SEM,'o','Color',ColorChoice(5,:),'CapSize',0,'MarkerFaceColor',ColorChoice(5,:))
plot(APaxis, FC_model_fit, 'LineWidth',1.5,'Color',ColorChoice(5,:))
% plot(APaxis, FC_competition, 'LineWidth',1.5,'Color',ColorChoice(2,:))
% plot(APaxis, FC_quenching, 'LineWidth',1.5,'Color',ColorChoice(1,:))
% plot(APaxis, FC_direct, 'LineWidth',1.5,'Color',ColorChoice(4,:))


xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1 1.2])

xlabel('embryo length')
ylabel('fold-change')

box on
StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC_fit','.tif']); 
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC_fit','.pdf']); 

%% Part 2 : Fit the Fold-change first, then use the parameters to recover 
% the rate profiles for the Runt WT and Runt nulls.

construct = 5; % 5th element in the DataTypes nomenclature structure.

Rate = compiledData{construct+1,9};
Rate_null = compiledData{construct+1+8,9};

FC_data = Rate./Rate_null;

APpos1 = 20;% [% of embryo length]
APpos2 = 45;% [% of embryo length]

APbin1 = APpos1/2.5 + 1;
APbin2 = APpos2/2.5 + 1;

fitRange = APbin1:APbin2;

% set the parameter bounds and initial value for the query
% params =[Kb, Kr, w_a, w_ap, w_rp, p];
lb = [1 0.1 1 1 0 0];
ub = [100 100 3 3 1 1];
%options.Algorithm = 'levenberg-marquardt';

params0 = [5, 5, 1.5, 1.5, 0.5, 0.001]; % example

% Fit for the Runt null
fun = @(params)model_FC_6A1R_direct_repression_V1(Bcd(fitRange), Runt(fitRange),...
                        params)  - FC_data(fitRange);
params_fit_FC = lsqnonlin(fun, params0, lb, ub)%, options)

%% Check the fitting result
% generate the predicted rate profile (over AP) based on the fitted
% parameters, for both Runt WT and Runt nulls.

%params_fit_FC = [1.0000    0.1000    3.0000    1.0016    0.7866    0.0301];

% Runt WT
FC_FittedParams = model_FC_6A1R_direct_repression_V1(Bcd, Runt,...
                        params_fit_FC);    

hold on
% Runt WT
plot(APaxis, FC_data)
plot(APaxis, FC_FittedParams)

xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1 1.2])

xlabel('embryo length')
ylabel('fold-change')

box on
StandardFigure(gcf,gca)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part3. Fitting the Runt null and Runt WT simultaneously 

%% Fit with the real data 
% Import the real data from the compiledData
% Then, pick one dataset as our starting point.(let's start with one Run site, for example, [001])                

construct = 6; % 5th element in the DataTypes nomenclature structure.

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
% params = [Kb, Kr, w_a, w_ap, w_rp, p, R_max];
lb = [0.1 0.1 1 1 0 0 50];
ub = [100 100 10 10 1 1 1000];
% options.Algorithm = 'levenberg-marquardt';
options = optimoptions('lsqnonlin','Display','iter', 'Algorithm', 'trust-region-reflective')

params0 = [1, 5, 2, 2, 0.2, 0.001, 200]; % example

% Fit for the Runt null
fun = @(params)abs(model_6A1R_direct_repression_V1(Bcd(fitRange), RuntNull(fitRange),...
                        params)  - Rate_null(fitRange)) + ...
               abs(model_6A1R_direct_repression_V1(Bcd(fitRange), Runt(fitRange),...
                        params)  - Rate(fitRange));
% find the set of parameters that minimizes the lsq sum of two functions
% with shared parameters
params_fit = lsqnonlin(fun, params0, lb, ub, options)




%% Calculate the fit using the fitted parameters
% Runt null
Rate_null_FittedParams = model_6A1R_direct_repression_V1(Bcd, RuntNull,...
                        params_fit);
               
% Runt WT
Rate_FittedParams = model_6A1R_direct_repression_V1(Bcd, Runt,...
                        params_fit);    

%% First, Runt nulls
figure
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


%% Second, Runt WT with its model fit         
figure
hold on
% Runt WT
errorbar(APaxis, Rate, Rate_SEM,'o','CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
plot(APaxis, Rate_FittedParams, 'LineWidth',1.5,'Color',ColorChoice(1,:))


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

% calculate the fold-change, FC
FC = Rate./Rate_null;
% calculate the joint error first, using the fractional error
fracError1 = Rate_SEM./Rate;
fracError2 = Rate_null_SEM./Rate_null;
FC_SEM = sqrt(fracError1.^2 + fracError2.^2).*FC;

% FC from the model-fit
FC_model_fit = Rate_FittedParams./Rate_null_FittedParams;

hold on
errorbar(APaxis, FC, FC_SEM,'o','Color',ColorChoice(5,:),'CapSize',0,'MarkerFaceColor',ColorChoice(5,:))
plot(APaxis, FC_model_fit, 'LineWidth',1.5,'Color',ColorChoice(5,:))
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


%% Part 4. Using the nlinfit/lsqcurvefit



end