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
%% Generate Predictions using different models

%% Model-type1 : direct repression + competitive
% 2 interactions : Runt-Bcd (binding), Runt-RNAP (direct repression).

% Define the Runt null as zeros for the Runt protein input.
RuntNull = zeros(41,1);

% parameters 
params =[Kb, Kr, w_a, w_ap, w_ar, w_rp, p];

params1 = [5,10,5,5,1,0.2,0.001]; % example

[P_bound] = model_6A1R_direct_repression_V1(Bcd, Runt,...
                        params1) 

[P_bound_Runtnull] = model_6A1R_direct_repression_V1(Bcd, RuntNull,...
                        params1) 

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
% params =[Kb, Kr, w_a, w_ap, w_ar, w_rp, p, R_max];
lb = [0 0 1 1 0 0 0 0];
ub = [10 10 100 100 1 1 100 1000];
options.Algorithm = 'levenberg-marquardt';

params0 = [10, 5, 2, 2, 1, 0.2, 0.001, 200]; % example

% Fit for the Runt null
fun = @(params)model_6A1R_direct_repression_V1(Bcd(fitRange), RuntNull(fitRange),...
                        params)  - Rate_null(fitRange);
params_fit_null = lsqnonlin(fun,params0, lb, ub, options)

% Fit for the Runt WT (with more parameters)
fun = @(params)model_6A1R_direct_repression_V1(Bcd(fitRange), Runt(fitRange),...
                        params)  - Rate(fitRange);
params_fit = lsqnonlin(fun,params0, lb, ub, options)

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
% params =[Kb, Kr, w_a, w_ap, w_ar, w_rp, p, R_max];
lb = [0 0 1 1 1 0 0 0];
ub = [100 100 10 10 1 1 100 1000];
options.Algorithm = 'levenberg-marquardt';

params0 = [10, 20, 1.5, 1.5, 1, 0.5, 10, 200]; % example

% Fit for the Runt null
fun = @(params)model_6A1R_direct_repression_V1(Bcd(fitRange), RuntNull(fitRange),...
                        params)  - Rate_null(fitRange);
params_fit_null = lsqnonlin(fun,params0, lb, ub, options)

% set the parameter bounds and initial value for the query
% params =[Kb, Kr, w_a, w_ap, w_ar, w_rp, p, R_max];
lb = [params_fit_null(1) 0 params_fit_null(3) params_fit_null(4) 0 0 params_fit_null(7) params_fit_null(8)];
ub = [params_fit_null(1) 100 params_fit_null(3) params_fit_null(4) 1 1 params_fit_null(7) params_fit_null(8)];
options.Algorithm = 'levenberg-marquardt';

% recycle the parameters for the activation, then start with some range of
% paramter for an initial conditon.
params1 = [params_fit_null(1) 20 params_fit_null(3) params_fit_null(4) 0.5 0.5 params_fit_null(7) params_fit_null(8)];

% Fit for the Runt WT (with more parameters)
fun = @(params)model_6A1R_direct_repression_V1(Bcd(fitRange), Runt(fitRange),...
                        params1)  - Rate(fitRange);
params_fit = lsqnonlin(fun,params0, lb, ub, options)
%% Check the fitting result
% generate the predicted rate profile (over AP) based on the fitted
% parameters, for both Runt WT and Runt nulls.

% Generate the model prediction with a set of parameters fitted above.

%params_fit = [92.9275   1.0000    1.0000    8.0408    1.0000    0.2000    0.0008  327.8060];
Rate_FittedParams = model_6A1R_direct_repression_V1(Bcd, Runt,...
                        params_fit);    
% Runt null
Rate_null_FittedParams = model_6A1R_direct_repression_V1(Bcd, RuntNull,...
                        params_fit_null);

%% First, Runt WT with its model fit         
hold on
% Runt WT
errorbar(APaxis, Rate, Rate_SEM,'o','CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
plot(APaxis, Rate_FittedParams, 'LineWidth',1.5,'Color',ColorChoice(1,:))


xlim([0.1 0.7])
xticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})

StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_WT_fit','.tif']); 
% saveas(gcf,[FigPath,filesep,constructNames{construct},'_WT_fit','.pdf']); 
%% Second, Runt nulls
hold on

errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
plot(APaxis, Rate_null_FittedParams, 'LineWidth',1.5,'Color',ColorChoice(4,:))


xlim([0.1 0.7])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,constructNames{construct},'_null_fit','.tif']); 
saveas(gcf,[FigPath,filesep,constructNames{construct},'_null_fit','.pdf']); 


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

hold on
errorbar(APaxis, FC, FC_SEM,'o','Color',ColorChoice(5,:),'CapSize',0,'MarkerFaceColor',ColorChoice(5,:))
plot(APaxis, FC_model_fit, 'LineWidth',1.5,'Color',ColorChoice(5,:))


xlim([0.1 0.7])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1 1.2])

xlabel('embryo length')
ylabel('fold-change')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC_fit','.tif']); 
saveas(gcf,[FigPath,filesep,constructNames{construct},'_FC_fit','.pdf']); 

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



















%% Step2. Start with a simple model for r0 (6 Bcd binding sites)
% Note that this could be expanded to 11 or more Bcd sites.

%% Case0. hb P2, without any Runt binding sites
% free parameters : K_A, p,w_ap, r (rate when RNAP is bound to the
% promoter)
% Use rate_r0_Hill for predicting the initial rate of RNAP loading

X = 0.2:0.025:0.6;

% initial rate of RNAP loading 
initialSlopes = load([filePath, filesep, 'InitialSlopes_ONnuclei_AllConstructs.mat'])
initialSlopes = initialSlopes.initialSlopes_ONnuclei;

initialSlope_r0 = cell2mat(initialSlopes(2,5));
initialSlope_r0_SEM = cell2mat(initialSlopes(2,6));

initialSlope_r0 = initialSlope_r0(:,3); % NC14
initialSlope_r0_SEM = initialSlope_r0_SEM(:,3); % NC14

fitRange = 9:18; % 20-42.5% of AP axis as a fit range

% nonlinear least square fitting to the model
x0 = [20 0.1 10 200]; %[K_a p wAP r]
% Note that the r is constrained by the maximum of initial rate of RNAP
% loading of r0 (or the maximum value of rate from all constructs, probably
% r1-close).

lb = [0 0 0 0];
ub = [1000 100 10000 10000];
options.Algorithm = 'levenberg-marquardt';
%lsqOptions=optimset('Display','none');

fun = @(x)rate_r0_Hill(x,Bcd(fitRange)) - initialSlope_r0(fitRange);
x_r0 = lsqnonlin(fun,x0,lb, ub, options);


% x_r0 = lsqnonlin(@rate_r0_Hill, x0, lb, ub, options, Bcd(fitRange),...
%         initialSlope_r0(fitRange))%
% x = lsqnonlin(fun,x0,[],[],[]) 

%% Check the parameter sensitivity
% for two sets of parameters
Prediction1 = rate_r0_Hill(x_r0_1,Bcd)
Prediction2 = rate_r0_Hill(x_r0_2,Bcd)

hold on
errorbar(APaxis, initialSlope_r0, initialSlope_r0_SEM)
% fitted by chi2
plot(APaxis, Prediction1)
plot(APaxis, Prediction2)

% legend('data',num2str(x_r0_1),num2str(x_r0_2))
legend('data','parameter1','parameter2')
xlabel('AP axis')
ylabel('initial rate of RNAP loading (AU)')

StandardFigure(gcf,gca)
saveas(gcf,[figPath,filesep,'initialRate_chi2_r0_parameter_sensitivityTest' , '.tif']); 
saveas(gcf,[figPath,filesep,'initialRate_chi2_r0_parameter_sensitivityTest' , '.tif']); 
%% Plot the data and the fit
hold on
% data (initial fit)
errorbar(APaxis, initialSlope_r0, initialSlope_r0_SEM)
% fitted by chi2
plot(APaxis, rate_r0_Hill(x_r0,Bcd))

xlabel('AP axis')
ylabel('initial rate of RNAP loading (AU)')
StandardFigure(gcf,gca)
% save the plot
figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_hbP2_nRuntsites';
saveas(gcf,[figPath,filesep,'initialRate_chi2_r0' , '.tif']); 
saveas(gcf,[figPath,filesep,'initialRate_chi2_r0' , '.pdf']); 

%% r1

% Data extract
initialSlope_r1 = cell2mat(initialSlopes(3,5));
initialSlope_r1_SEM = cell2mat(initialSlopes(3,6));

initialSlope_r1 = initialSlope_r1(:,3); % NC14
initialSlope_r1_SEM = initialSlope_r1_SEM(:,3); % NC14

% Use the parameters fitted from the r0, then just adjust additional
% parameters.
x_r0; %[K_a p wAP r]
y0 = [5 0.01]; % K_r, wRP

lb = [5 0]; 
ub = [100 10];
options.Algorithm = 'levenberg-marquardt';
%lsqOptions=optimset('Display','none');

fun1 = @(y)rate_r1_Hill(x_r0,y,Bcd(fitRange),Runt(fitRange)) - initialSlope_r1(fitRange);
y_r1 = lsqnonlin(fun1,y0,lb, ub, options);

%% Test the parameter sensitivity, K_r
K_r = 10.^[-5 -4 -3 -2 -1 0 1 2 3 4 5];

hold on
errorbar(APaxis, initialSlope_r1, initialSlope_r1_SEM)
for i=1:length(K_r)
    Prediction = rate_r1_Hill(x_r0, [K_r(i), y_r1(2)], Bcd, Runt)
    plot(APaxis, Prediction)
end
%% Plot to check r1 and r1-from chi2
Prediction = rate_r1_Hill(x_r0, y_r1, Bcd, Runt)

hold on 
errorbar(APaxis, initialSlope_r1, initialSlope_r1_SEM)
plot(APaxis, Prediction)

xlabel('AP axis')
ylabel('initial rate of RNAP loading (AU)')
StandardFigure(gcf,gca)
% save the plot
figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_hbP2_nRuntsites';
saveas(gcf,[figPath,filesep,'initialRate_chi2_r1' , '.tif']); 
saveas(gcf,[figPath,filesep,'initialRate_chi2_r1' , '.pdf']); 
%% r1-close
% Data extract
initialSlope_r1_close = cell2mat(initialSlopes(6,5));
initialSlope_r1_close_SEM = cell2mat(initialSlopes(6,6));

initialSlope_r1_close = initialSlope_r1_close(:,3); % NC14
initialSlope_r1_close_SEM = initialSlope_r1_close_SEM(:,3); % NC14

% Use the parameters fitted from the r0, then just adjust additional
% parameters.
x_r0; %[K_a p wAP r]
y0 = [5 0.01];% K_r, wRP

% constrain the K_r as obtained in r1 case
lb = [y_r1(1) 0]; 
ub = [y_r1(1) 10];
% lb = [1 0]; 
% ub = [100 10];
options.Algorithm = 'levenberg-marquardt';
%lsqOptions=optimset('Display','none');

fun1 = @(y)rate_r1_Hill(x_r0,y,Bcd(fitRange),Runt(fitRange)) - initialSlope_r1_close(fitRange);
y_r1_close = lsqnonlin(fun1,y0,lb, ub, options);

%% Plot to check r1-close and r1-close-from chi2
Prediction = rate_r1_Hill(x_r0, y_r1_close, Bcd, Runt)

hold on 
errorbar(APaxis, initialSlope_r1_close, initialSlope_r1_close_SEM)
plot(APaxis, Prediction)

xlabel('AP axis')
ylabel('initial rate of RNAP loading (AU)')
StandardFigure(gcf,gca)
% save the plot
figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_hbP2_nRuntsites';
saveas(gcf,[figPath,filesep,'initialRate_chi2_r1_close' , '.tif']); 
saveas(gcf,[figPath,filesep,'initialRate_chi2_r1_close' , '.pdf']);

%% Predict the 2 binding sites ([1,0,0] + [0,0,1])
Prediction_r2_far = rate_r2_Hill_multiRep(x0,[y_r1(1), y_r1(2), y_r1_close(2)],Bcd, Runt);

% actual r2([1,1,0]) initial slope
initialSlope_r2_far = cell2mat(initialSlopes(9,5));
initialSlope_r2_far_SEM = cell2mat(initialSlopes(9,6));

initialSlope_r2_far = initialSlope_r2_far(:,3); % NC14
initialSlope_r2_far_SEM = initialSlope_r2_far_SEM(:,3); % NC14

%% plot r2 prediciton and data
hold on
plot(APaxis, Prediction_r2_far)
errorbar(APaxis, initialSlope_r2_far, initialSlope_r2_far_SEM)
end