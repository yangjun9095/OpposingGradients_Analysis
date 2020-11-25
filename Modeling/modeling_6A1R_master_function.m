function modeling_6A1R_master_function(varargin)
% Description : 
% 
% Last updated : 11/25/2020, by YJK

% Assumptions : 

%% Step0. Set up the directories

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

%% If the data pre-processing steps are done, then just import the output
TFData = load([FilePath, filesep, 'TFinput.mat']);
TFdata = TFData.TFinput;

Bcd = TFdata(:,1);
Run = TFdata(:,2);
RunNull = TFdata(:,3);

% Make a matrix whose each column is each TF
TF(:,1) = Bcd;
TF(:,2) = Run;
TF(:,3) = RunNull;

%% Model-type1 : direct repression

%% Pre-process the data (both input TF and output rate(slope)) such that we can fit to a given model.

% Import the real data from the compiledData
% Then, pick one dataset as our starting point.(let's start with one Run site, for example, [001])                
% Define the AP range
APaxis = 0:0.025:1;

% Pick the construct for fitting
% N-th element in the DataTypes nomenclature structure.
data_model_input = struct;
k=1; %counter
for construct = [2,5,6]
    % initialize the variables
    data = [];
    input_combined =[];
    data_interp = [];
    input_combined_interp = [];

    Rate = compiledData{construct+1,9};
    Rate_SEM = compiledData{construct+1,10};

    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};

    % Define the range of fitting
    APaxis = 0:0.025:1;
    
    if construct == 2
        APpos1 = 20;% [% of embryo length]
        APpos2 = 42.5;% [% of embryo length]
    else
        APpos1 = 20;% [% of embryo length]
        APpos2 = 45;% [% of embryo length]
    end

    APbin1 = APpos1/2.5 + 1;
    APbin2 = APpos2/2.5 + 1;

    fitRange = APbin1:APbin2;

    % Trim the data of WT and Null into one vector
    data = [Rate_null(fitRange); Rate(fitRange)];

    input_combined(:,1) = [Bcd(fitRange); Bcd(fitRange)];
    input_combined(:,2) = [RunNull(fitRange); Run(fitRange)];
    

    % save into a structure
    data_model_input(k).construtName = DataTypes(construct)
    data_model_input(k).APbins = APaxis(fitRange);
    data_model_input(k).input = input_combined;
    data_model_input(k).output = data;
    k=k+1;

end   

%% (Optional) : Interpolation that could be done in the above for loop)
%     % interplate to see if this makes our inference of parameters better
%     AP_interp = APaxis(fitRange(1)):0.025/10:APaxis(fitRange(end));
%     AP_interp = AP_interp'; % transpose
%     APbins_interp = [AP_interp ; AP_interp];
%     
%     % interpolate the input TF data
%     Bcd_interp = interp1(APaxis(fitRange), Bcd(fitRange), AP_interp);
%     RunNull_interp= interp1(APaxis(fitRange), RunNull(fitRange), AP_interp);
%     Run_interp = interp1(APaxis(fitRange), Run(fitRange), AP_interp);
%     
%     % combine those interpolated inputs into a matrix
%     input_combined_interp(:,1) = [Bcd_interp ; Bcd_interp];
%     input_combined_interp(:,2) = [RunNull_interp ; Run_interp];
%     
%     % interpolate the output data
%     Rate_null_interp = interp1(APaxis(fitRange), Rate_null(fitRange), AP_interp);
%     Rate_interp = interp1(APaxis(fitRange), Rate(fitRange), AP_interp);
%     data_interp = [Rate_null_interp ; Rate_interp];


%     data_model_input(construct).APbins =  AP_interp;
%     data_model_input(construct).input = input_combined_interp;
%     data_model_input(construct).output = data_interp;
    
%% Fitting process using the global_fit_construct_generalThermoV2.m
% Inputs : xdata, ydata, and model for the fitting

% Pick a model from different modes
% Define anonymous function from the given model input. (either direct,
% competition, quenching, etc.)
mdl0 = @model_6A1R_direct_repression_V2;
mdl = @(params, TF) model_6A1R_direct_repression_V2(params, TF);

% Define the parameter bounds for the fitting
% params = [Kb, Kr, w_a, w_ap, w_rp, p, R_max];
lb = [10^(-2) 10^(-2) 1 1 0 10^(-4) 50];
ub = [10^5 10^5 10^2 10^2 1 10 1000];

% options = optimoptions('lsqnonlin','Display','iter', 'Algorithm', 'trust-region-reflective');

% initial input of parameter guess
params0 = [1, 5, 2, 2, 0.2, 0.001, 300]; % example

% Initialize the FitResult structure
FitResult = struct;
k=1; % counter

for construct = [2,5,6]
    % initialize the variables
    params_fit = [];
    Res = [];
    Jacobian = [];
    CI = [];
    STD = [];
    Ypred = [];
    delta = [];
    
    % define the inputs for the fitting
    inputTF = data_model_input(k).input;
    data = data_model_input(k).output;
    
    [params_fit, Res, Jacobian, CI, STD, Ypred, delta] = ...
                global_fit_construct_generalThermoV2(inputTF, data, mdl0, mdl, params0, lb, ub);
            
    % save the fit results (parameters and Jacobian, etc.) into the FitResult structure.        
    FitResult(k).constructName = DataTypes(construct);
    FitResult(k).params_fit = params_fit;
    FitResult(k).Res = Res;
    FitResult(k).Jacobian = Jacobian;
    FitResult(k).CI = CI;
    FitResult(k).STD = STD;
    FitResult(k).Ypred = Ypred;
    FitResult(k).delta = delta;

    k=k+1;
end

%% Save the fitting results into structures
save('S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\generalThermoFit_direct_V3.mat',...
        'data_model_input', 'FitResult')
%% Covariance of parameters
% CovB = inv(J'*J).*MSE;
% covfig = figure;
% cv = @(x, y) sqrt(abs(x)) ./ sqrt((y'*y));
% imagesc(cv(CovB, b));
% colorbar;
% ylabel('parameter 1')
% xlabel('parameter 2')
% title('Covariance matrix of the parameters');


%% Part2. Competition
%% Fitting process using the global_fit_construct_generalThermoV2.m

% Pick a model from different modes
% Define anonymous function from the given model input. (either direct,
% competition, quenching, etc.)
mdl0 = @model_6A1R_competition_V2;
mdl = @(params, TF) model_6A1R_competition_V2(params, TF);

% Define the parameter bounds for the fitting
% params = [Kb, Kr, w_a, w_ap, w_br, p, R_max];
lb = [10^(-2) 10^(-2) 1 1 0 10^(-4) 50];
ub = [10^3 10^3 10^2 10^2 1 10 1000];

% options = optimoptions('lsqnonlin','Display','iter', 'Algorithm', 'trust-region-reflective');

% initial input of parameter guess
params0 = [1, 5, 2, 2, 0.2, 0.001, 300]; % example

% Initialize the FitResult structure
FitResult = struct;
k=1; % counter

for construct = [2,5,6]
    % initialize the variables
    params_fit = [];
    Res = [];
    Jacobian = [];
    CI = [];
    STD = [];
    Ypred = [];
    delta = [];
    
    % define the inputs for the fitting
    inputTF = data_model_input(k).input;
    data = data_model_input(k).output;
    
    [params_fit, Res, Jacobian, CI, STD, Ypred, delta] = ...
                global_fit_construct_generalThermoV2(inputTF, data, mdl0, mdl, params0, lb, ub);
            
    % save the fit results (parameters and Jacobian, etc.) into the FitResult structure.        
    FitResult(k).constructName = DataTypes(construct);
    FitResult(k).params_fit = params_fit;
    FitResult(k).Res = Res;
    FitResult(k).Jacobian = Jacobian;
    FitResult(k).CI = CI;
    FitResult(k).STD = STD;
    FitResult(k).Ypred = Ypred;
    FitResult(k).delta = delta;

    k=k+1;
end

%% plot the fitting results
APaxis = 0:0.025:1;

% counter
k=1;
for construct = [2,5,6]
    clf
    
    % initialize the variables
    Rate = [];
    Rate_SEM = [];
    Rate_null = [];
    Rate_null_SEM = [];
    
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

%   %Save the plot
    saveas(gcf,[FigPath,filesep,'fit_competition_V3_',constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,'fit_competition_V3_',constructNames{construct},'.pdf']); 
end
%% Plot the parameters

% initialize the indexing
k=1;

hold on
for construct = [2,5,6]
    errorbar(FitResult(k).params_fit, FitResult(k).STD,'o','Color', ColorChoice(construct,:),'MarkerFaceColor', ColorChoice(construct,:))
    k=k+1;
end

% errorbar(params_fit_100, STD_100,'o','Color', ColorChoice(2,:),'MarkerFaceColor', ColorChoice(2,:))
% errorbar(params_fit_010, STD_010,'o','Color', ColorChoice(6,:),'MarkerFaceColor', ColorChoice(6,:))
% errorbar(params_fit_001, STD_001,'o','Color', ColorChoice(5,:),'MarkerFaceColor', ColorChoice(5,:))

set(gca, 'YScale','log')
legend('100','001','010')

xlim([0 8])
xticklabels({'','K_{b}','K_{r}','\omega_{b}','\omega_{bp}','\omega_{br}','p','R_{max}',''})
xlabel('parameters')
ylabel('fitted values')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'fittedParams_competition','.tif']); 
saveas(gcf,[FigPath,filesep,'fittedParams_competition','.pdf']); 

%% Save the fitting results into structures
save('S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\generalThermoFit_direct_V3.mat',...
        'data_model_input', 'FitResult')
%% Part3. Quenching

end