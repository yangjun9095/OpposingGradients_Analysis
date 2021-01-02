function modeling_6A1R_master_function(varargin)
% Description : 
% 
% Last updated : 11/25/2020, by YJK

% Assumptions : 

%% Step0. Set up the directories
% file to save
SavePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';
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

Bicoid = TFdata(:,1);
Runt = TFdata(:,2);
RuntNull = TFdata(:,3);

% Make a matrix whose each column is each TF
TF(:,1) = Bicoid;
TF(:,2) = Runt;
TF(:,3) = RuntNull;

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

    input_combined(:,1) = [Bicoid(fitRange); Bicoid(fitRange)];
    input_combined(:,2) = [RuntNull(fitRange); Runt(fitRange)];
    

    % save into a structure
    data_model_input(k).construtName = DataTypes(construct);
    data_model_input(k).APbins = APaxis(fitRange);
    data_model_input(k).input = input_combined;
    data_model_input(k).output = data;
    k=k+1;

end   

    
%% Fitting process using the global_fit_construct_generalThermoV2.m
% Model : Direct repression
% Inputs : xdata, ydata, and model for the fitting

% Pick a model from different modes
% Define anonymous function from the given model input. (either direct,
% competition, quenching, etc.)
mdl0 = @model_6A1R_HillModel_V3;
mdl = @(params, TF) model_6A1R_HillModel_V3(params, TF);

% Define the parameter bounds for the fitting
% params = [Kb, Kr, w_a, w_ap, w_rp, p, R_max];
lb = [1, 1, 1, 0.001, 0, 50];
ub = [10^2, 10^2, 10^3, 1, 1, 350];

% initial input of parameter guess
params0 = [10, 10, 10, 0.2, 0.01, 300]; % example

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
    
%% check the fitted parameters by pluggin back into the model,
% to see whether the fit matches well with the data

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
    pause
%   %Save the plot
%     saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.pdf']); 
end

%% Dissect each term in the function of HillV3.
% input_combined(:,1)
% P_bound = p*(1+b.^6*w_bp + r*w_rp + b.^6.*r*w_bp*w_rp)./...
%             (1 + p + b.^6 + r + b.^6.*r + b.^6*p*w_bp + r*p*w_rp + b.^6.*r*p*w_bp*w_rp);

end