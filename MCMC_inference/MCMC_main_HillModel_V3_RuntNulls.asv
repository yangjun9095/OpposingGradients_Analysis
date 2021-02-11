function MCMC_main_HillModel_V3_RuntNulls
%% Description
% This script is doing MCMC fit for all Runt null data one by one, to
% compare the parameters for the Hill.V3 model, to see what varies across
% constructs.

%% Variables

%% Import data for the MCMC inference
% From the "preprocess_data_for_MCMC.m" script
% xdata(TFinputs) and ydata(initial rate), note that we will do
% a simultaneous fitting for the Runt WT and Runt nulls.

% We need another separate script to process the data for inputs in this
% script. : This is now done in the "preprocess_data_for_MCMC.m" script.
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';

% Load the output data
preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
data = preprocessedData.data;

% Load the input TF data
load([FilePath, filesep, 'TFinput.mat'])

Bicoid = TFinput(:,1);
Runt = TFinput(:,2);
RuntNull = TFinput(:,3);

%% Define the FilePath, FigPath to save the results
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls';

%% Set up the MCMC inference
% Note that we don't have any parameters regarding the repressor, thus we
% will just use the simpler form of the model, "model_6A0R_HillModel_V3.m"

%% MCMC Default settings (also make it optional)
% fileDir = pwd;
% saveLoc = pwd;
numParPools = 8;
% n_burn = 5000;
n_steps = 2*10^4;
n_simu = n_steps;
% ratePriorWidth = 50;
AP_start = 20; % [% of embryo length]
AP_end = 50;   % [% of embryo length]
% loadPrevious = false;
% globalFit = 1;

%% Pick a model for the MCMC inference
mdl = @(TF, params) model_6A0R_HillModel_V3(params, TF);
model.modelfun = mdl;

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, InputData) sum((InputData.ydata-mdl(InputData.xdata, params)).^2);

%% Loop over all constructs to perform the MCMC inference on the Bcd-dependent parameters
% [Kb, w_bp, p, R_max];
for construct=1:length(data)
    %% Pull the construct of our interest (for the parameter inference)
    % Choose a construct 
    % Pick the dataset from the data.mat
    Data = data(construct);

    % MCMC analysis on the initial slope (averaged over embryos)
    % initialize the AP bins
    APaxis = Data.APbins;

    %Truncate APbins to user-specified range (input was optional argument in
    %this function.
    NoNaN_index_null = ~isnan(Data.Rate_null);
    NoNaN_index_WT = ~isnan(Data.Rate_WT);
    % calculate the AP bins that are not NaNs in both WT and Null datasets
    NoNaN_index = NoNaN_index_null.*NoNaN_index_WT;

    NoNaNindices = find(NoNaN_index);

    % Range that is set as an initial guess. We will get a common set of
    % APbins that does not have NaN values in these AP bins.
    APbin_start = 20/2.5 + 1;
    APbin_end = 45/2.5 + 1;

    APbinRange = (APbin_start:APbin_end)';

    % find the common elements of AP bins between Not-NaNs, and also the
    % pre-set range (20-45%)
    APbins_fit = intersect(NoNaNindices, APbinRange);

    % initialize the initial rate (slope)
    Rate_WT = Data.Rate_WT;
    Rate_null = Data.Rate_null;

    % Truncate the vectors using the range of AP bins
    APbins = APaxis(APbins_fit);
    Rate_WT = Rate_WT(APbins_fit);
    Rate_null = Rate_null(APbins_fit);

    Bcd = Bicoid(APbins_fit);
    Run = Runt(APbins_fit);
    RunNull = RuntNull(APbins_fit);

    % Decide whether we want to fit the Runt null/WT data together or not.
    % depending on this, we will set the xdata and ydata for the fitting.
    MCMCdata = struct;
    MCMCdata.APdata = [APbins'];
    MCMCdata.ydata = [Rate_null];
    % input TF
    MCMCdata.Bcd = [Bcd];
    MCMCdata.Run = [RunNull];
    MCMCdata.xdata = [MCMCdata.Bcd];

    MCMCdata.R_max = max(MCMCdata.ydata);
    MCMCdata.R_min = min(MCMCdata.ydata);
    
    %% Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}', 'w_bp', 'p', 'R_{max}'};
    params = cell(1, length(names));

    % Initialize MCMC parameters.
    Kb0 = 10;   % 100*rand;
    w_bp0 = 40;
    p0 = 0.1;
    R_max0 = MCMCdata.R_max;
    % R_min0 = MCMCdata.R_min;

    params0 = [Kb0, w_bp0, p0, R_max0];

    % Bounds of the parameters
    LB = [0.1, 1, 0, 50];
    UB = [10^2, 10^2, 1, 500];


    for i = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?

        if i==4
            pri_mu = R_max0;
            pri_sig = 100;
    %         targetflag = 0; % Fix this parameter
    %     elseif i==2
    %         pri_mu = NaN;
    %         pri_sig = Inf;
    %         localflag = 1; % keep this parameter consistent across batches (different constructs of 1 Run site)
        else
            pri_mu = NaN;
            pri_sig = Inf;
        end
        % give a pretty narrow prior
    %     pri_mu = params0(i);
    %     pri_sig = pri_mu*0.1;
        params{1, i} = {names{i}, params0(i), LB(i), UB(i), pri_mu, pri_sig, targetflag, localflag};

    end

    %This is the variance in the parameter proposal distribution. Change these
    %to change the proposal acceptance rate, or if the convergence doesn't look
    %good.

    Kb_step = 0.1;
    w_bp_step = 0.1;
    p_step = 0.01;
    R_max_step = 1;


    % Initialize the covariance matrix
    J0 = diag([Kb_step, w_bp_step, p_step, R_max_step]);


    %% MCMC - Options
    options = [];
    n_steps = 2*10^4;
    options.nsimu = n_steps; %n_steps; %Number of steps
    options.updatesigma = 1; %Update error variance
    options.qcov = J0; %Initial covariance
    % options.burnintime = 0.5*n_steps; %Burn in time
    options.adaptint = 100;
    options.method = 'dram';
    options.verbosity = 0; %Decrease text output

    %% Run the MCMC (this whole block can be inserted inside the for loop above,
    % to run the MCMC for different constructs.

    %we're gonna run this three times and use the initial results of one
    %run as conditions for the next. this is an alternative when common least
    %squares gives results too poor to initialize with


    results = [];
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    %[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,chain,s2chain,~]=mcmcrun(model,MCMCdata,params,options,results);
    
    
    % Extract chain results into individual parameters
    params_inferred = []; 
    params_inferred_sigma = [];
    n_burn = 0.5*n_steps;

    for k=1:length(names)
        params_inferred(1,k) = mean(chain(n_burn+1:end,k));
        params_inferred_std(1,k) = std(chain(n_burn+1:end,k));
    end
    
    %% Save the MCMC result into a structure for future usage
    MCMC_6A0R_RuntNulls(construct).results = results;
    MCMC_6A0R_RuntNulls(construct).chain = chain;
    MCMC_6A0R_RuntNulls(construct).s2chain = s2chain;
    
%     % inferred parameters
    MCMC_6A0R_RuntNulls(construct).params_inferred = params_inferred;
    MCMC_6A0R_RuntNulls(construct).params_inferred_std = params_inferred_std;
    
end   

%% Post-processing for inferred parameters (mean and std)
% for construct=1:length(data)
%     % extract the MCMC chain from each dataset
%     chain_temp = MCMC_6A0R_RuntNulls(construct).chain;
%     
%     params_inferred = []; 
%     params_inferred_std = [];
%     for k=1:length(names)
%         params_inferred(1,k) = mean(chain(n_burn+1:end,k));
%         params_inferred_std(1,k) = std(chain(n_burn+1:end,k));
%     end
%     
%     MCMC_6A0R_RuntNulls(construct).params_inferred = params_inferred;
%     MCMC_6A0R_RuntNulls(construct).params_inferred_std = params_inferred_std;
% end

%% generate plots of inferred parameters

end