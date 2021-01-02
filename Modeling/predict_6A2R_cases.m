%% predict_6A2R_cases
function predict_6A2R_cases
%% Description
% This script is to predict the level of hbP2 + 2Run sites

% input : 
% model type (parameters)
% dataset (constructs)
% 
% 
%% FilePaths
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\Predict_HillV3';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\Predict_HillV3';

%% Load the data (Bcd, Run, and output)

% First, output : [001], [010], [100], load all of the data and parameters inferred
preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
data = preprocessedData.data;

% Second, the input TF data : Bcd, Run
load([FilePath, filesep, 'TFinput.mat'])

Bicoid = TFinput(:,1);
Runt = TFinput(:,2);
RuntNull = TFinput(:,3);

%% Fit the Runt null to extract the Bcd dependent parameters
% Kb, w_bp, p, R_max

% Default settings for the MCMC
% MCMC - Options
options = [];
n_steps = 2*10^4;
options.nsimu = n_steps; %n_steps; %Number of steps
options.updatesigma = 1; %Update error variance
options.qcov = J0; %Initial covariance
% options.burnintime = 0.5*n_steps; %Burn in time
options.adaptint = 100;
options.method = 'dram';
options.verbosity = 0; %Decrease text output

% Define a structure to save the MCMC results for the parameters inferred
% from the Runt null cases.
MCMC_2RunSites = struct;
k=1; % counter

for construct = [3, 7, 8] % indices from the data.mat
    % perform the MCMC for each construct for a Hill.V3 model (as the Runt
    % dependent part is all zero, thus the Runt dependent parameters would
    % not be inferred.
    % Pull the construct of our interest (for the parameter inference)

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
    MCMCdata.APdata = [APbins' ; APbins'];
    MCMCdata.ydata = [Rate_null ; Rate_WT];
    % input TF
    MCMCdata.Bcd = [Bcd ; Bcd];
    MCMCdata.Run = [RunNull ; Run];
    MCMCdata.xdata = [MCMCdata.Bcd, MCMCdata.Run];

    MCMCdata.R_max = max(MCMCdata.ydata);
    MCMCdata.R_min = min(MCMCdata.ydata);

    % Pick a model for the fitting
    mdl = @(TF, params) model_6A1R_HillModel_V3(params, TF);

    %leaving this here in case it'll be useful in the future
    model.ssfun = @(params, data) sum((data.ydata-mdl(data.xdata, params)).^2);

    model.modelfun = mdl;  %use mcmcrun generated ssfun 

    % Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}', 'K_{r}', 'w_bp', 'w_rp', 'p', 'R_{max}'};
    params = cell(1, length(names));

    % Initialize MCMC parameters.
    Kb0 = 10;   % 100*rand;
    Kr0 = 5;     % 100*rand;
    w_bp0 = 40;
    w_rp0 = 0.2;
    p0 = 0.1;
    R_max0 = MCMCdata.R_max; 
    % R_min0 = MCMCdata.R_min;

    params0 = [Kb0, Kr0, w_bp0, w_rp0, p0, R_max0];

    % Bounds of the parameters
    LB = [0.1, 0.1, 1, 0, 0, 50];
    UB = [10^2, 10^2, 10^2, 1.2, 1, 500];


    for i = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?

        if i==6
            pri_mu = MCMCdata.R_max;
            pri_sig = 20;
    %         targetflag = 0; % Fix this parameter
        elseif i==2 || i==4
            pri_mu = NaN;
            pri_sig = Inf;
            targetflag = 0; % Fix this parameter
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
    Kr_step = 0.01;
    w_bp_step = 0.1;
    w_rp_step = 0.01;
    p_step = 0.01;
    R_max_step = 1;


    % Initialize the covariance matrix
    J0 = diag([Kb_step, Kr_step, w_bp_step, w_rp_step, p_step, R_max_step]);

    % Run the MCMC (this whole block can be inserted inside the for loop above,
    % to run the MCMC for different constructs.

    %we're gonna run this three times and use the initial results of one
    %run as conditions for the next. this is an alternative when common least
    %squares gives results too poor to initialize with
    results = [];
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    %[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,chain,s2chain,~]=mcmcrun(model,MCMCdata,params,options,results);
    
    % save the result into the structure
    MCMC_2RunSites(k).results = results;
    MCMC_2RunSites(k).chain = chain;
    MCMC_2RunSites(k).construct = data(construct).constructName;
    % update the counter
    k = k+1; 
end

%% Check the MCMC result for the Bcd dependent parameters
datanum = 1;
chain = MCMC_2RunSites(datanum).chain;
results = MCMC_2RunSites(datanum).results;
%% Diagnose the MCMC result
% stats = chainstats(chain,results);
n_burn = 0.5*n_steps;
chainstats(chain(n_burn+1:end,:),results);

%% generate corner plots
%n_burn = 1;%0.5*n_steps;% 0.5*10^4;
m = [chain(n_burn:end,1), chain(n_burn:end,2), chain(n_burn:end,3), chain(n_burn:end,4)];
corner = figure;
names = {'K_{b}','\omega_{bp}', 'p','R_{max}'};
ecornerplot(m,'names',names);


% Save the plot
% FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';
% saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct}  ,'.tif']); 
% saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']);  

%% Grab the Runt dependent parameters (chain) from the MCMC

%% Predict the level for the Runt WT

end