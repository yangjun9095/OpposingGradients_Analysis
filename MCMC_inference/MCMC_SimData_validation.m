%% MCMC inference validation
% Last updated : 
% This script is for validating the MCMC inference protocol, meaning that
% we perform the MCMC inference on simulated datasets, to see whether we can
% retrieve the original parameters.

%% File directories
FilePath
DataPath
FigPath

%% Load the real data (TF input) to generate the simulated data (SimData)
% preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
% data = preprocessedData.data;

% Load the input TF data
load([FilePath, filesep, 'TFinput.mat'])

Bicoid = TFinput(:,1);
Runt = TFinput(:,2);
RuntNull = TFinput(:,3);

% Define the AP range to extract the data
APrange = 9:19; % 20-45% of the embryo length.

% input TF
MCMCdata.Bcd = [Bicoid(APrange) ; Bicoid(APrange)];
MCMCdata.Run = [RuntNull(APrange) ; Runt(APrange)];
MCMCdata.xdata = [MCMCdata.Bcd, MCMCdata.Run];

TF = MCMCdata.xdata;
%% Define the model to simulate the data
mdl = @(TF, params) model_6A1R_direct_repression_V2(params, TF);

%% Import the inferred parameter distribution from a dataset ([001])
MCMCresult = load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_inference_ThermoModelV2\MCMCresults\MCMCresult_2001.mat');

% extract the chain
MCMCchain = MCMCresult.chain;

% extract
n_simu = MCMCresult.results.nsimu;
n_burn = 0.5*n_simu;

% get the mean and std from the inferred parameter dist.
names = {'Kb','Kr','w_b','w_{bp}','w_{rp}','p','R_{max}'};

mean_Kb = mean(MCMCchain(n_burn+1:end,1));
sigma_Kb = std(MCMCchain(n_burn+1:end,1));

mean_Kr = mean(MCMCchain(n_burn+1:end,2));
sigma_Kr = std(MCMCchain(n_burn+1:end,2));

mean_w_b = mean(MCMCchain(n_burn+1:end,3));
sigma_w_b = std(MCMCchain(n_burn+1:end,3));

mean_w_bp = mean(MCMCchain(n_burn+1:end,4));
sigma_w_bp = std(MCMCchain(n_burn+1:end,4));

mean_w_rp = mean(MCMCchain(n_burn+1:end,5));
sigma_w_rp = std(MCMCchain(n_burn+1:end,5));

mean_p = mean(MCMCchain(n_burn+1:end,6));
sigma_p = std(MCMCchain(n_burn+1:end,6));

mean_R_max = mean(MCMCchain(n_burn+1:end,7));
sigma_R_max = std(MCMCchain(n_burn+1:end,7));

params_inferred = [mean_Kb, mean_Kr, mean_w_b, mean_w_bp, mean_w_rp, mean_p, mean_R_max];
params_inferred_std = [sigma_Kb, sigma_Kr, sigma_w_b, sigma_w_bp, sigma_w_rp, sigma_p, sigma_R_max];
%% generate simulated data by sampling from the normal dist. of parameters with mu and sigma defined above (from the inferred dist.)
n_sim = 100; % number of simulation

MCMC_SimData = struct;

for i=1:n_sim
    % sample the parameters from normal dist. with mu and sigma as above
    Kb = normrnd(mean_Kb, sigma_Kb);
    Kr = normrnd(mean_Kr, sigma_Kr);
    w_b = normrnd(mean_w_b, sigma_w_b);
    w_bp = normrnd(mean_w_bp, sigma_w_bp);
    w_rp = normrnd(mean_w_rp, sigma_w_rp);
    p = normrnd(mean_p, sigma_p);
    R_max = normrnd(mean_R_max, sigma_R_max);

    params_set = [Kb, Kr, w_b, w_bp, w_rp, p, R_max];
    SimData = model_6A1R_direct_repression_V2(params_set, TF);
    
    % define the MCMC data as an input for the model
    MCMCdata.SimData = SimData;
    MCMCdata.TF = TF;
    
    % Perform the MCMC inference on the simulated data (SimData)
    %  Pick a model for the fitting
    mdl = @(TF, params) model_6A1R_direct_repression_V2(params, TF);

    %leaving this here in case it'll be useful in the future
    model.ssfun = @(params, MCMCdata) sum((MCMCdata.SimData-mdl(MCMCdata.TF, params)).^2);

    model.modelfun = mdl;  %use mcmcrun generated ssfun 
    model.N = length(SimData);

    % Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}','K_{r}','w_{b}','w_{bp}','w_{rp}','p','R_{max}'};
    params = cell(1, length(names));

    % Initialize MCMC parameters.
    Kb0 = 1000;   % 100*rand;
    Kr0 = 10;     % 100*rand;
    w_b0 = 2;    % 10*rand;
    w_bp0 = 5;   % 10*rand;
    % repression (0< w <1)
    w_rp0 = 0.1;
    p0 = 0.004; %
    R_max0 = R_max; %500*rand;


    params0 = [Kb0, Kr0, w_b0, w_bp0, w_rp0, p0, R_max0];

    % Bounds of the parameters
    LB = [0.1, 0.1, 1, 1, 10^(-6), 0, 50];
    UB = [5*10^3, 10^2, 10^6, 10^6, 1, 10, 10^3];

    % pri_mu = NaN; %default prior gaussian mean
    % pri_sig = Inf; %default prior gaussian variance
    localflag = 0; %is this local to this dataset or shared amongst batches?

    for j = 1:length(names)

        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.

        if j==7
            pri_mu(j) = R_max;
            pri_sig(j) = 20;
    %         targetflag = 0; % Fix this parameter
        else
            pri_mu(j) = NaN;
            pri_sig(j) = Inf;
        end
        params{1, j} = {names{j},params0(j), LB(j), UB(j), pri_mu(j), pri_sig(j), targetflag, localflag};

    end

    %This is the variance in the parameter proposal distribution. Change these
    %to change the proposal acceptance rate, or if the convergence doesn't look
    %good.

    Kb_step = 10;
    Kr_step = 1;
    w_b_step = 0.1;
    w_bp_step = 0.1;
    p_step = 0.001;
    R_max_step = 1;

    % repression
    w_rp_step = 10^-(5);


    % Initialize the covariance matrix
    J0 = diag([Kb_step, Kr_step, w_b_step, w_bp_step,...
        w_rp_step, p_step, R_max_step]);


    % MCMC - Options
    options = [];
    n_steps = 10^4;
    options.nsimu = n_steps; %n_steps; %Number of steps
    options.updatesigma = 1; %Update error variance
    options.qcov = J0; %Initial covariance
    % options.burnintime = 0.5*n_steps; %Burn in time
    options.adaptint = 100;
    options.method = 'dram';
    options.verbosity = 0; %Decrease text output

    % Run the MCMC (this whole block can be inserted inside the for loop above,
    % to run the MCMC for different constructs.

    %we're gonna run this three times and use the initial results of one
    %run as conditions for the next. this is an alternative when common least
    %squares gives results too poor to initialize with
    results = [];
%     [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
%     [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    %[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,chain,s2chain,~]=mcmcrun(model,MCMCdata,params,options,results);
    
    % save the result
    MCMC_SimData(i).params_set = params_set;
%     MCMC_SimData(i).SimData = SimData;
    MCMC_SimData(i).chain = chain;
    MCMC_SimData(i).results = results;
    MCMC_SimData(i).params_inferred = mean(chain(5001:end,:));
    
end

%% Extract the MCMC results, then get the distribution

% A = extractfield(MCMC_SimData, 'params_inferred')
% reshape(A, n_sim, 7)

for i=1:n_sim
    params_compiled(i,:) = MCMC_SimData(i).params_inferred;
end

%% generate histograms
figure
% hold on
for i=1:length(params_inferred)-1
    subplot(2,3,i)

    if i<6
        hold on
        nbins = 10;
        histogram(params_compiled(:,i), nbins)
        xline(params_inferred(i),'--')
    elseif i==6
        hold on
        edges = [0.001:0.0005:0.01];
        histogram(params_compiled(:,i),edges)
        xline(params_inferred(i),'--')
    end

    xlabel(names{i})
    box on
    pause
%     StandardFigure(gcf,gca)
end