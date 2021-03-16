function MCMC_main_HillV3_direct_6A2R_RunWT_usingDist_6AnR_RuntNulls
%% Description
% This script is doing MCMC fit for 6A2R, Run WT data one by one, to
% compare the Run-dependent parameters for the Hill.V3 model, to see what varies across
% constructs.

% One important caveat here is that we'll use the "distribution of posteriors" 
% or the info about mean and std of posteriors for each parameter from
% the previous MCMC inference using 6A0R_HillModelV3 for each construct
% with Runt null (w/o Run protein).

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
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference';

%% Load the 6AnR Run null MCMC result (to plug into [Kb, w_bp, p, R_max]
% tempPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\w_bp_0-100_WeakPrior_w_bp';
tempPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\w_bp_0-100';
load([tempPath, filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat'],'MCMC_6A0R_RuntNulls')

%% Load the 6A1R MCMC inference result for the Run-dependent parameters (Make sure to import the 6A1R, Hill.V3 result)
tempPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\Posterior_dist_Bcd_params_fromRunNulls';
A = load([tempPath, filesep, 'MCMC_6A1R_RuntWT_params.mat'])
MCMC_6A1R_RuntWT = A.MCMC_6A1R_RuntWT;
%% Set up the MCMC inference
% Note that we don't have any parameters regarding the repressor, thus we
% will just use the simpler form of the model, "model_6A0R_HillModel_V3.m"

%% Pick a model for the MCMC inference
mdl = @(TF, params) model_6A2R_HillModel_V3_direct(params, TF);
model.modelfun = mdl;

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, InputData) sum((InputData.ydata-mdl(InputData.xdata, params)).^2);

%% MCMC - Options
options = [];
n_steps = 2*10^5;
options.nsimu = n_steps; %n_steps; %Number of steps
options.updatesigma = 1; %Update error variance
%     options.qcov = J0; %Initial covariance
% options.burnintime = 0.5*n_steps; %Burn in time
options.adaptint = 100;
options.method = 'dram';
options.verbosity = 0; %Decrease text output
%% Loop over all constructs to perform the MCMC inference
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params;

% initialize the structure to save the MCMC result
MCMC_6A2R_RuntWT = struct;
m = 1; % counter to save into the structure, as we're only dealing with 3 constructs (6A1R)

% r2-new : [011] : 3rd in the MCMC_6A0R_RuntNulls
% r2-close : [110] : 7th
% r2-far : [101] : 8th

for construct= [3,7,8] 
    
    %% Pull the construct of our interest (for the parameter inference)
    % Choose a construct 
    % Pick the dataset from the data.mat
    Data = data(construct);

    % MCMC analysis on the initial slope (averaged over embryos)
    % initialize the AP bins
    APaxis = Data.APbins;

    %Truncate APbins to user-specified range (input was optional argument in
    %this function.
%     NoNaN_index_null = ~isnan(Data.Rate_null);
    NoNaN_index_WT = ~isnan(Data.Rate_WT);
    % calculate the AP bins that are not NaNs in both WT and Null datasets
    NoNaN_index = NoNaN_index_WT; %NoNaN_index_null;

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
    Rate_WT_forFit = Rate_WT(APbins_fit);
    %Rate_null_forFit = Rate_null(APbins_fit);

    Bcd = Bicoid(APbins_fit);
    Run = Runt(APbins_fit);
%     RunNull = RuntNull(APbins_fit);

    % Decide whether we want to fit the Runt null/WT data together or not.
    % depending on this, we will set the xdata and ydata for the fitting.
    MCMCdata = struct;
    MCMCdata.APdata = [APbins'];
    MCMCdata.ydata = [Rate_WT_forFit];
    % input TF
    MCMCdata.Bcd = [Bcd];
    MCMCdata.Run = [Run];
    MCMCdata.xdata = [MCMCdata.Bcd , MCMCdata.Run];

    MCMCdata.R_max = max(MCMCdata.ydata);
    MCMCdata.R_min = min(MCMCdata.ydata);
    
    %% Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}', 'w_{bp}', 'p', 'R_{max}','K_{r1}','K_{r2}','w_{rp1}','w_{rp2}'};
    params = cell(1, length(names));
        
    % Pull the fitted parameters from the MCMC inference on the Runt
    % null datasets for one-by-one.
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_Bcd_std = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
    
    % Initialize MCMC parameters.
    Kb0 = params_Bcd(1);
    w_bp0 = params_Bcd(2);
    p0 = params_Bcd(3);
    R_max0 = params_Bcd(4);
    
    % Initial condition for the Run-dependent parameters
    % we will utilize the posteriors from the 6A1R constructs as the prior.
    % There's a nail in the coffin that we're assuming the prior as
    % Gaussian, let's assume that it's a fair assumption for now.
    
    if i==1 % r2-new : [011]
        params_temp1 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
        params_temp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
        params_std_temp1 =  MCMC_6A1R_RuntWT(2).params_inferred_sigma; %[001]
        params_std_temp2 =  MCMC_6A1R_RuntWT(3).params_inferred_sigma; %[010]
    elseif i==2 % r2-close : [110]
        params_temp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        params_temp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
        params_std_temp1 =  MCMC_6A1R_RuntWT(1).params_inferred_sigma; %[100]
        params_std_temp2 =  MCMC_6A1R_RuntWT(3).params_inferred_sigma; %[010]
    elseif i==3 % r2-far : [101]
        params_temp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        params_temp2 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
        params_std_temp1 =  MCMC_6A1R_RuntWT(1).params_inferred_sigma; %[100]
        params_std_temp2 =  MCMC_6A1R_RuntWT(2).params_inferred_sigma; %[001]
    end
    
    K_r1_0 = params_temp1(2);
    K_r1_0_std = params_std_temp1(2);
    K_r2_0 = params_temp2(2);
    K_r2_0_std = params_std_temp2(2);
    
    w_rp1_0 = params_temp1(4);
    w_rp1_0_std = params_std_temp1(4);
    w_rp2_0 = params_temp2(2);
    w_rp2_0_std = params_std_temp2(2);
    
    

    params0 = [Kb0, w_bp0, p0, R_max0, K_r1_0, K_r2_0, w_rp1_0, w_rp2_0];
    params0_std = [params_Bcd_std(1), params_Bcd_std(2), params_Bcd_std(3), params_Bcd_std(4),...
                    K_r1_0_std, K_r2_0_std, w_rp1_0_std, w_rp2_0_std];

    % Bounds of the parameters
    LB = [0.1, 1, 0, 50, 0.1, 0.1, 0, 0];
    UB = [10^2, 10^2, 1, 500, 100, 100, 1.2, 1.2];


    for i = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?
        
        % For all parameters, define the priors from the posteriors from
        % the previous rounds of MCMC (either on Runt nulls or 6A1R cases)
        pri_mu = params0(i);
        pri_sig = params0_std(i);

        
    %     elseif i==2
    %         pri_mu = NaN;
    %         pri_sig = Inf;
    %         localflag = 1; % keep this parameter consistent across batches (different constructs of 1 Run site)
        % give a pretty narrow prior
    %     pri_mu = params0(i);
    %     pri_sig = pri_mu*0.1;
    
        params{1, i} = {names{i}, params0(i), LB(i), UB(i), pri_mu, pri_sig, targetflag, localflag};

    end

    %This is the variance in the parameter proposal distribution. Change these
    %to change the proposal acceptance rate, or if the convergence doesn't look
    %good.

%     Kb_step = 0.1;
%     Kr_step = 0.1;
%     w_bp_step = 1;
%     w_rp_step = 0.01;
%     p_step = 0.01;
%     R_max_step = 1;
% 
% 
%     % Initialize the covariance matrix
%     J0 = diag([Kb_step, Kr_step, w_bp_step, w_rp_step, p_step, R_max_step]);


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
        params_inferred_sigma(1,k) = std(chain(n_burn+1:end,k));
    end
    
    %% Save the MCMC result into a structure for future usage
    MCMC_6A2R_RuntWT(m).name = constructNames{construct};
    MCMC_6A2R_RuntWT(m).results = results;
    MCMC_6A2R_RuntWT(m).chain = chain;
    MCMC_6A2R_RuntWT(m).s2chain = s2chain;
    
%     % inferred parameters
    MCMC_6A2R_RuntWT(m).params_inferred = params_inferred;
    MCMC_6A2R_RuntWT(m).params_inferred_sigma = params_inferred_sigma;
    
    % update the counter
    m=m+1;
end   


%% Part2. Treating the two Run sites the same : meaning the same parameters
% This is because the two binding sites are "symmetric" in the sense that 

%% Pick a model for the MCMC inference
mdl = @(TF, params) model_6A2R_HillModel_V3_direct_sameRunParams(params, TF);
model.modelfun = mdl;

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, InputData) sum((InputData.ydata-mdl(InputData.xdata, params)).^2);

%% MCMC - Options
options = [];
n_steps = 2*10^4;
options.nsimu = n_steps; %n_steps; %Number of steps
options.updatesigma = 1; %Update error variance
%     options.qcov = J0; %Initial covariance
% options.burnintime = 0.5*n_steps; %Burn in time
options.adaptint = 100;
options.method = 'dram';
options.verbosity = 0; %Decrease text output
%% Loop over all constructs to perform the MCMC inference
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params;

% initialize the structure to save the MCMC result
MCMC_6A2R_RuntWT = struct;
m = 1; % counter to save into the structure, as we're only dealing with 3 constructs (6A1R)

% r2-new : [011] : 3rd in the MCMC_6A0R_RuntNulls
% r2-close : [110] : 7th
% r2-far : [101] : 8th

for construct= [3,7,8] 
    
    %% Pull the construct of our interest (for the parameter inference)
    % Choose a construct 
    % Pick the dataset from the data.mat
    Data = data(construct);

    % MCMC analysis on the initial slope (averaged over embryos)
    % initialize the AP bins
    APaxis = Data.APbins;

    %Truncate APbins to user-specified range (input was optional argument in
    %this function.
%     NoNaN_index_null = ~isnan(Data.Rate_null);
    NoNaN_index_WT = ~isnan(Data.Rate_WT);
    % calculate the AP bins that are not NaNs in both WT and Null datasets
    NoNaN_index = NoNaN_index_WT; %NoNaN_index_null;

    NoNaNindices = find(NoNaN_index);

    % Range that is set as an initial guess. We will get a common set of
    % APbins that does not have NaN values in these AP bins.
    APbin_start = 20/2.5 + 1;
    APbin_end = 50/2.5 + 1;

    APbinRange = (APbin_start:APbin_end)';

    % find the common elements of AP bins between Not-NaNs, and also the
    % pre-set range (20-45%)
    APbins_fit = intersect(NoNaNindices, APbinRange);

    % initialize the initial rate (slope)
    Rate_WT = Data.Rate_WT;
    Rate_null = Data.Rate_null;

    % Truncate the vectors using the range of AP bins
    APbins = APaxis(APbins_fit);
    Rate_WT_forFit = Rate_WT(APbins_fit);
    %Rate_null_forFit = Rate_null(APbins_fit);

    Bcd = Bicoid(APbins_fit);
    Run = Runt(APbins_fit);
%     RunNull = RuntNull(APbins_fit);

    % Decide whether we want to fit the Runt null/WT data together or not.
    % depending on this, we will set the xdata and ydata for the fitting.
    MCMCdata = struct;
    MCMCdata.APdata = [APbins'];
    MCMCdata.ydata = [Rate_WT_forFit];
    % input TF
    MCMCdata.Bcd = [Bcd];
    MCMCdata.Run = [Run];
    MCMCdata.xdata = [MCMCdata.Bcd , MCMCdata.Run];

    MCMCdata.R_max = max(MCMCdata.ydata);
    MCMCdata.R_min = min(MCMCdata.ydata);
    
    %% Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}', 'w_{bp}', 'p', 'R_{max}','K_{r1}','w_{rp1}'};
    params = cell(1, length(names));
        
    % Pull the fitted parameters from the MCMC inference on the Runt
    % null datasets for one-by-one.
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_Bcd_std = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
    
    % Initialize MCMC parameters.
    Kb0 = params_Bcd(1);
    w_bp0 = params_Bcd(2);
    p0 = params_Bcd(3);
    R_max0 = params_Bcd(4);
    
    % Initial condition for the Run-dependent parameters
    % we will utilize the posteriors from the 6A1R constructs as the prior.
    % There's a nail in the coffin that we're assuming the prior as
    % Gaussian, let's assume that it's a fair assumption for now.
    
    if i==1 % r2-new : [011]
        params_temp1 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
        params_temp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
        params_std_temp1 =  MCMC_6A1R_RuntWT(2).params_inferred_sigma; %[001]
        params_std_temp2 =  MCMC_6A1R_RuntWT(3).params_inferred_sigma; %[010]
    elseif i==2 % r2-close : [110]
        params_temp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        params_temp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
        params_std_temp1 =  MCMC_6A1R_RuntWT(1).params_inferred_sigma; %[100]
        params_std_temp2 =  MCMC_6A1R_RuntWT(3).params_inferred_sigma; %[010]
    elseif i==3 % r2-far : [101]
        params_temp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        params_temp2 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
        params_std_temp1 =  MCMC_6A1R_RuntWT(1).params_inferred_sigma; %[100]
        params_std_temp2 =  MCMC_6A1R_RuntWT(2).params_inferred_sigma; %[001]
    end
    
    K_r1_0 = params_temp1(2);
    K_r1_0_std = params_std_temp1(2);
    K_r2_0 = params_temp2(2);
    K_r2_0_std = params_std_temp2(2);
    
    w_rp1_0 = params_temp1(4);
    w_rp1_0_std = params_std_temp1(4);
    w_rp2_0 = params_temp2(2);
    w_rp2_0_std = params_std_temp2(2);
    
    

    params0 = [Kb0, w_bp0, p0, R_max0, K_r1_0, w_rp1_0];
    params0_std = [params_Bcd_std(1), params_Bcd_std(2), params_Bcd_std(3), params_Bcd_std(4),...
                    K_r1_0_std, w_rp1_0_std];

    % Bounds of the parameters
    LB = [0.1, 1, 0, 50, 0.1, 0];
    UB = [10^2, 10^2, 1, 500, 100, 1.2];

   
    for i = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?
        
        % For all parameters, define the priors from the posteriors from
        % the previous rounds of MCMC (either on Runt nulls or 6A1R cases)
        if i<5
            pri_mu = params0(i);
            pri_sig = params0_std(i);
            targetflag = 0; % fix the parameter
        else    
            pri_mu = params0(i);
            pri_sig = params0_std(i);
        end
        
    %     elseif i==2
    %         pri_mu = NaN;
    %         pri_sig = Inf;
    %         localflag = 1; % keep this parameter consistent across batches (different constructs of 1 Run site)
        % give a pretty narrow prior
    %     pri_mu = params0(i);
    %     pri_sig = pri_mu*0.1;
    
        params{1, i} = {names{i}, params0(i), LB(i), UB(i), pri_mu, pri_sig, targetflag, localflag};

    end

    %This is the variance in the parameter proposal distribution. Change these
    %to change the proposal acceptance rate, or if the convergence doesn't look
    %good.

%     Kb_step = 0.1;
%     Kr_step = 0.1;
%     w_bp_step = 1;
%     w_rp_step = 0.01;
%     p_step = 0.01;
%     R_max_step = 1;
% 
% 
%     % Initialize the covariance matrix
%     J0 = diag([Kb_step, Kr_step, w_bp_step, w_rp_step, p_step, R_max_step]);


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

    for k=1:2%length(names)
        params_inferred(1,k) = mean(chain(n_burn+1:end,k));
        params_inferred_sigma(1,k) = std(chain(n_burn+1:end,k));
    end
    
    %% Save the MCMC result into a structure for future usage
    MCMC_6A2R_RuntWT(m).name = constructNames{construct};
    MCMC_6A2R_RuntWT(m).results = results;
    MCMC_6A2R_RuntWT(m).chain = chain;
    MCMC_6A2R_RuntWT(m).s2chain = s2chain;
    
%     % inferred parameters
    MCMC_6A2R_RuntWT(m).params_inferred = params_inferred;
    MCMC_6A2R_RuntWT(m).params_inferred_sigma = params_inferred_sigma;
    
    % update the counter
    m=m+1;
end   


%% generate raw fits 
APaxis = 0:0.025:1;
% APbin_start;
% APbin_end;
Bcd = Bicoid(APbin_start:APbin_end); 
Run = Runt(APbin_start:APbin_end); 
RunNull = RuntNull(APbin_start:APbin_end);

TF = [Bcd, Run];
TF_null = [Bcd, RunNull];

for index = 1:3
    % define the construct index (which is consistent with the way it's
    % defiend in the compiledData.mat)
    constructIndex = [3,7,8];
    construct = constructIndex(index);
    
    % extract parameters from the MCMC results (we're pulling the Bcd
    % dependent parameters and Run dependent parameters respectively from
    % their own inferences. (Bcd parameters from the Run null, and Run
    % parameters from the Run WT).
%     params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
%     params_Run = MCMC_6A1R_RuntWT(index).params_inferred;
    
    % extract individual parameters to construct the input parameter
%     Kb = params_Bcd(1);
%     w_bp = params_Bcd(2);
%     p = params_Bcd(3);
%     R_max = params_Bcd(4);
%     Kr = params_Run(1);
%     w_rp = params_Run(2);
%   
    % point estimates from the Runt nulls
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    % MCMC inferred parameters for the Runt WT
    params_temp = MCMC_6A2R_RuntWT(index).params_inferred;
    % combine the two sets of parameters
    params_MCMC = [params_Bcd, params_temp];
    
    output = model_6A2R_HillModel_V3_direct_sameRunParams(params_MCMC, TF);
    % Runt nulls
    fit_nulls = model_6A2R_HillModel_V3_direct_sameRunParams(params_MCMC, TF_null);
    
    % data (Runt nulls)
    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};
    
    % data (WT)
    Rate_WT = compiledData{construct+1,9};
    Rate_WT_SEM = compiledData{construct+1,10};
    
    clf
    hold on
    errorbar(APaxis, Rate_null, Rate_null_SEM, 'o', 'Color', ColorChoice(4,:),'LineWidth', 1)    
    errorbar(APaxis, Rate_WT, Rate_WT_SEM, 'o', 'Color', ColorChoice(1,:),'LineWidth', 1)
    
    plot(APaxis(APbin_start:APbin_end), fit_nulls, 'Color', ColorChoice(4,:),'LineWidth', 2)
    plot(APaxis(APbin_start:APbin_end), output, 'Color', ColorChoice(1,:),'LineWidth', 2)
    xlim([0.2 0.5])
    xticks([0.2 0.3 0.4 0.5])
    ylim([0 400])
    yticks([0 100 200 300 400])
    
    xlabel('embryo length')
    ylabel('initial rate (AU/min)')
    
    box on
    legend('data(null)','data(WT)','Fit (null)', 'Fit (WT)')
    StandardFigure(gcf,gca)
    pause
    saveas(gcf,[FigPath,filesep,'raw_fits_null_WT_', constructNames{construct}  ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'raw_fits_null_WT_', constructNames{construct} ,'.pdf']); 
end

%% generate corner plots
for index = 1:3
%     clf
    construct = constructIndex(index);
    chain = MCMC_6A1R_RuntWT(index).chain;
    n_burn = 0.5*n_steps;
    m = [chain(n_burn+1:end,1), chain(n_burn+1:end,2)];
    corner = figure;
%     names = {'K_{b}','\omega_{bp}','p','R_{max}'};
    names = {'K_{r}','\omega_{rp}'};
    ecornerplot(m,'names',names);
    
    saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct}  ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']); 
end

%% generate plots of inferred parameters
hold on
for index = 1:3
    construct = constructIndex(index);
    params_inferred = MCMC_6A2R_RuntWT(index).params_inferred;
    params_inferred_std = MCMC_6A2R_RuntWT(index).params_inferred_sigma;
    errorbar(1:2, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 3])
xticks([1,2])
xticklabels({'K_{r}', '\omega_{rp}'})
% yticks([0 100 200 300 400])
xlabel('inferred parameters')
ylabel('inferred parameters')
legend('100','001','010','Location', 'NorthWest')

box on
StandardFigure(gcf,gca)

% % save the plots
saveas(gcf,[FigPath,filesep, 'MCMCfit_6A2R_RuntWT_pt_estimate_BcdParams','.tif']);
saveas(gcf,[FigPath,filesep, 'MCMCfit_6A2R_RuntWT_pt_estimate_BcdParams','.pdf']);

%% generate the plot of inferred parameters (log scale)
hold on
for index = 1:3
    construct = constructIndex(index);
    params_inferred = MCMC_6A2R_RuntWT(index).params_inferred;
    params_inferred_std = MCMC_6A2R_RuntWT(index).params_inferred_sigma;
    errorbar(1:2, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
end
 
xlim([0 3])
xticks([1,2])
xticklabels({'K_{r}', '\omega_{rp}'})
% yticks([0 100 200 300 400])
xlabel('inferred parameters')
ylabel('inferred parameters')
legend('100','001','010','Location', 'NorthWest')

set(gca,'YScale','Log')
box on
StandardFigure(gcf,gca)

% % save the plots
saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_post_dist_BcdParams_LogScale','.tif']);
saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_post_dist_BcdParams_LogScale','.pdf']);



%% generate the plot for C.V. of inferred parameters

% calculate the Coefficient of Variation (C.V.) of inferred parameters
% across constructs (different enhancers).

% First, initialize the matrix
params_inferred_all = [];
params_inferred_error_all = [];

constructIndex = [3,7,8];

for index = 1:3
    construct = constructIndex(index);
    % params from the Runt null
    params_inferred_null(index,:) = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_inferred_std_null(index,:) = MCMC_6A0R_RuntNulls(construct).params_inferred;
    % params from the Runt WT
    params_inferred_all(index,:) = MCMC_6A2R_RuntWT(index).params_inferred;
    params_inferred_error_all(index,:) = MCMC_6A2R_RuntWT(index).params_inferred_sigma;
end

% calculate the mean and std of parameters "over constructs".
params_null_mean = mean(params_inferred_null);
params_null_std = std(params_inferred_std_null);

params_mean = mean(params_inferred_all);
params_std = std(params_inferred_all);
% params_error = sqrt(sum(params_inferred_error_all.^2)/length(data))
n_boots = 100;
params_null_std_boostrap = bootstrp(n_boots, @std, params_inferred_null);
parms_std_boostrap = bootstrp(n_boots, @std, params_inferred_all);

params_null_SEM = std(params_null_std_boostrap)./sqrt(n_boots-1);
params_SEM = std(parms_std_boostrap)./sqrt(n_boots-1);

params_null_CV = params_null_std./params_null_mean;
params_null_CV_error = params_null_SEM./params_null_mean;

params_CV = params_std./params_mean;
params_CV_error = params_SEM./params_mean;

% plot
hold on
errorbar(5:6, params_CV, params_CV_error,'o','LineWidth',2)
errorbar(1:4, params_null_CV, params_null_CV_error, 'o','LineWidth',2)

xlim([0 7])
xticks([1,2,3,4,5,6])
xticklabels({'K_{b}','\omega_{bp}','p','R_{max}','K_{r}','\omega_{rp}'})
% yticks([0 100 200 300 400])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2])
% set(gca,'YScale','log')
% xlabel('')
ylabel('Coefficient of Variation')
legend('MCMC (Run WT)','MCMC (Run null)')

box on
StandardFigure(gcf,gca)

% save the plots
saveas(gcf,[FigPath,filesep, 'CV_parameters','.tif']);
saveas(gcf,[FigPath,filesep, 'CV_parameters','.pdf']);
%% test the sensitivity of certain parameters (for the case of [011] or [111], where the raw fits don't quite look great.
% construct = 4;
% params = MCMC_6A0R_RuntNulls(construct).params_inferred;
% output = model_6A0R_HillModel_V3(params, Bcd);
% 
% % data
% Rate_null = compiledData{construct+9,9};
% Rate_null_SEM = compiledData{construct+9,10};
% 
% figure
% hold on
% errorbar(APaxis, Rate_null, Rate_null_SEM)
% plot(APaxis(9:19), output)

%% save the result into mat files (MCMC_6A0R_RuntNulls)
FilePath = FigPath;
save([FilePath, filesep, 'MCMC_6A2R_RuntWT_params.mat'],'MCMC_6A2R_RuntWT')

end