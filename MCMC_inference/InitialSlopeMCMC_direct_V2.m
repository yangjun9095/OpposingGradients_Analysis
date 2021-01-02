function InitialSlopeMCMC_direct_V2(varargin)

%% Description
% Fits the initial slope profile over the AP axis for Runt WT/Null (or both simultaneously) 
% using Markov Chain Monte Carlo, implemented in the MCMCstat package. 
% Data specification

% The code can load multiple datasets, and possesses user-specified options
% as variable arguments. The details of these inputs are described in the
% readme and in the manuscript, and are initialized with default values if
% not specified by the user.

%Variable arguments:
%   'fileDir', fileDir: directory of dataset files (default = root
%   directory)
%   'saveLoc', saveLoc: directory of saved MCMC results (default = root
%   directory)
%   'numParPools', numParPools: number of parallel workers to use (default
%   = 8)
%   'n_burn', n_burn: number of burn-in steps in MCMC algorithm (default = 10000)
%   'n_steps', n_steps: number of steps in MCMC algorithm including burn-in (default =
%   20000)
%   'construct', construct: option to specify which enhancer construct to
%   fit.

%% Variable inputs
%Default settings
FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';

fileDir = pwd;
saveLoc = pwd;
numParPools = 16;
% n_burn = 5000;
n_steps = 10^6;
% ratePriorWidth = 50;
AP_start = 20; % [% of embryo length]
AP_end = 45;   % [% of embryo length]
loadPrevious = false;
globalFit = 1;

% Define the model to be used for the parameter inference. This will be
% realized in utilizing different forms and input parameters of
% SumOfSquares function below. % competition, direct repression, quenching,
% and combination are the options.
% Model = 'default'; 
% 
% for i=1:length(varargin)
%     if strcmpi(varargin{i},'fileDir')
%         fileDir = varargin{i+1};
%     end
%     if strcmpi(varargin{i},'saveLoc')
%         saveLoc = varargin{i+1};
%     end
%     if strcmpi(varargin{i},'numParPools')
%         numParPools = varargin{i+1};
%     end
%     if strcmpi(varargin{i},'n_burn')
%         n_burn = varargin{i+1};
%     end
%     if strcmpi(varargin{i},'n_steps')
%         n_steps = varargin{i+1};
%     end

%     if strcmpi(varargin{i},'Model')
%          % one of the categories : "direct", "competition", "quenching", or
%          % combination of all three (or more).
%         Model = varargin{i+1};
%     end
%     if strcmpi(varargin{i},'construct')
%          % one of the enhancer constructs 
%         construct = varargin{i+1};
%     end
%     if strcmpi(varargin{i},'globalFit')
%          % global fit (fitting Runt WT and null simultaneously) 
%         globalFit = 1;
%     end
% end

%% Load dataset
% We need another separate script to process the data for inputs in this
% script. : This is now done in the "preprocess_data_for_MCMC.m" script.
% Ideas : compiledData.mat -> extract the initial slope (Rate) as mean or
% from individual embryos 

preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
data = preprocessedData.data;
% this load a structure named "data"

% Load the input TF data
load([FilePath, filesep, 'TFinput.mat'])

Bicoid = TFinput(:,1);
Runt = TFinput(:,2);
RuntNull = TFinput(:,3);
%% Extracting specific constructs, redefine as "Data"

% Loop through all constructs (data), to compare the constructName with our
% input, "construct", then extract the cells in that row.
for i=1:length(data)
    if strcmpi(data(i).constructName, construct)
        Data = data(i);
    end
end

% APbinRanges = {{[20,42.5]},{20,45},{20,45}}
% k=1; % counter
% for construct = [2,5,6]

%% Pull the construct of our interest (for the parameter inference)

% Choose a construct 
construct = 5;
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

%% Set up data to be saved in structure
MCMCchain = struct('Kb_chain',{},'Kr_chain',{},'w_b_chain',{},'w_bp_chain',{},...
    'w_rp_chain',{},'p_chain',{},'R_max_chain',{});
MCMCresults = struct('mean_Kb',{},'sigma_Kb',{},'mean_Kr',{},'sigma_Kr',{},...
    'mean_w_b',{},'sigma_w_b',{},'mean_w_bp',{},'sigma_w_bp',{},...
    'mean_w_rp',{},'sigma_w_rp',{},...
    'mean_p',{},'sigma_p',{},'mean_R_max',{},'sigma_R_max',{});
MCMCplot = struct('Rate_plot',{},'simRate',{});

% DatasetName = data_all(k).data(1).name; %Assuming all the cells come from the same dataset


%% MCMC analysis loop through single embryos

%% Setup MCMC Fit

% Define the anonymous functions for the sum of squares function and model.
% We will define different SumOfSquared functions for each 
if strcmpi(Model, 'direct')    
    ssfun = @(params,MCMCdata) SumofSquaresFunction_DirectRepression_V1(MCMCdata, params);
elseif strcmpi(Model, 'competition')   
    ssfun = @(params,MCMCdata) SumofSquaresFunction_Competition_V1(MCMCdata, params);
elseif strcmpi(Model, 'quenching')   
    ssfun = @(params,MCMCdata) SumofSquaresFunction_Quenching_AllBcdBound_V1(MCMCdata, params);
else
    ssfun = @(params,MCMCdata) SumofSquaresFunction_InitialSlopeMCMC(MCMCdata, params);
end


% Initialize MCMC parameters.

Kb0 = 100;   % 100*rand;
Kr0 = 10;     % 100*rand;
w_b0 = 1;    % 10*rand;
w_bp0 = 5;   % 10*rand;
p0 = 0.001; %
R_max0 = 500; %500*rand;

% repression (0< w <1)
% w_br0 = rand;
w_rp0 = 0.2;
% w_brp0 = rand;

params0 = [Kb0, Kr0, w_b0, w_bp0, w_rp0, p0, R_max0, ];

sigma2_0 = 1; %Initial guess for error variance

%This is the variance in the parameter proposal distribution. Change these
%to change the proposal acceptance rate, or if the convergence doesn't look
%good.

Kb_step = 1;
Kr_step = 1;
w_b_step = 0.1;
w_bp_step = 0.1;
p_step = 0.0001;
R_max_step = 1;

% repression
% w_br_step = 0.01;
w_rp_step = 0.01;
% w_brp_step = 0.01;

% Initialize the covariance matrix
J0 = diag([Kb_step, Kr_step, w_b_step, w_bp_step,...
    p_step, R_max_step, w_rp_step]);

%% Setup MCMC parameters and options
% put the initial parameters and bounds in a form that the mcmc function
% accepts

% pri_mu = NaN; %default prior gaussian mean
% pri_sig = Inf; %default prior gaussian variance
% localflag = 0; %is this local to this dataset or shared amongst batches?

% PARAMS  theta structure
%   {  {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      {'par2',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      ... }

% params with a weak prior
% params = {  {'Kb', params0(1), 10^(-2), 10^5, 10^3, 10^4}
%             {'Kr', params0(2),  10^(-2), 10^5, 10^2, 10^3}
%             {'w_b', params0(3),  1, 10^2, 5, 10}
%             {'w_bp', params0(4), 1, 10^2, 5, 10}
%             {'p', params0(5), 0, 0.1}
%             {'R_max', params0(6), 0, 1000, 500, 200}
%             {'w_rp', params0(7), 0, 1}
%             };

% params with a uniform prior
params = {  {'Kb', params0(1), 10^(-1), 10^4}
            {'Kr', params0(2),  10^(-1), 10^4}
            {'w_b', params0(3),  1, 10^3}
            {'w_bp', params0(4), 1, 10^3}
            {'w_rp', params0(5), 0, 1}
            {'R_max', params0(6), 0, 1000}
            {'p', params0(7), 0, 0.1}
            };


model = [];
model.ssfun = ssfun;
model.sigma2 = sigma2_0;
model.N = length(MCMCdata.ydata);

%MCMC options
options = [];
options.nsimu = n_steps; %Number of steps
options.updatesigma = 1; %Update error variance
options.qcov = J0; %Initial covariance
% options.burnintime = n_burn; %Burn in time
options.adaptint = 100;
options.method = 'dram';
options.verbosity = 0; %Decrease text output

% AR options
options.drscale = 5; % a high value (5) is important for multimodal parameter spaces
% options.waitbar = wb; %the waitbar is rate limiting sometimes
% options.nsimu = nSimu; %should be between 1E3 and 1E6
% options.updatesigma = 1; %honestly don't know what this does



%Run the MCMC
[results,chain,s2chain] = mcmcrun(model,MCMCdata,params,options);



%% Extract chain results into individual parameters
n_burn = 1;
Kb_chain = chain(n_burn:end,1);
Kr_chain = chain(n_burn:end,2);
w_b_chain = chain(n_burn:end,3);
w_bp_chain = chain(n_burn:end,4);
w_rp_chain = chain(n_burn:end,5);
R_max_chain = chain(n_burn:end,6);
p_chain = chain(n_burn:end,7);


%Mean/STD results
mean_Kb = mean(Kb_chain);
sigma_Kb = std(Kb_chain,1);
mean_Kr = mean(Kr_chain);
sigma_Kr = std(Kr_chain,1);
mean_w_b = mean(w_b_chain);
sigma_w_b = std(w_b_chain,1);
mean_w_bp = mean(w_bp_chain);
sigma_w_bp = std(w_bp_chain,1);
mean_p = mean(p_chain);
sigma_p = std(p_chain,1);
mean_R_max = mean(R_max_chain);
sigma_R_max = std(R_max_chain,1);
mean_w_rp = mean(w_rp_chain);
sigma_w_rp = std(w_rp_chain,1);

mean_sigma = sqrt(mean(s2chain));
sigma_sigma = std(sqrt(s2chain),1);



%% Part2. Factor out for a separate script : 

%% Simulated fluorescences of best fit
% define the parameters based on which model we ended up using
if strcmpi(Model, 'direct')    
    fit_model = @(Bcd, Run, params) model_6A1R_direct_repression_V1(Bcd, Run, params);
    params_fitted = [mean_Kb mean_Kr mean_w_b mean_w_bp mean_w_rp mean_p mean_R_max];
elseif strcmpi(Model, 'competition')   
    fit_model = @(Bcd, Run, params) model_6A1R_competition_V1(Bcd, Run, params);
elseif strcmpi(Model, 'quenching')   
    fit_model = @(Bcd, Run, params) model_6A1R_quenching_V1(Bcd, Run, params);
else
    fit_model = @(Bcd, Run, params) model_6A1R_combination_all_V1(Bcd, Run, params);
end

Rate_fit = fit_model([Bcd Bcd], [RunNull Run], params_fitted)
% 
% MCMCplot = Rate_fit;

% Plot to check
figure
hold on
plot(Data.APbins, Data.Rate_null)
plot(Data.APbins, Data.Rate_WT)

plot(APaxis, Rate_fit(:,1))
plot(APaxis, Rate_fit(:,2))


%% Diagnose if the inference was done good
chainstats(chain,results)

%% generate plots with error bars of the output (similar to the nlpredci)
mdl = @(params, TF) model_6A1R_direct_repression_V2;
xdata = [MCMCdata.Bcd, MCMCdata.Run];
out = mcmcpred(results,chain,[],xdata, mdl);

%% generate cornerplots
mcmcplot(chain,[],results,'pairs');

%% generate plots of inferred parameters
params_inferred = [mean_Kb mean_Kr mean_w_b mean_w_bp mean_w_rp mean_p mean_R_max];

params_inferred_STD = [sigma_Kb sigma_Kr sigma_w_b sigma_w_bp sigma_w_rp sigma_p sigma_R_max];

hold on
construct = 6;
errorbar(params_inferred,params_inferred_STD,'o','Color', ColorChoice(construct,:),'MarkerFaceColor', ColorChoice(construct,:))

% errorbar(params_fit_100, STD_100,'o','Color', ColorChoice(2,:),'MarkerFaceColor', ColorChoice(2,:))
% errorbar(params_fit_010, STD_010,'o','Color', ColorChoice(6,:),'MarkerFaceColor', ColorChoice(6,:))
% errorbar(params_fit_001, STD_001,'o','Color', ColorChoice(5,:),'MarkerFaceColor', ColorChoice(5,:))

set(gca, 'YScale','log')
% legend('100','001','010')
legend('001')

xlim([0 8])
xticklabels({'','K_{b}','K_{r}','\omega_{b}','\omega_{bp}','\omega_{rp}','p','R_{max}',''})
xlabel('parameters')
ylabel('fitted values')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'MCMC_Params_direct_001','.tif']); 
saveas(gcf,[FigPath,filesep,'MCMC_Params_direct_001','.pdf']); 

%% (Optional) generate Pearson's correleation coefficient plots
rho = @(x, y) x ./ (y'*y); %pearson's correlation coefficient

imagesc(rho(results.cov, sqrt(diag(results.cov))));
colorbar;
ylabel('parameter 1')
xlabel('parameter 2')
title('Correlation coefficient');
colormap(viridis);
    
%% Save the data
MCMCchain = struct;
MCMCchain.Kb_chain = Kb_chain;
MCMCchain.Kr_chain = Kr_chain;
MCMCchain.w_b_chain = w_b_chain;
MCMCchain.w_bp_chain = w_bp_chain;
MCMCchain.p_chain = p_chain;
MCMCchain.R_max_chain = R_max_chain;
% MCMCchain.w_br_chain = w_br_chain;
MCMCchain.w_rp_chain = w_rp_chain;
% MCMCchain.w_brp_chain = w_brp_chain;

MCMCresults = struct;
% Mean/STD values of inferred parameters
MCMCresults.mean_Kb = mean_Kb;
MCMCresults.sigma_Kb = sigma_Kb;
MCMCresults.mean_Kr = mean_Kr;
MCMCresults.sigma_Kr = sigma_Kr;
MCMCresults.mean_w_b = mean_w_b;
MCMCresults.sigma_w_b = sigma_w_b;
MCMCresults.mean_w_bp = mean_w_bp;
MCMCresults.sigma_w_bp = sigma_w_bp;
MCMCresults.mean_p = mean_p;
MCMCresults.sigma_p = sigma_p;
MCMCresults.mean_R_max = mean_R_max;
MCMCresults.sigma_R_max = sigma_R_max;
% MCMCresults.mean_w_br = mean_w_br;
% MCMCresults.sigma_w_br = sigma_w_br;
MCMCresults.mean_w_rp = mean_w_rp;
MCMCresults.sigma_w_rp = sigma_w_rp;
% MCMCresults.mean_w_brp = mean_w_brp;
% MCMCresults.sigma_w_brp = sigma_w_brp;

%% Save data into .mat structure
%MCMC results and plots
filename = [date,'-',DatasetName];
save([saveLoc,'\',filename,'.mat'],'MCMCresults','MCMCplot','DatasetName');

%MCMC raw chains
filename = [date,'-',DatasetName,'_RawChain'];
save([saveLoc,'\',filename,'.mat'],'MCMCchain');
% end

close(w);

disp(['MCMC analysis complete. Information stored in: ',saveLoc]);
end