function Fig7_MCMC_HillV3_direct_6A3R_RunWT_fixed_Kr_withCoop_V1
%% Description
% This script is doing MCMC fit for 6A3R, Run WT data one by one, to
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
% load the MCMC results for the 6A1R constructs for 1) fixed K_r (global, shared) 
% and 2)w_rp (local) for each constructs.
% load : "MCMC_6A1R_RuntWT.mat"
% Note that the chain is constructed as "K_r(shared) w_rp1 w_rp2 w_rp3"
tempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A1R_Bcd_params_fromRuntNulls_fixed_Kr';
load([tempPath2, filesep, 'MCMC_6A1R_RuntWT_params.mat'])

% extract the Runt-dependent parameters
K_r = MCMC_6A1R_RuntWT.params_inferred(1);
w_rp1 = MCMC_6A1R_RuntWT.params_inferred(2); %[100]
w_rp2 = MCMC_6A1R_RuntWT.params_inferred(3); %[001]
w_rp3 = MCMC_6A1R_RuntWT.params_inferred(4); %[010]

% extract the std for the future usage
params_Run = MCMC_6A1R_RuntWT.params_inferred;
params_Run_std = MCMC_6A1R_RuntWT.params_inferred_sigma;


% Second, the inferred parameters from the hbP2 + 2 Runt sites
% load : MCMC_6A2R_RuntWT.mat
% tempPath3 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\both_Run_Run_higher_order\10 limit';
% load([tempPath3, filesep, 'MCMC_6A2R_RuntWT_params.mat'])

% extract the Runt-Runt cooperativity, and higher-order cooperativity
% w_rr1 = MCMC_6A2R_RuntWT(1).params_inferred(1); % [011]
% w_ho1 = MCMC_6A2R_RuntWT(1).params_inferred(2); % [011]
% w_rr2 = MCMC_6A2R_RuntWT(2).params_inferred(1); % [110]
% w_ho2 = MCMC_6A2R_RuntWT(2).params_inferred(2); % [110]
% w_rr3 = MCMC_6A2R_RuntWT(3).params_inferred(1); % [101]
% w_ho3 = MCMC_6A2R_RuntWT(3).params_inferred(2); % [101]


% Third, hybrid form of the parameters, w_ho only for [011], [101], w_rr & w_ho for [110]
tempPath4 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\hybrid_w_ho_only_w_rr_w_ho_only_[110]';
load([tempPath4, filesep, 'MCMC_6A2R_hybrid_w_rr_w_ho.mat']); 
% this will load two matrices, "MCMC_6A2R_params", and
% "MCMC_6A2R_params_std".


% extract the Runt-Runt cooperativity, and higher-order cooperativity
w_rr1 = MCMC_6A2R_params(1,1); % [011]
w_ho1 = MCMC_6A2R_params(1,2); % [011]
w_rr2 = MCMC_6A2R_params(2,1); % [110]
w_ho2 = MCMC_6A2R_params(2,2); % [110]
w_rr3 = MCMC_6A2R_params(3,1); % [101]
w_ho3 = MCMC_6A2R_params(3,2); % [101]



%% First, let's try to see how the prediction looks like with only pair-wise parameters
% Using the model_6A3R_HillModel_V3_direct_HigherCoop
% higher-order coop when all three Runt molecules are bound as well as the
% RNAP : w_hoho
w_hoho = 1;

% take the Bcd/RNAP params for the [111] construct
params_Bcd = MCMC_6A0R_RuntNulls(4).params_inferred;
params_Run = MCMC_6A1R_RuntWT.params_inferred;
params_temp = [params_Bcd, params_Run, w_rr1, w_rr2, w_rr3, w_ho1, w_ho2, w_ho3, w_hoho];
% params_temp = [params_Bcd, params_Run, 1, 6.88, 1, 8.61, 0.0468, 0.1964, w_hoho];

output = model_6A3R_HillModel_V3_direct_HigherCoop(params_temp, TFinput);
fit_nulls = model_6A0R_HillModel_V3(params_Bcd, TF_null);

        
% generate raw fits 
APaxis = 0:0.025:1;
APbin_start = 9;
APbin_end = 21;

Bcd = Bicoid(APbin_start:APbin_end); 
Run = Runt(APbin_start:APbin_end); 
RunNull = RuntNull(APbin_start:APbin_end);

TF = [Bcd, Run];
TF_null = [Bcd, RunNull];


% define the construct index (which is consistent with the way it's
% defiend in the compiledData.mat)

construct = 4;
    
% Runt nulls
fit_nulls = model_6A0R_HillModel_V3(params_Bcd, TF_null);

% data (Runt nulls)
Rate_null = compiledData{construct+1+8,9};
Rate_null_SEM = compiledData{construct+1+8,10};

% data (WT)
Rate_WT = compiledData{construct+1,9};
Rate_WT_SEM = compiledData{construct+1,10};

% find the APbins that have the data points
NoNaN_index_WT = ~isnan(Rate_WT);
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
APbins_fit = 9:14; %APbinRange;

clf
hold on
errorbar(APaxis, Rate_null, Rate_null_SEM, 'o', 'Color', ColorChoice(4,:),'LineWidth', 1, 'MarkerFaceColor',ColorChoice(4,:) )%, 'MarkerSize',7)    
errorbar(APaxis, Rate_WT, Rate_WT_SEM, 'o', 'Color', ColorChoice(1,:),'LineWidth', 1, 'MarkerFaceColor',ColorChoice(1,:)) %, 'MarkerSize',7)

% plot(APaxis(APbin_start:APbin_end), fit_nulls, 'Color', ColorChoice(4,:),'LineWidth', 2)
% plot(APaxis(APbin_start:APbin_end), output, 'Color', ColorChoice(1,:),'LineWidth', 2)

plot(APaxis(APbin_start:APbin_end), fit_nulls, 'Color', ColorChoice(4,:),'LineWidth', 2)
% plot(APaxis(APbin_start:APbin_end), output, 'Color', ColorChoice(1,:),'LineWidth', 2)
plot(APaxis(APbins_fit), output(APbins_fit), 'Color', ColorChoice(1,:),'LineWidth', 2)

xlim([0.2 0.5])
xticks([0.2 0.3 0.4 0.5])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel('initial rate (AU/min)')

box on
legend('data(null)','data(WT)','Fit (null)', 'Fit (WT)')
StandardFigure(gcf,gca)
% FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\6A3R_prediction_w_rr_w_ho_pairwiseOnly';
% saveas(gcf,[FigPath,filesep,'raw_fits_null_WT_w_hoho_300', constructNames{construct}  ,'.tif']); 
% saveas(gcf,[FigPath,filesep,'raw_fits_null_WT_w_hoho_300', constructNames{construct} ,'.pdf']); 


%% Set up the MCMC inference
% Note that we don't have any parameters regarding the repressor, thus we
% will just use the simpler form of the model, "model_6A0R_HillModel_V3.m"

% index for the constructs
constructIndex = [4];
% AP axis
APaxis = 0:0.025:1;
APbin_start=9; % 20% of embryo length
APbin_end=21; % 50% of embryo length
Bcd = Bicoid(APbin_start:APbin_end); 
Run = Runt(APbin_start:APbin_end); 
RunNull = RuntNull(APbin_start:APbin_end);

TF = [Bcd, Run];
TF_null = [Bcd, RunNull];

%% Pick a model for the MCMC inference
% take the Bcd/RNAP params for the [111] construct
params_Bcd = MCMC_6A0R_RuntNulls(4).params_inferred;
params_Run = MCMC_6A1R_RuntWT.params_inferred;
params_temp = [params_Bcd, params_Run, w_rr1, w_rr2, w_rr3, w_ho1, w_ho2, w_ho3, w_hoho];
% params_temp = [params_Bcd, params_Run, 1, 6.88, 1, 8.61, 0.0468, 0.1964, w_hoho];
params_fixed = [params_Bcd, params_Run, w_rr1, w_rr2, w_rr3, w_ho1, w_ho2, w_ho3];

% output = model_6A3R_HillModel_V3_direct_HigherCoop(params_temp, TFinput);

mdl = @(TF, w_hoho, params_fixed) model_6A3R_HillModel_V3_direct_HigherCoop([params_fixed, w_hoho], TF);
model.modelfun = mdl;

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, InputData) sum((InputData.ydata-mdl(InputData.xdata, params, params_fixed)).^2);
% model.ssfun = @model_6A2R_HillV3_direct_batch_SS;
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

% params_fixed = [Kb, w_bp, p, R_max, K_r, w_rp1, w_rp2];
% params = [omega_rr, omega_ho];

% initialize the structure to save the MCMC result
MCMC_6A3R_RuntWT = struct;
m = 1; % counter to save into the structure, as we're only dealing with 3 constructs (6A1R)

% index for the constructs
constructIndex = 4;

for index=1 % [111]
    %% Pull the construct of our interest (for the parameter inference)
    % Choose a construct 
    construct = constructIndex(index);
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
    MCMCdata.APdata = [APbins'];
    MCMCdata.ydata = [Rate_WT_forFit];
    % input TF
    MCMCdata.Bcd = [Bcd];
    MCMCdata.Run = [Run];
    MCMCdata.xdata = [Bcd , Run];

%     MCMCdata{index}.R_max = max(Rate_WT_forFit);
%     MCMCdata{index}.R_min = min(Rate_WT_forFit);
    
    


    MCMCdata.params_fixed = [params_Bcd, params_Run, w_rr1, w_rr2, w_rr3, w_ho1, w_ho2, w_ho3];
    
    %% Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'w_{hoho}'};
    params = cell(1, length(names));


    params0 = [1,1];
    % params0_std

    % Bounds of the parameters
    LB = [0];
%     UB = [10^(1), 10^(2)];
    UB = [1000];

    for i = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?

        % For all parameters, define the priors from the posteriors from
        % the previous rounds of MCMC (either on Runt nulls or 6A1R cases)
%         if i==2 % w_ho
%             targetflag = 0; % fix the parameter
%         end
%         if i==1 % w_rr
%             targetflag = 0; % fix the parameter
%         end
        
        params{1, i} = {names{i}, params0(i), LB(i), UB(i), pri_mu, pri_sig, targetflag, localflag};

    end
    
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

    %% Extract chain results into individual parameters
    params_inferred = []; 
    params_inferred_sigma = [];
    n_burn = 0.5*n_steps;

    for k=1
        params_inferred(1,k) = mean(chain(n_burn+1:end,k));
        params_inferred_sigma(1,k) = std(chain(n_burn+1:end,k));
    end

    %% Save the MCMC result into a structure for future usage
    MCMC_6A3R_RuntWT(index).name = 'w_hoho';
    MCMC_6A3R_RuntWT(index).results = results;
    MCMC_6A3R_RuntWT(index).chain = chain;
    MCMC_6A3R_RuntWT(index).s2chain = s2chain;

    % inferred parameters
    MCMC_6A3R_RuntWT(index).params_inferred = params_inferred;
    MCMC_6A3R_RuntWT(index).params_inferred_sigma = params_inferred_sigma;

end 
    

%% Calculate the confidence interval for MCMC fitting
data_pred = MCMCdata; % not used
%params_Run = [Kr, w_rp];
params_fixed;
% TF was defined above (outside of this for loop)
model_mcmc = @(w_hoho) model_6A3R_HillModel_V3_direct_HigherCoop([params_fixed, w_hoho], TF);

MCMCresults = results;
MCMCresults.parind = [1];
MCMCresults.local = [1];
MCMCresults.theta = [results.theta(1)];
MCMCresults.nbatch = 1;
n_burn= 10000; % number of burn-in steps in MCMC (50% of the total chain)
predout = mcmcpred_V2(MCMCresults, chain(n_burn+1:end,:), [], data_pred, model_mcmc);

%% generate raw fits 
APaxis = 0:0.025:1;
APbin_start = 9;
APbin_end = 21;

Bcd = Bicoid(APbin_start:APbin_end); 
Run = Runt(APbin_start:APbin_end); 
RunNull = RuntNull(APbin_start:APbin_end);

TF = [Bcd, Run];
TF_null = [Bcd, RunNull];

FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A3R_MCMC_fits_w_hoho';
for index = 1
    % define the construct index (which is consistent with the way it's
    % defiend in the compiledData.mat)
    constructIndex = [3,7,8];
    construct = constructIndex(index);

    % grab the inferred parameters from this MCMC
    w_hoho = MCMC_6A3R_RuntWT.params_inferred;

    output = model_6A3R_HillModel_V3_direct_HigherCoop([params_fixed, w_hoho], TF);
    % Runt nulls
    fit_nulls = model_6A0R_HillModel_V3(params_Bcd, TF);
    
    % Rate for both Runt null and Runt WT
    Rate_null = data(construct).Rate_null;
    [~,num_samples] = size(data(construct).Rate_null_individual);
    Rate_null_SEM = zeros(41,1);
    Rate_null_SEM = nanstd(data(construct).Rate_null_individual,0,2)./sqrt(num_samples);

    Rate_WT = data(construct).Rate_WT;
    [~,num_samples] = size(data(construct).Rate_WT_individual);
    Rate_WT_SEM = zeros(41,1);
    Rate_WT_SEM = nanstd(data(construct).Rate_WT_individual,0,2)./sqrt(num_samples);
    
%     % data (Runt nulls)
%     Rate_null = compiledData{construct+1+8,9};
%     Rate_null_SEM = compiledData{construct+1+8,10};
%     
%     % data (WT)
%     Rate_WT = compiledData{construct+1,9};
%     Rate_WT_SEM = compiledData{construct+1,10};
    
    clf
    hold on
    errorbar(APaxis, Rate_null, Rate_null_SEM, 'o', 'Color', ColorChoice(4,:),'LineWidth', 1, 'MarkerFaceColor', ColorChoice(4,:))    
    errorbar(APaxis, Rate_WT, Rate_WT_SEM, 'o', 'Color', ColorChoice(1,:),'LineWidth', 1, 'MarkerFaceColor', ColorChoice(1,:))
    
    plot(APaxis(APbin_start:APbin_end), fit_nulls, 'Color', ColorChoice(4,:),'LineWidth', 2)
%     plot(APaxis(APbin_start:APbin_end), output, 'Color', ColorChoice(1,:),'LineWidth', 2)
    % Run WT MCMC fit with 95% confidence interval
    out = predout;
    nn = (size(out.predlims{1}{1},1) + 1) / 2;
    plimi = out.predlims{1}{1};
    yl = plimi(3,:);
    yu = plimi(2*nn-3,:);
%     plot(APbins, output)
    %plot(APaxis(APbin_start:APbin_end), output)
    shadedErrorBar(APaxis(APbin_start:APbin_end), output, [yu-transpose(output);transpose(output)-yl],'lineProps','-k')%['color', ColorChoice(1,:),'LineWidth', 2])

    xlim([0.2 0.5])
    xticks([0.2 0.3 0.4 0.5])
    ylim([0 400])
    yticks([0 100 200 300 400])
    
    xlabel('embryo length')
    ylabel('initial rate (AU/min)')
    
    box on
    legend('data(null)','data(WT)','Fit (null)', 'Fit (WT)')
    StandardFigure(gcf,gca)
    %pause
    saveas(gcf,[FigPath,filesep,'raw_fits_null_WT_95%CI_', constructNames{construct}  ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'raw_fits_null_WT_95%CI_', constructNames{construct} ,'.pdf']); 
end

%% generate the raw fits for only Runt WT, also three constructs combined into one plot

APaxis = 0:0.025:1;
APbin_start = 9;
APbin_end = 21;

Bcd = Bicoid(APbin_start:APbin_end); 
Run = Runt(APbin_start:APbin_end); 
RunNull = RuntNull(APbin_start:APbin_end);

TF = [Bcd, Run];
TF_null = [Bcd, RunNull];

hold on
for index = 1
    % define the construct index (which is consistent with the way it's
    % defiend in the compiledData.mat)
    construct = constructIndex(index);
    
    % extract parameters from the MCMC results (we're pulling the Bcd
    % dependent parameters and Run dependent parameters respectively from
    % their own inferences. (Bcd parameters from the Run null, and Run
    % parameters from the Run WT).
%     params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
%     params_Run = MCMC_6A1R_RuntWT(index).params_inferred;
    
    
    % grab the inferred parameters from this MCMC
    w_hoho = MCMC_6A3R_RuntWT(index).params_inferred;
   
    output = model_6A3R_HillModel_V3_direct_HigherCoop(params_MCMC, TF, params_fixed);

    % data (WT)
    Rate_WT = compiledData{construct+1,9};
    Rate_WT_SEM = compiledData{construct+1,10};
    
    errorbar(APaxis, Rate_WT, Rate_WT_SEM, 'o', 'Color', ColorChoice(construct,:),'LineWidth', 1)
    
    plot(APaxis(APbin_start:APbin_end), output, 'Color', ColorChoice(construct,:),'LineWidth', 2)

    
end

% figure format
xlim([0.2 0.5])
xticks([0.2 0.3 0.4 0.5])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel('initial rate (AU/min)')

box on
legend('011','fit','110','fit','101','fit')
StandardFigure(gcf,gca)


saveas(gcf,[FigPath,filesep,'raw_fits_RuntWT_6A2R_compiled','.tif']); 
saveas(gcf,[FigPath,filesep,'raw_fits_RuntWT_6A2R_compiled','.pdf']); 
%% generate corner plots

for index = 1:3
    construct = constructIndex(index);
    chain = MCMC_6A3R_RuntWT(index).chain;
    n_burn = 0.5*n_steps;
    m = [chain(n_burn+1:end,1), chain(n_burn+1:end,2)];
    corner = figure;
    %     names = {'K_{b}','\omega_{bp}','p','R_{max}'};
    names = {'\omega_{rr}', '\omega_{ho}'};
    ecornerplot(m,'names',names);
%     histogram(chain(n_burn+1:end),50,'Normalization','probability')

%     saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.tif']); 
%     saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']);
    exportgraphics(gcf,[FigPath, filesep,'Corner_plot_highres_', constructNames{construct},'.pdf'],'ContentType','vector')
end

end