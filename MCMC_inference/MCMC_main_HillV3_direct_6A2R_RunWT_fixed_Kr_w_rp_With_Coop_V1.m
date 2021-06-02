function MCMC_main_HillV3_direct_6A2R_RunWT_fixed_Kr_w_rp_With_Coop_V1
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
% load the MCMC results for the 6A1R constructs for 1) fixed K_r (global, shared) 
% and 2)w_rp (local) for each constructs.
% loads : "MCMC_6A1R_RuntWT"
% Note that the chain is constructed as "K_r(shared) w_rp1 w_rp2 w_rp3"
tempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\Point_estimate_Bcd_params_fromRunNulls_fixed_K_r';
load([tempPath2, filesep, 'MCMC_6A1R_RuntWT_params.mat'])

% extract the Runt-dependent parameters
K_r = MCMC_6A1R_RuntWT.params_inferred(1);
w_rp1 = MCMC_6A1R_RuntWT.params_inferred(2); %[100]
w_rp2 = MCMC_6A1R_RuntWT.params_inferred(3); %[001]
w_rp3 = MCMC_6A1R_RuntWT.params_inferred(4); %[010]

% std
K_r_std = MCMC_6A1R_RuntWT.params_inferred_sigma(1);
w_rp1_std = MCMC_6A1R_RuntWT.params_inferred_sigma(2); %[100]
w_rp2_std = MCMC_6A1R_RuntWT.params_inferred_sigma(3); %[001]
w_rp3_std = MCMC_6A1R_RuntWT.params_inferred_sigma(4); %[010]

% extract the std for the future usage
params_Run = MCMC_6A1R_RuntWT.params_inferred;
params_Run_std = MCMC_6A1R_RuntWT.params_inferred_sigma;

%% Set up the MCMC inference
% Note that we don't have any parameters regarding the repressor, thus we
% will just use the simpler form of the model, "model_6A0R_HillModel_V3.m"

%% Pick a model for the MCMC inference
% params_fixed = [Kb, w_bp, p, R_max, K_r, w_rp1, w_rp2];
% params = [omega_rr, omega_ho];
mdl = @(TF, params, params_fixed) model_6A2R_HillModel_V3_direct_fixed_Kr_Bcd_RNAP_w_rp_params(params, TF, params_fixed);
model.modelfun = mdl;

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, InputData) sum((InputData.ydata-mdl(InputData.xdata, params, InputData.params_fixed)).^2);
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
MCMC_6A2R_RuntWT = struct;
m = 1; % counter to save into the structure, as we're only dealing with 3 constructs (6A1R)
    
% r2-new : [011] : 3rd in the MCMC_6A0R_RuntNulls
% r2-close : [110] : 7th
% r2-far : [101] : 8th

% index for the constructs
constructIndex = [3,7,8];

for index=1:3 % [100, 001, 010]
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
    
    % parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    % Runt parameters
    if index==1 % r2-new : [011]
        params_Run = [K_r, w_rp2, w_rp3];
    elseif index==2 % r2-close : [110]
        params_Run = [K_r, w_rp1, w_rp3];
    elseif index==3 % r2-far : [101]
        params_Run = [K_r, w_rp1, w_rp2];
    end

    MCMCdata.params_fixed = [params_Bcd, params_Run];
    
    %% Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'w_{rr}', 'w_{ho}'};
    params = cell(1, length(names));


    params0 = [1,1];
    % params0_std

    % Bounds of the parameters
    LB = [0,0];
%     UB = [10^(1), 10^(2)];
    UB = [3, 20];

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

    for k=1:2
        params_inferred(1,k) = mean(chain(n_burn+1:end,k));
        params_inferred_sigma(1,k) = std(chain(n_burn+1:end,k));
    end

    %% Save the MCMC result into a structure for future usage
    MCMC_6A2R_RuntWT(index).name = 'fixed_Kr_w_rp_coop'
    MCMC_6A2R_RuntWT(index).results = results;
    MCMC_6A2R_RuntWT(index).chain = chain;
    MCMC_6A2R_RuntWT(index).s2chain = s2chain;

    % inferred parameters
    MCMC_6A2R_RuntWT(index).params_inferred = params_inferred;
    MCMC_6A2R_RuntWT(index).params_inferred_sigma = params_inferred_sigma;

end 
    
%% generate raw fits 
APaxis = 0:0.025:1;
APbin_start = 9;
APbin_end = 21;

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
    
    % point estimates from the Runt nulls
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    % MCMC inferred parameters for the Runt WT 
    % : [K_r, w_rp1, w_rp2, w_rp3]; 
    % w_rp1 : [100]
    % w_rp2 : [001]
    % w_rp3 : [010]
    
    % Runt parameters
    if index==1 % r2-new : [011]
        params_Run = [K_r, w_rp2, w_rp3];
    elseif index==2 % r2-close : [110]
        params_Run = [K_r, w_rp1, w_rp3];
    elseif index==3 % r2-far : [101]
        params_Run = [K_r, w_rp1, w_rp2];
    end
    params_fixed = [params_Bcd, params_Run];
    
    % grab the inferred parameters from this MCMC
    params_MCMC = MCMC_6A2R_RuntWT(index).params_inferred;
    
%     % fixing the w_rr = 1;
%     params_MCMC(1) = 1;
    
    % when the w_rr or w_ho is set to be 1 using the targetflag
    if length(params_MCMC)==1 &&  params{1,2}{1,7}==0 % targetflag 
        params_MCMC = [params_MCMC, 1];
    elseif length(params_MCMC)==1 &&  params{1,1}{1,7} ==0
        params_MCMC = [1, params_MCMC];
    end
    
    output = model_6A2R_HillModel_V3_direct_fixed_Kr_Bcd_RNAP_w_rp_params(params_MCMC, TF, params_fixed);
    % Runt nulls
    fit_nulls = model_6A0R_HillModel_V3(params_Bcd, TF);
    
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
%     saveas(gcf,[FigPath,filesep,'raw_fits_null_WT_', constructNames{construct}  ,'.tif']); 
%     saveas(gcf,[FigPath,filesep,'raw_fits_null_WT_', constructNames{construct} ,'.pdf']); 
end

%% generate corner plots

for index = 1:3
    construct = constructIndex(index);
    chain = MCMC_6A2R_RuntWT(index).chain;
    n_burn = 0.5*n_steps;
    m = [chain(n_burn+1:end,1), chain(n_burn+1:end,2)];
    corner = figure;
    %     names = {'K_{b}','\omega_{bp}','p','R_{max}'};
    names = {'\omega_{rr}', '\omega_{ho}'};
    ecornerplot(m,'names',names);
%     histogram(chain(n_burn+1:end),50,'Normalization','probability')

%     saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.tif']); 
%     saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']);
%     exportgraphics(gcf,[FigPath, filesep,'Corner_plot_highres_', constructNames{construct},'.pdf'],'ContentType','vector')
end

%% Error bar using MCMCPred
% out=mcmcpred(results,chain,s2chain,data,modelfun,nsample,varargin)
%% generate histogram for w_rr ONLY case (w_ho = 1)

for index = 1:3
    construct = constructIndex(index);
    chain = MCMC_6A2R_RuntWT(index).chain;
    n_burn = 0.5*n_steps;
    
    chain_temp = chain(n_burn+1:end,1);
    mean_posterior = mean(chain_temp);
    std_posterior = std(chain_temp);
    
    % histogram of posterior
    histogram(chain_temp,50,'Normalization','probability')
    % mean and std of posterior
    xline(mean_posterior,'LineWidth',2)
    xline(mean_posterior-std_posterior,'--','LineWidth',2)
    xline(mean_posterior+std_posterior,'--','LineWidth',2)
    
    ylim([0 0.04])
    
    xlabel('\omega_{rr}')
    ylabel('frequency')
    legend('postrior','mean')

    box on
    StandardFigure(gcf,gca)

    saveas(gcf,[FigPath,filesep,'histogram_', constructNames{construct} ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'histogram_', constructNames{construct} ,'.pdf']); 
end

%% 
%% generate histogram for w_ho ONLY case (w_rr = 1)

for index = 1:3
    construct = constructIndex(index);
    chain = MCMC_6A2R_RuntWT(index).chain;
    n_burn = 0.5*n_steps;
    
    chain_temp = chain(n_burn+1:end,1);
    mean_posterior = mean(chain_temp);
    std_posterior = std(chain_temp);
    
    chain_temp = log10(chain_temp);
    mean_posterior = log10(mean_posterior);
    std_posterior = log10(std_posterior);
    
    % histogram of posterior
    histogram(chain_temp,50,'Normalization','probability')
    % mean and std of posterior
    xline(mean_posterior,'LineWidth',2)
    xline(mean_posterior-std_posterior,'--','LineWidth',2)
    xline(mean_posterior+std_posterior,'--','LineWidth',2)
    
%     ylim([0 0.04])
    
    xlabel('log_{10} \omega_{ho}')
    ylabel('frequency')
    legend('postrior','mean','Location','NorthWest')

    box on
    StandardFigure(gcf,gca)
%     pause

    saveas(gcf,[FigPath,filesep,'histogram_logscale_', constructNames{construct} ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'histogram_logscale_', constructNames{construct} ,'.pdf']); 
end

%% generate plots of inferred parameters (w_rr and w_ho) for all 6A2R constructs (log scale)

% parse the parameters
for i=1:3
    params_inferred(i,:) = MCMC_6A2R_RuntWT(i).params_inferred;
    params_inferred_std(i,:) = MCMC_6A2R_RuntWT(i).params_inferred_sigma;
end
hold on
errorbar([1;2], params_inferred(1,:), params_inferred_std(1,:),'o','LineWidth',2,'Color',ColorChoice(3,:))
errorbar([1;2], params_inferred(2,:), params_inferred_std(2,:),'o','LineWidth',2,'Color',ColorChoice(7,:))
errorbar([1;2], params_inferred(3,:), params_inferred_std(3,:),'o','LineWidth',2,'Color',ColorChoice(8,:))

xlim([0 3])
ylim([0 100])
xticks([1,2])
xticklabels({'\omega_{rr}', '\omega_{ho}'})
% yticks([0 100 200 300 400])
xlabel('parameters')
ylabel('inferred values')
legend('011','110','101','Location', 'NorthWest')

set(gca,'YScale','Log')
box on
StandardFigure(gcf,gca)

% save the plot
% FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\both_Run_Run_higher_order\10 limit';
% saveas(gcf,[FigPath,filesep, 'inferred_params_w_rr_w_ho_logScale','.tif']);
% saveas(gcf,[FigPath,filesep, 'inferred_params_w_rr_w_ho_logScale','.pdf']);


%% generate plots of inferred parameters (compare the w_rp values from 6A1R or 6A2R)
hold on

% Runt-dependent parameters from the MCMC-6A1R-HillV3, direct
params_Run;
params_Run_std;
errorbar(1:3, params_Run(2:4), params_inferred_std(2:4),'o','LineWidth',2,'Color',ColorChoice(1,:))

% params inferred from 6A2R, HillV3, direct : [w_rp1, w_rp2, w_rp3]
params_inferred = MCMC_6A2R_RuntWT.params_inferred;
params_inferred_std = MCMC_6A2R_RuntWT.params_inferred_sigma;
errorbar(1:3, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(4,:))


xlim([0 4])
xticks([1,2,3])
% xticklabels({'\omega_{rp1}', '\omega_{rp2}', '\omega_{rp3}'})
xticklabels({'[100]','[001]','[010]'})
% yticks([0 100 200 300 400])
xlabel('Runt binding site')
ylabel('inferred \omega_{rp}')
% legend('100','001','010','Location', 'NorthWest')
legend('1 Runt site','2 Runt sites')

box on
StandardFigure(gcf,gca)

% % save the plots
% saveas(gcf,[FigPath,filesep, 'inferred_params_comparison_6A1R_6A2R','.tif']);
% saveas(gcf,[FigPath,filesep, 'inferred_params_comparison_6A1R_6A2R','.pdf']);

%% generate plots of inferred parameters (w_rr and w_ho) for all 6A2R constructs

% parse the parameters
for i=1:3
    params_inferred(i,:) = MCMC_6A2R_RuntWT(i).params_inferred;
    params_inferred_std(i,:) = MCMC_6A2R_RuntWT(i).params_inferred_sigma;
end
hold on
errorbar([1;2], params_inferred(1,:), params_inferred_std(1,:),'o','LineWidth',2,'Color',ColorChoice(3,:))
errorbar([1;2], params_inferred(2,:), params_inferred_std(2,:),'o','LineWidth',2,'Color',ColorChoice(7,:))
errorbar([1;2], params_inferred(3,:), params_inferred_std(3,:),'o','LineWidth',2,'Color',ColorChoice(8,:))

xlim([0 3])
ylim([0 15])
xticks([1,2])
yticks([0 5 10 15 ])
xticklabels({'\omega_{rr}', '\omega_{ho}'})
% yticks([0 100 200 300 400])
xlabel('parameters')
ylabel('inferred values')
legend('011','110','101','Location', 'NorthWest')

% set(gca,'YScale','Log')
box on
StandardFigure(gcf,gca)

% save the plot
% FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\both_Run_Run_higher_order\10 limit';
% saveas(gcf,[FigPath,filesep, 'inferred_params_w_rr_w_ho','.tif']);
% saveas(gcf,[FigPath,filesep, 'inferred_params_w_rr_w_ho','.pdf']);

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

%% generate the plot of inferred parameters

%% save the result into mat files (MCMC_6A0R_RuntNulls)
FilePath = FigPath;
save([FilePath, filesep, 'MCMC_6A2R_RuntWT_params.mat'],'MCMC_6A2R_RuntWT')

end