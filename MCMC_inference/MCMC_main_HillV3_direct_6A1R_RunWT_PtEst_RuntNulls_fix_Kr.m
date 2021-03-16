function MCMC_main_HillV3_direct_6A1R_RunWT_PtEst_RuntNulls_fix_Kr
%% Description
% This script is doing MCMC fit for 6A1R, Run WT data one by one, to
% compare the Run-dependent parameters for the Hill.V3 model, to see what varies across
% constructs.

% One important caveat here is that we'll use the "point estimates" from
% the previous MCMC inference using 6A0R_HillModelV3 for each construct
% with Runt null (w/o Run protein).

% Here, we will fix the K_r value across constructs 

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
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct';

%% Load the 6AnR Run null MCMC result (to plug into [Kb, w_bp, p, R_max]
% tempPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\w_bp_0-100_WeakPrior_w_bp';
tempPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\w_bp_0-100';
load([tempPath, filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat'],'MCMC_6A0R_RuntNulls')
%% Set up the MCMC inference
% Note that we don't have any parameters regarding the repressor, thus we
% will just use the simpler form of the model, "model_6A0R_HillModel_V3.m"

%% MCMC Default settings (also make it optional)
options = [];
n_steps = 2*10^4;
options.nsimu = n_steps; %n_steps; %Number of steps
options.updatesigma = 1; %Update error variance
% options.qcov = J0; %Initial covariance
% options.burnintime = 0.5*n_steps; %Burn in time
options.adaptint = 100;
options.method = 'dram';
options.verbosity = 0; %Decrease text output

%% Pick a model for the MCMC inference
% modified model with pre-fixed parameters for Bcd/RNAP
mdl = @(TF, params, params_fixed) model_6A1R_HillModel_V3_direct_fixedBcdRNAPparams(params, TF, params_fixed);
model.modelfun = mdl;

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, InputData) sum((InputData.ydata-mdl(InputData.xdata, params, InputData.params_fixed)).^2);
%% Loop over all constructs to perform the MCMC inference on the Bcd-dependent parameters
% [Kb, w_bp, p, R_max];

% initialize a structure for compiled data as an input for MCMC
% MCMCdata = struct;
% initialize a structure to save the MCMC result
MCMC_6A1R_RuntWT = struct;

% index for the constructs
constructIndex = [2,5,6];

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
    MCMCdata{index}.APdata = [APbins'];
    MCMCdata{index}.ydata = [Rate_WT_forFit];
    % input TF
    MCMCdata{index}.Bcd = [Bcd];
    MCMCdata{index}.Run = [Run];
    MCMCdata{index}.xdata = [Bcd , Run];

    MCMCdata{index}.R_max = max(Rate_WT_forFit);
    MCMCdata{index}.R_min = min(Rate_WT_forFit);
    
    % parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    MCMCdata{index}.params_fixed = params_Bcd;
end 

%% MCMC settings
%% Define the parameters for the model
% put the initial parameters and bounds in a form that the mcmc function
% accepts
names = {'K_{r}', 'w_rp'};
params = cell(1, length(names));

% Pull the fitted parameters from the MCMC inference on the Runt
% null datasets for one-by-one.
% params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
% params_Bcd_std = MCMC_6A0R_RuntNulls(construct).params_inferred_std;

% Initialize MCMC parameters.
% Kb0 = params_Bcd(1);
% w_bp0 = params_Bcd(2);
% p0 = params_Bcd(3);
% R_max0 = params_Bcd(4);

% Initial condition for the Run-dependent parameters
Kr0 = 5;
w_rp0 = 0.2;

params0 = [Kr0, w_rp0];
params0_std = [100, 10];

% Bounds of the parameters
LB = [0.1, 0];
UB = [10^3, 1.2];


for i = 1:length(names)
    % default values
    pri_mu = NaN; %default prior gaussian mean
    pri_sig = Inf; %default prior gaussian variance
    targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
    localflag = 0; %is this local to this dataset or shared amongst batches?

    if i==1 % K_r
%             pri_mu = 5;
%             pri_sig = 20;
        targetflag = 1;
        localflag = 0;
    elseif i==4 % w_rp
%             pri_mu = rand;
%             pri_sig = 10;
        targetflag = 1;
        localflag = 1;
%     else
%         pri_mu = params0(i);
%         pri_sig = params0_std(i);
%         targetflag = 0; % Fix this parameter
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
% 
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

for k=1:2 %length(names)
    params_inferred(1,k) = mean(chain(n_burn+1:end,k));
    params_inferred_sigma(1,k) = std(chain(n_burn+1:end,k));
end

%% Save the MCMC result into a structure for future usage
MCMC_6A1R_RuntWT(m).name = constructNames{construct};
MCMC_6A1R_RuntWT(m).results = results;
MCMC_6A1R_RuntWT(m).chain = chain;
MCMC_6A1R_RuntWT(m).s2chain = s2chain;

%     % inferred parameters
MCMC_6A1R_RuntWT(m).params_inferred = params_inferred;
MCMC_6A1R_RuntWT(m).params_inferred_sigma = params_inferred_sigma;
    

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
    constructIndex = [2,5,6];
    construct = constructIndex(index);
    
    % extract parameters from the MCMC results (we're pulling the Bcd
    % dependent parameters and Run dependent parameters respectively from
    % their own inferences. (Bcd parameters from the Run null, and Run
    % parameters from the Run WT).
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_Run = MCMC_6A1R_RuntWT(index).params_inferred;
    
    % extract individual parameters to construct the input parameter
    Kb = params_Bcd(1);
    w_bp = params_Bcd(2);
    p = params_Bcd(3);
    R_max = params_Bcd(4);
    Kr = params_Run(1);
    w_rp = params_Run(2);
    
    params_MCMC = [Kb, Kr, w_bp, w_rp, p, R_max];
    output = model_6A1R_HillModel_V3_direct(params_MCMC, TF);
    fit_nulls = model_6A1R_HillModel_V3_direct(params_MCMC, TF_null);
    
    % data (WT)
    Rate_WT = compiledData{construct+1,9};
    Rate_WT_SEM = compiledData{construct+1,10};
    % data (null)
    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};
    
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
%     pause
%     saveas(gcf,[FigPath,filesep,'raw_fits_WT_null_', constructNames{construct}  ,'.tif']); 
%     saveas(gcf,[FigPath,filesep,'raw_fits_WT_null_', constructNames{construct} ,'.pdf']); 
end

%% generate corner plots
for index = 1:3
%     clf
    construct = constructIndex(index);
    chain = MCMC_6A1R_RuntWT(index).chain;
    n_burn = 0.5*n_steps;
    m = [log10(chain(n_burn+1:end,1)), chain(n_burn+1:end,2)];
    corner = figure;
%     names = {'K_{b}','\omega_{bp}','p','R_{max}'};
    names = {'log(K_{r})','\omega_{rp}'};
    ecornerplot(m,'names',names);
    
    saveas(gcf,[FigPath,filesep,'Corner_plot_K_r_1000_logscale_', constructNames{construct} ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'Corner_plot_K_r_1000_logscale_', constructNames{construct} ,'.pdf']); 
end

%% generate plots of inferred parameters
hold on
for index = 1:3
    construct = constructIndex(index);
    params_inferred = MCMC_6A1R_RuntWT(index).params_inferred;
    params_inferred_std = MCMC_6A1R_RuntWT(index).params_inferred_sigma;
    errorbar(1:2, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 4])
xticks([1,2])
xticklabels({'K_{r}','\omega_{rp}'})
yticks([0 20 40 60 80 100])
xlabel('inferred parameters')
ylabel('inferred parameters')
legend('100','001','010','Location', 'NorthEast')

box on
StandardFigure(gcf,gca)

% % save the plots
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_point_estimate_BcdParams','.tif']);
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_point_estimate_BcdParams','.pdf']);

%% generate the plot of inferred parameters (log scale)
hold on
for index = 1:3
    construct = constructIndex(index);
    params_inferred = MCMC_6A1R_RuntWT(index).params_inferred;
    params_inferred_std = MCMC_6A1R_RuntWT(index).params_inferred_sigma;
    errorbar(1:2, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 4])
xticks([1,2])
xticklabels({'K_{r}','\omega_{rp}'})
% yticks([0 20 40 60 80 100])
xlabel('inferred parameters')
ylabel('inferred parameters')
legend('100','001','010','Location', 'NorthEast')

set(gca,'YScale','log')
box on
StandardFigure(gcf,gca)

% % save the plots
saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_point_estimate_BcdParams_LogScale','.tif']);
saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_point_estimate_BcdParams_LogScale','.pdf']);

%% generate the plot for C.V. of inferred parameters

% calculate the Coefficient of Variation (C.V.) of inferred parameters
% across constructs (different enhancers).

% First, initialize the matrix
params_inferred_all = [];
params_inferred_error_all = [];

for index = 1:3
    construct = constructIndex(index);
    % params from the Runt null
    params_inferred_null(index,:) = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_inferred_std_null(index,:) = MCMC_6A0R_RuntNulls(construct).params_inferred;
    % params from the Runt WT
    params_inferred_all(index,:) = MCMC_6A1R_RuntWT(index).params_inferred;
    params_inferred_error_all(index,:) = MCMC_6A1R_RuntWT(index).params_inferred_sigma;
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
errorbar([2,4], params_CV, params_CV_error,'o','LineWidth',2)
errorbar([1,3,5,6], params_null_CV, params_null_CV_error, 'o','LineWidth',2)

xlim([0 7])
xticks([1,2,3,4,5,6])
xticklabels({'K_{b}','K_{r}','\omega_{bp}','\omega_{rp}','p','R_{max}'})
% yticks([0 100 200 300 400])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2])
% set(gca,'YScale','log')
% xlabel('')
ylabel('Coefficient of Variation')
legend('MCMC (Run WT)','MCMC (Run null)','Location','NorthWest')

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
save([FilePath, filesep, 'MCMC_6A1R_RuntWT_params.mat'],'MCMC_6A1R_RuntWT')

end