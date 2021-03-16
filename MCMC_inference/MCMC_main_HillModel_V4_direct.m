%% MCMC_initial_slope_DirectRepression_test_Bounds_Priors


% file path
FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';

%% Default settings (also make it optional)
% fileDir = pwd;
% saveLoc = pwd;
numParPools = 8;
% n_burn = 5000;
n_steps = 2*10^4;
n_simu = n_steps;
% ratePriorWidth = 50;
AP_start = 20; % [% of embryo length]
AP_end = 45;   % [% of embryo length]
% loadPrevious = false;
% globalFit = 1;

%% Import data for the MCMC inference
% From the "preprocess_data_for_MCMC.m" script
% xdata(TFinputs) and ydata(initial rate), note that we will do
% a simultaneous fitting for the Runt WT and Runt nulls.

% We need another separate script to process the data for inputs in this
% script. : This is now done in the "preprocess_data_for_MCMC.m" script.

preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
data = preprocessedData.data;

% Load the input TF data
load([FilePath, filesep, 'TFinput.mat'])

Bicoid = TFinput(:,1);
Runt = TFinput(:,2);
RuntNull = TFinput(:,3);

%% Load the pre-processed MCMCdata that has incorporated [100], [001], and [101] datasets
% MCMCdata_1RunSite.mat has MCMCdata structure for this.
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';
load([FilePath, filesep,'MCMCdata_1RunSite.mat'],'MCMCdata')
%%  Pick a model for the fitting
mdl = @(TF, params) model_6A1R_HillModel_V4_direct(params, TF);
model.modelfun = mdl;

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, InputData) sum((InputData.ydata-mdl(InputData.xdata, params)).^2);

%% MCMC inference for all 6A1R constructs : 

% initialize the structure : MCMC_6A1R_global
MCMC_6A1R_global = struct;
% initialize the counter : m
m=1;

for construct = [2,5,6]
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


    %% Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}', 'K_{r}', 'w_bp', 'w_rp', 'p', 'R_{max}'};
    params = cell(1, length(names));

    % Initialize MCMC parameters.
    Kb0 = 60;   % 100*rand;
    Kr0 = 5;     % 100*rand;
    w_bp0 = 40;
    w_rp0 = 0.2;
    p0 = 0.1;
    R_max0 = MCMCdata.R_max;
    % R_max0 = MCMCdata(2).R_max; 
    % R_min0 = MCMCdata.R_min;

    params0 = [Kb0, Kr0, w_bp0, w_rp0, p0, R_max0];

    % Bounds of the parameters
    LB = [0.01, 0.01, 0, 0, 0, 50];
    UB = [10^2, 100, 2*10^2, 2, 1, 500];


    for i = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?

        if i==6
            pri_mu = R_max0;
            pri_sig = 100;
    %         targetflag = 0; % Fix this parameter
        elseif i==2
            pri_mu = Kr0;
            pri_sig = Inf;
            targetflag = 0; % Fix this parameter
%             localflag = 1; % keep this parameter consistent across batches (different constructs of 1 Run site)
%         elseif i==2
%             pri_mu = 10;
%             pri_sig = 100;
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
    Kr_step = 0.1;
    w_bp_step = 0.1;
    w_rp_step = 0.01;
    p_step = 0.01;
    R_max_step = 1;


    % Initialize the covariance matrix
    J0 = diag([Kb_step, Kr_step, w_bp_step, w_rp_step, p_step, R_max_step]);


    %% MCMC - Options
    options = [];
    n_steps = 2*10^5;
    n_burn = 0.5*n_steps;
    options.nsimu = n_steps; %n_steps; %Number of steps
    options.updatesigma = 1; %Update error variance
%     options.qcov = J0; %Initial covariance
    % options.burnintime = 0.5*n_steps; %Burn in time
    options.adaptint = 100;
    options.method = 'dram';
    options.verbosity = 0; %Decrease text output

    %% Run the MCMC (this whole block can be inserted inside the for loop above,
    % to run the MCMC for different constructs.

    %we're gonna run this three times and use the initial results of one
    %run as conditions for the next. this is an alternative when common least
    %squares gives results too poor to initialize with
    % index = find([2,5,6]==construct);

    results = [];
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    %[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,chain,s2chain,~]=mcmcrun(model,MCMCdata,params,options,results);
    
    %% calculate the mean and std of the inferred parameters (note that I'm
    % taking out the n_burn steps at the beginning).
    params_inferred = [];
    params_inferred_sigma = [];
    
    for k=1:length(names)-1
        params_inferred(1,k) = mean(chain(n_burn+1:end,k));
        params_inferred_sigma(1,k) = std(chain(n_burn+1:end,k));
    end

    MCMC_6A1R_global(m).name = constructNames{construct};
    MCMC_6A1R_global(m).results = results;
    MCMC_6A1R_global(m).chain = chain;
    MCMC_6A1R_global(m).s2chain = s2chain;
    
    MCMC_6A1R_global(m).params_inferred = params_inferred;
    MCMC_6A1R_global(m).params_inferred_std = params_inferred_sigma;
    
    % update the counter 
    m=m+1;
    
end

%% post-processing
% for i=1:3
%     chain = MCMC_6A1R_global(i).chain;
%     for k=1:length(names)
%         params_inferred(1,k) = mean(chain(n_burn+1:end,k));
%         params_inferred_sigma(1,k) = std(chain(n_burn+1:end,k));
%     end
%     MCMC_6A1R_global(i).params_inferred = params_inferred;
%     MCMC_6A1R_global(i).params_inferred_std = params_inferred_sigma;
% end
%% Diagnose the MCMC result
% stats = chainstats(chain,results);
n_burn = 0.5*n_steps;
chainstats(chain(n_burn+1:end,:),results);

%% generate corner plots
%n_burn = 1;%0.5*n_steps;% 0.5*10^4;
constructIndex = [2,5,6];

for index = 1:3
    construct = constructIndex(index);
    chain = [];
    chain = MCMC_6A1R_global(index).chain;
    m = [chain(n_burn:end,1), chain(n_burn:end,2), chain(n_burn:end,3), chain(n_burn:end,4), chain(n_burn:end,5)];
    corner = figure;
    names = {'K_{b}','\omega_{bp}', '\omega_{rp}','p','R_{max}'};
    ecornerplot(m,'names',names);


    % Save the plot
    % FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\direct';
    saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct}  ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']);  
end


%% Extract MAP estimates from the posterior chains

%% generate raw fits 
APaxis = 0:0.025:1;
% APbin_start;
% APbin_end;
Bcd = Bicoid(APbin_start:APbin_end); 
Run = Runt(APbin_start:APbin_end); 

TF = [Bcd, Run];

for index = 1:3
    % define the construct index (which is consistent with the way it's
    % defiend in the compiledData.mat)
    constructIndex = [2,5,6];
    construct = constructIndex(index);
    
    % extract parameters from the MCMC results (we're pulling the Bcd
    % dependent parameters and Run dependent parameters respectively from
    % their own inferences. (Bcd parameters from the Run null, and Run
    % parameters from the Run WT).
    params_temp = MCMC_6A1R_global(index).params_inferred;
    
    params_MCMC = [params_temp(1), Kr0, params_temp(2:end)];
    output = model_6A1R_HillModel_V4_direct(params_MCMC, TF);
    
    
    params_Bcd = [params_temp(1:2), params_temp(4:5)];
    output_null = model_6A0R_HillModel_V3(params_Bcd, TF);
    
    % data (WT)
    Rate_WT = compiledData{construct+1,9};
    Rate_WT_SEM = compiledData{construct+1,10};
    
    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};
    

    
    clf
    hold on
    % data (Runt null and Runt WT)
    errorbar(APaxis, Rate_null, Rate_null_SEM, 'o', 'Color', ColorChoice(4,:),'LineWidth', 1)    
    errorbar(APaxis, Rate_WT, Rate_WT_SEM, 'o', 'Color', ColorChoice(1,:),'LineWidth', 1)
    
    % MCMC fit
    plot(APaxis(APbin_start:APbin_end), output_null, 'Color', ColorChoice(4,:),'LineWidth', 2)
    plot(APaxis(APbin_start:APbin_end), output, 'Color', ColorChoice(1,:),'LineWidth', 2)
    xlim([0.2 0.5])
    xticks([0.2 0.3 0.4 0.5])
    ylim([0 400])
    yticks([0 100 200 300 400])
    
    xlabel('embryo length')
    ylabel('initial rate (AU/min)')
    
    box on
    legend('data','MCMC fit')
    legend('Runt null','Runt WT','Fit (null)', 'Fit (WT)')

    StandardFigure(gcf,gca)
%     pause
    saveas(gcf,[FigPath,filesep,'raw_fits_', constructNames{construct}  ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'raw_fits_', constructNames{construct} ,'.pdf']); 
end

%% Save the MCMC results
modelName = 'HillV3';

FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';

save([FilePath,filesep,'MCMCresult_',modelName,'_',constructNames{construct},'.mat'],'results','chain','s2chain');

%% generate plots of the inferred parameters
% initialize the indexing
k=1;

hold on
for construct = [2,5,6]
    errorbar(1:6, params_MCMC(k,:), params_MCMC_std(k,:),'o','Color', ColorChoice(construct,:), 'MarkerFaceColor', ColorChoice(construct,:))% 'MarkerFaceColor', ColorChoice(construct,:)
    k=k+1;
end

set(gca, 'YScale','log')
legend('100','001','010')

xlim([0 8])
xticklabels({'','K_{b}','K_{r}','\omega_{bp}','\omega_{rp}','p','R_{max}',''})
xlabel('parameters')
ylabel('inferred values')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'params_MCMC_1RunSite','.tif']); 
saveas(gcf,[FigPath,filesep,'params_MCMC_1RunSite','.pdf']); 
 
%% save the result into mat files (MCMC_6A1R_global)
FilePath = FigPath;
save([FilePath, filesep, 'MCMC_6A1R_global_params.mat'],'MCMC_6A1R_global')