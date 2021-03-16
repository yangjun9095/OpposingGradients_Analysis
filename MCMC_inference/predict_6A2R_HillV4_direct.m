%% predict_6A2R_cases
function predict_6A2R_HillV4_direct
%% Description
% This script is to predict the level of hbP2 + 2Run sites
% Note that it's for Hill.V4     

% input : 
% model type (parameters)
% dataset (constructs)
% 
% 
%% FilePaths
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';

%% Load the data (Bcd, Run, and output)

% First, output : [001], [010], [100], load all of the data and parameters inferred
preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
data = preprocessedData.data;

% Second, the input TF data : Bcd, Run
load([FilePath, filesep, 'TFinput.mat'])

Bicoid = TFinput(:,1);
Runt = TFinput(:,2);
RuntNull = TFinput(:,3);

%% (Updated on 2/12/2021)
%% We can also pull the result from the MCMC_6A0R_RuntNulls_BcdParams.mat 
% path : 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls'
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';
MCMC_6A0R_RuntNulls = load([FilePath, filesep, 'RuntNulls', filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat'],'MCMC_6A0R_RuntNulls');

%% Extract the Bcd-dependent parameters : 

%% Load the Runt dependent parameters (chain) from the MCMC - 1Run site cases
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\HillV4_direct\HillV4_pt_estimate_BcdParams';
load([FilePath,filesep,'MCMC_6A1R_RuntWT_params.mat']);

% Note that this loads a structure, MCMC_6A0R_RuntNulls, containing fields
% like chain, results, params_inferred, etc.

% Use the params_MCMC for the K_r and w_rp
%% Predict the level for the Runt WT
% Use a custom-function for the model : model_6A2R_HillV3
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params

% n_simu = length(MCMC_2RunSites(1).chain);
% n_burn = 0.5*n_simu;

% r2-new : [011] : 3rd in the MCMC_6A0R_RuntNulls
% r2-close : [110] : 7th
% r2-far : [101] : 8th
ConstructIndex = [3,7,8];

% Predict the cases of 2 binding sites
for i=1:3
    % step1. Bcd-dependent parameters
%     chain_temp = MCMC_2RunSites(i).chain;
%     params_Bcd = mean(chain_temp(n_burn+1:end,:)); % [Kb, w_bp, p, R_max];
    construct = ConstructIndex(i);
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
%     params_Run = MCMC_6A1R_RuntWT.params_inferred;
    
    % step2. Run-dependent parameters
    if i==1 % r2-new : [011]
        w_rp1 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
        w_rp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
    elseif i==2 % r2-close : [110]
        w_rp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        w_rp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
    elseif i==3 % r2-far : [101]
        w_rp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        w_rp2 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
    end
    
    params_6A2R = [params_Bcd, w_rp1, w_rp2];
    
    Prediction(i,:) = model_6A2R_HillModel_V4_direct(params_6A2R, TF);
    fit_RuntNull(i,:) = model_6A0R_HillModel_V3(params_Bcd, TF);
end

%% generate plots of Prediction versus Data (2 Run sites)

for i=1:3
    construct = ConstructIndex(i);
    
    clf
   % initialize the variables
    Rate = [];
    Rate_SEM = [];
    Rate_null = [];
    
    % Pull the relevant data from the compiledData
    Rate = compiledData{construct+1,9};
    Rate_SEM = compiledData{construct+1,10};

    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};
    
    hold on
    % Runt nulls
    errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
    % Runt WT
    errorbar(APaxis, Rate, Rate_SEM,'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
    
    % Prediction
    % Runt null
    plot(APaxis(9:21), fit_RuntNull(i,:),'Color',ColorChoice(4,:), 'LineWidth', 2)
    % Runt WT
    plot(APaxis(9:21), Prediction(i,:),'Color',ColorChoice(1,:), 'LineWidth', 2)

    xlim([0.2 0.6])
    xticks([0.2 0.3 0.4 0.5 0.6])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
%     ylabel({'initial RNAP', 'loading rate (AU/min)'})
    ylabel('initial rate (AU/min)')
    legend('data(null)','data(WT)','MCMC fit','prediction')
    box on
    StandardFigure(gcf,gca)
%   %Save the plot
    FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\HillV4_direct\HillV4_6A2R_prediction';
    saveas(gcf,[FigPath,filesep,'Prediction_HillV4_',constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,'Prediction_HillV4_',constructNames{construct},'.pdf']);
    pause
end



%% Part2. Perform MCMC inference for the 6A2R cases, using the posteriors from 
% the previous round of MCMC : 6AnR Runt nulls for Bicoid/RNAP parameters,
% and 6A1R for the Run-parameters. Then, use those posteriors as priors.

%% Pick a model for the MCMC inference
mdl = @(TF, params) model_6A1R_HillModel_V4_direct(params, TF);
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
%% Loop over all constructs to perform the MCMC inference on the Bcd-dependent parameters
% [Kb, w_bp, p, R_max];

% initialize the structure to save the MCMC result
MCMC_6A2R_RuntWT = struct;
m = 1; % counter to save into the structure, as we're only dealing with 3 constructs (6A1R)

% r2-new : [011] : 3rd in the MCMC_6A0R_RuntNulls
% r2-close : [110] : 7th
% r2-far : [101] : 8th

for i=1:3
    
    %% Pull the construct of our interest (for the parameter inference)
    % Choose a construct 
    % Pick the dataset from the data.mat
    constructIndex = [3,7,8]; % [011, , 110, 101]
    construct = constructIndex(i);
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
    
    %% Define the parameters for the model (Note the order of the params) 
    % this should match with that (order) of the referenced model above.
    
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}', 'w_bp', 'p', 'R_{max}', 'w_rp1', 'w_rp2'};
    params = cell(1, length(names));
        
    % Pull the fitted parameters from the MCMC inference on the Runt
    % null datasets for one-by-one.
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_Bcd_std = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
    
    % Pull the fitted pramaeters from the MCMC inference on the 6A1R cases
    % for the Run-dependent parameters
    if i==1 % r2-new : [011]
        w_rp1 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
        w_rp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
        w_rp1_std =  MCMC_6A1R_RuntWT(2).params_inferred_sigma; %[001]
        w_rp2_std =  MCMC_6A1R_RuntWT(3).params_inferred_sigma; %[010]
    elseif i==2 % r2-close : [110]
        w_rp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        w_rp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
        w_rp1_std =  MCMC_6A1R_RuntWT(1).params_inferred_sigma; %[100]
        w_rp2_std =  MCMC_6A1R_RuntWT(3).params_inferred_sigma; %[010]
    elseif i==3 % r2-far : [101]
        w_rp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        w_rp2 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
        w_rp1_std =  MCMC_6A1R_RuntWT(1).params_inferred_sigma; %[100]
        w_rp2_std =  MCMC_6A1R_RuntWT(2).params_inferred_sigma; %[001]
    end
    % Initialize MCMC parameters.
    Kb0 = params_Bcd(1);
    w_bp0 = params_Bcd(2);
    p0 = params_Bcd(3);
    R_max0 = params_Bcd(4);
    
    % Initial condition for the Run-dependent parameters
    % get the multiplication of the two w_rp as those two are joint
    % parameters in the end, get the std appropriately with error
    % propagation.
    w_rp0 = w_rp1*w_rp2;
    w_rp0_std = sqrt((w_rp1_std/w_rp1)^2 + (w_rp2_std/w_rp2)^2);
    
    
    % In 6A2R, HillV4 model, the w_rp1 and w_rp2 are degenerate. Thus, we
    % will use the trick of fixing the w_rp2 as 1, then infer the w_rp1
    % only, basically one joint parameter \omega.
    params0 = [Kb0, w_bp0, p0, R_max0, w_rp0, 1];
    params0_std = [params_Bcd_std(1), params_Bcd_std(2), params_Bcd_std(3), params_Bcd_std(4), w_rp0_std, 10];

    % Bounds of the parameters
    LB = [0.1, 1, 0, 50, 0, 0];
    UB = [10^2, 10^2, 1, 500, 1.2, 1.2];


    for j = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?

        if j==6 % w_rp2
%             pri_mu = 0.2;
%             pri_sig = 1;
            targetflag = 0; % Fix this parameter
        else
            pri_mu = params0(j);
            pri_sig = params0_std(j);
%             targetflag = 0; % Fix this parameter
        end
        
    %     elseif i==2
    %         pri_mu = NaN;
    %         pri_sig = Inf;
    %         localflag = 1; % keep this parameter consistent across batches (different constructs of 1 Run site)
        % give a pretty narrow prior
    %     pri_mu = params0(i);
    %     pri_sig = pri_mu*0.1;
    
        params{1, j} = {names{j}, params0(j), LB(j), UB(j), pri_mu, pri_sig, targetflag, localflag};

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
    
    
    % Extract chain results into individual parameters
    params_inferred = []; 
    params_inferred_sigma = [];
    n_burn = 0.5*n_steps;

    for k=1:length(names)-1
        params_inferred(1,k) = mean(chain(n_burn+1:end,k));
        params_inferred_sigma(1,k) = std(chain(n_burn+1:end,k));
    end
    
    %% Save the MCMC result into a structure for future usage
    MCMC_6A2R_RuntWT(i).name = constructNames{construct};
    MCMC_6A2R_RuntWT(i).results = results;
    MCMC_6A2R_RuntWT(i).chain = chain;
    MCMC_6A2R_RuntWT(i).s2chain = s2chain;
    
%     % inferred parameters
    MCMC_6A2R_RuntWT(i).params_inferred = params_inferred;
    MCMC_6A2R_RuntWT(i).params_inferred_sigma = params_inferred_sigma;
    
    % update the counter
%     m=m+1;
end   

%% generate plots of raw fits, MCMC corner plots, etc.

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
%     params_inferred = [Kb, w_bp, p, R_max, w_rp1];
    params_temp = MCMC_6A2R_RuntWT(index).params_inferred;
    params_MCMC = [params_temp, 1]; % the last element 1 is for the w_rp2
    output = model_6A2R_HillModel_V4_direct(params_MCMC, TF);
    
    % generate Runt null fit (using the 6A0R_model_HillV3
    params_Bcd = params_temp(1:4);
    fit_RuntNull = model_6A0R_HillModel_V3(params_Bcd, TF);
    
    % data (Runt null)
    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};
    
    % data (WT)
    Rate_WT = compiledData{construct+1,9};
    Rate_WT_SEM = compiledData{construct+1,10};
    
    clf
    hold on
    % Runt nulls
    errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
    % Runt WT
    errorbar(APaxis, Rate_WT, Rate_WT_SEM, 'o', 'Color', ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
    
    % MCMC
    plot(APaxis(APbin_start:APbin_end), fit_RuntNull, 'Color', ColorChoice(4,:),'LineWidth', 2)
    plot(APaxis(APbin_start:APbin_end), output, 'Color', ColorChoice(1,:),'LineWidth', 2)
    xlim([0.2 0.5])
    xticks([0.2 0.3 0.4 0.5])
    ylim([0 400])
    yticks([0 100 200 300 400])
    
    xlabel('embryo length')
    ylabel('initial rate (AU/min)')
    
    box on
    legend('data(null)','data(WT)','MCMC fit','MCMC fit')
    StandardFigure(gcf,gca)
%     pause
    saveas(gcf,[FigPath,filesep,'raw_fits_', constructNames{construct}  ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'raw_fits_', constructNames{construct} ,'.pdf']); 
end

%% generate corner plots
for index = 1:3
%     clf
    construct = constructIndex(index);
    chain = MCMC_6A2R_RuntWT(index).chain;
    n_burn = 0.5*n_steps;
    m = [chain(n_burn+1:end,1), chain(n_burn+1:end,2), chain(n_burn+1:end,3), chain(n_burn+1:end,4), chain(n_burn+1:end,5)];
    corner = figure;
    names = {'K_{b}','\omega_{bp}','p','R_{max}','w_{rp-joint}'};
    ecornerplot(m,'names',names);
    
    saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct}  ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']); 
end

%% save the result into mat files , MCMC_6A2R_RuntWT
%% save the result into mat files (MCMC_6A0R_RuntNulls)
FilePath = FigPath;
save([FilePath, filesep, 'MCMC_6A2R_RuntWT_params.mat'],'MCMC_6A2R_RuntWT')

%% generate plots of inferred parameters
% hold on
% for index = 1:3
%     construct = constructIndex(index);
%     params_inferred = MCMC_6A1R_RuntWT(index).params_inferred;
%     params_inferred_std = MCMC_6A1R_RuntWT(index).params_inferred_sigma;
%     errorbar(1:5, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
% end
% 
% xlim([0 6])
% xticks([1,2,3,4,5])
% xticklabels({'K_{b}','\omega_{bp}','\omega_{rp}','p','R_{max}'})
% yticks([0 100 200 300 400])
% xlabel('inferred parameters')
% ylabel('inferred parameters')
% legend('100','001','010','Location', 'NorthWest')
% 
% box on
% StandardFigure(gcf,gca)
% 
% % % save the plots
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_post_dist_BcdParams','.tif']);
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_post_dist_BcdParams','.pdf']);
% 
% %% generate the plot of inferred parameters (log scale)
% hold on
% for index = 1:3
%     construct = constructIndex(index);
%     params_inferred = MCMC_6A1R_RuntWT(index).params_inferred;
%     params_inferred_std = MCMC_6A1R_RuntWT(index).params_inferred_sigma;
%     errorbar(1:5, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
% end
% 
% xlim([0 7])
% xticks([1,2,3,4,5])
% xticklabels({'K_{b}','\omega_{bp}','\omega_{rp}','p','R_{max}'})
% % yticks([0 100 200 300 400])
% xlabel('inferred parameters')
% ylabel('inferred parameters')
% legend('100','001','010','Location', 'NorthEast')
% 
% set(gca,'YScale','Log')
% box on
% StandardFigure(gcf,gca)
% 
% % % save the plots
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_post_dist_BcdParams_LogScale','.tif']);
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A1R_RuntWT_post_dist_BcdParams_LogScale','.pdf']);
% 
% 
% 
% %% generate the plot for C.V. of inferred parameters
% 
% % calculate the Coefficient of Variation (C.V.) of inferred parameters
% % across constructs (different enhancers).
% 
% % First, initialize the matrix
% params_inferred_all = [];
% params_inferred_error_all = [];
% 
% for index = 1:3
%     construct = constructIndex(index);
%     params_inferred_all(index,:) = MCMC_6A1R_RuntWT(index).params_inferred;
%     params_inferred_error_all(index,:) = MCMC_6A1R_RuntWT(index).params_inferred_sigma;
% end
% 
% % calculate the mean and std of parameters "over constructs".
% params_mean = mean(params_inferred_all);
% params_std = std(params_inferred_all);
% % params_error = sqrt(sum(params_inferred_error_all.^2)/length(data))
% n_boots = 100;
% parms_std_boostrap = bootstrp(n_boots, @std, params_inferred_all);
% params_SEM = std(parms_std_boostrap)./sqrt(n_boots-1);
% 
% params_CV = params_std./params_mean;
% params_CV_error = params_SEM./params_mean;
% 
% errorbar(1:5, params_CV, params_CV_error,'o','LineWidth',2)
% 
% xlim([0 6])
% xticks([1,2,3,4,5])
% xticklabels({'K_{b}','\omega_{bp}','\omega_{rp}','p','R_{max}'})
% % yticks([0 100 200 300 400])
% ylim([0 1.2])
% yticks([0 0.2 0.4 0.6 0.8 1 1.2])
% % set(gca,'YScale','log')
% % xlabel('')
% ylabel('Coefficient of Variation')
% 
% box on
% StandardFigure(gcf,gca)
% 
% % save the plots
% saveas(gcf,[FigPath,filesep, 'CV_parameters','.tif']);
% saveas(gcf,[FigPath,filesep, 'CV_parameters','.pdf']);






end
