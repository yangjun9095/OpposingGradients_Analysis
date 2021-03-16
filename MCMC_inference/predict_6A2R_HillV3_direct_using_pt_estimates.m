%% predict_6A2R_cases
function predict_6A2R_HillV3_competition
%% Description
% This script is to predict the level of hbP2 + 2Run sites

% input : 
% model type (parameters)
% dataset (constructs)
% 
% 
%% FilePaths
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Competition';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Competition';

%% Load the data (Bcd, Run, and output)
% Bcd, Run, TF, TF_null
% APbin_start, APbin_end

%% Load the MCMC results

%% Load the MCMC results 
% 1) Bcd-dependent parameters from 6A0R_RuntNulls.mat
% 2) Run-dependent parameters from 6A1R_RuntWT.mat
tempPath1 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls';
load([tempPath1, filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat']);
tempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\Point_estimate_Bcd_params_fromRunNulls';
load([tempPath2, filesep, 'MCMC_6A1R_RuntWT_params.mat']);

%% Predict the level for the Runt WT
% Use a custom-function for the model : model_6A2R_HillV3_competition
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_br1, w_br2] = params

% r2-new : [011]
% r2-close : [110]
% r2-far : [101]
constructIndex = [3,7,8];

% initialize the Prediction matrix
Prediction = [];
fit_nulls = [];

% Predict the cases of 2 binding sites
for i=1:3
    
    construct = constructIndex(i);
    % extract the Bcd-parameters from 6AnR, Runt null datasets
    params_Bcd = [];
    params_temp1 = [];
    params_temp2 = [];
    params_6A2R = [];
    
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred; %[Kb, w_bp, p, R_max]
    
    if i==1 % r2-new : [011]
        params_temp1 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
        params_temp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
    elseif i==2 % r2-close : [110]
        params_temp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        params_temp2 = MCMC_6A1R_RuntWT(3).params_inferred; %[010]
    elseif i==3 % r2-far : [101]
        params_temp1 = MCMC_6A1R_RuntWT(1).params_inferred; %[100]
        params_temp2 = MCMC_6A1R_RuntWT(2).params_inferred; %[001]
    end
    Kr1 = params_temp1(1);
    w_rp1 = params_temp1(2);
    Kr2 = params_temp2(1);
    w_rp2 = params_temp2(2);
    
    params_6A2R = [params_Bcd, Kr1, Kr2, w_rp1, w_rp2];
    
    Prediction(i,:) = model_6A2R_HillModel_V3_direct(params_6A2R, TF);
    fit_nulls(i,:) = model_6A0R_HillModel_V3(params_6A2R, TF);
end

%% generate plots of Prediction versus Data (2 Run sites)

for i=1:3
    construct = constructIndex(i);
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
    plot(APaxis(APbin_start:APbin_end), fit_nulls(i,:),'Color',ColorChoice(4,:), 'LineWidth', 2)
    % Runt WT
    plot(APaxis(APbin_start:APbin_end), Prediction(i,:),'Color',ColorChoice(1,:), 'LineWidth', 2)

    xlim([0.2 0.5])
    xticks([0.2 0.3 0.4 0.5])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
    ylabel({'initial RNAP', 'loading rate (AU/min)'})
    legend('data (null)','data (WT)', 'MCMC fit (null)', 'prediction')
    
    box on
    StandardFigure(gcf,gca)
%   %Save the plot
    saveas(gcf,[FigPath,filesep,'Prediction_HillV3_direct_',constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,'Prediction_HillV3_direct_',constructNames{construct},'.pdf']);
    pause
end

%% Predict the Rates with Run-Run cooperativity
% mnodel_6A2R_HillModelV3_RunCoop.m
% w_rr > 1 : stronger repression than independent Runt binding.
w_rr = [10^(-6), 1, 10^6, 10^12, 10^18, 10^24];

% Predict the level for the Runt WT
% Use a custom-function for the model : model_6A2R_HillV3
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params

n_simu = length(MCMC_2RunSites(1).chain);
n_burn = 0.5*n_simu;

% r2-new : [011]
% r2-close : [110]
% r2-far : [101]
% Predict the cases of 2 binding sites
Prediction_Run_Run_coop_compete = zeros(3, length(w_rr), length(TF_combined));

% Loop over constructs (2 Run sites)
for i=1:3
    % step1. Bcd-dependent parameters
    construct = ConstructIndex(i);
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    if i==1 % r2-new : [011]
        params_temp1 = params_MCMC(2,:); %[001]
        params_temp2 = params_MCMC(3,:); %[010]
    elseif i==2 % r2-close : [110]
        params_temp1 = params_MCMC(1,:); %[100]
        params_temp2 = params_MCMC(3,:); %[010]
    elseif i==3 % r2-far : [101]
        params_temp1 = params_MCMC(1,:); %[100]
        params_temp2 = params_MCMC(2,:); %[001]
    end
    Kr1 = params_temp1(2);
    w_rp1 = params_temp1(4);
    Kr2 = params_temp2(2);
    w_rp2 = params_temp2(4);
    
    for j=1:length(w_rr)
        omega_rr = w_rr(j);

        params_6A2R = [params_Bcd, Kr1, Kr2, w_rp1, w_rp2, omega_rr];

        Prediction_Run_Run_coop_compete(i,j,:) = model_6A2R_HillModel_V3_competition_RunCoop(params_6A2R, TF_combined);
    end
end

%% 


for i=1:3
    % construct define
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
    
    % Run-Run interaction strength
    for j=1:length(w_rr)
        cIndex = floor((log10(w_rr(j))-(-6))/30*256);
        if cIndex ==0
            cIndex = cIndex + 1;
        end
    
        color = cmap(cIndex,:);
        % Prediction
        % Runt null
        plot(APaxis(9:19), squeeze(Prediction_Run_Run_coop_compete(i,j,1:11)),'Color',ColorChoice(4,:), 'LineWidth', 2)
        % Runt WT
        plot(APaxis(9:19), squeeze(Prediction_Run_Run_coop_compete(i,j,12:end)),'Color',color, 'LineWidth', 2)
%         pause
    end

    xlim([0.2 0.6])
    xticks([0.2 0.3 0.4 0.5 0.6])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
    ylabel({'initial RNAP', 'loading rate (AU/min)'})
    
    box on
    colorbar
    StandardFigure(gcf,gca)
    pause
%   %Save the plot
    FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Competition';
    saveas(gcf,[FigPath,filesep,'Prediction_HillV3_',constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,'Prediction_HillV3_',constructNames{construct},'.pdf']);
%     pause
end


%% 
%% Modification 2: Higher-order cooperativity between Run's action
% w term for when the two Run molecules are bound, and exerting some effect to
% the RNAP.

%% Predict the Rates with higher order cooperativity between Run-RNAP
% mnodel_6A2R_HillModelV3_RunCoop.m
% w_rr < 1 : stronger repression than independent Runt effect.
w_ho = 10.^([-10 -1 -0.1 0 1 10]);

% Predict the level for the Runt WT
% Use a custom-function for the model : model_6A2R_HillV3
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params

n_simu = length(MCMC_2RunSites(1).chain);
n_burn = 0.5*n_simu;

% r2-new : [011]
% r2-close : [110]
% r2-far : [101]
% Predict the cases of 2 binding sites
Prediction_higher_coop = zeros(3, length(w_ho), length(TF_combined));
for i=1:3
    construct = ConstructIndex(i);
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    if i==1 % r2-new : [011]
        params_temp1 = params_MCMC(2,:); %[001]
        params_temp2 = params_MCMC(3,:); %[010]
    elseif i==2 % r2-close : [110]
        params_temp1 = params_MCMC(1,:); %[100]
        params_temp2 = params_MCMC(3,:); %[010]
    elseif i==3 % r2-far : [101]
        params_temp1 = params_MCMC(1,:); %[100]
        params_temp2 = params_MCMC(2,:); %[001]
    end
    Kr1 = params_temp1(2);
    w_rp1 = params_temp1(4);
    Kr2 = params_temp2(2);
    w_rp2 = params_temp2(4);
    
    for j=1:length(w_ho)
        omega_ho = w_ho(j);

        params_6A2R = [params_Bcd, Kr1, Kr2, w_rp1, w_rp2, omega_ho];

        Prediction_higher_coop(i,j,:) = model_6A2R_HillModel_V3_competition_HigherCoop(params_6A2R, TF_combined);
    end
end



%% generate plots showing the effect of higher order cooperativity between Run-RNAP
cmap = colormap('viridis');
caxis([-10 10])

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
    
%     cbar = linspace(0, 256, 22);

    

    % Run-Run interaction strength
    for j=1:length(w_ho)
        cIndex = floor((log10(w_ho(j))-(-10))/20*256);
        if cIndex ==0
            cIndex = cIndex + 1;
        end
    
        color = cmap(cIndex,:);
        % Prediction
        % Runt null
        plot(APaxis(9:19), squeeze(Prediction_higher_coop(i,j,1:length(TF_combined)/2)),'Color',ColorChoice(4,:), 'LineWidth', 2)
        % Runt WT
        plot(APaxis(9:19), squeeze(Prediction_higher_coop(i,j,length(TF_combined)/2+1:end)),'Color',color, 'LineWidth', 2)

    end

    xlim([0.2 0.6])
    xticks([0.2 0.3 0.4 0.5 0.6])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
    ylabel({'initial RNAP', 'loading rate (AU/min)'})
    
    colorbar
    box on
    StandardFigure(gcf,gca)

%   %Save the plot
    saveas(gcf,[FigPath,filesep,'Prediction_HillV3_HighCoop_',constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep, 'Prediction_HillV3_HighCoop_',constructNames{construct},'.pdf']);
    pause
end
%% Approach 2 : fixing the Kr,w_rp from the 1 Run site cases, fit for the Bcd dependent parameters
% Use a custom-function for the model : model_6A2R_HillV3
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params

n_simu = length(MCMC_2RunSites(1).chain);
n_burn = 0.5*n_simu;

% r2-new : [011]
% r2-close : [110]
% r2-far : [101]

% MCMC inference for the cases of 2 binding sites
for i=1:3
    chain_temp = MCMC_2RunSites(i).chain;
    params_Bcd = mean(chain_temp(n_burn+1:end,:)); % [Kb, w_bp, p, R_max];
    
    if i==1 % r2-new : [011]
        construct = 3;
        params_temp1 = params_MCMC(2,:); %[001]
        params_temp2 = params_MCMC(3,:); %[010]
    elseif i==2 % r2-close : [110]
        construct = 7;
        params_temp1 = params_MCMC(1,:); %[100]
        params_temp2 = params_MCMC(3,:); %[010]
    elseif i==3 % r2-far : [101]
        construct = 8;
        params_temp1 = params_MCMC(1,:); %[100]
        params_temp2 = params_MCMC(2,:); %[001]
    end
    Kr1 = params_temp1(2);
    w_rp1 = params_temp1(4);
    Kr2 = params_temp2(2);
    w_rp2 = params_temp2(4);
    
%     params_6A2R = [params_Bcd, Kr1, Kr2, w_rp1, w_rp2];
    
    % Define the datasets to plug into the MCMC protocol
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



    % Perform the MCMC inference on the Kb, w_bp, p, R_max
    % Pick a model for the fitting
    mdl = @(TF, params) model_6A2R_HillModel_V3(params, TF);
    %leaving this here in case it'll be useful in the future
    model.ssfun = @(params, data) sum((data.ydata-mdl(data.xdata, params)).^2);
    model.modelfun = mdl;  %use mcmcrun generated ssfun 

    % Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}', 'w_bp', 'p', 'R_{max}', 'K_{r1}', 'K_{r2}', 'w_{rp1}', 'w_{rp2}'};
    params = cell(1, length(names));

    % Initialize MCMC parameters.

    params0 = [params_Bcd, Kr1, Kr2, w_rp1, w_rp2];

    % Bounds of the parameters
    LB = [0.1, 1, 0, 50, 0.1 0.1 0 0 ];
    UB = [10^2, 10^2, 0.5, 500, 10^2, 10^2, 1, 1];


    for j = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?

        if j==4
            pri_mu = MCMCdata.R_max;
            pri_sig = 20;
    %         targetflag = 0; % Fix this parameter
        elseif j>4
            targetflag = 0; % Fix these parameters (for repressors)
        else
            pri_mu = NaN;
            pri_sig = Inf;
        end
        % give a pretty narrow prior
    %     pri_mu = params0(i);
    %     pri_sig = pri_mu*0.1;
        params{1, j} = {names{j}, params0(j), LB(j), UB(j), pri_mu, pri_sig, targetflag, localflag};

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
    J0 = diag([Kb_step, w_bp_step, p_step, R_max_step, ...
                Kr_step, Kr_step, w_rp_step, w_rp_step]);


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
    
    MCMC_2RunSites_global(i).results = results;
    MCMC_2RunSites_global(i).chain = chain;
end

%% Diagnose the MCMC result
% stats = chainstats(chain,results);
n_burn = 0.5*n_steps;
chainstats(chain(n_burn+1:end,:),results);

%% generate corner plots
%n_burn = 1;%0.5*n_steps;% 0.5*10^4;
m = [chain(n_burn:end,1), chain(n_burn:end,2), chain(n_burn:end,3), chain(n_burn:end,4), chain(n_burn:end,5),  chain(n_burn:end,6)];
corner = figure;
names = {'K_{b}','K_{r}','\omega_{bp}', '\omega_{rp}','p','R_{max}'};
ecornerplot(m,'names',names);


% Save the plot
% FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';
% saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct}  ,'.tif']); 
% saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']); 

%% Let's use the inferred [Kb, w_bp, p, R_max] to see if the models with fixed repressor parameters can predict the 2 Run sites
clear Prediction
for i=1:3
    params_Bcd_temp = mean(MCMC_2RunSites_global(i).chain(n_burn+1:end,:));
        if i==1 % r2-new : [011]
        params_temp1 = params_MCMC(2,:); %[001]
        params_temp2 = params_MCMC(3,:); %[010]
    elseif i==2 % r2-close : [110]
        params_temp1 = params_MCMC(1,:); %[100]
        params_temp2 = params_MCMC(3,:); %[010]
    elseif i==3 % r2-far : [101]
        params_temp1 = params_MCMC(1,:); %[100]
        params_temp2 = params_MCMC(2,:); %[001]
        end
    
    Kr1 = params_temp1(2);
    w_rp1 = params_temp1(4);
    Kr2 = params_temp2(2);
    w_rp2 = params_temp2(4);
    
    params_6A2R = [params_Bcd_temp, Kr1, Kr2, w_rp1, w_rp2];
    
    Prediction(i,:) = model_6A2R_HillModel_V3(params_6A2R, TF_combined);
end
%% generate plots of Prediction versus Data (2 Run sites)

for i=1:3
    if i==1
        construct = 3;
    elseif i==2
        construct = 7;
    elseif i==3
        construct = 8;
    end
    
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
    plot(APaxis, Prediction(i,1:41),'Color',ColorChoice(4,:), 'LineWidth', 2)
    % Runt WT
    plot(APaxis, Prediction(i,42:end),'Color',ColorChoice(1,:), 'LineWidth', 2)

    xlim([0.2 0.6])
    xticks([0.2 0.3 0.4 0.5 0.6])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
    ylabel({'initial RNAP', 'loading rate (AU/min)'})

    StandardFigure(gcf,gca)

%   %Save the plot
%     saveas(gcf,[FigPath,filesep,'Prediction_HillV3_',constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,'Prediction_HillV3_',constructNames{construct},'.pdf']);
    pause
end

