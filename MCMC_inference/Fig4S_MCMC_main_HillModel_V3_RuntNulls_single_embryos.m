function Fig4S_MCMC_main_HillModel_V3_RuntNulls_single_embryos
%% Description
% This script is doing MCMC fit for all Runt null data one by one, to
% compare the parameters for the Hill.V3 model, to see what varies across
% constructs.

% This expended script is for fitting to individual embryos. As an example,
% we will use the [010] construct, which seems to have the highest
% embryo-to-embryo variability.

%% Variables

%% Import data for the MCMC inference
% From the "preprocess_data_for_MCMC.m" script
% xdata(TFinputs) and ydata(initial rate), note that we will do
% a simultaneous fitting for the Runt WT and Runt nulls.

% We need another separate script to process the data for inputs in this
% script. : This is now done in the "preprocess_data_for_MCMC.m" script.
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';

% Load the output  
preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
data = preprocessedData.data;

% Load the input TF data
%load([FilePath, filesep, 'TFinput.mat'])

%Bicoid = TFinput(:,1);
%Runt = TFinput(:,2);
%RuntNull = TFinput(:,3);

compiledData = load('S:\YangJoon\Dropbox\OpposingGradientsFigures\mat files\compiledData.mat')

%% import necessary modules, such as color schemes
ColorChoice = [];
%% Define the FilePath, FigPath to save the results
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls';

%% If for regenerating plots, just load the result file
A = load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\MCMC_6A0R_RuntNulls_BcdParams.mat')

% with weakly informed prior
% FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\w_bp_0-100_WeakPrior_w_bp'
% load([FilePath, filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat'])
%% Set up the MCMC inference
% Note that we don't have any parameters regarding the repressor, thus we
% will just use the simpler form of the model, "model_6A0R_HillModel_V3.m"

%% MCMC Default settings (also make it optional)
% fileDir = pwd;
% saveLoc = pwd;
numParPools = 8;
% n_burn = 5000;
n_steps = 2*10^4;
n_simu = n_steps;
% ratePriorWidth = 50;
% AP_start = 20; % [% of embryo length]
% AP_end = 50;   % [% of embryo length]
% loadPrevious = false;
% globalFit = 1;

%% Pick a model for the MCMC inference
mdl = @(TF, params) model_6A0R_HillModel_V3(params, TF);
model.modelfun = mdl;

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, InputData) sum((InputData.ydata-mdl(InputData.xdata, params)).^2);


%% MCMC - Options
options = [];
n_steps = 2*10^4;
n_burn = 0.5*n_steps;
options.nsimu = n_steps; %n_steps; %Number of steps
options.updatesigma = 1; %Update error variance
%     options.qcov = J0; %Initial covariance
% options.burnintime = 0.5*n_steps; %Burn in time
options.adaptint = 1;
options.method = 'dram';
options.verbosity = 0; %Decrease text output
%% Loop over all constructs to perform the MCMC inference on the Bcd-dependent parameters
% [Kb, w_bp, p, R_max];

construct = 6; % [010] (r1-mid)

Data = data(construct);

[~, num_embryos] = size(Data.Rate_null_individual);

for embryo = 1:num_embryos
    %% Pull the construct of our interest (for the parameter inference)
    % Choose a construct 
    % Pick the dataset from the data.mat
    Data = data(construct);

    % MCMC analysis on the initial slope (averaged over embryos)
    % initialize the AP bins
    APaxis = Data.APbins;

    %Truncate APbins to user-specified range (input was optional argument in
    %this function.
    NoNaN_index_null = ~isnan(Data.Rate_null_individual(:,embryo));
%     NoNaN_index_WT = ~isnan(Data.Rate_WT);
    % calculate the AP bins that are not NaNs in both WT and Null datasets
    NoNaN_index = NoNaN_index_null;%.*NoNaN_index_WT;

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
    Rate_null = Data.Rate_null_individual(:,embryo);

    % Truncate the vectors using the range of AP bins
    APbins = APaxis(APbins_fit);
%     Rate_WT = Rate_WT(APbins_fit);
    Rate_null_forFit = Rate_null(APbins_fit);

    Bcd = Bicoid(APbins_fit);
%     Run = Runt(APbins_fit);
%     RunNull = RuntNull(APbins_fit);

    % Decide whether we want to fit the Runt null/WT data together or not.
    % depending on this, we will set the xdata and ydata for the fitting.
    MCMCdata = struct;
    MCMCdata.APdata = [APbins'];
    MCMCdata.ydata = [Rate_null_forFit];
    % input TF
    MCMCdata.Bcd = [Bcd];
%     MCMCdata.Run = [RunNull];
    MCMCdata.xdata = [MCMCdata.Bcd];

    MCMCdata.R_max = max(MCMCdata.ydata);
    MCMCdata.R_min = min(MCMCdata.ydata);
    
    %% Define the parameters for the model
    % put the initial parameters and bounds in a form that the mcmc function
    % accepts
    names = {'K_{b}', 'w_bp', 'p', 'R_{max}'};
    params = cell(1, length(names));

    % Initialize MCMC parameters.
    Kb0 = 10;   % 100*rand;
    w_bp0 = 50;
    p0 = 0.1;
    R_max0 = MCMCdata.R_max;
    % R_min0 = MCMCdata.R_min;

    params0 = [Kb0, w_bp0, p0, R_max0];

    % Bounds of the parameters
    LB = [0.1, 0, 0.0001, 50];
    UB = [100, 100, 1, R_max0+100];


    for i = 1:length(names)
        % default values
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.    
        localflag = 0; %is this local to this dataset or shared amongst batches?

        if i==4
            pri_mu = R_max0;
            pri_sig = 100;
    %         targetflag = 0; % Fix this parameter
    %     elseif i==2
    %         pri_mu = NaN;
    %         pri_sig = Inf;
    %         localflag = 1; % keep this parameter consistent across batches (different constructs of 1 Run site)
%         elseif i==1
%             pri_mu = 10;
%             pri_sig = 100;
%         elseif i==2
%             pri_mu = 50;
%             pri_sig = 1000;
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
    w_bp_step = 1;
    p_step = 0.01;
    R_max_step = 1;


    % Initialize the covariance matrix
%     J0 = diag([Kb_step, w_bp_step, p_step, R_max_step]);




    %% Run the MCMC (this whole block can be inserted inside the for loop above,
    % to run the MCMC for different constructs.

    %we're gonna run this three times and use the initial results of one
    %run as conditions for the next. this is an alternative when common least
    %squares gives results too poor to initialize with


    results = [];
    %[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,chain,s2chain,~]=mcmcrun(model,MCMCdata,params,options,results);
    
    % Extract chain results into individual parameters
    params_inferred = []; 
    params_inferred_sigma = [];
    n_burn = 0.5*n_steps;

    for k=1:length(names)
        params_inferred(1,k) = mean(chain(n_burn+1:end,k));
        params_inferred_std(1,k) = std(chain(n_burn+1:end,k));
    end
    
    %% Save the MCMC result into a structure for future usage
    MCMC_6A0R_RuntNulls_010(embryo).results = results;
    MCMC_6A0R_RuntNulls_010(embryo).chain = chain;
    MCMC_6A0R_RuntNulls_010(embryo).s2chain = s2chain;
    
%     % inferred parameters
    MCMC_6A0R_RuntNulls_010(embryo).params_inferred = params_inferred;
    MCMC_6A0R_RuntNulls_010(embryo).params_inferred_std = params_inferred_std;

    %MCMC_6A0R_RuntNulls(construct).MCMCpred_out = out;
    
end   


%% post-processing for inferred parameters (from individual embryos)

params_inferred = zeros(num_embryos,4);
params_inferred_std = zeros(num_embryos,4);

for i=1:num_embryos
    params_inferred(i,:) = MCMC_6A0R_RuntNulls_010(i).params_inferred;
    params_inferred_std(i,:) = MCMC_6A0R_RuntNulls_010(i).params_inferred_std;
end

params_inferred_mean = mean(params_inferred);
params_inferred_sem = std(params_inferred)./sqrt(num_embryos);

%% plot the mean and std of inferred parameters
% 1) inferred from individual embryos, 2) inferred from "averaged embryos"

hold on
% inferred from "averaged embryo"
errorbar([1,3,5,7], MCMC_6A0R_RuntNulls(6).params_inferred, ...
            MCMC_6A0R_RuntNulls(6).params_inferred_std,'o','CapSize',0)

% inferred from individual embryos
errorbar([2,4,6,8], params_inferred_mean, params_inferred_sem,'o','CapSize',0)


%% plot the mean and std of inferred parameters 
% one parameter at a time, as the y-axis scale is different

yscale = {[0, 100], [0,100], [0, 1], [0, 400]};
ytickslabels = {[0,20,40,60,80,100],[0,20,40,60,80,100],[0, 0.2 0.4 0.6 0.8 1], [0,100,200,300,400] };
params_names = {'K_{b}', '\omega_{bp}','p','R'};

params_names_save = {'K_b','w_bp','p','R'};


for i= 1:4 % number of parameters
    clf
    hold on
    errorbar(0.5, MCMC_6A0R_RuntNulls(6).params_inferred(i),...
                MCMC_6A0R_RuntNulls(6).params_inferred_std(i),'o','MarkerFaceColor','auto')
    errorbar(1.5, params_inferred_mean(i), params_inferred_sem(i),'o','MarkerFaceColor','auto')
    
    legend('averaged','individual')

    xlim([0 2])
    ylim(yscale{i})
    yticks(ytickslabels{i})
    xticks([1])
    xticklabels(params_names{i})
    box on

    StandardFigure(gcf,gca)

    %Save the plot
    saveas(gcf,[FigPath,filesep,'params_inferred_individual_vs_averaged_',params_names_save{i},'_',constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,'params_inferred_individual_vs_averaged_',params_names_save{i},'_',constructNames{construct},'.pdf']);
    pause
end


%% Post-processing for inferred parameters (mean and std)
% for construct=1:length(data)
%     % extract the MCMC chain from each dataset
%     chain_temp = MCMC_6A0R_RuntNulls(construct).chain;
%     
%     params_inferred = []; 
%     params_inferred_std = [];
%     for k=1:length(names)
%         params_inferred(1,k) = mean(chain(n_burn+1:end,k));
%         params_inferred_std(1,k) = std(chain(n_burn+1:end,k));
%     end
%     
%     MCMC_6A0R_RuntNulls(construct).params_inferred = params_inferred;
%     MCMC_6A0R_RuntNulls(construct).params_inferred_std = params_inferred_std;
% end

%% generate raw fits 

% define the fitting range along the AP axis
APaxis = 0:0.025:1;
APbin_start=9;
APbin_end=21;
Bcd = Bicoid(APbin_start:APbin_end); 

% define Bicoid as global variable to streamline the function
%Bicoid_global = Bcd;
%global Bicoid_global


%% generate plots of MCMC fitting for all constructs (Runt null)
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls'

for construct = 1:length(data)

    % Pick the dataset from the data.mat
    Data = data(construct);

    % MCMC analysis on the initial slope (averaged over embryos)
    % initialize the AP bins
    APaxis = Data.APbins;

    %Truncate APbins to user-specified range (input was optional argument in
    %this function.
    NoNaN_index_null = ~isnan(Data.Rate_null);
    % calculate the AP bins that are not NaNs in both WT and Null datasets
    NoNaN_index = NoNaN_index_null;%.*NoNaN_index_WT;

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
%     Rate_WT = Rate_WT(APbins_fit);
    Rate_null_forFit = Rate_null(APbins_fit);

    Bcd = Bicoid(APbins_fit);
%     Run = Runt(APbins_fit);
%     RunNull = RuntNull(APbins_fit);

    % Decide whether we want to fit the Runt null/WT data together or not.
    % depending on this, we will set the xdata and ydata for the fitting.
    MCMCdata = struct;
    MCMCdata.APdata = [APbins'];
    MCMCdata.ydata = [Rate_null_forFit];
    % input TF
    MCMCdata.Bcd = [Bcd];
%     MCMCdata.Run = [RunNull];
    MCMCdata.xdata = [MCMCdata.Bcd];

    MCMCdata.R_max = max(MCMCdata.ydata);
    MCMCdata.R_min = min(MCMCdata.ydata);
    Data = [MCMCdata.APdata; MCMCdata.ydata];

    % MCMCpred (model and data?)
    data_pred = MCMCdata.APdata;
    model_mcmc = @(params) model_6A0R_HillModel_V3(params,Bcd);

    % MCMC fit result
    params = MCMC_6A0R_RuntNulls(construct).params_inferred;
    output = model_6A0R_HillModel_V3(params, Bcd);

    % MCMC estimation error (MCMCpred)
    results = MCMC_6A0R_RuntNulls(construct).results;
    chain = MCMC_6A0R_RuntNulls(construct).chain;
    s2chain = MCMC_6A0R_RuntNulls(construct).s2chain;

    out = mcmcpred(results, chain(n_burn+1:end,:), [], data_pred, model_mcmc);

    
    % data
    Rate_null = compiledData{construct+9,9};
    Rate_null_SEM = compiledData{construct+9,10};
    
    clf
    hold on
    % data
    errorbar(APaxis, Rate_null, Rate_null_SEM)
    % MCMCpred
    nn = (size(out.predlims{1}{1},1)+1)/2;
    plimi = out.predlims{1};
    yl = plimi{1}(3,:);
    yu = plimi{1}(2*nn-3,:);
    %plot(APaxis(APbin_start:APbin_end), output)
    % MCMC pred plot
    %shadedErrorBar(APbins, output, [yu;yl],'lineProps','-r')
    %h=mcmcpredplot(out,MCMCdata.APdata)
    xx = MCMCdata.APdata; % x-axis
    nn = (size(out.predlims{1}{1},1) + 1) / 2;
    plimi = out.predlims{1};
    yl = plimi{1}(3,:); % y-lower limit
    yu = plimi{1}(2*nn-3,:); % y-upper limit
    yf = plimi{1}(nn,:); % y-fit
    yy=yf;


    fillyy(xx,yl,yu,[0.9 0.9 0.9])
    plot(xx,yy,'-k')
    
    xlim([0.2 0.5])
    xticks([0.2 0.3 0.4 0.5])
    ylim([0 400])
    yticks([0 100 200 300 400])
    
    xlabel('embryo length')
    ylabel('initial rate (AU/min)')
    
    box on
    legend('data','MCMC fit')
    StandardFigure(gcf,gca)
    %pause
%     
    saveas(gcf,[FigPath,filesep,'raw_fits_yLimFree_95%CI_', constructNames{construct}  ,'.tif']); 
    saveas(gcf,[FigPath,filesep,'raw_fits_yLimFree_95%CI_', constructNames{construct} ,'.pdf']); 
end

%% filter out some weird posteriors
% this is for the visualization purpose only.
% Fro example, the R_max never goes above 200, for the most cases, thus we
% will limit the maximum bounds to be 200.


%% generate corner plots
for construct = 5 %1:length(data)
%     clf
    chain = MCMC_6A0R_RuntNulls(construct).chain;
    n_burn = 0.5*n_steps;
    m = [chain(n_burn+1:end,1), chain(n_burn+1:end,2), chain(n_burn+1:end,3), chain(n_burn+1:end,4)];
    corner = figure;
    names = {'K_{b}','\omega_{bp}','p','R'}; 
    ecornerplot(m,'names',names);
    
    
%     saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct}  ,'.tif']); 
%     saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']); 
%   exportgraphics
%     exportgraphics(gcf,[FigPath, filesep,'Corner_plot_highres_w_bp_0-200_', constructNames{construct},'.pdf'],'ContentType','vector')
end

%% generate histogram for posterior chains from each parameter (with point estimates)
for construct = 5 %1:length(data)
%     clf
    chain = MCMC_6A0R_RuntNulls(construct).chain;
    n_burn = 0.5*n_steps;
%     m = [chain(n_burn+1:end,1), chain(n_burn+1:end,2), chain(n_burn+1:end,3), chain(n_burn+1:end,4)];
%     corner = figure;
    names = {'K_{b}','\omega_{bp}','p','R_{max}'};
    filenames = {'K_b', 'w_bp', 'p', 'R_max'};
    for i=1:4
%         chain_temp = [];
        chain_temp = chain(n_burn+1:end,i);
        mean_posterior = mean(chain_temp);
        std_posterior = std(chain_temp);
        figure
        hold on
        % histogram of posterior
        histogram(chain_temp,50,'Normalization','probability')
        % mean and std of posterior
        xline(mean_posterior,'LineWidth',2)
        xline(mean_posterior-std_posterior,'--','LineWidth',2)
        xline(mean_posterior+std_posterior,'--','LineWidth',2)
        xlabel(names(i))
        ylabel('frequency')
        legend('postrior','mean')
        
        box on
        StandardFigure(gcf,gca)
        saveas(gcf,[FigPath,filesep, 'histogram_mean_',filenames{i},'.tif']);
        saveas(gcf,[FigPath,filesep, 'histogram_mean_',filenames{i},'.pdf']);
    end
    
    
    
%     saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct}  ,'.tif']); 
%     saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct} ,'.pdf']); 
end
%% generate plots of inferred parameters
hold on
for construct = 1:length(data)
    params_inferred = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_inferred_std = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
    errorbar(1:4, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 7])
xticks([1,2,3,4])
xticklabels({'K_{b}','\omega_{bp}','p','R_{max}'})
yticks([0 100 200 300 400])
% xlabel('')
ylabel('inferred parameters')
legend('000','100','011','111','001','010','110','101', 'Location', 'NorthEast')

box on
StandardFigure(gcf,gca)

% % save the plots
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A0R_RuntNull_AllConstructs','.tif']);
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A0R_RuntNull_AllConstructs','.pdf']);

%% generate the plot of inferred parameters (log scale)
hold on
for construct = 1:length(data)
    params_inferred = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_inferred_std = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
    errorbar(1:4, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 7])
xticks([1,2,3,4])
xticklabels({'K_{b}','\omega_{bp}','p','R_{max}'})
% yticks([0 100 200 300 400])
set(gca,'YScale','log')
% xlabel('')
ylabel('inferred parameters')
legend('000','100','011','111','001','010','110','101', 'Location', 'NorthEast')

box on
StandardFigure(gcf,gca)

% save the plots
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A0R_RuntNull_AllConstructs_LogScale','.tif']);
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A0R_RuntNull_AllConstructs_LogScale','.pdf']);

%% generate plots of inferred parameters (as individual panels)

% extract the parameters across constructs (8x4 matrice)
for construct = 1:8%length(data)
    params_inferred(construct,:) = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_inferred_std(construct,:) = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
%     errorbar(1:4, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
end

params_names = {'K_{b}','w_{bp}','p','R'};

% generate the plot for each inferred parameter as individual panels
for params=1:4
    xindex = [0.85 0.9 0.95 1 1.05 1.1 1.15 1.2];
    figure(params)
    hold on
    for j=1:8
        errorbar(xindex(j), params_inferred(j,params),...
                    params_inferred_std(j,params),...
                    'o','LineWidth',2,'Color',ColorChoice(j,:),'CapSize',0)
    end
        xlim([0.7 1.3])
        xticks([1])
        xticklabels(params_names(params))
        
        % y limit
        if params==3 % p
            ylim([0 1])
%             yticks([0 0.2 0.4 0.6 0.8 1])
        elseif params==4 % R
            ylim([0 400]) 
%             yticks([0 100 200 300 400])
        else % K_b and/or w_bp
            ylim([0 100])
%             yticks([0 20 40 60 80 100])
            
        end

        % xlabel('')
        ylabel('inferred parameters')
%         legend('000','100','011','111','001','010','110','101', 'Location', 'NorthEast')
        
        %set(gca,'YScale','log')
        
        box on
        StandardFigure(gcf,gca)

        
        % save the plots
        FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls'
        saveas(gcf,[FigPath,filesep, 'inferred_params_',params_names{params},'.tif']);
        saveas(gcf,[FigPath,filesep, 'inferred_params_',params_names{params},'.pdf']);
end





%% generate the plot for C.V. of inferred parameters (bootstrap error)
% calculate the Coefficient of Variation (C.V.) of inferred parameters
% across constructs (different enhancers).
for construct = 1:length(data)
    params_inferred_all(construct,:) = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_inferred_error_all(construct,:) = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
end

% calculate the mean and std of parameters "over constructs".
params_mean = mean(params_inferred_all);
params_std = std(params_inferred_all);
% params_error = sqrt(sum(params_inferred_error_all.^2)/length(data))
n_boots = 100;
parms_std_boostrap = bootstrp(n_boots, @std, params_inferred_all);
params_SEM = std(parms_std_boostrap)./sqrt(n_boots-1);

params_CV = params_std./params_mean;
params_CV_error = params_SEM./params_mean;

errorbar(1:4, params_CV, params_CV_error,'o','LineWidth',2)

xlim([0 5]) 
xticks([1,2,3,4])
xticklabels({'K_{b}','\omega_{bp}','p','R_{max}'})
% yticks([0 100 200 300 400])
ylim([0 1])
% set(gca,'YScale','log')
% xlabel('')
ylabel('Coefficient of Variation')

box on
StandardFigure(gcf,gca)

% save the plots
% saveas(gcf,[FigPath,filesep, 'CV_parameters','.tif']);
% saveas(gcf,[FigPath,filesep, 'CV_parameters','.pdf']);

%% generate the plot for C.V. of inferred parameters (propagated error)
% calculate the Coefficient of Variation (C.V.) of inferred parameters
% across constructs (different enhancers).
for construct = 1:length(data)
    params_inferred_all(construct,:) = MCMC_6A0R_RuntNulls(construct).params_inferred;
    params_inferred_error_all(construct,:) = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
end

% calculate the mean and std of parameters "over constructs".
params_mean = mean(params_inferred_all);
params_std = std(params_inferred_all);
% params_error = sqrt(sum(params_inferred_error_all.^2)/length(data))

N_data = length(data); % number of constructs

% propagate the error from individual measurements
error_ind = params_inferred_error_all./params_inferred_all;
params_mean_error_propagated = sqrt(sum(error_ind.^2))/ N_data;
parms_std_error_propagated =  sqrt(sum(error_ind.^2*4))/ N_data;

params_CV = params_std./params_mean;
params_CV_error = params_CV.*sqrt(params_mean_error_propagated.^2 + parms_std_error_propagated.^2);

errorbar(1:4, params_CV, params_CV_error,'o','LineWidth',2)

xlim([0 5]) 
xticks([1,2,3,4])
xticklabels({'K_{b}','\omega_{bp}','p','R'})
% yticks([0 100 200 300 400])
ylim([0 1])
% set(gca,'YScale','log')
% xlabel('')
ylabel('Coefficient of Variation')

box on
StandardFigure(gcf,gca)

% save the plots
saveas(gcf,[FigPath,filesep, 'CV_parameters_propagated_error','.tif']);
saveas(gcf,[FigPath,filesep, 'CV_parameters_propagated_error','.pdf']);

%% generate a plot of Rate/Rate_max over AP axis for all construts : to see if the plots are collapsed
hold on
for construct = 1:length(data)
    Rate_null = compiledData{construct+9,9};
    Rate_null_SEM = compiledData{construct+9,10};
    
    % extract the R_max parameter from the MCMC inference
    R_max = MCMC_6A0R_RuntNulls(construct).params_inferred(4);
    
    Rate_null_norm = Rate_null./R_max;
    Rate_null_SEM_norm = Rate_null_SEM./R_max;
    
    errorbar(APaxis, Rate_null_norm, Rate_null_SEM_norm, 'o','LineWidth',2, 'Color',ColorChoice(construct,:))
end

xlim([0.2 0.5])
xticks([0.2 0.3 0.4 0.5])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1 1.2])

xlabel('embryo length')
ylabel('normalized rate')
legend('000','100','011','111','001','010','110','101', 'Location', 'NorthEast')

box on
StandardFigure(gcf,gca)

saveas(gcf,[FigPath,filesep,'normalized_raw_fits' ,'.tif']); 
saveas(gcf,[FigPath,filesep,'normalized_raw_fits' ,'.pdf']); 

%% generate a plot of Rate over AP axis for all construts 
hold on
for construct = 1:length(data)
    Rate_null = compiledData{construct+9,9};
    Rate_null_SEM = compiledData{construct+9,10};
    
%     errorbar(APaxis, Rate_null, Rate_null_SEM, 'o',)
    errorbar(APaxis, Rate_null, Rate_null_SEM, ...
                    'LineWidth',2, 'Color',ColorChoice(construct,:))
end

xlim([0.2 0.5])
xticks([0.2 0.3 0.4 0.5])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel('initial rate (AU/min)')
legend('000','100','011','111','001','010','110','101', 'Location', 'NorthEast')

box on
StandardFigure(gcf,gca)

saveas(gcf,[FigPath,filesep,'raw_initial_slope_errorbar','.tif']); 
saveas(gcf,[FigPath,filesep,'raw_initial_slope_errorbar','.pdf']); 

%% generate a plot of Rate over AP axis for all construts 
hold on
for construct = 1:length(data)
    Rate_null = compiledData{construct+9,9};
    Rate_null_SEM = compiledData{construct+9,10};
    
%     errorbar(APaxis, Rate_null, Rate_null_SEM, 'o','LineWidth',2, 'Color',ColorChoice(construct,:))
    shadedErrorBar(APaxis, Rate_null, Rate_null_SEM,...
                    'lineprops',{'color',ColorChoice(construct,:),'linewidth',2})
%     pause
end

xlim([0.2 0.5])
xticks([0.2 0.3 0.4 0.5])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel('initial rate (AU/min)')
legend('000','100','011','111','001','010','110','101', 'Location', 'NorthEast')

box on
StandardFigure(gcf,gca)

saveas(gcf,[FigPath,filesep,'raw_initial_slope_shadedError','.tif']); 
saveas(gcf,[FigPath,filesep,'raw_initial_slope_shadedError','.pdf']); 
    
%% save the result into mat files (MCMC_6A0R_RuntNulls)
FilePath = FigPath;
save([FilePath, filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat'],'MCMC_6A0R_RuntNulls')

end