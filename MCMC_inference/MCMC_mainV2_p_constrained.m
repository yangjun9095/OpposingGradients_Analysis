%% MCMC_DirectRepression_test_Bounds_Priors

%% Default settings (also make it optional)
% file path
FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';

% fileDir = pwd;
% saveLoc = pwd;
numParPools = 16;
% n_burn = 5000;
n_steps = 10^5;      
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

%% Pull the construct
% Pick the dataset from the data.mat
construct = 5;
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

% Fine the common indices of APbinRange that are not NaNs for both WT and
% Null
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
MCMCdata.Bcd = [Bcd ; Bcd];
MCMCdata.Run = [RunNull ; Run];
MCMCdata.xdata = [MCMCdata.Bcd, MCMCdata.Run];
MCMCdata.R_max = max(MCMCdata.ydata);
MCMCdata.R_min = min(MCMCdata.ydata);

%% construct a structure to save the inference output
% Here, we will save the fields of chain, results, bounds, prior, burn-in
% time, construct, etc.
MCMC_output = struct;

%%  Define the model and SS for the MCMC
mdl = @(TF, params) model_6A1R_direct_repression_V3_R_p_constrained(params, TF);

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, data) sum( (data.ydata-mdl(data.xdata, params)) .^2);

model.modelfun   = mdl;  %use mcmcrun generated ssfun 

%% Define the parameters for the model
% put the initial parameters and bounds in a form that the mcmc function
% accepts
names = {'K_{b}','K_{r}','w_{b}','w_{bp}','w_{rp}','R_{max}','R_{min}'};
params = cell(1, length(names));

% Initialize MCMC parameters.
Kb0 = 10;   % 100*rand;
Kr0 = 10;     % 100*rand;
w_b0 = 2;    % 10*rand;
w_bp0 = 5;   % 10*rand;
% repression (0< w <1)
w_rp0 = 0.2;
% p0 = 0.001; %
R_max0 = MCMCdata.R_max;
R_min0 = MCMCdata.R_min;


params0 = [Kb0, Kr0, w_b0, w_bp0, w_rp0, R_max0, R_min0];

% Bounds of the parameters
LB = [1, 1, 1, 1, 10^(-3), 0, 0];
UB = [10^2, 10^2, 10^3, 10^3, 1, 10^3, 10^3];
% UB = [Inf, Inf, Inf, Inf, 1, R_max0, R_min0];

pri_mu = NaN; %default prior gaussian mean
pri_sig = Inf; %default prior gaussian variance
localflag = 0; %is this local to this dataset or shared amongst batches?

for i = 1:length(names)
    
    targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.
    
%     if ~isnan(fixedKD) && contains(names(i), "KD")
%         k0(i) = fixedKD;
%         targetflag = 0;
%     end
    if i==1
%         pri_mu = 50;
%         pri_sig = 100;
    elseif i==2
%         pri_mu = 50;
%         pri_sig = 100;
    elseif i==5
%         pri_mu = 0.5;
%         pri_sig = 2;
%         disp(['putting a weak prior on w_rp',' \mu =', pri_mu,' \sigma =', pri_sig])
    elseif i==6 || i==7
        targetflag = 0; % fix the parameters to a constant value
    end
    params{1, i} = {names{i}, params0(i), LB(i), UB(i), pri_mu, pri_sig, targetflag, localflag};
    
end

%This is the variance in the parameter proposal distribution. Change these
%to change the proposal acceptance rate, or if the convergence doesn't look
%good.

Kb_step = 0.1;
Kr_step = 0.1;
w_b_step = 0.1;
w_bp_step = 0.1;
% p_step = 0.0001;
R_max_step = 0.1;
R_min_step = 0.1;

% repression
w_rp_step = 0.01;


% Initialize the covariance matrix
J0 = diag([Kb_step, Kr_step, w_b_step, w_bp_step,...
    w_rp_step, R_max_step, R_min_step]);


%% MCMC - Options
options = [];
n_steps = 2*10^4; %n_steps; %Number of steps
options.nsimu = n_steps;
options.updatesigma = 1; %Update error variance
options.qcov = J0; %Initial covariance
% options.burnintime = n_burn; %Burn in time
options.adaptint = 100;
options.method = 'dram';
options.verbosity = 0; %Decrease text output

% AR options
% options.drscale = 5; % a high value (5) is important for multimodal parameter spaces
% options.waitbar = wb; %the waitbar is rate limiting sometimes
% options.nsimu = nSimu; %should be between 1E3 and 1E6
% options.updatesigma = 1; %honestly don't know what this does

%% Run the MCMC (this whole block can be inserted inside the for loop above,
% to run the MCMC for different constructs.

% initialize a structure to save the data
% MCMC_result = struct;

%we're gonna run this three times and use the initial results of one
%run as conditions for the next. this is an alternative when common least
%squares gives results too poor to initialize with 
results = [];
[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
%[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
[results,chain,s2chain,~]=mcmcrun(model,MCMCdata,params,options,results);

% Set up the structure to save the results
% MCMC_result(k).results = results;
% MCMC_result(k).chain = chain;
% MCMC_result(k).s2chain = s2chain;
%% Diagnose the MCMC result
chainstats(chain,results)

%% generate corner plots
n_burn = 0.5*n_steps;% 0.5*10^4;
m = [chain(n_burn:end,1), chain(n_burn:end,2), chain(n_burn:end,3), chain(n_burn:end,4), chain(n_burn:end,5)];
corner = figure;
names = {'Kb','Kr','w_b','w_{bp}','w_{rp}'};
ecornerplot(m,'names',names);

%% Optional (trial and error step to settle down the initial conditions, etc.)
setNum = setNum+1;
MCMCTemp(setNum).chain = chain;
MCMCTemp(setNum).results = results;
MCMCTemp(setNum).params0 = params0;
MCMCTemp(setNum).LB = LB;
MCMCTemp(setNum).UB = UB;
MCMCTemp(setNum).J0 = J0;

%% Plot for checking the fit
%% Calculate the model output using the mean parameter values
% This is assuming a unimodal distribution of parameters.

%% Calculate the output range using the mcmcpred
% MCMCdata.xdata
% R_max_chain = R_max0*ones(n_steps-n_burn, 1);
% R_min_chain = R_min0*ones(n_steps-n_burn, 1);
% 
% MCMCchain = chain(n_burn+1:end,:);
% 
% MCMCchain(:,6) = R_max_chain;
% MCMCchain(:,7) = R_min_chain;

out = mcmcpred(results,chain(n_burn:end,:),[],MCMCdata.xdata, mdl);

%% MCMC prediction plot

% plothandle = mcmcpredplot(out, MCMCdata.ydata, [])

% mcmcpredplot(out);
nn = (size(out.predlims{1}{1},1) + 1) / 2;
plimi = out.predlims{1};
yl = plimi{1}(1,:);
yu = plimi{1}(2*nn-1,:);

% km = mean(MCMCchain);
% ks = std(MCMCchain);
    

yf = plimi{1}(nn,:);
    
yy = yf;
yyl = yl;
yyu = yu;

index1 = 1:11;
index2 = 12:22;

figure
hold on
% Runt null
errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10}, 'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
% shadedErrorBar(APaxis(APbins_fit), yf(1,index1), [yu(1,index1);  yl(1,index1)],'lineProps',{'markerfacecolor',ColorChoice(4,:)})
plot(APaxis(APbins_fit), yf(1,index1), 'Color',ColorChoice(4,:))
plot(APaxis(APbins_fit), yl(1,index1),'Color',(ColorChoice(4,:)+[1 1 1])/2)
plot(APaxis(APbins_fit), yu(1,index1),'Color',(ColorChoice(4,:)+[1 1 1])/2)

% Runt WT
errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10}, 'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
% shadedErrorBar(APaxis(APbins_fit), yf(1,index2), [yu(1,index2);  yl(1,index2)],'lineProps',{'markerfacecolor',ColorChoice(1,:)})
plot(APaxis(APbins_fit), yf(1,index2),'Color',ColorChoice(1,:))
plot(APaxis(APbins_fit), yl(1,index2),'Color',(ColorChoice(1,:)+[1 1 1])/2)
plot(APaxis(APbins_fit), yu(1,index2),'Color',(ColorChoice(1,:)+[1 1 1])/2)                                 

xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})

StandardFigure(gcf,gca)

%Save the plot
% saveas(gcf,[FigPath,filesep,'MCMC_fit_direct_001','.tif']); 
% saveas(gcf,[FigPath,filesep,'MCMC_fit_direct_001','.pdf']); 
%% Set up data to be saved in structure
MCMCchain = struct('Kb_chain',{},'Kr_chain',{},'w_b_chain',{},'w_bp_chain',{},...
    'w_rp_chain',{},'p_chain',{},'R_max_chain',{});
MCMCresults = struct('mean_Kb',{},'sigma_Kb',{},'mean_Kr',{},'sigma_Kr',{},...
    'mean_w_b',{},'sigma_w_b',{},'mean_w_bp',{},'sigma_w_bp',{},...
    'mean_w_rp',{},'sigma_w_rp',{},'mean_p',{},'sigma_p',{},'mean_R_max',{},'sigma_R_max',{});
MCMCplot = struct('Rate_plot',{},'simRate',{});


%% Extract chain results into individual parameters
n_burn = n_simu*0.25;
Kb_chain = chain(n_burn:end,1);
Kr_chain = chain(n_burn:end,2);
w_b_chain = chain(n_burn:end,3);
w_bp_chain = chain(n_burn:end,4);
w_rp_chain = chain(n_burn:end,5);
p_chain = chain(n_burn:end,6);
R_max_chain = chain(n_burn:end,7);

%Mean/STD results
mean_Kb = mean(Kb_chain);
sigma_Kb = std(Kb_chain,1);

mean_Kr = mean(Kr_chain);
sigma_Kr = std(Kr_chain,1);

mean_w_b = mean(w_b_chain);
sigma_w_b = std(w_b_chain,1);

mean_w_bp = mean(w_bp_chain);
sigma_w_bp = std(w_bp_chain,1);

mean_w_rp = mean(w_rp_chain);
sigma_w_rp = std(w_rp_chain,1);

mean_p = mean(p_chain);
sigma_p = std(p_chain,1);

mean_R_max = mean(R_max_chain);
sigma_R_max = std(R_max_chain,1);


mean_sigma = sqrt(mean(s2chain));
sigma_sigma = std(sqrt(s2chain),1);

%% 
%% generate plots of inferred parameters
params_inferred = [mean_Kb mean_Kr mean_w_b mean_w_bp mean_w_rp mean_p mean_R_max];

params_inferred_STD = [sigma_Kb sigma_Kr sigma_w_b sigma_w_bp sigma_w_rp sigma_p sigma_R_max];

% params_comparison = figure;
hold on
construct = 5;
errorbar(params_inferred,params_inferred_STD,'o','Color', ColorChoice(construct,:),'MarkerFaceColor', ColorChoice(construct,:))
errorbar(params_fit_001, params_fit_SD_001,'o','Color', ColorChoice(2,:),'MarkerFaceColor', ColorChoice(2,:))
% errorbar(params_fit_100, STD_100,'o','Color', ColorChoice(2,:),'MarkerFaceColor', ColorChoice(2,:))
% errorbar(params_fit_010, STD_010,'o','Color', ColorChoice(6,:),'MarkerFaceColor', ColorChoice(6,:))
% errorbar(params_fit_001, STD_001,'o','Color', ColorChoice(5,:),'MarkerFaceColor', ColorChoice(5,:))


% legend('100','001','010')
legend('nonlinear least squared', 'MCMC')


xlim([0 8])
xticklabels({'','K_{b}','K_{r}','\omega_{b}','\omega_{bp}','\omega_{rp}','p','R_{max}',''})
xlabel('parameters')
ylabel('inferred values')
set(gca, 'YScale','log')

StandardFigure(gcf,gca)

% Save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_inference_ThermoModelV2';
saveas(gcf,[FigPath,filesep,'MCMC_lsq_Params_direct_001','.tif']); 
saveas(gcf,[FigPath,filesep,'MCMC_lsq_Params_direct_001','.pdf']); 

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
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_inference_ThermoModelV2\MCMCresults';
%MCMC results and plots
% DatasetName = constructName(construct);
DatasetName = '001'
% filename = [date,'-',DatasetName];
save([FilePath,filesep,'MCMCresult_',DatasetName,'.mat'],'MCMCresults','MCMCplot','DatasetName');

%MCMC raw chains
filename = [date,'-',DatasetName,'_RawChain'];
save([FilePath,filesep,'MCMCchain_',DatasetName,'.mat'],'MCMCchain');

% disp(['MCMC analysis complete. Information stored in: ',saveLoc]);