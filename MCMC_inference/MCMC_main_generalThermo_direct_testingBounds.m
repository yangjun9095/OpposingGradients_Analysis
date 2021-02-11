%% MCMC_initial_slope_DirectRepression_test_Bounds_Priors


% file path
FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';

%% Default settings (also make it optional)
% fileDir = pwd;
% saveLoc = pwd;
numParPools = 8;
% n_burn = 5000;
n_steps = 2*10^5;
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


% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};
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

%%  Pick a model for the fitting
mdl = @(TF, params) model_6A1R_direct_repression_V2(params, TF);

%leaving this here in case it'll be useful in the future
model.ssfun = @(params, data) sum((data.ydata-mdl(data.xdata, params)).^2);

model.modelfun = mdl;  %use mcmcrun generated ssfun 

%% Define the parameters for the model
% put the initial parameters and bounds in a form that the mcmc function
% accepts
names = {'K_{b}','K_{r}','w_{b}','w_{bp}','w_{rp}','p','R_{max}'};
params = cell(1, length(names));

% Initialize the MCMC parameters.
Kb0 = 20;   % 100*rand;
Kr0 = 5;     % 100*rand;
w_b0 = 20;    % 10*rand;
w_bp0 = 2;   % 10*rand;
% repression (0< w <1)
w_rp0 = 0.5;
p0 = 0.1; %
R_max0 = MCMCdata.R_max; %500*rand;


params0 = [Kb0, Kr0, w_b0, w_bp0, w_rp0, p0, R_max0];

% Define the vector for the upper/lower limits of the parameters, for now,
% we will only vary the K_{b} and K_{r}.
maxBound = [10^2, 10^3, 10^4, 10^5, 10^6, 10^7];

for j=1:length(maxBound)
    % Bounds of the parameters
    LB = [0.1, 0.1, 1, 1, 10^(-6), 0, 50];
    UB = [maxBound(j), 10^2, 10^3, 10^3, 1, 10, 10^3];

    for i = 1:length(names)
        pri_mu = NaN; %default prior gaussian mean
        pri_sig = Inf; %default prior gaussian variance
        localflag = 0; %is this local to this dataset or shared amongst batches?
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.

        if i==7
            pri_mu = MCMCdata.R_max;
            pri_sig = 20;
    %         targetflag = 0; % Fix this parameter
        else
    %         pri_mu = NaN;
    %         pri_sig = Inf;
        end
        params{1, i} = {names{i},params0(i), LB(i), UB(i), pri_mu, pri_sig, targetflag, localflag};

    end

    %This is the variance in the parameter proposal distribution. Change these
    %to change the proposal acceptance rate, or if the convergence doesn't look
    %good.

    Kb_step = 10;
    Kr_step = 1;
    w_b_step = 0.1;
    w_bp_step = 0.1;
    p_step = 0.001;
    R_max_step = 1;

    % repression
    w_rp_step = 0.01;


    % Initialize the covariance matrix
    J0 = diag([Kb_step, Kr_step, w_b_step, w_bp_step,...
        w_rp_step, p_step, R_max_step]);


    %% MCMC - Options
    options = [];
    n_steps = 2*10^5;
    options.nsimu = n_steps; %n_steps; %Number of steps
    options.updatesigma = 1; %Update error variance
    options.qcov = J0; %Initial covariance
    % options.burnintime = 0.5*n_steps; %Burn in time
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

    %we're gonna run this three times and use the initial results of one
    %run as conditions for the next. this is an alternative when common least
    %squares gives results too poor to initialize with
    results = [];
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    %[results,~,~,~]=mcmcrun(model,MCMCdata,params,options,results);
    [results,chain,s2chain,~]=mcmcrun(model,MCMCdata,params,options,results);
    
    MCMCResults_testBounds(j).results = results;
    MCMCResults_testBounds(j).chain = chain;
    MCMCResults_testBounds(j).maxBound = maxBound(j);
end

%% Check the MCMC result
figure
for j=3:length(maxBound)

    chain = MCMCResults_testBounds(j).chain;
    results = MCMCResults_testBounds(j).results;
    %% Diagnose the MCMC result
    % stats = chainstats(chain,results);
    n_burn = 0.5*n_steps;
    chainstats(chain(n_burn+1:end,:),results);

    %% generate corner plots
    clf
    n_burn = 0.5*n_steps; %0.5*n_steps;% 0.5*10^4;
        m = [chain(n_burn:end,1), chain(n_burn:end,2), chain(n_burn:end,3), chain(n_burn:end,4), chain(n_burn:end,5), chain(n_burn:end,6)];%, chain(n_burn:end,7)];
        corner = figure;
        names = {'Kb','Kr','w_b','w_{bp}','w_{rp}','p','R_{max}'};
        ecornerplot(m,'names',names);

    % Save the plot
        FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_ThermoModelV2';
        saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct},'max_bound_K_b_' num2str(maxBound(j))  ,'.tif']); 
        saveas(gcf,[FigPath,filesep,'Corner_plot_', constructNames{construct},'max_bound_K_b_' num2str(maxBound(j))  ,'.pdf']); 
end

%% Extract the parameters from different bounds
for j=1:length(maxBound)
    %% pull the chain and results from the main cell/structure
    chain = MCMCResults_testBounds(j).chain;
    results = MCMCResults_testBounds(j).results;
   
    
    %% Extract chain results into individual parameters
    n_burn = n_steps*0.5;
    params_fit(j,:) = mean(chain(n_burn+1:end,:));
    params_fit_std(j,:) = std(chain(n_burn+1:end,:));

end

%% Color module
% This is defining the line color
% We have 8 distinct datasets, with or without Runt protein.
% I think selecting 8 distinguishable color sets, then changing the
% brightness by either adding/subtracting white would be a better idea than
% selecting 16 different color sets.

colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.purple = [171,133,172]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

%colorDict.magenta = [208,109,171]/255;
%colorDict.lightBlue = [115,142,193]/255;
colorDict.lightgreen = [205,214,209]/255;
colorDict.pink = [232,177,157]/255;
colorDict.thickpink = [132,27,69]/255;

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.darkgreen; colorDict.thickpink]; 

% For now, I'll add white (color+[1 1 1])/2 to make thinner color (for the
% Runt nulls)
%% generate MCMC fits for the raw initial slope data

figure
hold on
% 1) Data
construct = 5;
% Runt Null
errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10}, 'o', 'Color',ColorChoice(4,:),'LineWidth', 1)
% Runt WT
errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10}, 'o', 'Color',ColorChoice(1,:),'LineWidth', 1)

% 2) MCMC fits
for j=1:length(maxBound)

    fit = model_6A1R_direct_repression_V2(params_fit(j,:), MCMCdata.xdata);
    
    % define the color index, note that we have to avoid 1 or 4 as those
    % are used for the data plots.
    if j==1
        cindex = 7;
    elseif j==4
        cindex = 8;
    else
        cindex = j;
    end
    % Runt Null
    plot(APaxis(APbinRange), fit(1:length(APbinRange)),'Color',ColorChoice(cindex,:),'LineWidth',2)
    % Runt WT
    plot(APaxis(APbinRange), fit(1+length(APbinRange):end),'Color',ColorChoice(cindex,:),'LineWidth',2)

    % plot(APaxis(APbinRange), MCMCdata.ydata(1:length(APbinRange)), 'Color',ColorChoice(4,:))
    % plot(APaxis(APbinRange), MCMCdata.ydata(1+length(APbinRange):end),'Color',ColorChoice(1,:))

    xlim([0.2 0.6])
    xticks([0.2 0.3 0.4 0.5 0.6])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
    ylabel({'initial RNAP', 'loading rate (AU/min)'})
    legend('Runt null','Runt WT')

    box on
    StandardFigure(gcf,gca)
    pause
    % Save the plot
    % saveas(gcf,[FigPath,filesep,'raw_data_slope_fits_', construct ,'.tif']); 
    % saveas(gcf,[FigPath,filesep,'raw_data_slope_fits_', construct ,'.pdf']);  

end
%% Optional - hierarchical running of the MCMC
% In here, we use the fact that most of the parameters converge easily
% except the "p".
% We have two ways : One is constraining the p further using the relation
% of R_max and R_min, the other is using the hierarchical inference to
% give a prior from the previous inference on the other parameters ,then,
% re-run the MCMC just for the "p". In here, we will do the second
% approach.

% Inferred parameters
% params_inferred = [mean_Kb, mean_Kr, mean_w_b, mean_w_bp, mean_w_rp, mean_p, mean_R_max];
params_inferred = [mean_Kb, mean_Kr, mean_w_b, mean_w_bp, mean_w_rp, mean_p, mean_R_max];

for i = 1:length(names)
    
    targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.
    
%     if ~isnan(fixedKD) && contains(names(i), "KD")
%         k0(i) = fixedKD;
%     end
    if ~i==6
        pri_mu(i) = mean(n_burn:end,i);
        pri_sig(i) = 0.1*pri_mu(i);%std(n_burn:end,i);
%         targetflag(i) = 1;
    elseif i==6
        pri_mu(i) = mean_p;
        pri_sig(i) = 0.1*mean_p;
%         targetflag(i) = 0;
    end
%     if i==7
%         targetflag = 0;
%     end
    params{1, i} = {names{i},params_inferred(i), LB(i), UB(i), pri_mu(i), pri_sig(i), targetflag, localflag};
    
end


Kb_step = 0.1*sigma_Kb;
Kr_step = 0.1*sigma_Kr;
w_b_step = 0.1*sigma_w_b;
w_bp_step = 0.1*sigma_w_bp;
p_step = 0.1*sigma_p;
R_max_step = 0.1;

% repression
w_rp_step = 0.1*sigma_w_rp;


% Initialize the covariance matrix
J0 = diag([Kb_step, Kr_step, w_b_step, w_bp_step,...
    w_rp_step, p_step, R_max_step]);
%% calculate the mean and std of the inferred parameters (using the chains)
%% generate plots of fits
% mean_R_max = R_max0;

params_inferred = [mean_Kb, mean_Kr, mean_w_b, mean_w_bp, mean_w_rp, mean_p, mean_R_max];
% params_lsq = FitResult(2).params_fit;
output = model_6A1R_direct_repression_V2(params_inferred, MCMCdata.xdata);

figure
hold on

% Data
construct = 6;
% Runt Null
errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10}, 'o', 'Color',ColorChoice(4,:),'LineWidth', 1)
% Runt WT
errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10}, 'o', 'Color',ColorChoice(1,:),'LineWidth', 1)

% Runt Null
plot(APaxis(APbinRange), output(1:length(APbinRange)),'Color',ColorChoice(4,:),'LineWidth',2)
% Runt WT
plot(APaxis(APbinRange), output(1+length(APbinRange):end),'Color',ColorChoice(1,:),'LineWidth',2)

% plot(APaxis(APbinRange), MCMCdata.ydata(1:length(APbinRange)), 'Color',ColorChoice(4,:))
% plot(APaxis(APbinRange), MCMCdata.ydata(1+length(APbinRange):end),'Color',ColorChoice(1,:))

xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})
legend('Runt null','Runt WT')

box on
StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,'raw_data_slope_fits_', construct ,'.tif']); 
% saveas(gcf,[FigPath,filesep,'raw_data_slope_fits_', construct ,'.pdf']);  

%% Calculate the output range using the mcmcpred
% MCMCdata.xdata
out = mcmcpred(results,chain(n_burn:end,:),[],MCMCdata.xdata, mdl);

%% MCMC prediction plot

% % plothandle = mcmcpredplot(out, MCMCdata.ydata, [])
% 
% % mcmcpredplot(out);
% nn = (size(out.predlims{1}{1},1) + 1) / 2;
% plimi = out.predlims{1};
% yl = plimi{1}(1,:);
% yu = plimi{1}(2*nn-1,:);
% 
% % km = mean(MCMCchain);
% % ks = std(MCMCchain);
%     
% 
% yf = plimi{1}(nn,:);
%     
% yy = yf;
% yyl = yl;
% yyu = yu;
% 
% index1 = 1:11;
% index2 = 12:22;
% 
% figure
% hold on
% % Runt null
% errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10}, 'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
% % shadedErrorBar(APaxis(APbins_fit), yf(1,index1), [yu(1,index1);  yl(1,index1)],'lineProps',{'markerfacecolor',ColorChoice(4,:)})
% plot(APaxis(APbins_fit), yf(1,index1), 'Color',ColorChoice(4,:))
% plot(APaxis(APbins_fit), yl(1,index1),'Color',(ColorChoice(4,:)+[1 1 1])/2)
% plot(APaxis(APbins_fit), yu(1,index1),'Color',(ColorChoice(4,:)+[1 1 1])/2)
% 
% % Runt WT
% errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10}, 'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
% % shadedErrorBar(APaxis(APbins_fit), yf(1,index2), [yu(1,index2);  yl(1,index2)],'lineProps',{'markerfacecolor',ColorChoice(1,:)})
% plot(APaxis(APbins_fit), yf(1,index2),'Color',ColorChoice(1,:))
% plot(APaxis(APbins_fit), yl(1,index2),'Color',(ColorChoice(1,:)+[1 1 1])/2)
% plot(APaxis(APbins_fit), yu(1,index2),'Color',(ColorChoice(1,:)+[1 1 1])/2)                                 
% 
% xlim([0.2 0.6])
% xticks([0.2 0.3 0.4 0.5 0.6])
% ylim([0 400])
% yticks([0 100 200 300 400])
% 
% xlabel('embryo length')
% ylabel({'initial RNAP', 'loading rate (AU/min)'})
% 
% StandardFigure(gcf,gca)
% 
% %Save the plot
% % saveas(gcf,[FigPath,filesep,'MCMC_fit_direct_001','.tif']); 
% % saveas(gcf,[FigPath,filesep,'MCMC_fit_direct_001','.pdf']); 

%% Optional (trial and error step to settle down the initial conditions, etc.)
setNum = setNum+1;
MCMCTemp(setNum).chain = chain;
MCMCTemp(setNum).results = results;
MCMCTemp(setNum).params0 = params0;
MCMCTemp(setNum).LB = LB;
MCMCTemp(setNum).UB = UB;
MCMCTemp(setNum).J0 = J0;

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
DatasetName = '001';
% filename = [date,'-',DatasetName];
save([FilePath,filesep,'MCMCresult_2',DatasetName,'.mat'],'results','chain','MCMCdata');
% save([FilePath,filesep,'MCMCresult_',DatasetName,'.mat'],'MCMCresults','MCMCplot','DatasetName');

%MCMC raw chains
% filename = [date,'-',DatasetName,'_RawChain'];
% save([FilePath,filesep,'MCMCchain_',DatasetName,'.mat'],'MCMCchain');

% disp(['MCMC analysis complete. Information stored in: ',saveLoc]);