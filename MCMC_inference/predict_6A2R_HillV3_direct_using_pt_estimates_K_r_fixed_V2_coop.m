%% predict_6A2R_cases
function predict_6A2R_HillV3_direct_using_pt_estimates_K_r_fixed_V2_coop
%% Description
% This script is to predict the level of hbP2 + 2Run sites

% input : 
% model type (parameters)
% dataset (constructs)
% 
% 
%% FilePaths
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\Point_estimate_Bcd_params_fromRunNulls_fixed_K_r';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr';

%% Load the data (Bcd, Run, and output)
% Bcd, Run, TF, TF_null
% APbin_start, APbin_end

%% Load the MCMC results
% 1) Bcd-dependent parameters from MCMC_6A0R_RuntNulls.mat
% 2) Run-dependent parameters from MCMC_6A1R_RuntWT.mat

% load the MCMC results for all construts with Runt null. (point estimates
% for the Bicoid-dependent params.)
% loads : "MCMC_6A0R_RuntNulls"
tempPath1 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls';
load([tempPath1, filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat']);

% load the MCMC results for the 6A1R constructs for 1) fixed K_r (global, shared) 
% and 2)w_rp (local) for each constructs.
% loads : "MCMC_6A1R_RuntWT"
% Note that the chain is constructed as "K_r(shared) w_rp1 w_rp2 w_rp3"
tempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\Point_estimate_Bcd_params_fromRunNulls_fixed_K_r';
load([tempPath2, filesep, 'MCMC_6A1R_RuntWT_params.mat'])

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

% extract the Runt-dependent parameters
K_r = MCMC_6A1R_RuntWT.params_inferred(1);
w_rp1 = MCMC_6A1R_RuntWT.params_inferred(2); %[100]
w_rp2 = MCMC_6A1R_RuntWT.params_inferred(3); %[001]
w_rp3 = MCMC_6A1R_RuntWT.params_inferred(4); %[010]

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
        params_6A2R = [params_Bcd, K_r, K_r, w_rp2, w_rp3];
    elseif i==2 % r2-close : [110]
        params_6A2R = [params_Bcd, K_r, K_r, w_rp1, w_rp3];
    elseif i==3 % r2-far : [101]
        params_6A2R = [params_Bcd, K_r, K_r, w_rp1, w_rp2];
    end
    
    Prediction(i,:) = model_6A2R_HillModel_V3_direct(params_6A2R, TF);
    fit_nulls(i,:) = model_6A0R_HillModel_V3(params_6A2R, TF);
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
for index = 1:3
    % define the construct index (which is consistent with the way it's
    % defiend in the compiledData.mat)
    constructIndex = [3,7,8];
    construct = constructIndex(index);

    % data (WT)
    Rate_WT = compiledData{construct+1,9};
    Rate_WT_SEM = compiledData{construct+1,10};
    
    errorbar(APaxis, Rate_WT, Rate_WT_SEM, 'o', 'Color', ColorChoice(construct,:),'LineWidth', 1)
    
    plot(APaxis(APbin_start:APbin_end), Prediction(index,:), 'Color', ColorChoice(construct,:),'LineWidth', 2)

    
end

% figure format
xlim([0.2 0.5])
xticks([0.2 0.3 0.4 0.5])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel('initial rate (AU/min)')

box on
legend('011','','110','','101','')
StandardFigure(gcf,gca)

FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr';
saveas(gcf,[FigPath,filesep,'raw_fits_RuntWT_6A2R_compiled','.tif']); 
saveas(gcf,[FigPath,filesep,'raw_fits_RuntWT_6A2R_compiled','.pdf']); 

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

%% (optional) calculate the RMSE values across constructs

%% Predict the Rates with Run-Run cooperativity
% model_6A2R_HillModelV3_RunCoop.m
% w_rr > 1 : stronger repression than independent Runt binding.
w_rr = [10^(-12), 10^(-6), 1, 10^6, 10^12, 10^24, 10^(48)];

% Predict the level for the Runt WT
% Use a custom-function for the model : model_6A2R_HillV3
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params

n_simu = length(MCMC_6A1R_RuntWT.chain);
n_burn = 0.5*n_simu;

% r2-new : [011]
% r2-close : [110]
% r2-far : [101]
% Predict the cases of 2 binding sites
Prediction_Run_Run_coop_direct = zeros(3, length(w_rr), length(TF));

% Loop over constructs (2 Run sites)
for i=1:3
    % step1. Bcd-dependent parameters
    construct = constructIndex(i);
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    if i==1 % r2-new : [011]
        params_6A2R = [params_Bcd, K_r, K_r, w_rp2, w_rp3];
    elseif i==2 % r2-close : [110]
        params_6A2R = [params_Bcd, K_r, K_r, w_rp1, w_rp3];
    elseif i==3 % r2-far : [101]
        params_6A2R = [params_Bcd, K_r, K_r, w_rp1, w_rp2];
    end
    
    for j=1:length(w_rr)
        omega_rr = w_rr(j);

        params_6A2R_Run_coop = [params_6A2R, omega_rr];

        Prediction_Run_Run_coop_direct(i,j,:) = model_6A2R_HillModel_V3_competition_RunCoop(params_6A2R_Run_coop, TF);
    end
    
    % Runt null data (MCMC fit)
%     fit_nulls(i,:) = model_6A0R_HillModel_V3(params_6A2R, TF);
end

%% generate plots for Runt-Runt cooperativity

% pick a colormap that is linearly scaled
cmap = colormap('viridis');

% w_rr : take logarithm
log_w_rr = log10(w_rr);
w_rr_range = max(log_w_rr) - min(log_w_rr);

for i=1:3
    % construct define
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
    
    % color scale
    caxis = [min(log_w_rr) max(log_w_rr)];
   
    % Run-Run interaction strength
    for j=1:length(w_rr)
        cIndex = floor((log10(w_rr(j))-(min(log_w_rr)))/w_rr_range*256);
        if cIndex ==0
            cIndex = cIndex + 1;
        end
    
        color = cmap(cIndex,:);
        % Prediction
        % Runt null
        plot(APaxis(9:21),fit_nulls(i,:),'Color',ColorChoice(4,:), 'LineWidth', 2)
        % Runt WT
        plot(APaxis(9:21), squeeze(Prediction_Run_Run_coop_direct(i,j,:)),'Color',color, 'LineWidth', 2)
%         pause
    end

    xlim([0.2 0.5])
    xticks([0.2 0.3 0.4 0.5])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
    ylabel({'initial RNAP', 'loading rate (AU/min)'})
    
    box on
    colorbar
%     set(gca,'colorscale','log')
%     cbh = colorbar('XTickLabel',{'10^(-12)','10^(-6)','1','10^6','10^12','10^24','10^(48)'},...
%                     'XTick', [-12, -6, 0, 6, 12, 24, 48]);
    
    StandardFigure(gcf,gca)

%   %Save the plot
%     saveas(gcf,[FigPath,filesep,'Prediction_HillV3_',constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,'Prediction_HillV3_',constructNames{construct},'.pdf']);
    pause
end


%% 
%% Modification 2: Higher-order cooperativity between Run's action
% w term for when the two Run molecules are bound, and exerting some effect to
% the RNAP.

%% Predict the Rates with higher order cooperativity between Run-RNAP
% mnodel_6A2R_HillModelV3_RunCoop.m
% w_rr < 1 : stronger repression than independent Runt effect.
w_ho = 10.^([-10 -5  0 5 10]);

% Predict the level for the Runt WT
% Use a custom-function for the model : model_6A2R_HillV3
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params

n_simu = 20000;
n_burn = 0.5*n_simu;

% r2-new : [011]
% r2-close : [110]
% r2-far : [101]
% Predict the cases of 2 binding sites
Prediction_higher_coop = zeros(3, length(w_ho), length(TF));

for i=1:3
    construct = constructIndex(i);
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    if i==1 % r2-new : [011]
        params_6A2R = [params_Bcd, K_r, K_r, w_rp2, w_rp3];
    elseif i==2 % r2-close : [110]
        params_6A2R = [params_Bcd, K_r, K_r, w_rp1, w_rp3];
    elseif i==3 % r2-far : [101]
        params_6A2R = [params_Bcd, K_r, K_r, w_rp1, w_rp2];
    end
    
    for j=1:length(w_ho)
        omega_ho = w_ho(j);

        params_fit = [params_6A2R, omega_ho];

        Prediction_higher_coop(i,j,:) = model_6A2R_HillModel_V3_direct_HigherCoop(params_fit, TF);
    end
end



%% generate plots showing the effect of higher order cooperativity between Run-RNAP
cmap = colormap('viridis');
% caxis([-10 10])
w_ho_range = max(log10(w_ho)) - min(log10(w_ho));

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
    
%     cbar = linspace(0, 256, 22);

    

    % Run-Run interaction strength
    for j=1:length(w_ho)
        cIndex = floor((log10(w_ho(j))-(min(log10(w_ho))))/w_ho_range*256);
        if cIndex ==0
            cIndex = cIndex + 1;
        end
    
        color = cmap(cIndex,:);
        % Prediction
        % Runt null
        plot(APaxis(9:21), fit_nulls(i,:),'Color',ColorChoice(4,:), 'LineWidth', 2)
        % Runt WT
        plot(APaxis(9:21), squeeze(Prediction_higher_coop(i,j,:)),'Color',color, 'LineWidth', 2)

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
%     saveas(gcf,[FigPath,filesep,'Prediction_HillV3_HighCoop_',constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep, 'Prediction_HillV3_HighCoop_',constructNames{construct},'.pdf']);
    pause
end
end