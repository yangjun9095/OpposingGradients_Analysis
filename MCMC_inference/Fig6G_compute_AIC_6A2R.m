function Fig6G_compute_AIC_6A2R
%% Description
% This script is to compute the Akaike Information Criterion (AIC) from
% 6A2R constructs with different sets of models with different number,
% identity of parameters.

% models: Hill.v3 with with zero free parameter, with w_rr, with w_ho, or with both.

% data: 6A2R initial rate of transcription (Run WT) datasets from [101],
% [011], and [110]
% r2-new : [011] : 3rd in the MCMC_6A0R_RuntNulls
% r2-close : [110] : 7th
% r2-far : [101] : 8th

%% load the Bicoid and Runt datasets
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';
% Load the output data
preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
data = preprocessedData.data;
% Load the input TF data
load([FilePath, filesep, 'TFinput.mat'])

Bicoid = TFinput(:,1);
Runt = TFinput(:,2);
RuntNull = TFinput(:,3);

%% import the data and data_sigma 
% (this is consistent regardless of the choice of the models)

% index for the constructs
constructIndex = [3,7,8];
% AP axis
APaxis = 0:0.025:1;
APbin_start=9; % 20% of embryo length
APbin_end=21; % 50% of embryo length
Bcd = Bicoid(APbin_start:APbin_end); 
Run = Runt(APbin_start:APbin_end); 
RunNull = RuntNull(APbin_start:APbin_end);

% we will concatenate the vectors on the order of [011],[110], and [101].

% initialize the vectors for 1)rate, 2)rate_std
rate_vec = [];
rate_std_vec = [];

% loop over all 6A2R constructs
for i=1:3
    index = constructIndex(i);
    rate_wt = data(index).Rate_WT;
    rate_wt_std = nanstd(data(index).Rate_WT_individual,0,2);

    % truncate for 20-50% of embryo length
    rate_wt = rate_wt(APbin_start:APbin_end);
    rate_wt_std = rate_wt_std(APbin_start:APbin_end);

    % concatenate to one vector
    rate_vec = [rate_vec; rate_wt];
    rate_std_vec = [rate_std_vec; rate_wt_std];
end

% for the posterior AP bins where there is either small variability, or
% lack of data points for some embryos, we might have std=0. In those case,
% we don't want those values to be over-represented, thus we will just put
% nan here.
rate_std_vec(rate_std_vec==0)=nan; % we cannot have the nominator to be 0

%% import the parameters (6AnR_RuntNulls - params_Bcd, 6A1R - params_Run)
% import K_b, w_bp, p, R as well as K_r, w_rp1, w_rp2
tempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A1R_Bcd_params_fromRuntNulls_fixed_Kr';
load([tempPath2, filesep, 'MCMC_6A1R_RuntWT_params.mat'])

% extract the Runt-dependent parameters
K_r = MCMC_6A1R_RuntWT.params_inferred(1);
w_rp1 = MCMC_6A1R_RuntWT.params_inferred(2); %[100]
w_rp2 = MCMC_6A1R_RuntWT.params_inferred(3); %[001]
w_rp3 = MCMC_6A1R_RuntWT.params_inferred(4); %[010]
%% 1. 6A2R with zero free parameter prediction using K_r, w_rp1, w_rp2
% generate the model's expected values using a set of parameters
model_expected = [];

% loop over all 6A2R constructs
% we will concatenate the vectors (model_expected) on the order of [011],[110], and [101].
for index=1:3
    construct = constructIndex(index);
    % import the relevant parameters
    % : 6AnR_RuntNulls - params_Bcd, 6A1R - params_Run
    % 1) parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    % Runt parameters
    if index==1 % r2-new : [011]
        params_Run = [K_r, K_r, w_rp2, w_rp3];
    elseif index==2 % r2-close : [110]
        params_Run = [K_r, K_r, w_rp1, w_rp3];
    elseif index==3 % r2-far : [101]
        params_Run = [K_r, K_r, w_rp1, w_rp2];
    end

    % We will use the Bcd and Run to generate the expected data for 20-50% of embryo length
    rate_expected = model_6A2R_HillModel_V3_direct([params_Bcd, params_Run], TF);

    % concatenate to one vector
    model_expected = [model_expected; rate_expected];
end

% remove the elements that are NaNs in either rate_vec or rate_vec_std (and
% their matching indices in model_expected).

% noNaNindices = ~isnan(rate_vec) & ~isnan(rate_std_vec);
% rate_vec = rate_vec(noNaNindices);
% rate_std_vec = rate_std_vec(noNaNindices);
% model_expected = model_expected(noNaNindices);
rate_std_vec(rate_std_vec==0)=1; % we cannot have the nominator to be 0


AIC_params_free = compute_AIC(model_expected, rate_vec, rate_std_vec, 0);
%% 2. 6A2R with w_rr as the free parameter

% load the 6A2R MCMC result (the loaded varaiable name is MCMC_6A2R_RuntWT)
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\Run_Run_coop_only\MCMC_6A2R_RuntWT_params.mat');

% initialize the model output vector
model_expected = [];
% loop over all 6A2R constructs
% we will concatenate the vectors (model_expected) on the order of [011],[110], and [101].
for index=1:3
    construct = constructIndex(index);
    % import the relevant parameters
    % : 6AnR_RuntNulls - params_Bcd, 6A1R - params_Run
    % 1) parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    w_rr = MCMC_6A2R_RuntWT(index).params_inferred;
    % Runt parameters
    if index==1 % r2-new : [011]
        params_Run = [K_r, w_rp2, w_rp3];
    elseif index==2 % r2-close : [110]
        params_Run = [K_r, w_rp1, w_rp3];
    elseif index==3 % r2-far : [101]
        params_Run = [K_r, w_rp1, w_rp2];
    end

    % We will use the Bcd and Run to generate the expected data for 20-50% of embryo length
    rate_expected = model_6A2R_HillModel_V3_direct_fixed_Kr_Bcd_RNAP_w_rp_params([w_rr, 1], TF, [params_Bcd,params_Run]);

    % concatenate to one vector
    model_expected = [model_expected; rate_expected];
end

% compute the AIC
AIC_w_rr = compute_AIC(model_expected, rate_vec, rate_std_vec, 1);
%% 3. 6A2R with w_ho as the free parameter
% load the 6A2R MCMC result (the loaded varaiable name is MCMC_6A2R_RuntWT)
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\higher_order_only\MCMC_6A2R_RuntWT_params.mat');

% initialize the model output vector
model_expected = [];
% loop over all 6A2R constructs
% we will concatenate the vectors (model_expected) on the order of [011],[110], and [101].
for index=1:3
    construct = constructIndex(index);
    % import the relevant parameters
    % : 6AnR_RuntNulls - params_Bcd, 6A1R - params_Run
    % 1) parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    w_ho = MCMC_6A2R_RuntWT(index).params_inferred;
    % Runt parameters
    if index==1 % r2-new : [011]
        params_Run = [K_r, w_rp2, w_rp3];
    elseif index==2 % r2-close : [110]
        params_Run = [K_r, w_rp1, w_rp3];
    elseif index==3 % r2-far : [101]
        params_Run = [K_r, w_rp1, w_rp2];
    end

    % We will use the Bcd and Run to generate the expected data for 20-50% of embryo length
    rate_expected = model_6A2R_HillModel_V3_direct_fixed_Kr_Bcd_RNAP_w_rp_params([1, w_ho], TF, [params_Bcd,params_Run]);

    % concatenate to one vector
    model_expected = [model_expected; rate_expected];
end

% compute the AIC
AIC_w_ho = compute_AIC(model_expected, rate_vec, rate_std_vec, 1);
%% 4. 6A2R with w_rr and w_ho as the free parameters
% load the 6A2R MCMC result (the loaded varaiable name is MCMC_6A2R_RuntWT)
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\both_Run_Run_higher_order\MCMC_6A2R_RuntWT_params.mat');

% initialize the model output vector
model_expected = [];
% loop over all 6A2R constructs
% we will concatenate the vectors (model_expected) on the order of [011],[110], and [101].
for index=1:3
    construct = constructIndex(index);
    % import the relevant parameters
    % : 6AnR_RuntNulls - params_Bcd, 6A1R - params_Run
    % 1) parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    % extract the MCMC run from 6A2R (w_rr and w_ho)
    w_rr = MCMC_6A2R_RuntWT(index).params_inferred(1);
    w_ho = MCMC_6A2R_RuntWT(index).params_inferred(2);
    % Runt parameters
    if index==1 % r2-new : [011]
        params_Run = [K_r, w_rp2, w_rp3];
    elseif index==2 % r2-close : [110]
        params_Run = [K_r, w_rp1, w_rp3];
    elseif index==3 % r2-far : [101]
        params_Run = [K_r, w_rp1, w_rp2];
    end

    % We will use the Bcd and Run to generate the expected data for 20-50% of embryo length
    rate_expected = model_6A2R_HillModel_V3_direct_fixed_Kr_Bcd_RNAP_w_rp_params([w_rr, w_ho], TF, [params_Bcd,params_Run]);

    % concatenate to one vector
    model_expected = [model_expected; rate_expected];
end

% compute the AIC
AIC_w_rr_w_ho = compute_AIC(model_expected, rate_vec, rate_std_vec, 2);


%% generate plots for AIC
plot([1,2,3,4], [AIC_params_free, AIC_w_rr, AIC_w_ho, AIC_w_rr_w_ho],'o','MarkerFaceColor',ColorChoice(1,:))
xlim([0 5])
ylim([0 600])
yticks([0 200 400 600])
xticks([1 2 3 4])
ylabel('AIC')
xticklabels({"zero free parameter","\omega_{rr}","\omega_{ho}","\omega_{rr} & \omega_{ho}"})
box on
StandardFigure(gcf,gca)
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop';
saveas(gcf, [FigPath, filesep, 'AIC_all_scenarios.pdf'])
end