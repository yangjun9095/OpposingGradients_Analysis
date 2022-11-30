function FigS_compute_AIC_6A1R_alternative_assumptions
%% Description
% This script is to compute the Akaike Information Criterion (AIC) from
% 6A1R constructs with different assumptions on K_r and w_rp

% models: Hill.v3 with with zero free parameter, with K_r, with w_rp, or with both.

% data: 6A1R initial rate of transcription (Run WT) datasets from [100],
% [001], and [010]
% r1-new : [100] : 2nd in the MCMC_6A0R_RuntNulls
% r1-close : [001] : 5th
% r1-far : [010] : 6th

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
constructIndex = [2,5,6];
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

% loop over all 6A1R constructs
for i=1:3
    index = constructIndex(i);
    rate_wt = compiledData{index+1,9};
    rate_wt_std = compiledData{index+1,10};
%     rate_wt = data(index).Rate_WT;
%     rate_wt_std = nanstd(data(index).Rate_WT_individual,0,2);

    % truncate for 20-50% of embryo length
    rate_wt = rate_wt(APbin_start:APbin_end);
    rate_wt_std = rate_wt_std(APbin_start:APbin_end);

    if i==1
        rate_wt = [rate_wt(1:8);nan;nan;nan;nan;nan];
        rate_wt_std = [rate_wt_std(1:8);nan;nan;nan;nan;nan]
    end

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
% % import K_b, w_bp, p, R as well as K_r, w_rp1, w_rp2
% tempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A1R_Bcd_params_fromRuntNulls_fixed_Kr';
% load([tempPath2, filesep, 'MCMC_6A1R_RuntWT_params.mat'])
% 
% % extract the Runt-dependent parameters
% K_r = MCMC_6A1R_RuntWT.params_inferred(1);
% w_rp1 = MCMC_6A1R_RuntWT.params_inferred(2); %[100]
% w_rp2 = MCMC_6A1R_RuntWT.params_inferred(3); %[001]
% w_rp3 = MCMC_6A1R_RuntWT.params_inferred(4); %[010]
%% 1. 6A1R with global K_r and global w_rp
% generate the model's expected values using a set of parameters
model_expected = [];

% load the 6A1R MCMC result (the loaded varaiable name is MCMC_6A1R_RuntWT)
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A1R_Bcd_parms_fromRuntNulls_fixed_Kr_fixed_w_rp\MCMC_6A1R_RuntWT_params.mat');

% loop over all 6A1R constructs
% we will concatenate the vectors (model_expected) on the order of [100],[001], and [010].
for index=1:3
    construct = constructIndex(index);
    % import the relevant parameters
    % : 6AnR_RuntNulls - params_Bcd, 6A1R - params_Run
    % 1) parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
   
    % Runt parameters (global K_r and global w_rp)
    params_Run = MCMC_6A1R_RuntWT.params_inferred;

    % rearrange the parameters for both Bicoid and Runt
    params_MCMC = [params_Bcd(1), params_Run(1), params_Bcd(2), params_Run(2), params_Bcd(3), params_Bcd(4)];

    % We will use the Bcd and Run to generate the expected data for 20-50% of embryo length
    rate_expected = model_6A1R_HillModel_V3_direct(params_MCMC, TF);

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


AIC_global_K_r_w_rp = compute_AIC(model_expected, rate_vec, rate_std_vec, 2);

hold on
plot(rate_vec)
plot(model_expected)
%% 2. 6A1R with global K_r and local w_rp

% load the 6A2R MCMC result (the loaded varaiable name is MCMC_6A2R_RuntWT)
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A1R_Bcd_params_fromRuntNulls_fixed_Kr\MCMC_6A1R_RuntWT_params.mat');

% initialize the model output vector
model_expected = [];

% loop over all 6A1R constructs
% we will concatenate the vectors (model_expected) on the order of [011],[110], and [101].
for index=1:3
    construct = constructIndex(index);
    % import the relevant parameters
    % : 6AnR_RuntNulls - params_Bcd, 6A1R - params_Run
    % 1) parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    K_r = MCMC_6A1R_RuntWT.params_inferred(1);
    % Runt parameters
    if index==1 % r1-new : [100]
        w_rp = MCMC_6A1R_RuntWT.params_inferred(2);
    elseif index==2 % r1-close : [001]
        w_rp = MCMC_6A1R_RuntWT.params_inferred(3);
    elseif index==3 % r1-far : [010]
        w_rp = MCMC_6A1R_RuntWT.params_inferred(4);
    end

    % rearrange the parameters for both Bicoid and Runt
    params_MCMC = [params_Bcd(1), K_r, params_Bcd(2), w_rp, params_Bcd(3), params_Bcd(4)];

    % We will use the Bcd and Run to generate the expected data for 20-50% of embryo length
    rate_expected = model_6A1R_HillModel_V3_direct(params_MCMC, TF);
    % concatenate to one vector
    model_expected = [model_expected; rate_expected];
end

% compute the AIC
AIC_global_K_r_local_w_rp = compute_AIC(model_expected, rate_vec, rate_std_vec, 4);

hold on
plot(rate_vec)
plot(model_expected)
%% 3. 6A1R with local K_r and global w_rp
% load the 6A1R MCMC result (the loaded varaiable name is MCMC_6A2R_RuntWT)
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A1R_Bcd_params_fromRuntNulls_fixed_w_rp\MCMC_6A1R_RuntWT_params.mat');

% initialize the model output vector
model_expected = [];

% loop over all 6A1R constructs
% we will concatenate the vectors (model_expected) on the order of [011],[110], and [101].
for index=1:3
    construct = constructIndex(index);
    % import the relevant parameters
    % : 6AnR_RuntNulls - params_Bcd, 6A1R - params_Run
    % 1) parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    w_rp = MCMC_6A1R_RuntWT.params_inferred(4);
    % Runt parameters
    if index==1 % r1-new : [100]
        K_r = MCMC_6A1R_RuntWT.params_inferred(1);
    elseif index==2 % r1-close : [001]
        K_r = MCMC_6A1R_RuntWT.params_inferred(2);
    elseif index==3 % r1-far : [010]
        K_r = MCMC_6A1R_RuntWT.params_inferred(3);
    end

    % rearrange the parameters for both Bicoid and Runt
    params_MCMC = [params_Bcd(1), K_r, params_Bcd(2), w_rp, params_Bcd(3), params_Bcd(4)];

    % We will use the Bcd and Run to generate the expected data for 20-50% of embryo length
    rate_expected = model_6A1R_HillModel_V3_direct(params_MCMC, TF);
    % concatenate to one vector
    model_expected = [model_expected; rate_expected];
end

% compute the AIC
AIC_local_K_r_global_w_rp = compute_AIC(model_expected, rate_vec, rate_std_vec, 4);

hold on
plot(rate_vec)
plot(model_expected)
%% 4. 6A1R with local K_r and local w_rp
% load the 6A1R MCMC result (the loaded varaiable name is MCMC_6A2R_RuntWT)
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A1R_Bcd_params_fromRuntNulls_local_Kr_local_w_rp\MCMC_6A1R_RuntWT_params.mat');

% initialize the model output vector
model_expected = [];

% loop over all 6A1R constructs
% we will concatenate the vectors (model_expected) on the order of [011],[110], and [101].
for index=1:3
    construct = constructIndex(index);
    % import the relevant parameters
    % : 6AnR_RuntNulls - params_Bcd, 6A1R - params_Run
    % 1) parameters for the Bicoid and RNAP (inferred from the Runt null
    % datasets for the same construct.)
    params_Bcd = MCMC_6A0R_RuntNulls(construct).params_inferred;
    
    % extract the MCMC run from 6A1R (K_r and w_rp)
    % Runt parameters
    if index==1 % r1-new : [100]
        K_r = MCMC_6A1R_RuntWT.params_inferred(1);
        w_rp = MCMC_6A1R_RuntWT.params_inferred(4);
    elseif index==2 % r1-close : [001]
        K_r = MCMC_6A1R_RuntWT.params_inferred(2);
        w_rp = MCMC_6A1R_RuntWT.params_inferred(5);
    elseif index==3 % r1-far : [010]
        K_r = MCMC_6A1R_RuntWT.params_inferred(3);
        w_rp = MCMC_6A1R_RuntWT.params_inferred(6);
    end

    % rearrange the parameters for both Bicoid and Runt
    params_MCMC = [params_Bcd(1), K_r, params_Bcd(2), w_rp, params_Bcd(3), params_Bcd(4)];

    % We will use the Bcd and Run to generate the expected data for 20-50% of embryo length
    rate_expected = model_6A1R_HillModel_V3_direct(params_MCMC, TF);

    % concatenate to one vector
    model_expected = [model_expected; rate_expected];
end

% compute the AIC
AIC_local_K_r_local_w_rp= compute_AIC(model_expected, rate_vec, rate_std_vec, 6);


%% generate plots for AIC
plot([1,2,3,4], [AIC_global_K_r_w_rp, AIC_global_K_r_local_w_rp, AIC_local_K_r_global_w_rp, AIC_local_K_r_local_w_rp],'o','MarkerFaceColor',ColorChoice(1,:))
xlim([0 5])
ylim([800 2100])
% yticks([300 350 2000])
xticks([1 2 3 4])
ylabel('AIC')
xticklabels({"global","\omega_{rp} local","K_{r} local","K_{r} & \omega_{rp} local"})
box on
StandardFigure(gcf,gca)
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference';
saveas(gcf, [FigPath, filesep, '6A1R_AIC_all_scenarios.pdf'])
end