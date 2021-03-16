function generate_plot_MCMC_inf_params_AztecPyramid
%% Desription
% this script is to generate plots for comparing the MCMC-inferred
% parameters from different constructs, 6A0R, 6A1R, 6A2R, etc.
% As we use the Aztec pyramid scheme for the Bcd-dependent parameters, we
% will only compare the Run-dependent parameters.

%% Load the MCMC results
TempPath1 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\Point_estimate_Bcd_params_fromRunNulls';
load([TempPath1, filesep, 'MCMC_6A1R_RuntWT_params.mat'])
% MCMC_6A1R_RuntWT structure
TempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\6A2R_MCMC_seq_inference_pt_estimate_fromNulls';
load([TempPath2, filesep, 'MCMC_6A2R_RuntWT_params.mat']);

%% generate a plot for K_{r} versus each Run site (over all constructs of 6A1R, 6A2R)

% extract the inferred parameters and std from 6A1R
params_6A1R = extractfield(MCMC_6A1R_RuntWT, 'params_inferred')
params_6A1R_std = extractfield(MCMC_6A1R_RuntWT, 'params_inferred_sigma')
% note that the order is [100], [010], [001]
Kr_6A1R = params_6A1R([1,5,3]);
w_rp_6A1R = params_6A1R([2,6,4]);

Kr_6A1R_std = params_6A1R_std([1,5,3]);
w_rp_6A1R_std = params_6A1R_std([2,6,4]);

% extract the inferred parameters and std from 6A2R
params_6A2R = extractfield(MCMC_6A2R_RuntWT, 'params_inferred')
params_6A2R_std = extractfield(MCMC_6A2R_RuntWT, 'params_inferred_sigma')
% note that the order is [011], [110], [101]
Kr_6A2R = params_6A2R([1,3,5]);
w_rp_6A2R = params_6A2R([2,4,6]);

Kr_6A2R_std = params_6A2R_std([1,3,5]);
w_rp_6A2R_std = params_6A2R_std([2,4,6]);

%% K_{r} for different Run site position
hold on
errorbar([1:3], Kr_6A1R, Kr_6A1R_std,'o','Color', ColorChoice(1,:),'LineWidth', 2)
errorbar([2,3]-0.05, [Kr_6A2R(1), Kr_6A2R(1)],[Kr_6A2R_std(1), Kr_6A2R_std(1)],'o','Color', ColorChoice(3,:),'LineWidth', 2)
errorbar([1,2]+0.05, [Kr_6A2R(2), Kr_6A2R(2)],[Kr_6A2R_std(2), Kr_6A2R_std(2)],'o','Color', ColorChoice(7,:),'LineWidth', 2)
errorbar([1,3]-0.1, [Kr_6A2R(3), Kr_6A2R(3)],[Kr_6A2R_std(3), Kr_6A2R_std(3)],'o','Color', ColorChoice(8,:),'LineWidth', 2)
legend('1 Run sites','[011]','[110]','[101]')

xlim([0 4])
xticks([1,2,3])
xticklabels({'[100]','[010]','[001]'})
yticks([0 20 40 60 80 100])
xlabel('binding site position')
ylabel('inferred K_{r}')
ylim([0 100])

box on
StandardFigure(gcf,gca)

% save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct'
saveas(gcf,[FigPath,filesep,'inferred_K_r_6A1R_6A2R' ,'.pdf']); 
saveas(gcf,[FigPath,filesep,'inferred_K_r_6A1R_6A2R' ,'.tif']); 

%% w_rp for different Run site position
hold on
errorbar([1:3], w_rp_6A1R, w_rp_6A1R_std,'o','Color', ColorChoice(1,:),'LineWidth', 2)
errorbar([2,3]-0.05, [w_rp_6A2R(1), w_rp_6A2R(1)],[w_rp_6A2R_std(1), w_rp_6A2R_std(1)],'o','Color', ColorChoice(3,:),'LineWidth', 2)
errorbar([1,2]+0.05, [w_rp_6A2R(2), w_rp_6A2R(2)],[w_rp_6A2R_std(2), w_rp_6A2R_std(2)],'o','Color', ColorChoice(7,:),'LineWidth', 2)
errorbar([1,3]-0.1, [w_rp_6A2R(3), w_rp_6A2R(3)],[w_rp_6A2R_std(3), w_rp_6A2R_std(3)],'o','Color', ColorChoice(8,:),'LineWidth', 2)
legend('1 Run sites','[011]','[110]','[101]')

xlim([0 4])
xticks([1,2,3])
xticklabels({'[100]','[010]','[001]'})
yticks([0 0.2 0.4 0.6 0.8 1 1.2])
xlabel('binding site position')
ylabel('inferred \omega_{rp}')
ylim([0 1.2])

box on
StandardFigure(gcf,gca)

% save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct'
saveas(gcf,[FigPath,filesep,'inferred_w_rp_6A1R_6A2R' ,'.pdf']); 
saveas(gcf,[FigPath,filesep,'inferred_w_rp_6A1R_6A2R' ,'.tif']); 
end