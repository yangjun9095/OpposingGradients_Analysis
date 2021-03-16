function calculate_CV_6A1R_inferredParameters
%% Description

%% Load the MCMC results from different "modes"
% FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct';
% load([FilePath, filesep, 'MCMCresult_HillV3_direct_params_1RunSite.mat'])
% It loads a structure called MCMCResult
load('S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Competition\MCMCResult_HillV3_competition.mat')

% Calculate the mean and std of the inferred parameters from the MCMC chain
n_steps = length(MCMCResult(1).chain);
n_burn = 0.5*n_steps;
for i=1:length(MCMCResult)
    params_6A1R_MCMC_mean(i,:) = mean(MCMCResult(i).chain(n_burn+1:end,:));
    params_6A1R_MCMC_std(i,:) = std(MCMCResult(i).chain(n_burn+1:end,:));
end

%% generate the plot of inferred parameters (log scale)
% hold on
% for construct = 1:length(data)
%     params_inferred = MCMC_6A0R_RuntNulls(construct).params_inferred;
%     params_inferred_std = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
%     errorbar(1:4, params_inferred, params_inferred_std,'o','LineWidth',2,'Color',ColorChoice(construct,:))
% end
% 
% xlim([0 7])
% xticks([1,2,3,4])
% xticklabels({'K_{b}','\omega_{bp}','p','R_{max}'})
% % yticks([0 100 200 300 400])
% set(gca,'YScale','log')
% % xlabel('')
% ylabel('inferred parameters')
% legend('000','100','011','111','001','010','110','101', 'Location', 'NorthEast')
% 
% box on
% StandardFigure(gcf,gca)
% 
% % save the plots
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A0R_RuntNull_AllConstructs_LogScale','.tif']);
% saveas(gcf,[FigPath,filesep, 'MCMCfit_6A0R_RuntNull_AllConstructs_LogScale','.pdf']);

%% generate the plot for C.V. of inferred parameters

% calculate the Coefficient of Variation (C.V.) of inferred parameters
% across constructs (different enhancers).
% for construct = 1:length(data)
%     params_inferred_all(construct,:) = MCMC_6A0R_RuntNulls(construct).params_inferred;
%     params_inferred_error_all(construct,:) = MCMC_6A0R_RuntNulls(construct).params_inferred_std;
% end

% calculate the mean and std of parameters "over constructs".
params_mean = mean(params_6A1R_MCMC_mean);
params_std = std(params_6A1R_MCMC_mean);
% params_error = sqrt(sum(params_inferred_error_all.^2)/length(data))

n_boots = 100;
parms_std_boostrap = bootstrp(n_boots, @std, params_6A1R_MCMC_mean);
params_SEM = std(parms_std_boostrap)./sqrt(n_boots-1);

params_CV = params_std./params_mean;
params_CV_error = params_SEM./params_mean;

errorbar(1:6, params_CV, params_CV_error,'o','LineWidth',2)

xlim([0 7])
xticks([1,2,3,4,5,6])
xticklabels({'K_{b}','K_{r}','\omega_{bp}','\omega_{rp}','p','R_{max}'})
% yticks([0 100 200 300 400])
ylim([0 1])
% set(gca,'YScale','log')
% xlabel('')
ylabel('Coefficient of Variation')

box on
StandardFigure(gcf,gca)

% save the plots
FigPath = FilePath;
saveas(gcf,[FigPath,filesep, 'CV_parameters','.tif']);
saveas(gcf,[FigPath,filesep, 'CV_parameters','.pdf']);
end