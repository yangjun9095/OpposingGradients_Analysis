function MCMC_assess_different_priors_HillV3_RuntNulls
%% Description 
% This document is for assessing the effect of different priors.
% From [Ref.]

% Here, we will use the MCMC results saved in different directories,
% basically load the posterior "chain", and use the posterior distribution
% to generate histograms for 1) prior, 2) posterior, 3) corner plot



%% Load the MCMC result
% FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\w_bp_0-100_WeakPrior_w_bp';
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\w_bp_0-1000';
A = load([FilePath, filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat'])
MCMC_6A0R_RuntNulls = A.MCMC_6A0R_RuntNulls;

%% Let's generate the prior/posterior histograms for the w_bp, for [100] construct
% There's no particular reason that it has to be [100] construct, but it's
% more of picking one construct to generate a representative example.

construct = 2;
% extract useful fields
chain = MCMC_6A0R_RuntNulls(construct).chain;

% results structure has info on prior, parameter limits, etc.
results =  MCMC_6A0R_RuntNulls(construct).results;
prior = results.prior;
priorFunction = results.priorfun;
theta = results.theta;
par_lim = results.limits;
par_lim = par_lim(2,:);

% generate the prior distribution
xvec = linspace(par_lim(1), par_lim(2),100);
yvec = normpdf(xvec,prior(2,1),prior(2,2));
yvec = unifpdf(xvec, par_lim(1), par_lim(2));


%% generate histogram of prior and posterior distribution of w_bp
hold on
plot(xvec,yvec,'LineWidth',2)
histogram(chain(:,2),50,'Normalization','probability')

% ticks
xticks([0 20 40 60 80 100])
yticks([0.01 0.02 0.03])
% label
xlabel('\omega_{bp}')
ylabel('prob. density')
legend('prior','posterior')

box on
StandardFigure(gcf,gca)

% save the plots
FigPath = FilePath;
saveas(gcf,[FigPath,filesep, 'histogram_prior_posterior_w_bp','.tif']);
saveas(gcf,[FigPath,filesep, 'histogram_prior_posterior_w_bp','.pdf']);

end