% generate plots comparing parameters inferred using different schemes
function gen_plots_inf_params_diff_schemes_MCMC
%% Description
% This script is to generate plots for the comparison of inferred
% parameters from two different schemes of MCMC.
% (1) sequential MCMC : doing MCMC for Run WT using the point-estimates from
% the Run nulls.
% (2) global MCMC : doing MCMC for Run WT and Run nulls simultaneously.

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
%% Load the MCMC result from the (1) sequential MCMC inference
tempPath1 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\RuntNulls\w_bp_0-100'
MCMC_seq = load([tempPath1, filesep, 'MCMC_6A0R_RuntNulls_BcdParams.mat']);
%% Load the MCMC result from the (2) global MCMC inference
tempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\global_MCMC_inference'
MCMC_global = load([tempPath2, filesep, 'MCMCresult_HillV3_direct_params_1RunSite.mat']);

%% pull the inferred parameters and their std into separate matrices
% Note that we will only treat the 1 Run site constructs.
m=1;
MCMC_seq_params = [];
MCMC_seq_params_std = [];

for i=[2,5,6]
    MCMC_seq_params(m,:) = MCMC_seq.MCMC_6A0R_RuntNulls(i).params_inferred;
    MCMC_seq_params_std(m,:) = MCMC_seq.MCMC_6A0R_RuntNulls(i).params_inferred_std;
    m=m+1;
end


MCMC_global_params = MCMC_global.params_MCMC(:,[1,3,5,6]);
MCMC_global_params_std =  MCMC_global.params_MCMC_std(:,[1,3,5,6]);
%% generate pairwise plots of inferred parameters : K_b, w_bp, p, R_max

% define the construct index
constructIndex = [2,5,6]; % for color

%% 1) K_{b}
Kb = [MCMC_seq_params(:,1), MCMC_global_params(:,1)];
Kb_std = [MCMC_seq_params_std(:,1), MCMC_global_params_std(:,1)];

figure(1)
hold on
for i=1:3
    construct = constructIndex(i);
    errorbar([1,2],Kb(i,:),Kb_std(i,:),'LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 3])
ylim([0 120])
yticks([0 20 40 60 80 100 120])
xticks([1 2])
xticklabels({'sequential','global'})
xlabel('K_{b}')
ylabel('inferred parameter')
legend('100','001','010')

box on
StandardFigure(gcf,gca)

% % save the plots
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\MCMC_seq_global_comparison'
saveas(gcf,[FigPath,filesep, 'K_b','.tif']);
saveas(gcf,[FigPath,filesep, 'K_b','.pdf']);

%% 2) w_bp
w_bp = [MCMC_seq_params(:,2), MCMC_global_params(:,2)];
w_bp_std = [MCMC_seq_params_std(:,2), MCMC_global_params_std(:,2)];

figure
hold on
for i=1:3
    construct = constructIndex(i);
    errorbar([1,2],w_bp(i,:),w_bp_std(i,:),'LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 3])
ylim([0 120])
yticks([0 20 40 60 80 100 120])
xticks([1 2])
xticklabels({'sequential','global'})
xlabel('\omega_{bp}')
ylabel('inferred parameter')
legend('100','001','010')

box on
StandardFigure(gcf,gca)

% % save the plots
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\MCMC_seq_global_comparison'
saveas(gcf,[FigPath,filesep, 'w_bp','.tif']);
saveas(gcf,[FigPath,filesep, 'w_bp','.pdf']);

%% 3) p
p = [MCMC_seq_params(:,3), MCMC_global_params(:,3)];
p_std = [MCMC_seq_params_std(:,3), MCMC_global_params_std(:,3)];

figure
hold on
for i=1:3
    construct = constructIndex(i);
    errorbar([1,2],p(i,:),p_std(i,:),'LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 3])
ylim([0 1])
yticks([0 0.2 0.4 0.6 0.8 1])
xticks([1 2])
xticklabels({'sequential','global'})
xlabel('\omega_{bp}')
ylabel('inferred parameter')
legend('100','001','010')

box on
StandardFigure(gcf,gca)

% % save the plots
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\MCMC_seq_global_comparison'
saveas(gcf,[FigPath,filesep, 'p','.tif']);
saveas(gcf,[FigPath,filesep, 'p','.pdf']);

%% 4) R_{max}
R_max = [MCMC_seq_params(:,4), MCMC_global_params(:,4)];
R_max_std = [MCMC_seq_params_std(:,4), MCMC_global_params_std(:,4)];

figure
hold on
for i=1:3
    construct = constructIndex(i);
    errorbar([1,2],R_max(i,:),R_max_std(i,:),'LineWidth',2,'Color',ColorChoice(construct,:))
end

xlim([0 3])
ylim([0 400])
yticks([0 100 200 300 400])
xticks([1 2])
xticklabels({'sequential','global'})
xlabel('\omega_{bp}')
ylabel('inferred parameter')
legend('100','001','010')

box on
StandardFigure(gcf,gca)

% % save the plots
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\MCMC_seq_global_comparison'
saveas(gcf,[FigPath,filesep, 'R_max','.tif']);
saveas(gcf,[FigPath,filesep, 'R_max','.pdf']);
end