% Compare the cooperativities, w_rr and w_ho for 1) w_ho only case, 
% 2) Both w_rr and w_ho case.

%% Load the settings for generating plots
% load ColorChoice, etc.

%% Load the datasets
% 1) load the MCMC results from the w_ho only case.
tempPath1 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\higher_order_only';
A = load([tempPath1, filesep, 'MCMC_6A2R_RuntWT_params.mat']);
MCMC_6A2R_w_ho_only = A.MCMC_6A2R_RuntWT;

tempPath2 = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\both_Run_Run_higher_order\10 limit';
B = load([tempPath2, filesep, 'MCMC_6A2R_RuntWT_params.mat']);
MCMC_6A2R_w_rr_w_ho = B.MCMC_6A2R_RuntWT;

%% plot the parameters

% parse the parameters
for i=1:3
    % we fixed the w_rr = 1, so std = 0
    params_inferred_w_ho(i,1) = 1;
    params_inferred_w_ho_std(i,1) = 0;
    % extract the w_ho values (w_ho only)
    params_inferred_w_ho(i,2) = MCMC_6A2R_w_ho_only(i).params_inferred;
    params_inferred_w_ho_std(i,2) = MCMC_6A2R_w_ho_only(i).params_inferred_sigma;
    
    % extract the w_rr and w_ho values
    params_inferred_w_rr_w_ho(i,:) = MCMC_6A2R_w_rr_w_ho(i).params_inferred;
    params_inferred_w_rr_w_ho_std(i,:) = MCMC_6A2R_w_rr_w_ho(i).params_inferred_sigma;
end




hold on
errorbar([1;2], params_inferred_w_ho(1,:), params_inferred_w_ho_std(1,:),'o','LineWidth',2,'Color',ColorChoice(3,:))
% errorbar([1;2], params_inferred_w_ho(2,:), params_inferred_w_ho_std(2,:),'o','LineWidth',2,'Color',ColorChoice(7,:))
errorbar([1;2], params_inferred_w_ho(3,:), params_inferred_w_ho_std(3,:),'o','LineWidth',2,'Color',ColorChoice(8,:))

% errorbar([1;2], params_inferred_w_rr_w_ho(1,:), params_inferred_w_rr_w_ho_std(1,:),'o','LineWidth',2,'Color',ColorChoice(3,:))
errorbar([1;2], params_inferred_w_rr_w_ho(2,:), params_inferred_w_rr_w_ho_std(2,:),'o','LineWidth',2,'Color',ColorChoice(7,:))
% errorbar([1;2], params_inferred_w_rr_w_ho(3,:), params_inferred_w_rr_w_ho_std(3,:),'o','LineWidth',2,'Color',ColorChoice(8,:))



xlim([0 3])
ylim([0 100])
xticks([1,2])
xticklabels({'\omega_{rr}', '\omega_{ho}'})
% yticks([0 100 200 300 400])
xlabel('parameters')
ylabel('inferred values')
legend('011','110','101','Location', 'NorthWest')


set(gca,'YScale','Log')
box on
StandardFigure(gcf,gca)

% save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop\hybrid_w_ho_only_w_rr_w_ho_only_[110]';
saveas(gcf,[FigPath,filesep, 'inferred_params_w_rr_w_ho_logScale','.tif']);
saveas(gcf,[FigPath,filesep, 'inferred_params_w_rr_w_ho_logScale','.pdf']);


%% create a new matrix containing w_rr and w_ho values: For [011] and [101], we set w_rr=1
MCMC_6A2R_params = [params_inferred_w_ho(1,:); params_inferred_w_rr_w_ho(2,:); params_inferred_w_ho(3,:) ];
MCMC_6A2R_params_std = [params_inferred_w_ho_std(1,:); params_inferred_w_rr_w_ho_std(2,:); params_inferred_w_ho_std(3,:) ];

% save these for future usage. Note that the ordering of the datsets are
% [011], [110], and [101].

%% compare the w_ho values from 1) w_ho only, and 2) both w_rr and w_ho 
hold on
errorbar([1;2], [params_inferred_w_ho(1,2),  params_inferred_w_rr_w_ho(1,2)],...
                    [params_inferred_w_ho_std(1,2), params_inferred_w_rr_w_ho_std(1,2)],...
                    'o','LineWidth',2,'Color',ColorChoice(3,:))
errorbar([1;2], [params_inferred_w_ho(2,2),  params_inferred_w_rr_w_ho(2,2)],...
                    [params_inferred_w_ho_std(2,2), params_inferred_w_rr_w_ho_std(2,2)],...
                    'o','LineWidth',2,'Color',ColorChoice(7,:))                
errorbar([1;2], [params_inferred_w_ho(3,2),  params_inferred_w_rr_w_ho(3,2)],...
                    [params_inferred_w_ho_std(3,2), params_inferred_w_rr_w_ho_std(3,2)],...
                    'o','LineWidth',2,'Color',ColorChoice(8,:))                

xlim([0 3])
% ylim([0 100])
xticks([1,2])
xticklabels({'\omega_{ho} (only)', '\omega_{ho} (with \omega_{rr})'})
% yticks([0 100 200 300 400])
xlabel('parameters')
ylabel('inferred values')
legend('011','110','101','Location', 'NorthWest')

set(gca,'YScale','Log')
box on
StandardFigure(gcf,gca)

% save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3\Direct\seq_MCMC_inference\6A2R_prediction_from_pt_estimates_fixed_Kr\fixed_K_r_w_rp_withCoop';
saveas(gcf,[FigPath,filesep, 'comparison_inferred_params_w_ho_logScale','.tif']);
saveas(gcf,[FigPath,filesep, 'comparison_inferred_params_w_ho_logScale','.pdf']);
