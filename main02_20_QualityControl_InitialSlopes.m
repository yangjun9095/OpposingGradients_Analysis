function main02_20_QualityControl_InitialSlopes
%% DESCRIPTION
% This script is for a quality check of fitted initial slopes for all the
% constructs

%% Quality control is finicky, since we really don't know what's the standard for good/bad.
% For now, I can try with a couple different thresholds, and see how much
% do different thresholds affect the final output.
% Ideally, different thresholds of good/bad fits wouldn't change the end
% result that much.

%% First, let's use C.V to see how good the fitting was for each construct.

% For each MeanFitsAsymmetric.m, there's RateFit (fitted initial slope),
% and SDRateFit (error estimated using confidence interval of 68%)
% We will use the C.V (SDRateFit/RateFit) at each AP bin (only for NC14),
% to see whether we can distinguish good/bad fits.

%% Preceding steps
% either running main02_08_plot_InitialSlopes_AllConstructs.m or load
% the datasets again.

%% Load the datasets
% We can use the custom function Extract_Fields_MeanFits for this!
% Data_r0 = LoadMS2Sets('r0-new','dontCompare')
% [fittedRate_r0,fittedRateSD_r0,fittedTon_r0] = Extract_Fields_MeanFits(Data_r0,'Asymmetric');

%% Define an empty matrix filled with NaNs.
% Extract RateFit, SDRateFit in NC14 for all APbins, then replace NaNs.

% So, the matrix dimension should be (rows: embryo index) x (columns: APbins)
% Make two of these. Replace these with RateFits, and also SDRateFits
% Divide to get the Coefficient of Variation for each AP bin, for each
% embryo. Then, let's plot it across the AP axis.

CV_r0 = squeeze(fittedRateSD_r0(:,3,:) ./ fittedRate_r0(:,3,:));

CV_r1 = squeeze(fittedRateSD_r1(:,3,:) ./ fittedRate_r1(:,3,:));
CV_r1_close = squeeze(fittedRateSD_r1_close(:,3,:) ./ fittedRate_r1_close(:,3,:));
CV_r1_mid = squeeze(fittedRateSD_r1_mid(:,3,:) ./ fittedRate_r1_mid(:,3,:));

CV_r2 = squeeze(fittedRateSD_r2(:,3,:) ./ fittedRate_r2(:,3,:));
CV_r2_close = squeeze(fittedRateSD_r2_close(:,3,:) ./ fittedRate_r2_close(:,3,:));
CV_r2_far = squeeze(fittedRateSD_r2_far(:,3,:) ./ fittedRate_r2_far(:,3,:));

CV_r3 = squeeze(fittedRateSD_r3(:,3,:) ./ fittedRate_r3(:,3,:));

%plot(0:0.025:1, CV_r1,'o-')

%% plot to check the C.V across the AP axis.
APaxis = 0:0.025:1;

hold on

plot(APaxis, CV_r0,'o-','Color',ColorChoice(1,:))
plot(APaxis, CV_r1,'o-','Color',ColorChoice(2,:))
plot(APaxis, CV_r2,'o-','Color',ColorChoice(3,:))
plot(APaxis, CV_r3,'o-','Color',ColorChoice(4,:))

plot(APaxis, CV_r1_close,'o-','Color',ColorChoice(5,:))
plot(APaxis, CV_r1_mid,'o-','Color',ColorChoice(6,:))

plot(APaxis, CV_r2_close,'o-','Color',ColorChoice(7,:))
plot(APaxis, CV_r2_far,'o-','Color',ColorChoice(8,:))

%ylim([0 1])

title('C.V across AP for all constructs')
xlabel('AP axis (EL)')
ylabel('C.V')

% legend('r0','r1','r2','r3','r1-close','r1-mid','r2-close','r2-far')

StandardFigure(gcf,gca)

% Save the figure
saveas(gcf,[FigPath 'CV_over_APaxis_AllConstructs_NC14' , '.tif']); 
saveas(gcf,[FigPath 'CV_over_APaxis_AllConstructs_NC14' , '.pdf']); 

%% Use the C.V threshold to cut off the bad fits
CV_thresh = 0.2; % this can be changed

% Define the C.V
CVfilter_r0 = fittedRateSD_r0./fittedRate_r0;
% Define a matrix filled with true or false (not NaN, and also smaller than
% the threshold)
CVfilter_r0 = ~isnan(CVfilter_r0).*(CVfilter_r0 < CV_thresh);
% Convert the zeros to NaNs
CVfilter_r0(CVfilter_r0==0) = nan;
% filter out the fitted initial rate using the CV filter
fittedRate_CVfilterd_r0 = fittedRate_r0.*CVfilter_r0;

% Repeat this for all the other constructs
% r1
CVfilter_r1 = fittedRateSD_r1./fittedRate_r1;
CVfilter_r1 = ~isnan(CVfilter_r1).*(CVfilter_r1 < CV_thresh);
CVfilter_r1(CVfilter_r1==0) = nan;
fittedRate_CVfilterd_r1 = fittedRate_r1.*CVfilter_r1;
% r1_close
CVfilter_r1_close = fittedRateSD_r1_close./fittedRate_r1_close;
CVfilter_r1_close = ~isnan(CVfilter_r1_close).*(CVfilter_r1_close < CV_thresh);
CVfilter_r1_close(CVfilter_r1_close==0) = nan;
fittedRate_CVfilterd_r1_close = fittedRate_r1_close.*CVfilter_r1_close;
% r1_mid
CVfilter_r1_mid = fittedRateSD_r1_mid./fittedRate_r1_mid;
CVfilter_r1_mid = ~isnan(CVfilter_r1_mid).*(CVfilter_r1_mid < CV_thresh);
CVfilter_r1_mid(CVfilter_r1_mid==0) = nan;
fittedRate_CVfilterd_r1_mid = fittedRate_r1_mid.*CVfilter_r1_mid;

% r2
CVfilter_r2 = fittedRateSD_r2./fittedRate_r2;
CVfilter_r2 = ~isnan(CVfilter_r2).*(CVfilter_r2 < CV_thresh);
CVfilter_r2(CVfilter_r2==0) = nan;
fittedRate_CVfilterd_r2 = fittedRate_r2.*CVfilter_r2;
% r2_close
CVfilter_r2_close = fittedRateSD_r2_close./fittedRate_r2_close;
CVfilter_r2_close = ~isnan(CVfilter_r2_close).*(CVfilter_r2_close < CV_thresh);
CVfilter_r2_close(CVfilter_r2_close==0) = nan;
fittedRate_CVfilterd_r2_close = fittedRate_r2_close.*CVfilter_r2_close;
% r2_far
CVfilter_r2_far = fittedRateSD_r2_far./fittedRate_r2_far;
CVfilter_r2_far = ~isnan(CVfilter_r2_far).*(CVfilter_r2_far < CV_thresh);
CVfilter_r2_far(CVfilter_r2_far==0) = nan;
fittedRate_CVfilterd_r2_far = fittedRate_r2_far.*CVfilter_r2_far;

% r3
CVfilter_r3 = fittedRateSD_r3./fittedRate_r3;
CVfilter_r3 = ~isnan(CVfilter_r3).*(CVfilter_r3 < CV_thresh);
CVfilter_r3(CVfilter_r3==0) = nan;
fittedRate_CVfilterd_r3 = fittedRate_r3.*CVfilter_r3;

%% Average over embryos (NC14 only) using nanmean, nanstd
% r0
average_fittedRate_r0_CVfiltered = nanmean(fittedRate_CVfilterd_r0,3);
SEM_fittedRate_r0_CVfiltered = nanstd(fittedRate_CVfilterd_r0,0,3)/sqrt(length(Data_r0));

% r1
average_fittedRate_r1_CVfiltered = nanmean(fittedRate_CVfilterd_r1,3);
SEM_fittedRate_r1_CVfiltered = nanstd(fittedRate_CVfilterd_r1,0,3)/sqrt(length(Data_r1));
% r1_close
average_fittedRate_r1_close_CVfiltered = nanmean(fittedRate_CVfilterd_r1_close,3);
SEM_fittedRate_r1_close_CVfiltered = nanstd(fittedRate_CVfilterd_r1_close,0,3)/sqrt(length(Data_r1_close));
%r1_mid
average_fittedRate_r1_mid_CVfiltered = nanmean(fittedRate_CVfilterd_r1_mid,3);
SEM_fittedRate_r1_mid_CVfiltered = nanstd(fittedRate_CVfilterd_r1_mid,0,3)/sqrt(length(Data_r1_mid));

%r2
average_fittedRate_r2_CVfiltered = nanmean(fittedRate_CVfilterd_r2,3);
SEM_fittedRate_r2_CVfiltered = nanstd(fittedRate_CVfilterd_r2,0,3)/sqrt(length(Data_r2));
%r2_close
average_fittedRate_r2_close_CVfiltered = nanmean(fittedRate_CVfilterd_r2_close,3);
SEM_fittedRate_r2_close_CVfiltered = nanstd(fittedRate_CVfilterd_r2_close,0,3)/sqrt(length(Data_r2_close));
%r2_far
average_fittedRate_r2_far_CVfiltered = nanmean(fittedRate_CVfilterd_r2_far,3);
SEM_fittedRate_r2_far_CVfiltered = nanstd(fittedRate_CVfilterd_r2_far,0,3)/sqrt(length(Data_r2_far));
%r3
average_fittedRate_r3_CVfiltered = nanmean(fittedRate_CVfilterd_r3,3);
SEM_fittedRate_r3_CVfiltered = nanstd(fittedRate_CVfilterd_r3,0,3)/sqrt(length(Data_r3));

%% Plot the averaged initial slope (SEM) post CV-filtering step.
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\';
%% 1) r1 variants
APaxis = 0:0.025:1;
nc = 3; % NC14

hold on
errorbar(0:0.025:1,average_fittedRate_r0_CVfiltered(:,nc),SEM_fittedRate_r0_CVfiltered(:,nc),'Color',ColorChoice(1,:))

errorbar(0:0.025:1,average_fittedRate_r1_CVfiltered(:,nc),SEM_fittedRate_r1_CVfiltered(:,nc),'Color',ColorChoice(2,:))
errorbar(0:0.025:1,average_fittedRate_r1_close_CVfiltered(:,nc),SEM_fittedRate_r1_close_CVfiltered(:,nc),'Color',ColorChoice(5,:))
errorbar(0:0.025:1,average_fittedRate_r1_mid_CVfiltered(:,nc),SEM_fittedRate_r1_mid_CVfiltered(:,nc),'Color',ColorChoice(6,:))

errorbar(0:0.025:1,average_fittedRate_r3_CVfiltered(:,nc),SEM_fittedRate_r3_CVfiltered(:,nc),'Color',ColorChoice(4,:)) % red

xlim([0.15 0.6])
ylim([0 250])

legend('r0','r1','r1-close','r1-mid', 'r3')%,'r1-male','r2-male','r3-male') 
xlabel('AP Position')
ylabel('Initial rate (AU/min)')
title('Initial rate of RNAP loading along AP axis, at NC 14')
StandardFigure(gcf, gca)

saveas(gcf,[FigPath 'InitialRate_AsymmetricFit_r1_variants' , '_NC14_CV_0p2_filtered' , '.tif']); 
saveas(gcf,[FigPath 'InitialRate_AsymmetricFit_r1_variants' , '_NC14_CV_0p2_filtered' , '.pdf']); 

%% r2 (variants)
APaxis = 0:0.025:1;
nc = 3; % NC14

hold on
errorbar(0:0.025:1,average_fittedRate_r0_CVfiltered(:,nc),SEM_fittedRate_r0_CVfiltered(:,nc),'Color',ColorChoice(1,:))

errorbar(0:0.025:1,average_fittedRate_r2_CVfiltered(:,nc),SEM_fittedRate_r2_CVfiltered(:,nc),'Color',ColorChoice(3,:))
errorbar(0:0.025:1,average_fittedRate_r2_close_CVfiltered(:,nc),SEM_fittedRate_r2_close_CVfiltered(:,nc),'Color',ColorChoice(7,:))
errorbar(0:0.025:1,average_fittedRate_r2_far_CVfiltered(:,nc),SEM_fittedRate_r2_far_CVfiltered(:,nc),'Color',ColorChoice(8,:))

errorbar(0:0.025:1,average_fittedRate_r3_CVfiltered(:,nc),SEM_fittedRate_r3_CVfiltered(:,nc),'Color',ColorChoice(4,:)) % red

xlim([0.15 0.6])
ylim([0 250])

legend('r0','r2','r2-close','r2-far', 'r3')%,'r1-male','r2-male','r3-male') 
xlabel('AP Position')
ylabel('Initial rate (AU/min)')
title('Initial rate of RNAP loading along AP axis, at NC 14')
StandardFigure(gcf, gca)

saveas(gcf,[FigPath 'InitialRate_AsymmetricFit_r2_variants' , '_NC14_CV_0p2_filtered' , '.tif']); 
saveas(gcf,[FigPath 'InitialRate_AsymmetricFit_r2_variants' , '_NC14_CV_0p2_filtered' , '.pdf']); 


%% Save the fields.
end