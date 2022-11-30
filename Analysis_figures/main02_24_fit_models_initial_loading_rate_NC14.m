function main02_02_fit_models_initial_loading_rate_NC14
%% Fitting the initial loading rates(Asymmetric) from NC14  with theoretical models
% Here, I will use Asymmetric-trapezoidal fitted initial slopes for fitting with some
% expression for different constructs with different numbers of binding
% sites. 
%Note that the initial rate of RNAP loading fit from individual embryos 
% are averaged over more than 3.

% Load the Bcd, Runt data
% Bicoid : Time-averaging is needed, but we know the Bcd gradient shape
% (length constant) is conserved, and only the amplitude changes over time
% (from Paul's fitting with Exponential), thus we'll just use the Bcd
% profile at one time point. The K_a will scale things properly anyway.
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

% First, we should build a prediction matrix using handful of parameters.
% Use r0_Hill_initial_rate_fit.m for optimizing the parameters, 
% r_max, r_basal, and K_d using lsqnonlin.

% Second, we need to decide which nc we're going to model.
% Here, we'll do NC14.

% Now, I will start from the initial loading rate in NC14.
% Bcd and Runt should be time-averaged for 0-10 minutes in NC14 (the initial
% rise phase) : YJK on 8/8/2019 : This assumption could change...

% 1) Bcd :Take the Bicoid data from Liz&Jonathan
% time resolution : 30 sec.
BcdAnt = load([FilePath, filesep, 'BcdGFPAnt.mat']);

BcdNC14 = BcdAnt.DataBcd.nc14;
% Take temporal average over 0-10 minutes into NC13.
BcdFluo = nanmean(BcdAnt.DataBcd.MeanVectorAP(BcdNC14: BcdNC14 + 20,:));
BcdFluoSD = nanmean(BcdAnt.DataBcd.SDVectorAP(BcdNC14: BcdNC14 + 20,:));

Runt = load([FilePath, filesep, 'Runt_time_averaged_female.mat']);
% Runt is already time-averaged. For 0-10 minutes, it's in the 5th cell.
RuntFluo = Runt.Averaged_Fluo(5,:);
RuntFluoSE = Runt.SE_Fluo(5,:);

RuntFluo = movmean(RuntFluo,3); % smoothening with 3 AP bins

% Extrapoloate the Runt data in anterior bins (since I need from 20%)
RuntFluo_extrap = interp1([0.3:0.025:0.625], RuntFluo(13:26), [0.2:0.025:0.275], 'pchip', 'extrap')

RuntFluo_extrapolated = nan(41,1);
RuntFluo_extrapolated(9:12) = RuntFluo_extrap; % 20%-27.5%
RuntFluo_extrapolated(13:27) = RuntFluo(13:27); % 30-65%

% plot to check
APaxis = 0:0.025:1;
RuntScale = 5; % Scaling for plot

hold on
yyaxis left
errorbar(APaxis, BcdFluo, BcdFluoSD)
ylabel('Bicoid concentration (AU)')

yyaxis right
errorbar(APaxis, RuntFluo * RuntScale, RuntFluoSE)
plot(APaxis, RuntFluo_extrapolated * RuntScale)
ylabel('Runt concentration (AU)')

xlim([0.2 0.6])

title({'Bcd and Runt profile over AP';'(time-averaged, 0-10 min in NC14)'})
xlabel('AP (Embryo Length)')
legend('Bcd','Runt')

StandardFigure(gcf,gca)

% Save figures
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope_Asymmetric';
% saveas(gcf,[FigPath,filesep,'InputTF_Time_Averaged_NC14_0-10min_BcdRunt' , '_NC14' , '.tif']); 
% saveas(gcf,[FigPath,filesep,'InputTF_Time_Averaged_NC14_0-10min_BcdRunt' , '_NC14' , '.pdf']); 
%% Save the extrapolated Runt profile
Runt = load([FilePath,filesep,'Runt_time_averaged_female.mat'])
Runt.Extrapolated_Fluo_0_10min_NC14 = RuntFluo_extrapolated;

save([FilePath,filesep,'Runt_time_averaged_female.mat'],...
   '-struct', 'Runt','-v7.3');
%% Load the fitted initial slope (Data)
% Data calculated from the script :
% main02_plot_Initial_Loading_Rates_Asymmetric.m
% 
InitialSlope = load([FilePath, filesep, 'AveragedInitialRate_AsymmetricFit.mat']);

%% Extract useful fields 
% Extract the initial slope from the saved structure.
nc = 14;% NC14
NC= nc - 11;
InitialRate_r0_NC14 = InitialSlope.average_fittedRate_r0(:,NC);
InitialRate_r1_NC14 = InitialSlope.average_fittedRate_r1_female(:,NC);
InitialRate_r2_NC14 = InitialSlope.average_fittedRate_r2_female(:,NC);
InitialRate_r3_NC14 = InitialSlope.average_fittedRate_r3_female(:,NC);

InitialRate_SEM_r0_NC14 = InitialSlope.SEM_fittedRate_r0(:,NC);
InitialRate_SEM_r1_NC14 = InitialSlope.SEM_fittedRate_r1_female(:,NC);
InitialRate_SEM_r2_NC14 = InitialSlope.SEM_fittedRate_r2_female(:,NC);
InitialRate_SEM_r3_NC14 = InitialSlope.SEM_fittedRate_r3_female(:,NC);

%% r0 : Initial Rate extrapolation for the posterior bins
% APbin1 = 0.15./0.025 + 1;
% APbin2 = 0.5 / 0.025 + 1;
% %InitialRate_r0_NC13_extrap = interp1([0.15:0.025:0.5], InitialRate_r0_NC13(APbin1:APbin2), [0.525:0.025:0.625], 'pchip', 'extrap')
% 
% %InitialRate_r0_NC14(22:26) = min(InitialRate_r0_NC14);
% 
% plot(APaxis, InitialRate_r0_NC14)
%% limit the AP bins for fitting 
% (We're trying to fit 20% of the AP axis to the most posterior point 
%  where there's data point)

% 1) find the non-NaNs from the AP axis
NoNaNIndices = ~isnan(InitialRate_r0_NC14);
NoNaNIndex = find(NoNaNIndices);
APbinEnd = max(NoNaNIndex);

APbinstart = 9;
APbinend = APbinEnd; %25;

% 2) Initial condition
R_max = 300; % maximum
R_bas = 50; % basal, probably due to ectopic expression...
Kd = 10;
N = 6;
% fitting for r0 first,
p0=[R_bas,R_max,Kd,N];
%lb=[0,R_max,R_bas];
lb = [0, 0, 0, 1];
ub=[R_bas*10, R_max*10, Kd*10^2, N*10];
options = optimset('Display','iter');

fun= @(p)r0_Hill_Thermo_initial_rate_fit(p,BcdFluo(APbinstart:APbinend)) - InitialRate_r0_NC14(APbinstart:APbinend)';

P = lsqnonlin(fun, p0, lb, ub, options);

% Plot the fit to check if it's reasonable
r_basal = P(1);
r = P(2);
K_a = P(3);
N = P(4);
%N=6; % Fix the Hill coeff. to be 6, corresponding to the # of binding sites.
BcdData = BcdFluo;

% Get the fitted profile
Fit = nan(41,1);
%Fit = (r_basal*ones(size(BcdData)) + (BcdData./K_a).^N  * r) ./ (1 + (BcdData./K_a).^N)
Fit = r0_Hill_Thermo_initial_rate_fit(P,BcdData);

hold on
%plot(APaxis, (r_basal*ones(size(BcdData))))
%plot(APaxis, (BcdData./K_a).^6  * r)

fitRange = APbinstart:APbinend;
% fitRange =1:41;
plot(APaxis(fitRange), Fit(fitRange))
plot(APaxis, InitialRate_r0_NC14,'-o')
xlim([0.15 0.6])
title('Fit r0 model to initial RNAP loading rate (NC14)')
xlabel('AP axis (embryo length)')
ylabel('Initial loaind rate (AU)')
legend('Fit','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope_Asymmetric';
saveas(gcf,[FigPath,filesep,'r0_InitialSlope_Hill_Fits' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath,filesep,'r0_InitialSlope_Hill_Fits' , '_NC14' , '.pdf']); 

%% r1 fitting
% %% First, extrapolate the r1 at posterior positions
% APbin1 = 0.15./0.025 + 1;
% APbin2 = 0.5 / 0.025 + 1;
% %InitialRate_r0_NC13_extrap = interp1([0.15:0.025:0.5], InitialRate_r0_NC13(APbin1:APbin2), [0.525:0.025:0.625], 'pchip', 'extrap')
% 
% InitialRate_r1_NC13(21:25) = min(InitialRate_r1_NC13);
% 
% plot(APaxis, InitialRate_r1_NC13)

% %% Fitting (r1) - Thermodynamic model
% 
% % 1) find the non-NaNs from the AP axis
% NoNaNIndices = ~isnan(InitialRate_r1_NC14);
% NoNaNIndex = find(NoNaNIndices);
% APbinEnd = max(NoNaNIndex);
% 
% APbinstart = 9;
% APbinend = APbinEnd; %25;
% 
% % Recycle the fitted parameters from r0, which is P
% param_r0 = P;
% 
% % Initial condition
% r_R = 200; % Rate at which repressor is bound
% K_r = 10; % Dissociation of repressor to DNA
% 
% % fitting for r1 first,
% p0=[K_r, r_R];
% %lb=[0,R_max,R_bas];
% lb = [0, 0];
% ub=[K_r,r_R]*10;
% options = optimset('Display','iter');
% 
% fun= @(q)r1_Thermo_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - InitialRate_r1_NC14(APbinstart:APbinend)';
% 
% Q1_Thermo = lsqnonlin(fun, p0, lb, ub, options);
% 
% Fit_r1_Thermo = r1_Thermo_initial_rate_fit(Q1_Thermo, param_r0, BcdFluo,RuntFluo_extrapolated')

%% Hill fit (r1)

% 1) find the non-NaNs from the AP axis
NoNaNIndices = ~isnan(InitialRate_r1_NC14);
NoNaNIndex = find(NoNaNIndices);
APbinEnd = max(NoNaNIndex);

APbinstart = 9;
APbinend = APbinEnd; %25;

% 2) Parameters
% Recycle the fitted parameters from r0, which is P
param_r0 = P;
% Initial condition
K_r = 10; % Dissociation of repressor to DNA
r_R1 = 0;

% fitting for r1 first,
p0=[K_r r_R1];
%lb=[K_r r_R1];
lb = [1 0];
ub=[K_r r_R1]*100;
options = optimset('Display','iter');

fun= @(q)r1_Hill_Thermo_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - InitialRate_r1_NC14(APbinstart:APbinend)';

Q1 = lsqnonlin(fun, p0, lb, ub, options);

Fit_r1 = r1_Hill_Thermo_initial_rate_fit(Q1, param_r0, BcdFluo,RuntFluo_extrapolated')
% Plot to check (r1 fitting)

% Defined fit range
fitRange = APbinstart:APbinend;

hold on
plot(APaxis(fitRange), Fit_r1(fitRange))
plot(APaxis, InitialRate_r1_NC14,'-o')
xlim([0.15 0.6])

title('Fit r1 model to initial RNAP loading rate (NC14)')
xlabel('AP axis (embryo length)')
ylabel('Initial loaind rate (AU)')
legend('Fit-model','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope_Asymmetric';
saveas(gcf,[FigPath,filesep,'r1_InitialSlope_Hill_Fits' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath,filesep,'r1_InitialSlope_Hill_Fits' , '_NC14' , '.pdf']); 

%% r2 - Prediction (Hill)
% % Using the previously determined parameters, plugged into the r2_Hill_initial_rate_fit
% Prediction_r2_Hill_r1param = r2_Hill_initial_rate_fit(Q1, param_r0, BcdFluo,RuntFluo_extrapolated');
% 
% hold on
% plot(APaxis, Prediction_r2_Hill_r1param)
% plot(APaxis, InitialRate_r2_NC14)

%% r3 - Prediction (Hill)
% % Using the previously determined parameters, plugged into the r3_Hill_initial_rate_fit
% Prediction_r3_Hill_r1param = r3_Hill_initial_rate_fit(Q1, param_r0, BcdFluo,RuntFluo_extrapolated');
% 
% hold on
% plot(APaxis, Prediction_r3_Hill_r1param)
% plot(APaxis, InitialRate_r3_NC14)
% 
% % r3 doesn't match that well.

%% Fit parameters for r2 (Hill)
% r2 interpolation in the middle AP bins
X = [0:0.025:0.35,0.4,0.425];
indices = round(X/0.025 + 1);
% It's the 16th AP bin that we're interpolating.
InitialRate_r2_NC14(16) = interp1(X,InitialRate_r2_NC14(indices),0.375);

% % r2 extrapolation in posterior bins
% InitialRate_r2_NC13(21:25) = min(InitialRate_r2_NC14);
% plot(APaxis, InitialRate_r2_NC13)

% 1) find the non-NaNs from the AP axis
NoNaNIndices = ~isnan(InitialRate_r2_NC14);
NoNaNIndex = find(NoNaNIndices);
APbinEnd = max(NoNaNIndex);

APbinstart = 9;
APbinend = APbinEnd; %25;

% 2) Parameters
% Recycle the fitted parameters from r0, which is P
param_r0 = P;

% Initial condition
K_r = 10; % Dissociation of repressor to DNA
r_R1 = 0;
r_R2 = 0;
w_R = 1;

% fitting for r2
p0=[K_r, r_R1, r_R2, w_R];
%lb=[0,R_max,R_bas];
lb = [1 0 0 1];
ub=[100 0 0 100];
options = optimset('Display','iter');

fun= @(q)r2_Hill_Thermo_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - InitialRate_r2_NC14(APbinstart:APbinend)';

Q2 = lsqnonlin(fun, p0, lb, ub, options);

Fit_r2 = r2_Hill_Thermo_initial_rate_fit(Q2, param_r0, BcdFluo,RuntFluo_extrapolated')

% Plot to check (r2 - fitting)
fitRange = APbinstart:APbinend;

hold on
plot(APaxis(fitRange), Fit_r2(fitRange))
%plot(APaxis, Prediction_r2_Hill_r1param)
plot(APaxis, InitialRate_r2_NC14,'-o')
xlim([0.15 0.6])

title('Fit r2 model to initial RNAP loading rate (NC14)')
xlabel('AP axis (embryo length)')
ylabel('Initial loaind rate (AU)')
%legend('Fit-Hill','Prediction_{r1 parameters}','Data')
legend('Fit-model','Data')
% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope_Asymmetric';
saveas(gcf,[FigPath,filesep,'r2_InitialSlope_Hill_Fits' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath,filesep,'r2_InitialSlope_Hill_Fits' , '_NC14' , '.pdf']); 

%% Fit parameters for r3 (Hill-Thermo)
% % r3 extrapolation in posterior bins
%InitialRate_r3_NC14(14:17) = InitialRate_r3_NC14(18);%min(InitialRate_r3_NC14);
% plot(APaxis, InitialRate_r3_NC13)

% 1) find the non-NaNs from the AP axis
NoNaNIndices = ~isnan(InitialRate_r3_NC14);
NoNaNIndex = find(NoNaNIndices);
APbinEnd = max(NoNaNIndex(1:end-1)); %max(NoNaNIndex(1:end)); 
% I'm taking out the most posterior bin since it doesn't look like a 

APbinstart = 9;
APbinend = APbinEnd; %25;

% 2) Parameters
% Recycle the fitted parameters from r0, which is P
param_r0 = P;

% Initial condition
K_r = 10; % Dissociation of repressor to DNA
r_R1 = 0;
r_R2 = 0;
w_R = 10;
r_R3 = 0;

% fitting for r3
p0=[K_r, r_R1, r_R2, w_R, r_R3];
lb = [1 0 0 1 0];
ub=[10^6 0 0 10^6 0];
options = optimset('Display','iter');

fun= @(q)r3_Hill_Thermo_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - InitialRate_r3_NC14(APbinstart:APbinend)';

[Q3 ,~,residual_r3,~,output_r3] = lsqnonlin(fun, p0, lb, ub, options);

Fit_r3 = r3_Hill_Thermo_initial_rate_fit(Q3, param_r0, BcdFluo,RuntFluo_extrapolated')

% Plot to check (r3 - fitting)

fitRange = APbinstart:APbinend;

hold on
plot(APaxis(fitRange), Fit_r3(fitRange))
% plot(APaxis, Prediction_r3_Hill_r1param)
plot(APaxis, InitialRate_r3_NC14,'-o')
xlim([0.15 0.6])
%ylim([0 120])

title('Fit r3 model to initial RNAP loading rate (NC14)')
xlabel('AP axis (embryo length)')
ylabel('Initial loaind rate (AU)')
% legend('Fit-Hill','Prediction_{r1 parameters}','Data')
legend('Fit-model','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope_Asymmetric';
saveas(gcf,[FigPath,filesep,'r3_InitialSlope_Hill_Fits' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath,filesep,'r3_InitialSlope_Hill_Fits' , '_NC14' , '.pdf']); 

%% Plot the gMM coefficent as # of Runt binding sites

% I need proper error bar from the Hill fitting... How do I get an error
% bar from the fitting? like confidence intervals?
semilogy([1,2,3], [Q1(1),Q2(1),Q3(1)],'-o')
xlim([0 4])
ylim([1 10^4])

title('K_{R} vs Number of Runt sites')
xlabel('Number of Runt binding sites')
ylabel('K_{R} (AU)')
xticks([0 1 2 3 4])

StandardFigure(gcf,gca)
saveas(gcf,[FigPath,filesep,'K_R_vs_Num_Runt_sites' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath,filesep,'K_R_vs_Num_Runt_sites' , '_NC14' , '.pdf']); 

%% Plot the w_{R}, the cooperativity term for the Repressors
% I need proper error bar from the Hill fitting... How do I get an error
% bar from the fitting? like confidence intervals?
semilogy([1,2,3], [nan,Q2(4),Q3(4)],'-o')
xlim([0 4])
ylim([1 10^3])

title('w_{R} vs Number of Runt sites')
xlabel('Number of Runt binding sites')
ylabel('w_{R} (AU)')
xticks([0 1 2 3 4])

StandardFigure(gcf,gca)
saveas(gcf,[FigPath,filesep,'w_R_vs_Num_Runt_sites' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath,filesep,'w_R_vs_Num_Runt_sites' , '_NC14' , '.pdf']); 
%% Exploraration : Hill effect of # of Runt sites
% hold on
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend)))
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend).^2))
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend).^3))
% %set(gca, 'YScale', 'log')
% legend('1','2','3')
end