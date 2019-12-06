function main02_02_fit_models_initial_loading_rate_NC13
%% Part2. Fitting the initial loading rates (mean) with theoretical models
% Here, I will use Mean12, Mean13, and Mean14 for fitting with some
% expression for different constructs with different numbers of binding
% sites. Note that the Mean12,13, and 14 means that initial rate of RNAP
% loading fit from individual embryos, then averaged.

% Load the Bcd, Runt data
% Bicoid : Time-averaging is needed, but we know the Bcd gradient shape
% (length constant) is conserved, and only the amplitude changes over time
% (from Paul's fitting with Exponential), thus we'll just use the Bcd
% profile at one time point. The K_a will scale things properly anyway.
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

% Bicoid = load([FilePath, filesep, 'Bcd-Averaged.mat']);

% First, we should build a prediction matrix using handful of parameters.
% Use Rate_r0.m for optimizing the parameters, r_max, r_basal, and K_d
% using lsqnonlin.

% Second, we need to decide which nc we're going to model.
% Here, we'll start with NC13, and maybe think about NC14 later.

% Now, I will start from the initial loading rate in NC13.
% Bcd and Runt should be time-averaged for 3-5 minutes in NC13 (the initial
% rise phase).
% BcdTime = Bicoid.ElapsedTime;
% BcdTres = median(diff(BcdTime));
% BcdTstart = floor(3/BcdTres);   % frame index for 3 min into NC13
% BcdTend = floor(5/BcdTres);     % frame index for 5 min into NC13
% BcdNC13 = Bicoid.nc13;
% 
% BcdFluo = nanmean(Bicoid.MeanVectorAP(BcdNC13 + BcdTstart:BcdNC13 + BcdTend,:),1);
% BcdFluoSE = nanmean(Bicoid.SEVectorAP(BcdNC13 + BcdTstart:BcdNC13 + BcdTend,:),1);

% % Bcd extrapolation
% BcdFluo_extrap = interp1(10:27, BcdFluo(10:27), 0.2, 'pchip','extrap')
% BcdFluo(9) = BcdFluo_extrap;

% Optional : Take the Bicoid data from Liz&Jonathan
% time resolution : 30 sec.
BcdAnt = load([FilePath, filesep, 'BcdGFPAnt.mat']);

BcdNC13 = BcdAnt.DataBcd.nc13;
% Take temporal average over 3-5 minutes into NC13.
BcdFluo = nanmean(BcdAnt.DataBcd.MeanVectorAP(BcdNC13 + 6: BcdNC13 + 10,:));
BcdFluoSD = nanmean(BcdAnt.DataBcd.SDVectorAP(BcdNC13 + 6: BcdNC13 + 10,:));

Runt = load([FilePath, filesep, 'Runt_time_averaged_female.mat']);
% Runt is already time-averaged. For 3-5 minutes, it's in the 1st cell.
RuntFluo = Runt.Averaged_Fluo(1,:);
RuntFluoSE = Runt.SE_Fluo(1,:);

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

title({'Bcd and Runt profile over AP';'(time-averaged, 3-5 min in NC13)'})
xlabel('AP (Embryo Length)')
legend('Bcd','Runt')
StandardFigure(gcf,gca)

% Save figures
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope';
% saveas(gcf,[FigPath,filesep,'InputTF_Time_Averaged_NC13_3-5min_BcdRunt' , '_NC13' , '.tif']); 
% saveas(gcf,[FigPath,filesep,'InputTF_Time_Averaged_NC13_3-5min_BcdRunt' , '_NC13' , '.pdf']); 
%% Save the extrapolated Runt profile
Runt.Extrapolated_Fluo_3_5min_NC13 = RuntFluo_extrapolated;

save([FilePath,filesep,'Runt_time_averaged_female.mat'],...
   '-struct', 'Runt','-v7.3');
%% Load the fitted initial slope (Data)
% Data calculated from the script :
% main02_plot_Initial_Loading_Rates_Asymmetric.m
% 
InitialSlope = load([FilePath, filesep, 'AveragedInitialRate_AsymmetricFit.mat']);

%% Extract useful fields 
InitialRate_r0_NC13 = InitialSlope.average_fittedRate_r0(:,2);
InitialRate_r1_NC13 = InitialSlope.average_fittedRate_r1_female(:,2);
InitialRate_r2_NC13 = InitialSlope.average_fittedRate_r2_female(:,2);
InitialRate_r3_NC13 = InitialSlope.average_fittedRate_r3_female(:,2);

InitialRate_SEM_r0_NC13 = InitialSlope.SEM_fittedRate_r0(:,2);
InitialRate_SEM_r1_NC13 = InitialSlope.SEM_fittedRate_r1_female(:,2);
InitialRate_SEM_r2_NC13 = InitialSlope.SEM_fittedRate_r2_female(:,2);
InitialRate_SEM_r3_NC13 = InitialSlope.SEM_fittedRate_r3_female(:,2);

%% Initial Rate extrapolation for the posterior bins
APbin1 = 0.15./0.025 + 1;
APbin2 = 0.5 / 0.025 + 1;
%InitialRate_r0_NC13_extrap = interp1([0.15:0.025:0.5], InitialRate_r0_NC13(APbin1:APbin2), [0.525:0.025:0.625], 'pchip', 'extrap')

InitialRate_r0_NC13(22:26) = min(InitialRate_r0_NC13);

plot(APaxis, InitialRate_r0_NC13)
%% limit the AP bins for fitting
APbinstart = 9;
APbinend = 25;

% Initial condition
R_max = 300; % maximum
R_bas = 50; % basal, probably due to ectopic expression...
Kd = 10;
N = 6;
% fitting for r0 first,
p0=[R_bas,R_max,Kd];%,N];
%lb=[0,R_max,R_bas];
lb = [0, 0, 0];%, 0];
ub=[R_bas*10, R_max*10, Kd*10];%, N*10];
options = optimset('Display','iter');

fun= @(p)r0_Hill_initial_rate_fit(p,BcdFluo(APbinstart:APbinend)) - InitialRate_r0_NC13(APbinstart:APbinend)';

P = lsqnonlin(fun, p0, lb, ub, options);

% Plot the fit to check if it's reasonable
r_basal = P(1);
r = P(2);
K_a = P(3);
%N = P(4);
N=6; % Fix the Hill coeff. to be 6, corresponding to the # of binding sites.
BcdData = BcdFluo;

Prediction = nan(41,1);

Prediction = (r_basal*ones(size(BcdData)) + (BcdData./K_a).^N  * r) ./ (1 + (BcdData./K_a).^N)

hold on
%plot(APaxis, (r_basal*ones(size(BcdData))))
%plot(APaxis, (BcdData./K_a).^6  * r)

plot(APaxis, Prediction)
plot(APaxis, InitialRate_r0_NC13,'-o')
title('Fit r0 model to initial RNAP loading rate (NC13)')
xlabel('AP axis (embryo length)')
ylabel('Initial loaind rate (AU)')
legend('Fit','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope';
% saveas(gcf,[FigPath,filesep,'r0_InitialSlope_Hill_Fits' , '_NC13' , '.tif']); 
% saveas(gcf,[FigPath,filesep,'r0_InitialSlope_Hill_Fits' , '_NC13' , '.pdf']); 

%% r1 fitting
%% First, extrapolate the r1 at posterior positions
APbin1 = 0.15./0.025 + 1;
APbin2 = 0.5 / 0.025 + 1;
%InitialRate_r0_NC13_extrap = interp1([0.15:0.025:0.5], InitialRate_r0_NC13(APbin1:APbin2), [0.525:0.025:0.625], 'pchip', 'extrap')

InitialRate_r1_NC13(21:25) = min(InitialRate_r1_NC13);

plot(APaxis, InitialRate_r1_NC13)

%% Fitting (r1) - Thermodynamic model
% Recycle the fitted parameters from r0, which is P
param_r0 = P;
APbinstart = 9;
APbinend = 25;

% Initial condition
r_R = 200; % Rate at which repressor is bound
K_r = 10; % Dissociation of repressor to DNA

% fitting for r1 first,
p0=[K_r, r_R];
%lb=[0,R_max,R_bas];
lb = [0, 0];
ub=[K_r,r_R]*10;
options = optimset('Display','iter');

fun= @(q)r1_Thermo_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - InitialRate_r1_NC13(APbinstart:APbinend)';

Q = lsqnonlin(fun, p0, lb, ub, options);

Prediction_r1_Thermo = r1_Thermo_initial_rate_fit(Q, param_r0, BcdFluo,RuntFluo_extrapolated')

%% Hill fit (r1)
% Recycle the fitted parameters from r0, which is P
param_r0 = P;
APbinstart = 9;
APbinend = 25;

% Initial condition
K_r = 10; % Dissociation of repressor to DNA

% fitting for r1 first,
p0=[K_r];
%lb=[0,R_max,R_bas];
lb = [0];
ub=[K_r]*10;
options = optimset('Display','iter');

fun= @(q)r1_Hill_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - InitialRate_r1_NC13(APbinstart:APbinend)';

Q_Hill = lsqnonlin(fun, p0, lb, ub, options);

Prediction_r1_Hill = r1_Hill_initial_rate_fit(Q_Hill, param_r0, BcdFluo,RuntFluo_extrapolated')
%% Plot to check (r1 fitting)
hold on
%plot(APaxis, (r_basal*ones(size(BcdData))))
%plot(APaxis, (BcdData./K_a).^6  * r)

plot(APaxis, Prediction_r1_Thermo)
pause
plot(APaxis, Prediction_r1_Hill)
pause
plot(APaxis, InitialRate_r1_NC13,'-o')
title('Fit r1 model to initial RNAP loading rate (NC13)')
xlabel('AP axis (embryo length)')
ylabel('Initial loaind rate (AU)')
legend('Fit-Thermo','Fit-Hill','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope';
% saveas(gcf,[FigPath,filesep,'r1_InitialSlope_Hill_Fits' , '_NC13' , '.tif']); 
% saveas(gcf,[FigPath,filesep,'r1_InitialSlope_Hill_Fits' , '_NC13' , '.pdf']); 

%% r2 - Prediction (Hill)
% Using the previously determined parameters, plugged into the r2_Hill_initial_rate_fit
Prediction_r2_Hill_r1param = r2_Hill_initial_rate_fit(Q_Hill, param_r0, BcdFluo,RuntFluo_extrapolated');

hold on
plot(APaxis, Prediction_r2_Hill_r1param)
plot(APaxis, InitialRate_r2_NC13)

%% r3 - Prediction (Hill)
% Using the previously determined parameters, plugged into the r3_Hill_initial_rate_fit
Prediction_r3_Hill_r1param = r3_Hill_initial_rate_fit(Q_Hill, param_r0, BcdFluo,RuntFluo_extrapolated');

hold on
plot(APaxis, Prediction_r3_Hill_r1param)
plot(APaxis, InitialRate_r3_NC13)

% r3 doesn't match that well.

%% Fit parameters for r2 (Hill)
% r2 extrapolation in posterior bins
InitialRate_r2_NC13(21:25) = min(InitialRate_r2_NC13);

%plot(APaxis, InitialRate_r2_NC13)

% Recycle the fitted parameters from r0, which is P
param_r0 = P;
APbinstart = 9;
APbinend = 25;

% Initial condition
K_r = 10; % Dissociation of repressor to DNA

% fitting for r2 first,
p0=[K_r];
%lb=[0,R_max,R_bas];
lb = [0];
ub=[K_r]*10;
options = optimset('Display','iter');

fun= @(q)r2_Hill_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - InitialRate_r2_NC13(APbinstart:APbinend)';

Q_Hill2 = lsqnonlin(fun, p0, lb, ub, options);

Prediction_r2_Hill = r2_Hill_initial_rate_fit(Q_Hill2, param_r0, BcdFluo,RuntFluo_extrapolated')

%% Plot to check (r2 - fitting)

hold on

plot(APaxis, Prediction_r2_Hill)
plot(APaxis, Prediction_r2_Hill_r1param)
plot(APaxis, InitialRate_r2_NC13,'-o')
title('Fit r2 model to initial RNAP loading rate (NC13)')
xlabel('AP axis (embryo length)')
ylabel('Initial loaind rate (AU)')
legend('Fit-Hill','Prediction_{r1 parameters}','Data')

% Save figure
StandardFigure(gcf,gca)
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope';
% saveas(gcf,[FigPath,filesep,'r2_InitialSlope_Hill_Fits' , '_NC13' , '.tif']); 
% saveas(gcf,[FigPath,filesep,'r2_InitialSlope_Hill_Fits' , '_NC13' , '.pdf']); 

%% Fit parameters for r3 (Hill)
% r2 extrapolation in posterior bins
InitialRate_r3_NC13(21:25) = min(InitialRate_r3_NC13);

plot(APaxis, InitialRate_r3_NC13)

% Recycle the fitted parameters from r0, which is P
param_r0 = P;
APbinstart = 9;
APbinend = 25;

% Initial condition
K_r = 10; % Dissociation of repressor to DNA

% fitting for r2 first,
p0=[K_r];
%lb=[0,R_max,R_bas];
lb = [0];
ub=[K_r]*10;
options = optimset('Display','iter');

fun= @(q)r3_Hill_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - InitialRate_r3_NC13(APbinstart:APbinend)';

Q_Hill3 = lsqnonlin(fun, p0, lb, ub, options);

Prediction_r3_Hill = r3_Hill_initial_rate_fit(Q_Hill3, param_r0, BcdFluo,RuntFluo_extrapolated')

%% Plot to check (r3 - fitting)

hold on

plot(APaxis, Prediction_r3_Hill)
plot(APaxis, Prediction_r3_Hill_r1param)
plot(APaxis, InitialRate_r3_NC13,'-o')
title('Fit r3 model to initial RNAP loading rate (NC13)')
xlabel('AP axis (embryo length)')
ylabel('Initial loaind rate (AU)')
legend('Fit-Hill','Prediction_{r1 parameters}','Data')

% Save figure
StandardFigure(gcf,gca)
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope';
% saveas(gcf,[FigPath,filesep,'r3_InitialSlope_Hill_Fits' , '_NC13' , '.tif']); 
% saveas(gcf,[FigPath,filesep,'r3_InitialSlope_Hill_Fits' , '_NC13' , '.pdf']); 

%% Plot the Hill coefficient as # of Runt binding sites

% I need proper error bar from the Hill fitting... How do I get an error
% bar from the fitting? like confidence intervals?
plot([1,2,3], [Q_Hill,Q_Hill2,Q_Hill3],'-o')
xlim([0 4])

title('Hill Coefficient vs # of Runt sites')
xlabel('Number of Runt binding sites')
ylabel('K_{R} (AU)')
xticks([0 1 2 3 4])

StandardFigure(gcf,gca)
% saveas(gcf,[FigPath,filesep,'Hill_Coeff_K_R_#Runt_sites' , '_NC13' , '.tif']); 
% saveas(gcf,[FigPath,filesep,'Hill_Coeff_K_R_#Runt_sites' , '_NC13' , '.pdf']); 
%% Exploraration : Hill effect of # of Runt sites
% hold on
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend)))
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend).^2))
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend).^3))
% %set(gca, 'YScale', 'log')
% legend('1','2','3')
end