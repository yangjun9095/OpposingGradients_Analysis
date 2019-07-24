function main03_02_fit_models_AccumulatedmRNA
%% Part2. Fitting the Accumulated mRNA profile with theoretical models
% Here, I will use AccumulatedmRNA from NC13 (or NC14) for fitting with some
% expression for different constructs with different numbers of binding
% sites. 
% Note that I used IntegratemRNA.m in here.

% Load the Bcd, Runt data
% Bicoid : Time-averaging is needed, but we know the Bcd gradient shape
% (length constant) is conserved, and only the amplitude changes over time
% (from Paul's fitting with Exponential), thus we'll just use the Bcd
% profile at one time point. The K_a will scale things properly anyway.
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

% First, we should build a prediction matrix using handful of parameters.
% Use Rate_r0.m for optimizing the parameters, r_max, r_basal, and K_d
% using lsqnonlin.
% Second, we need to decide which nc we're going to model.
% Here, we'll start with NC13, and maybe think about NC14 later.

% Optional : Take the Bicoid data from Liz&Jonathan
BcdAnt = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\eGFP-Bcd-From-Liz-Jonathan\BcdGFPAnt.mat');
% time resolution : 30 sec.
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
RuntFluo_extrap = interp1([0.3:0.025:0.625], RuntFluo(13:26), [0.2:0.025:0.275], 'pchip', 'extrap');

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
%% Load the Accumulated mRNA (Data)
% Data calculated from the script :
% main12_compare_r0123.m, using IntegratemRNA.m script.
% Optional : how this accumulated mRNA is calculated 
% 1) Load the datasets
r0Data = LoadMS2Sets('r0','dontCompare')
r1Data = LoadMS2Sets('r1-new-female','dontCompare')
r2Data = LoadMS2Sets('r2-new-female','dontCompare')
r3Data = LoadMS2Sets('r3-new-female','dontCompare')

% 2) IntegratemRNA
[TotalProd_r0,TotalProdError_r0,TotalProdN_r0,...
    MeanTotalProd_r0,SDTotalProd_r0,SETotalProd_r0]=IntegratemRNA(r0Data,1,2);

[TotalProd_r1,TotalProdError_r1,TotalProdN_r1,...
    MeanTotalProd_r1,SDTotalProd_r1,SETotalProd_r1]=IntegratemRNA(r1Data,1,2);

[TotalProd_r2,TotalProdError_r2,TotalProdN_r2,...
    MeanTotalProd_r2,SDTotalProd_r2,SETotalProd_r2]=IntegratemRNA(r2Data,1,2);

[TotalProd_r3,TotalProdError_r3,TotalProdN_r3,...
    MeanTotalProd_r3,SDTotalProd_r3,SETotalProd_r3]=IntegratemRNA(r3Data,1,2);

%% Plot the Accumulated mRNA
APaxis = 0:0.025:1;
% Figure Path

% NC13
NC= 13;
AccumulatedmRNA_NC13_figure = figure
hold on
errorbar(APaxis, MeanTotalProd_r0(:,NC), SETotalProd_r0(:,NC))
errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
errorbar(APaxis, MeanTotalProd_r2(:,NC), SETotalProd_r2(:,NC))
errorbar(APaxis, MeanTotalProd_r3(:,NC), SETotalProd_r3(:,NC))
xlim([0.2 0.6])

title('Accumulated mRNA over AP @ NC13')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r0','r1','r2','r3')

StandardFigure(AccumulatedmRNA_NC13_figure,AccumulatedmRNA_NC13_figure.CurrentAxes)
%saveas(AccumulatedmRNA_NC13_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC13','_SE','.tif'])
%saveas(AccumulatedmRNA_NC13_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC13','_SE','.pdf'])

% NC14
NC = 14;
AccumulatedmRNA_NC14_figure = figure
hold on
errorbar(APaxis, MeanTotalProd_r0(:,NC), SETotalProd_r0(:,NC))
errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
errorbar(APaxis, MeanTotalProd_r2(:,NC), SETotalProd_r2(:,NC))
errorbar(APaxis, MeanTotalProd_r3(:,NC), SETotalProd_r3(:,NC))

xlim([0.2 0.5])

title('Accumulated mRNA over AP @ NC14')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r0','r1','r2','r3')

StandardFigure(AccumulatedmRNA_NC14_figure,AccumulatedmRNA_NC14_figure.CurrentAxes)
%saveas(AccumulatedmRNA_NC14_figure,[FigPath,filesep,'AccumulaatedmRNA_r0123',DataType(1:end-1),'_NC14','_SE','.tif'])
%saveas(AccumulatedmRNA_NC14_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC14','_SE','.pdf'])

%% Defining the fields that will be used for the fitting.
AccumulatedmRNA_r0_NC13 = MeanTotalProd_r0(:,13);
AccumulatedmRNA_r1_NC13 = MeanTotalProd_r1(:,13);
AccumulatedmRNA_r2_NC13 = MeanTotalProd_r2(:,13);
AccumulatedmRNA_r3_NC13 = MeanTotalProd_r3(:,13);

AccumulatedmRNA_SE_r0_NC13 = SETotalProd_r0(:,13);
AccumulatedmRNA_SE_r1_NC13 = SETotalProd_r1(:,13);
AccumulatedmRNA_SE_r2_NC13 = SETotalProd_r2(:,13);
AccumulatedmRNA_SE_r3_NC13 = SETotalProd_r3(:,13);

%% Extrapolate the posterior data points
% (Not sure if I need to, but will do anyway)
% Here, I'll make some assumptions based on how the raw data looks like.
% There are data point up to 55%, thus will assume 55.25- 60% values as
% the same as 55%.
AccumulatedmRNA_r0_NC13(23:25) = AccumulatedmRNA_r0_NC13(22);
AccumulatedmRNA_r1_NC13(23:25) = AccumulatedmRNA_r1_NC13(22);
AccumulatedmRNA_r2_NC13(23:25) = AccumulatedmRNA_r2_NC13(22);
AccumulatedmRNA_r3_NC13(23:25) = AccumulatedmRNA_r3_NC13(22);

%% limit the AP bins for fitting
APbinstart = 9;
APbinend = 25;

% Initial condition
R_max = 15000; % maximum
R_bas = 50; % basal, probably due to ectopic expression...
Kd = 10;
N = 6;
% fitting for r0 first,
p0=[R_bas,R_max,Kd];%,N];
%lb=[0,R_max,R_bas];
lb = [0, 0, 0];%, 0];
ub=[R_bas*10, R_max*10, Kd*10];%, N*10];
options = optimset('Display','iter');

fun= @(p)r0_Hill_initial_rate_fit(p,BcdFluo(APbinstart:APbinend)) - AccumulatedmRNA_r0_NC13(APbinstart:APbinend)';

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
plot(APaxis, AccumulatedmRNA_r0_NC13,'-o')
xlim([0.2 0.6])

title('Fit r0 model to Accumulated mRNA (NC13)')
xlabel('AP axis (embryo length)')
ylabel('Accumulated mRNA (AU)')
legend('Fit','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_AccumulatedmRNA';
saveas(gcf,[FigPath,filesep,'r0_AccumulatedmRNA_Hill_Fits' , '_NC13' , '.tif']); 
saveas(gcf,[FigPath,filesep,'r0_AccumulatedmRNA_Hill_Fits' , '_NC13' , '.pdf']); 

%% r1 fitting
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

fun= @(q)r1_Hill_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - AccumulatedmRNA_r1_NC13(APbinstart:APbinend)';

Q_Hill = lsqnonlin(fun, p0, lb, ub, options);

Prediction_r1_Hill = r1_Hill_initial_rate_fit(Q_Hill, param_r0, BcdFluo,RuntFluo_extrapolated')
%% Plot to check (r1 fitting)
hold on
%plot(APaxis, (r_basal*ones(size(BcdData))))
%plot(APaxis, (BcdData./K_a).^6  * r)

%plot(APaxis, Prediction_r1_Thermo)
%pause
plot(APaxis, Prediction_r1_Hill)
%pause
plot(APaxis, AccumulatedmRNA_r1_NC13,'-o')
xlim([0.2 0.6])
title('Fit r1 model to Accumulated mRNA (NC13)')
xlabel('AP axis (embryo length)')
ylabel('Accumulated mRNA (AU)')
legend('Fit-Hill','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_AccumulatedmRNA';
saveas(gcf,[FigPath,filesep,'r1_AccumulatedmRNA_Hill_Fits' , '_NC13' , '.tif']); 
saveas(gcf,[FigPath,filesep,'r1_AccumulatedmRNA_Hill_Fits' , '_NC13' , '.pdf']); 

%% r2 - Prediction (Hill)
% Using the previously determined parameters, plugged into the r2_Hill_initial_rate_fit
Prediction_r2_Hill_r1param = r2_Hill_initial_rate_fit(Q_Hill, param_r0, BcdFluo,RuntFluo_extrapolated');

hold on
plot(APaxis, Prediction_r2_Hill_r1param)
plot(APaxis, AccumulatedmRNA_r2_NC13)
xlim([0.2 0.6])
%% r3 - Prediction (Hill)
% Using the previously determined parameters, plugged into the r3_Hill_initial_rate_fit
Prediction_r3_Hill_r1param = r3_Hill_initial_rate_fit(Q_Hill, param_r0, BcdFluo,RuntFluo_extrapolated');

hold on
plot(APaxis, Prediction_r3_Hill_r1param)
plot(APaxis, AccumulatedmRNA_r3_NC13)
xlim([0.2 0.6])
% r3 doesn't match that well.

%% Fit parameters for r2 (Hill)

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

fun= @(q)r2_Hill_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - AccumulatedmRNA_r2_NC13(APbinstart:APbinend)';

Q_Hill2 = lsqnonlin(fun, p0, lb, ub, options);

Prediction_r2_Hill = r2_Hill_initial_rate_fit(Q_Hill2, param_r0, BcdFluo,RuntFluo_extrapolated')

%% Plot to check (r2 - fitting)

hold on

plot(APaxis, Prediction_r2_Hill)
plot(APaxis, Prediction_r2_Hill_r1param)
plot(APaxis, AccumulatedmRNA_r2_NC13,'-o')
xlim([0.2 0.6])
title('Fit r2 model to Accumulated mRNA (NC13)')
xlabel('AP axis (embryo length)')
ylabel('Accumulated mRNA (AU)')
legend('Fit-Hill','Prediction_{r1 parameters}','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_AccumulatedmRNA';
saveas(gcf,[FigPath,filesep,'r2_AccumulatedmRNA_Hill_Fits' , '_NC13' , '.tif']); 
saveas(gcf,[FigPath,filesep,'r2_AccumulatedmRNA_Hill_Fits' , '_NC13' , '.pdf']); 

%% Fit parameters for r3 (Hill)

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
ub=[K_r]*1000;
options = optimset('Display','iter');

fun= @(q)r3_Hill_initial_rate_fit(q,param_r0,BcdFluo(APbinstart:APbinend),RuntFluo_extrapolated(APbinstart:APbinend)') - AccumulatedmRNA_r3_NC13(APbinstart:APbinend)';

Q_Hill3 = lsqnonlin(fun, p0, lb, ub, options);

Prediction_r3_Hill = r3_Hill_initial_rate_fit(Q_Hill3, param_r0, BcdFluo,RuntFluo_extrapolated')

%% Plot to check (r3 - fitting)

hold on

plot(APaxis, Prediction_r3_Hill)
plot(APaxis, Prediction_r3_Hill_r1param)
plot(APaxis, AccumulatedmRNA_r3_NC13,'-o')
title('Fit r3 model to Accumulated mRNA (NC13)')
xlabel('AP axis (embryo length)')
ylabel('Accumulated mRNA (AU)')
legend('Fit-Hill','Prediction_{r1 parameters}','Data')

% Save figure
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_AccumulatedmRNA';
saveas(gcf,[FigPath,filesep,'r3_AccumulatedmRNA_Hill_Fits' , '_NC13' , '.tif']); 
saveas(gcf,[FigPath,filesep,'r3_AccumulatedmRNA_Hill_Fits' , '_NC13' , '.pdf']); 
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
saveas(gcf,[FigPath,filesep,'Hill_Coeff_K_R_#Runt_sites' , '_NC13' , '.tif']); 
saveas(gcf,[FigPath,filesep,'Hill_Coeff_K_R_#Runt_sites' , '_NC13' , '.pdf']); 
%% Exploraration : Hill effect of # of Runt sites
% hold on
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend)))
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend).^2))
% plot(APaxis(APbinstart:APbinend),1./(1+ RuntFluo_extrapolated(APbinstart:APbinend).^3))
% %set(gca, 'YScale', 'log')
% legend('1','2','3')
end