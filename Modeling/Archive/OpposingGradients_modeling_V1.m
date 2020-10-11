function OpposingGradients_modeling_V1
% hbP2 + N Runt binding sites constructs modeling with actual Bcd and Runt
% gradient.
% Last updated : Mar 2020, by YJK

% Assumptions : 
% 1) Bicoid binds cooperatively, and could be desribed using Hill function.
% 2) Runt is not known for its cooperativity, thus I'll not assume anything
% on this.

%% Step1. Get the actual input TF 
% Here, I'll use the time-averaged Bcd and Runt profiles
% Actually, the time-averaging "time window" doesn't matter that much,
% since the gradient scales nicely over time.

filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

%BcdData = load([filePath, filesep, 'Bcd-Averaged.mat']);
RuntData = load([filePath, filesep, 'Runt_TimeAveraged_mixedSex_NC14.mat']);

BcdData = load([filePath, filesep, 'BcdGFPAnt.mat']);
%% Step1-1 : Processing the inputs - Bcd(time-averaging)
% We need a better way of doing this : calculating the time-averaged
% profile for individual embryo, then get mean and std(SEM) over embryos.

BcdData = BcdData.DataBcd;
BcdTime = BcdData.ElapsedTime;
BcdFluo = BcdData.MeanVectorAP;
BcdNC13= BcdData.nc13;
BcdNC14 = BcdData.nc14;

% Time-average with 10 minute time-window.
tInitial = 0; %0 minutes
tWindow = 10; % 10 minutes
timestep = nanmean(diff(BcdTime));
n_initial = ceil(tInitial/timestep);
n_timesteps = ceil(tWindow/timestep);
tWindowIndices = BcdNC14+n_initial:BcdNC14+n_timesteps;

BcdFluo_tAveraged = nanmean(BcdFluo(tWindowIndices,:));
SDBcdFluo_tAveraged = nanstd(BcdFluo(tWindowIndices,:),0);
SEBcdFluo_tAveraged = SDBcdFluo_tAveraged./sqrt(n_timesteps); % divide by sqrt of number of data points, we need to change this to number of embryos at some point.

%% Runt
RuntFluo_tAveraged = RuntData.AveragedFluo_tAveraged_mixed(2,:); % 0-10 mins into NC14
SERuntFluo_tAveraged = RuntData.SEFluo_tAveraged_mixed(2,:);% 0-10 mins into NC14
%% Quick plot to check
hold on
APaxis = 0:0.025:1;
yyaxis left
errorbar(APaxis, BcdFluo_tAveraged/max(BcdFluo_tAveraged), SEBcdFluo_tAveraged/max(BcdFluo_tAveraged))
xlim([0.2 0.6])
ylim([0 1.2])
ylabel('Bcd concentration (Normalized)')
yyaxis right
errorbar(APaxis, RuntFluo_tAveraged/max(RuntFluo_tAveraged), SERuntFluo_tAveraged/max(RuntFluo_tAveraged))
ylim([0 1.2])
ylabel('Runt concentration (Normalized)')
xlabel('AP axis (EL)')

StandardFigure(gcf,gca)
% save the plot
figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\Bcd_Runt_together_Mar2020';
saveas(gcf,[figPath,filesep,'Bcd_Runt_AP_tAveraged_0_10minNC14_20_60%.tif'])
saveas(gcf,[figPath,filesep,'Bcd_Runt_AP_tAveraged_0_10minNC14_20_60%.pdf'])
%saveas(gcf,[figPath,filesep,'Bcd_Runt_AP_tAveraged_0_10minNC14.pdf'])

%% 
%% Step 1-2 : Define Bcd and Runt using these time-averaged vectors

Bcd = BcdFluo_tAveraged;
Runt = RuntFluo_tAveraged;

% Transpose to plug in as inputs
Bcd = Bcd';
Runt = Runt';

%% Step2. Start with a simple model for r0 (6 Bcd binding sites)
% Note that this could be expanded to 11 or more Bcd sites.

%% Case0. hb P2, without any Runt binding sites
% free parameters : K_A, p,w_ap, r (rate when RNAP is bound to the
% promoter)
% Use rate_r0_Hill for predicting the initial rate of RNAP loading

X = 0.2:0.025:0.6;

% initial rate of RNAP loading 
initialSlopes = load([filePath, filesep, 'InitialSlopes_ONnuclei_AllConstructs.mat'])
initialSlopes = initialSlopes.initialSlopes_ONnuclei;

initialSlope_r0 = cell2mat(initialSlopes(2,5));
initialSlope_r0_SEM = cell2mat(initialSlopes(2,6));

initialSlope_r0 = initialSlope_r0(:,3); % NC14
initialSlope_r0_SEM = initialSlope_r0_SEM(:,3); % NC14

fitRange = 9:18; % 20-42.5% of AP axis as a fit range

% nonlinear least square fitting to the model
x0 = [20 0.1 10 200]; %[K_a p wAP r]
% Note that the r is constrained by the maximum of initial rate of RNAP
% loading of r0 (or the maximum value of rate from all constructs, probably
% r1-close).

lb = [0 0 0 0];
ub = [1000 100 10000 10000];
options.Algorithm = 'levenberg-marquardt';
%lsqOptions=optimset('Display','none');

fun = @(x)rate_r0_Hill(x,Bcd(fitRange)) - initialSlope_r0(fitRange);
x_r0 = lsqnonlin(fun,x0,lb, ub, options);


% x_r0 = lsqnonlin(@rate_r0_Hill, x0, lb, ub, options, Bcd(fitRange),...
%         initialSlope_r0(fitRange))%
% x = lsqnonlin(fun,x0,[],[],[]) 

%% Check the parameter sensitivity
% for two sets of parameters
Prediction1 = rate_r0_Hill(x_r0_1,Bcd)
Prediction2 = rate_r0_Hill(x_r0_2,Bcd)

hold on
errorbar(APaxis, initialSlope_r0, initialSlope_r0_SEM)
% fitted by chi2
plot(APaxis, Prediction1)
plot(APaxis, Prediction2)

% legend('data',num2str(x_r0_1),num2str(x_r0_2))
legend('data','parameter1','parameter2')
xlabel('AP axis')
ylabel('initial rate of RNAP loading (AU)')

StandardFigure(gcf,gca)
saveas(gcf,[figPath,filesep,'initialRate_chi2_r0_parameter_sensitivityTest' , '.tif']); 
saveas(gcf,[figPath,filesep,'initialRate_chi2_r0_parameter_sensitivityTest' , '.tif']); 
%% Plot the data and the fit
hold on
% data (initial fit)
errorbar(APaxis, initialSlope_r0, initialSlope_r0_SEM)
% fitted by chi2
plot(APaxis, rate_r0_Hill(x_r0,Bcd))

xlabel('AP axis')
ylabel('initial rate of RNAP loading (AU)')
StandardFigure(gcf,gca)
% save the plot
figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_hbP2_nRuntsites';
saveas(gcf,[figPath,filesep,'initialRate_chi2_r0' , '.tif']); 
saveas(gcf,[figPath,filesep,'initialRate_chi2_r0' , '.pdf']); 

%% r1

% Data extract
initialSlope_r1 = cell2mat(initialSlopes(3,5));
initialSlope_r1_SEM = cell2mat(initialSlopes(3,6));

initialSlope_r1 = initialSlope_r1(:,3); % NC14
initialSlope_r1_SEM = initialSlope_r1_SEM(:,3); % NC14

% Use the parameters fitted from the r0, then just adjust additional
% parameters.
x_r0; %[K_a p wAP r]
y0 = [5 0.01]; % K_r, wRP

lb = [5 0]; 
ub = [100 10];
options.Algorithm = 'levenberg-marquardt';
%lsqOptions=optimset('Display','none');

fun1 = @(y)rate_r1_Hill(x_r0,y,Bcd(fitRange),Runt(fitRange)) - initialSlope_r1(fitRange);
y_r1 = lsqnonlin(fun1,y0,lb, ub, options);

%% Test the parameter sensitivity, K_r
K_r = 10.^[-5 -4 -3 -2 -1 0 1 2 3 4 5];

hold on
errorbar(APaxis, initialSlope_r1, initialSlope_r1_SEM)
for i=1:length(K_r)
    Prediction = rate_r1_Hill(x_r0, [K_r(i), y_r1(2)], Bcd, Runt)
    plot(APaxis, Prediction)
end
%% Plot to check r1 and r1-from chi2
Prediction = rate_r1_Hill(x_r0, y_r1, Bcd, Runt)

hold on 
errorbar(APaxis, initialSlope_r1, initialSlope_r1_SEM)
plot(APaxis, Prediction)

xlabel('AP axis')
ylabel('initial rate of RNAP loading (AU)')
StandardFigure(gcf,gca)
% save the plot
figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_hbP2_nRuntsites';
saveas(gcf,[figPath,filesep,'initialRate_chi2_r1' , '.tif']); 
saveas(gcf,[figPath,filesep,'initialRate_chi2_r1' , '.pdf']); 
%% r1-close
% Data extract
initialSlope_r1_close = cell2mat(initialSlopes(6,5));
initialSlope_r1_close_SEM = cell2mat(initialSlopes(6,6));

initialSlope_r1_close = initialSlope_r1_close(:,3); % NC14
initialSlope_r1_close_SEM = initialSlope_r1_close_SEM(:,3); % NC14

% Use the parameters fitted from the r0, then just adjust additional
% parameters.
x_r0; %[K_a p wAP r]
y0 = [5 0.01];% K_r, wRP

% constrain the K_r as obtained in r1 case
lb = [y_r1(1) 0]; 
ub = [y_r1(1) 10];
% lb = [1 0]; 
% ub = [100 10];
options.Algorithm = 'levenberg-marquardt';
%lsqOptions=optimset('Display','none');

fun1 = @(y)rate_r1_Hill(x_r0,y,Bcd(fitRange),Runt(fitRange)) - initialSlope_r1_close(fitRange);
y_r1_close = lsqnonlin(fun1,y0,lb, ub, options);

%% Plot to check r1-close and r1-close-from chi2
Prediction = rate_r1_Hill(x_r0, y_r1_close, Bcd, Runt)

hold on 
errorbar(APaxis, initialSlope_r1_close, initialSlope_r1_close_SEM)
plot(APaxis, Prediction)

xlabel('AP axis')
ylabel('initial rate of RNAP loading (AU)')
StandardFigure(gcf,gca)
% save the plot
figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_hbP2_nRuntsites';
saveas(gcf,[figPath,filesep,'initialRate_chi2_r1_close' , '.tif']); 
saveas(gcf,[figPath,filesep,'initialRate_chi2_r1_close' , '.pdf']);

%% Predict the 2 binding sites ([1,0,0] + [0,0,1])
Prediction_r2_far = rate_r2_Hill_multiRep(x0,[y_r1(1), y_r1(2), y_r1_close(2)],Bcd, Runt);

% actual r2([1,1,0]) initial slope
initialSlope_r2_far = cell2mat(initialSlopes(9,5));
initialSlope_r2_far_SEM = cell2mat(initialSlopes(9,6));

initialSlope_r2_far = initialSlope_r2_far(:,3); % NC14
initialSlope_r2_far_SEM = initialSlope_r2_far_SEM(:,3); % NC14

%% plot r2 prediciton and data
hold on
plot(APaxis, Prediction_r2_far)
errorbar(APaxis, initialSlope_r2_far, initialSlope_r2_far_SEM)
end