function 
% hbP2 + N Runt binding sites constructs modeling with actual Bcd and Runt
% gradient.

% Assumptions : 
% 1) Bicoid binds cooperatively, and could be desribed using Hill function.
% 2) Runt is not known for its cooperativity, thus I'll not assume anything
% on this.

%% Step1. Think about realistic morphogen gradients



%% Step2. Predicting the transcriptional output for different scenarios
% First, I will start with 4 scenarios, competition, quenching, direct,
% shut-down, with different forms of equations with stat.mech.

% First, let's think about the scaling of Bcd and Runt, which is from
% the dissociation constant.
% This scaling represents Kd for each protein
BcdScale = [10^(-3) 10^(-2) 10^(-1) 1 10 10^2 10^3 ];
RuntScale = [ 10^(-3) 10^(-2) 10^(-1) 1 10 10^2 10^3 ];

%% Case0. hb P2, without any Runt binding sites
% free parameters : K_A, r_basal, and r

X = 0:0.01:1;

% Real data (r0) from inital slope fitting (FitMeanAPAsymmetric)
% compiled in main02_plot_Initial_Loading_Rates_Asymmetric.m
% Note, here are 4 datasets that are averaged.
nc=2; % NC13;
%errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc))
% Fill in the values in the posterior, up to 80%
average_fittedRate_r0(22:33,nc) = 40; %average_fittedRate_r0(21,nc)

% interpolate the averaged_fittedRate_r0
initial_rate_r0 = interp1(0:0.025:1,average_fittedRate_r0(:,nc),0:0.01:1);

% fitRange
fitRange = 21:81;

% fun = @(x)x(1)*ones(size(BcdData(fitRange))) + (BcdData(fitRange)./x(3)).^6 * x(2) - initial_rate_r0(fitRange);
% Rate = r_basal + (Bcd/K_a).^6*r ;
% nonlinear least square fitting to the model
x0 = [40 1 0.24];%[1  0.1 0.1]; %[r_basal r K_a]
lb = [0, 0, 0];
ub = [60 200 200];
options.Algorithm = 'levenberg-marquardt';
%lsqOptions=optimset('Display','none');

x = lsqnonlin(@r0_Hill_initial_rate_fit, x0, lb, ub, options, Bcd_interp_nc13(fitRange),...
        initial_rate_r0(fitRange))%
% x = lsqnonlin(fun,x0,[],[],[])    
%% Plot the data and the fit
hold on
% data (initial fit)
plot(X,Bcd_interp_nc13*100)
plot(X(fitRange),initial_rate_r0(fitRange))
% fitted 
plot(X(fitRange), (x(1)*ones(size(Bcd_interp_nc13(fitRange))) + (Bcd_interp_nc13(fitRange)./x(3)).^6  * x(2))./(1 + (Bcd_interp_nc13(fitRange)./x(3)).^6))

% save the plot

%% r1
% Use the parameters fitted from the r0, then just adjust additional
% parameters.
r_basal = x(1);
r = x(2);
K_a = x(3);


% Real data (r0) from inital slope fitting (FitMeanAPAsymmetric)
% compiled in main02_plot_Initial_Loading_Rates_Asymmetric.m
% Note, here are 4 datasets that are averaged.
nc=2; % NC13;
%errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc))
% Fill in the values in the posterior, up to 80%
average_fittedRate_r1(21:33,nc) = 40; %average_fittedRate_r0(21,nc)

% interpolate the averaged_fittedRate_r0
initial_rate_r1 = interp1(0:0.025:1,average_fittedRate_r1(:,nc),0:0.01:1);

% fitRange
fitRange = 21:71;

% fun = @(x)x(1)*ones(size(BcdData(fitRange))) + (BcdData(fitRange)./x(3)).^6 * x(2) - initial_rate_r0(fitRange);
% Rate = r_basal + (Bcd/K_a).^6*r ;
% nonlinear least square fitting to the model
x0 = [0.1 0.1];%[1  0.1 0.1]; %[r_basal r K_a]
lb = [0, 0];
ub = [100 100];
options.Algorithm = 'levenberg-marquardt';
%lsqOptions=optimset('Display','none');

x1 = lsqnonlin(@r1_Hill_initial_rate_fit, x0, lb, ub, options, x, Bcd_interp_nc13(fitRange),...
       Runt_interp_nc13_BGsubtracted(fitRange) + 10,initial_rate_r1(fitRange))%
   
   

%% Plot the data and the fit (r1)
hold on
% data (initial fit)
plot(X,Bcd_interp_nc13*100)
plot(X(fitRange),initial_rate_r1(fitRange))
% fitted 
BcdData = Bcd_interp_nc13(fitRange);
RuntData = Runt_interp_nc13_BGsubtracted(fitRange) + 10;

plot(X(fitRange), (x(1)*ones(size(BcdData)) + (BcdData./x(3)).^6  * x(2) + (RuntData./x1(1)).*(BcdData./x(3)).^6*x1(2)) ./...
                (1 + (BcdData./x(3)).^6 + (RuntData./x1(1)).*(BcdData./x(3)).^6*x1(2) + (RuntData./x1(1)) ))


% save the plot



%% Now try to predict the r2 and r3
%% r2 data
nc=2; % NC13;
%errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc))
% Fill in the values in the posterior, up to 80%
average_fittedRate_r2(21:33,nc) = 30; %average_fittedRate_r0(21,nc)

% interpolate the averaged_fittedRate_r0
initial_rate_r2 = interp1(0:0.025:1,average_fittedRate_r1(:,nc),0:0.01:1);

% fitRange
fitRange = 21:71;
%% r2
r2_Prediction = (x(1)*ones(size(BcdData)) + (BcdData./x(3)).^6  * x(2) + ...
    (RuntData./x1(1)).*(BcdData./x(3)).^6*x1(2) + (RuntData./x1(1)).^2.*(BcdData./x(3)).^6*x1(2)) ./...
                (1 + (BcdData./x(3)).^6 + (RuntData./x1(1)).*(BcdData./x(3)).^6 + (RuntData./x1(1)) + (RuntData./x1(1)).^2 )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Calculate the gradient features from the real r0 dataset

% Load the initial rate of RNAP loading from the r0 datasets
%Data for r0 on Data(1,:). Change the N value, number of datasets
NC13_InitialRates = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\AveragedInitialRate_NC13.mat');
NC13_InitialRates_r0 = NC13_InitialRates.Mean13(1,:);

InitialRate_r0_normalized = NC13_InitialRates_r0 ./ max(NC13_InitialRates_r0);
[Max, Min, DynamicRange, BoundaryPosition, Slope] = getGradientFeatures (X,InitialRate_r0_normalized);

Maximuim_r0_data = Max;
DynamicRange_r0_data = DynamicRange;
BoundaryPosition_r0_data = BoundaryPosition;
Slope_r0_data = Slope;

%% scatter plot (hbP2-r0)

hold on
scatter3(Maximum_r0, BoundaryPosition_r0, Slope_r0,10)
plot3(Maximuim_r0_data,BoundaryPosition_r0_data,(-1)*Slope_r0_data,'or')
title('Combination of P_{bound} features')
xlabel('Maximum')
ylabel('Boundary Position')
zlabel('Slope')
set(gca,'fontSize',15)
grid on
view(-40,20)

% Save the plots
saveas(gcf,'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\ModelingEffort-hbP2\hbP2_r0_Prediction_Max_BP_Slope.fig')
saveas(gcf,'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\ModelingEffort-hbP2\hbP2_r0_Prediction_Max_BP_Slope.tif')
%% Dynamic range (r0)
hold on
grid on
scatter3(DynamicRange_r0, BoundaryPosition_r0, Slope_r0,10)
plot3(DynamicRange_r0_data,BoundaryPosition_r0_data,(-1)*Slope_r0_data,'or')
title('Combination of P_{bound} features')
xlabel('Dynamic range')
ylabel('Boundary Position')
zlabel('Slope')
set(gca,'fontSize',20)
view(-40,20)

% Save the plots
saveas(gcf,'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\ModelingEffort-hbP2\hbP2_r0_Prediction_dynamicRange_BP_Slope.fig')
saveas(gcf,'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\ModelingEffort-hbP2\hbP2_r0_Prediction_dynamicRange_BP_Slope.tif')
%% 2D collapse version

hold on
plot(BoundaryPosition_r0, Slope_r0./ DynamicRange_r0,'o','Color',[213,108,85]/255)
plot(BoundaryPosition_r0_data, (-1)*Slope_r0_data./DynamicRange_r0_data,'o','Color',[115,143,193]/255)

title('Combination of P_{bound} features')
xlabel('Boundary Position')
ylabel('Slope/Dynamic Range')
legend('r0-model','Data')
set(gca,'fontSize',20)
%% Actual hbP2-r0 data
% We can use the result here, to comapre with actual data, then constrain
% the parameter sets.
% Then, use the parameter sets to start with, for more complex cases.
Rate_NC13 = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\AveragedInitialRate_NC13.mat');
Rate_NC13_r0 = Rate_NC13.Mean13(1,:);
Rate_NC13_r1 = Rate_NC13.Mean13(2,:);
Rate_NC13_r2 = Rate_NC13.Mean13(3,:);
Rate_NC13_r3 = Rate_NC13.Mean13(4,:);

% Smoothen the Rate_NC13_r0
Rate_NC13_r0_smooth = movmean(Rate_NC13_r0,3);
Rate_NC13_r1_smooth = movmean(Rate_NC13_r1,3);
Rate_NC13_r2_smooth = movmean(Rate_NC13_r2,3);
Rate_NC13_r3_smooth = movmean(Rate_NC13_r3,3);

% Interpolate with 1% binning
Rate_NC13_r0_interp = interp1(0:0.025:1,Rate_NC13_r0_smooth,0:0.01:1);
Rate_NC13_r1_interp = interp1(0:0.025:1,Rate_NC13_r1_smooth,0:0.01:1);
Rate_NC13_r2_interp = interp1(0:0.025:1,Rate_NC13_r2_smooth,0:0.01:1);
Rate_NC13_r3_interp = interp1(0:0.025:1,Rate_NC13_r3_smooth,0:0.01:1);

%% r0
% Normalize the expression
Rate_NC13_r0_interp_normalized = Rate_NC13_r0_interp/ max(Rate_NC13_r0_interp);
[Max, Min, dynamicRange, boundaryPosition, slope] = getGradientFeatures (0:0.01:1,Rate_NC13_r0_interp_normalized);


%% Case1. competition

P = [0.1:0.1:1,2:10];
wAR = [0:0.1:1]; % competition, about binding
wAP = [1 exp(1) exp(2) exp(3) exp(4) exp(5)];% interaction between Bcd and Pol II
wRP = 1;
wAR_prime = 1; % repression on activator's action, espeically in quenching


l = 1;
for i=1:length(P)
    for j=1:length(wAR)
        for k=1:length(wAP)
            for m=1:length(BcdScale)
                for n=1:length(RuntScale)
                    % Let's break down each statistical weight, and think about what's
                    % dominating, to figure out a reasonable parameter range.
                    Bcd_bound = Bcd*BcdScale(m);
                    Runt_bound = Runt*RuntScale(n);
                    Bcd_Runt_bound = Bcd*BcdScale(m).*Runt*RuntScale(n)*wAR(j);
                    Pol_bound = P(i);
                    Bcd_Pol_bound = Bcd*BcdScale(m)*P(i)*wAP(k);
                    Runt_Pol_bound = Runt*RuntScale(n)*P(i)*wRP;
                    everything_bound = Bcd*BcdScale(m).*Runt*RuntScale(n)*...
                                        P(i)*wAR(j)*wAP(k)*wRP*wAR_prime;

                    % Define Partition function, and P_bound states
        %             PartitionFunction_competition(:,i,j,k) = Bcd_bound + Runt_bound + Bcd_Runt_bound +...
        %                                 Pol_bound + Bcd_Pol_bound + Runt_Pol_bound + everything_bound;
                    PartitionFunction = Bcd_bound + Runt_bound + Bcd_Runt_bound +...
                                        Pol_bound + Bcd_Pol_bound + Runt_Pol_bound + everything_bound;
                    RNAP_bound = Pol_bound + Bcd_Pol_bound + Runt_Pol_bound + everything_bound;

                    P_bound_temp = RNAP_bound./PartitionFunction;

                    % Now, let's get the features of gradients, Max, Min, Dynamic
                    % Range, Boundary position, and Slope.
                    [max, min, dynamicRange, boundaryPosition, slope] = getGradientFeatures (X,P_bound_temp);

        %             Maximum(i,j,k) = max;
        %             Minimum(i,j,k) = min;
        %             DynamicRange(i,j,k) = dynamicRange;
        %             BoundaryPosition(i,j,k) = boundaryPosition;
        %             Slope(i,j,k) = slope;
                    Maximum(l) = max;
                    Minimum(l) = min;
                    DynamicRange(l) = dynamicRange;
                    BoundaryPosition(l) = boundaryPosition;
                    Slope(l) = -1*slope; % Making the slope positive, just to compare the amplitude.
                    l=l+1;
                end
            end
        end
    end
end

%% scatter plot (competition)
Maximum_competition = Maximum;
BoundaryPosition_competition = BoundaryPosition;
Slope_competition = Slope;
DynamicRange_competition = DynamicRange;

scatter3(Maximum_competition, BoundaryPosition_competition, Slope_competition,1)
xlabel('Maximum')
ylabel('Boundary Position')
zlabel('Slope')