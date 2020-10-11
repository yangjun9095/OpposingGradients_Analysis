function DataForModeling_deprecated
% Here, I will just think about from 20~80% (or 20~60%) of the embryo, since 
% 1) all boundary features are determined in this regime
% 2) we don't have to think about terminal system
% 3) this is usually the region we're imaging

X = 0.2:0.01:0.6; % AP axis

% (1) Bcd gradient, note that the length constant is 0.2
% Bcd0 = 60; % 60nM in the most anterior from Gregor 2007a
% Bcd = Bcd0*exp(-5*X);
% Kb = 6; % Bcd dissociation constant. This is for scaling, should change

% or we can just use actual Bcd data
% The Bcd data below is interpolated(temporally) using pchip, I recall that
% it's from Jonathan.
Bcd = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Bicoid.mat')

Bcd_averaged_nc13 = nanmean(Bcd.pchbcd(1:120,:));
% Spatial interpolation 
Bcd_interp_nc13 = interp1(0:0.025:1,Bcd_averaged_nc13,0:0.01:1);

plot(0:0.01:1,Bcd_interp_nc13)
BcdScale = 60; % This shoud change...

% (2) Runt gradient, this should come from actual data
% For now, let's take an average over time, for nc13 and nc14.
% The dataset below is from averaging 3~4 embryos, with 1min of frame rate.
% Also, APbin of 2.5%, thus I need to smooth, and interpolate
Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-1min-200Hz-Female-Averaged.mat');
%Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Runt-1min-200Hz-Female-Averaged.mat');

% Runt profile (averaged over nc13)
timeRange_nc13 = Runt.nc13:Runt.nc14-3;

Runt_averaged_nc13 = nanmean(Runt.MeanVectorAP(timeRange_nc13,:));
Runt_averaged_nc13_error = sqrt(nansum(Runt.SDVectorAP(timeRange_nc13,:).^2)./length(Runt.SDVectorAP(timeRange_nc13,:)));

% Runt profile (averaged over the beginning of 20 minutes of nc14)
timeRange_nc14 = Runt.nc14:Runt.nc14+20; % for the first 20 minutes
Runt_averaged_nc14 = nanmean(Runt.MeanVectorAP(timeRange_nc14,:));
Runt_averaged_nc14_error = sqrt(nansum(Runt.SDVectorAP(timeRange_nc14,:).^2)./length(Runt.SDVectorAP(timeRange_nc14,:)));

% Quick check for the Runt gradient
hold on
errorbar(0:0.025:1,Runt_averaged_nc13,Runt_averaged_nc13_error)
errorbar(0:0.025:1,Runt_averaged_nc14,Runt_averaged_nc14_error)
title('Time-averaged Runt nuclear fluorescence over AP')
xlabel('AP axis (EL)')
ylabel('Runt nuclear fluorescence (AU)')
legend('nc13','nc14 (20min)')
hold off

% Interpolation and smoothing (or the other way)
Runt_smoothed_nc13 = movmean(Runt_averaged_nc13,3);
Runt_interp_nc13 = interp1(0:0.025:1,Runt_smoothed_nc13,0:0.01:1);
%Runt_error_interp_nc13 = interp1(0:0.025:1,Runt_averaged_nc13_error,0:0.01:1);

% Now, subtracat the background. This should also be revisited once we get
% a better sense about the background fluorescence, like free eGFP. For
% now, I'm just subtracting the minimum value, assuming that Runt is
% almost zero at the very anterior (or posterior)
RuntBG = 150; %min(Runt_interp_nc13);%150; % Here, I just chose this since I didn't have a value for 20% or more anterior...
Runt_interp_nc13_BGsubtracted = Runt_interp_nc13 - RuntBG ;

% Check both Bcd and Runt inputs
figure(2)
hold on
plot(0:0.01:1,Bcd_interp_nc13*BcdScale)
plot(0:0.01:1,Runt_interp_nc13_BGsubtracted)
%errorbar(0:0.01:1,Runt_interp_nc13,Runt_error_interp_nc13)

title('Bcd and Runt over AP (averaged over nc13)')
xlabel('AP axis (EL)')
ylabel('Nuclear fluorescence (AU)')
legend('Bcd','Runt')
hold off

% Only think about 20% to 60% of the embryo, and scale for Bcd and Runt.
EL = 0:0.01:1;
X = 0.2:0.01:0.6;

BcdScale = 20;
Bcd = Bcd_interp_nc13*BcdScale;
Bcd = Bcd(21:61); % getting 0.2~0.6 of EL

RuntScale = 1;
Runt = Runt_interp_nc13_BGsubtracted*RuntScale;
Runt = Runt(21:61); % getting 0.2~0.6 of EL




%% 
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
end