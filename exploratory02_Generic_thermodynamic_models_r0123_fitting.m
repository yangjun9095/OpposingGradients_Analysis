function exploratory02_effective_thermodynamic_models_r0123_fitting

% I'll start with the initial rate of RNAP loading data from r0,1,2,3
% constructs.

%% Step1. load the input/output datasets for fitting
%% Bcd and Runt
%Bcd = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Bcd-Averaged.mat')
Bcd = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\Bicoid.mat');

Bcd_averaged_nc13 = nanmean(Bcd.pchbcd(1:120,:));
Bcd_averaged_nc14 = nanmean(Bcd.pchbcd(140:end,:));

% Spatial interpolation 
Bcd_interp_nc13 = interp1(0:0.025:1,Bcd_averaged_nc13,0:0.01:1);

plot(0:0.01:1,Bcd_interp_nc13)
BcdScale = 60; % This shoud change...

% (2) Runt gradient, this should come from actual data
% For now, let's take an average over time, for nc13 and nc14.
% The dataset below is from averaging 3~4 embryos, with 1min of frame rate.
% Also, APbin of 2.5%, thus I need to smooth, and interpolate
Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Runt-1min-200Hz-Male-Averaged.mat');
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
RuntBG = 175; %min(Runt_interp_nc13);%150; % Here, I just chose this since I didn't have a value for 20% or more anterior...
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
standardizeFigure_YJK(gca,legend,[])
%saveas()
% Only think about 20% to 60% of the embryo, and scale for Bcd and Runt.
EL = 0:0.01:1;
X = 0.2:0.01:0.6;

BcdScale = 60;
Bcd = Bcd_interp_nc13*BcdScale;
Bcd = Bcd(21:61); % getting 0.2~0.6 of EL

RuntScale = 2;
Runt = Runt_interp_nc13_BGsubtracted*RuntScale;
Runt = Runt(21:61); % getting 0.2~0.6 of EL

%% NC14
% Load datasets (NC14)
Rate_NC14 = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedInitialRate_NC14.mat');
Rate_NC14_r0 = Rate_NC14.Mean14(1,:);
Rate_NC14_r1 = Rate_NC14.Mean14(2,:);
Rate_NC14_r2 = Rate_NC14.Mean14(3,:);
Rate_NC14_r3 = Rate_NC14.Mean14(4,:);

% Smoothen the Rate_NC13_r0
Rate_NC14_r0_smooth = movmean(Rate_NC14_r0,1);
Rate_NC14_r1_smooth = movmean(Rate_NC14_r1,1);
Rate_NC14_r2_smooth = movmean(Rate_NC14_r2,1);
Rate_NC14_r3_smooth = movmean(Rate_NC14_r3,1);

% Interpolate with 1% binning
Rate_NC14_r0_interp = interp1(0:0.025:1,Rate_NC14_r0_smooth,0:0.01:1);
Rate_NC14_r1_interp = interp1(0:0.025:1,Rate_NC14_r1_smooth,0:0.01:1);
Rate_NC14_r2_interp = interp1(0:0.025:1,Rate_NC14_r2_smooth,0:0.01:1);
Rate_NC14_r3_interp = interp1(0:0.025:1,Rate_NC14_r3_smooth,0:0.01:1);

end