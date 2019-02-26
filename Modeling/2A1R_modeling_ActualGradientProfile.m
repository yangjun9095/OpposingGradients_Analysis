% Modeling for 2A1R given the actual gradient profiles of Bcd and Runt

% Last updated : Feb.2019

%% Step1. Think about realistic morphogen gradients

% Here, I will just think about from 20~80% (or 20~60%) of the embryo, since 
% 1) all boundary features are determined in this regime
% 2) we don't have to think about terminal system
% 3) this is usually the region we're imaging

X = 0.2:0.01:0.7; % AP axis

% (1) Bcd gradient, note that the length constant is 0.2
Bcd0 = 60; % 60nM in the most anterior from Gregor 2007a
Bcd = Bcd0*exp(-5*X);
Kb = 6; % Bcd dissociation constant. This is for scaling, should change

% or we can just use actual Bcd data
% The Bcd data below is interpolated(temporally) using pchip, I recall that
% it's from Jonathan.
Bcd = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Bicoid.mat')

Bcd_averaged_nc13 = nanmean(Bcd.pchbcd(1:120,:));
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

% Only think about 20% to 60% of the embryo, and scale for Bcd and Runt.
EL = 0:0.01:1;
X = 0.2:0.01:0.6;

BcdScale = 60;
Bcd = Bcd_interp_nc13*BcdScale;
Bcd = Bcd(21:61); % getting 0.2~0.6 of EL

RuntScale = 2;
Runt = Runt_interp_nc13_BGsubtracted*RuntScale;
Runt = Runt(21:61); % getting 0.2~0.6 of EL

%% 
%% Step4. Let's explore a much broader parameter space, using getGradientFeatures.m
% including BcdScale and RuntScale.
% For now, I will fix the Bcd scaling and Runt scaling, and just tune the
% P, wAR, wAP, wRP, wAR_prime, to see whether different scenarios of
% repression would actually lead to different predictions.

% Scaling of Bcd and Runt
% This scaling represents Kd for each protein
BcdScale = [10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 1 10 10^2 10^3 10^4 10^5];
RuntScale = [10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 1 10 10^2 10^3 10^4 10^5];

%% 1. Competition
P = [0.1:0.1:1,2:10];
wAR = [0:0.1:1]; % competition, about binding
wA = exp([0 1 2 3 4 5]); % Cooperativity between Bcd molecules (binding)
wAP = exp([0 1 2 3 4 5]);%exp(5);
wRP = 1;
wAR_prime = 1; % repression on activator's action, espeically in quenching
wRAA = 1; % repressing the Bcd's cooperativity

%PartitionFunction_competition = nan(length(X),length(P),length(wAR),length(wAP));

l = 1;
for i=1:length(P)
    for j=1:length(wAR)
        for k=1:length(wAP)
            for m=1:length(BcdScale)
                for n=1:length(RuntScale)
                    for p=1:length(wA)
                        % Let's break down each statistical weight, and think about what's
                        % dominating, to figure out a reasonable parameter range.
                        Nothing_bound = 1;
                        Bcd_bound = 2*Bcd*BcdScale(m);
                        Two_Bcd_bound = (Bcd*BcdScale(m)).^2*wA(p);
                        Runt_bound = Runt*RuntScale(n);
                        Bcd_Runt_bound = 2*Bcd*BcdScale(m).*Runt*RuntScale(n)*wAR(j);
                        Two_Bcd_Runt_bound = (Bcd*BcdScale(m)).^2*wA(p)*wRAA.*Runt*RuntScale(n)*(wAR(j)^2);
                        
                        Pol_bound = P(i);
                        Bcd_Pol_bound = 2*Bcd*BcdScale(m)*P(i)*wAP(k);
                        Two_Bcd_Pol_bound = (Bcd*BcdScale(m)).^2*wA(p)*P(i)*(wAP(k)).^2;
                        Runt_Pol_bound = Runt*RuntScale(n)*P(i)*wRP;
                        Bcd_Runt_Pol_bound = 2*Bcd*BcdScale(m).*Runt*RuntScale(n)*...
                                            P(i)*wAR(j)*wAP(k)*wRP*wAR_prime;
                        everything_bound = (Bcd*BcdScale(m)).^2*wA(p)*wRAA.*Runt*RuntScale(n)*...
                            P(i)*(wAR(j)).^2*(wAP(k))^2*wRP*(wAR_prime)^2;
                                        

                        % Define Partition function, and P_bound states
            %             PartitionFunction_competition(:,i,j,k) = Bcd_bound + Runt_bound + Bcd_Runt_bound +...
            %                                 Pol_bound + Bcd_Pol_bound + Runt_Pol_bound + everything_bound;
                        PartitionFunction = Nothing_bound + Bcd_bound + Two_Bcd_bound + ... 
                                            Runt_bound + Bcd_Runt_bound + Two_Bcd_Runt_bound + ...
                                            Pol_bound + Bcd_Pol_bound + Two_Bcd_Pol_bound + ...
                                            Runt_Pol_bound + Bcd_Runt_Pol_bound + everything_bound;
                        RNAP_bound =        Pol_bound + Bcd_Pol_bound + Two_Bcd_Pol_bound + ...
                                            Runt_Pol_bound + Bcd_Runt_Pol_bound + everything_bound;

                        P_bound_temp = RNAP_bound./PartitionFunction;

                        % Now, let's get the features of gradients, Max, Min, Dynamic
                        % Range, Boundary position, and Slope.
                        [Max, Min, dynamicRange, boundaryPosition, slope] = getGradientFeatures (X,P_bound_temp);

            %             Maximum(i,j,k) = max;
            %             Minimum(i,j,k) = min;
            %             DynamicRange(i,j,k) = dynamicRange;
            %             BoundaryPosition(i,j,k) = boundaryPosition;
            %             Slope(i,j,k) = slope;
                        Maximum(l) = Max;
                        Minimum(l) = Min;
                        DynamicRange(l) = dynamicRange;
                        BoundaryPosition(l) = boundaryPosition;
                        Slope(l) = -1*slope; % Making the slope positive, just to compare the amplitude.
                        l=l+1;
                    end
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
                            
%% 1-2. Anti-cooperativity
P = [0.1:0.1:1,2:10];
wAR = 1; % competition, about binding
wA = exp([0 1 2 3 4 5]); % Cooperativity between Bcd molecules (binding)
wAP = exp([0 1 2 3 4 5]);%exp(5);
wRP = 1;
wAR_prime = 1; % repression on activator's action, espeically in quenching
wRAA = exp(-1*[0 1 2 3 4 5]); % repressing the Bcd's cooperativity

%PartitionFunction_competition = nan(length(X),length(P),length(wAR),length(wAP));

l = 1;
for i=1:length(P)
    for j=1:length(wRAA)
        for k=1:length(wAP)
            for m=1:length(BcdScale)
                for n=1:length(RuntScale)
                    for p=1:length(wA)
                        % Let's break down each statistical weight, and think about what's
                        % dominating, to figure out a reasonable parameter range.
                        Nothing_bound = 1;
                        Bcd_bound = 2*Bcd*BcdScale(m);
                        Two_Bcd_bound = (Bcd*BcdScale(m)).^2*wA(p);
                        Runt_bound = Runt*RuntScale(n);
                        Bcd_Runt_bound = 2*Bcd*BcdScale(m).*Runt*RuntScale(n)*wAR;
                        Two_Bcd_Runt_bound = (Bcd*BcdScale(m)).^2*wA(p)*wRAA(j).*Runt*RuntScale(n)*(wAR^2);
                        
                        Pol_bound = P(i);
                        Bcd_Pol_bound = 2*Bcd*BcdScale(m)*P(i)*wAP(k);
                        Two_Bcd_Pol_bound = (Bcd*BcdScale(m)).^2*wA(p)*P(i)*(wAP(k)).^2;
                        Runt_Pol_bound = Runt*RuntScale(n)*P(i)*wRP;
                        Bcd_Runt_Pol_bound = 2*Bcd*BcdScale(m).*Runt*RuntScale(n)*...
                                            P(i)*wAR*wAP(k)*wRP*wAR_prime;
                        everything_bound = (Bcd*BcdScale(m)).^2*wA(p)*wRAA(j).*Runt*RuntScale(n)*...
                            P(i)*(wAR).^2*(wAP(k))^2*wRP*(wAR_prime)^2;
                                        

                        % Define Partition function, and P_bound states
            %             PartitionFunction_competition(:,i,j,k) = Bcd_bound + Runt_bound + Bcd_Runt_bound +...
            %                                 Pol_bound + Bcd_Pol_bound + Runt_Pol_bound + everything_bound;
                        PartitionFunction = Nothing_bound + Bcd_bound + Two_Bcd_bound + ... 
                                            Runt_bound + Bcd_Runt_bound + Two_Bcd_Runt_bound + ...
                                            Pol_bound + Bcd_Pol_bound + Two_Bcd_Pol_bound + ...
                                            Runt_Pol_bound + Bcd_Runt_Pol_bound + everything_bound;
                        RNAP_bound =        Pol_bound + Bcd_Pol_bound + Two_Bcd_Pol_bound + ...
                                            Runt_Pol_bound + Bcd_Runt_Pol_bound + everything_bound;

                        P_bound_temp = RNAP_bound./PartitionFunction;

                        % Now, let's get the features of gradients, Max, Min, Dynamic
                        % Range, Boundary position, and Slope.
                        [Max, Min, dynamicRange, boundaryPosition, slope] = getGradientFeatures (X,P_bound_temp);

            %             Maximum(i,j,k) = max;
            %             Minimum(i,j,k) = min;
            %             DynamicRange(i,j,k) = dynamicRange;
            %             BoundaryPosition(i,j,k) = boundaryPosition;
            %             Slope(i,j,k) = slope;
                        Maximum(l) = Max;
                        Minimum(l) = Min;
                        DynamicRange(l) = dynamicRange;
                        BoundaryPosition(l) = boundaryPosition;
                        Slope(l) = -1*slope; % Making the slope positive, just to compare the amplitude.
                        l=l+1;
                    end
                end
            end
        end
    end
end

%% scatter plot (competition)
Maximum_anti_coop = Maximum;
BoundaryPosition_anti_coop = BoundaryPosition;
Slope_anti_coop = Slope;
DynamicRange_anti_coop = DynamicRange;

scatter3(Maximum_anti_coop, BoundaryPosition_anti_coop, Slope_anti_coop,1)
xlabel('Maximum')
ylabel('Boundary Position')
zlabel('Slope')
%% 2. Quenching
P = [0.1:0.1:1,2:10];
wAR = 1; % competition, about binding
wA = exp([0 1 2 3 4 5]); % Cooperativity between Bcd molecules (binding)
wAP = exp([0 1 2 3 4 5]);%exp(5);
wRP = 1;
wAR_prime = exp(-1*[0 1 2 3 4 5]); % repression on activator's action, espeically in quenching
wRAA = 1; %repressing the Bcd's cooperativity

%PartitionFunction_competition = nan(length(X),length(P),length(wAR),length(wAP));

l = 1;
for i=1:length(P)
    for j=1:length(wAR_prime)
        for k=1:length(wAP)
            for m=1:length(BcdScale)
                for n=1:length(RuntScale)
                    for p=1:length(wA)
                        % Let's break down each statistical weight, and think about what's
                        % dominating, to figure out a reasonable parameter range.
                        Nothing_bound = 1;
                        Bcd_bound = 2*Bcd*BcdScale(m);
                        Two_Bcd_bound = (Bcd*BcdScale(m)).^2*wA(p);
                        Runt_bound = Runt*RuntScale(n);
                        Bcd_Runt_bound = 2*Bcd*BcdScale(m).*Runt*RuntScale(n)*wAR;
                        Two_Bcd_Runt_bound = (Bcd*BcdScale(m)).^2*wA(p)*wRAA.*Runt*RuntScale(n)*(wAR^2);
                        
                        Pol_bound = P(i);
                        Bcd_Pol_bound = 2*Bcd*BcdScale(m)*P(i)*wAP(k);
                        Two_Bcd_Pol_bound = (Bcd*BcdScale(m)).^2*wA(p)*P(i)*(wAP(k)).^2;
                        Runt_Pol_bound = Runt*RuntScale(n)*P(i)*wRP;
                        Bcd_Runt_Pol_bound = 2*Bcd*BcdScale(m).*Runt*RuntScale(n)*...
                                            P(i)*wAR*wAP(k)*wRP*wAR_prime(j);
                        everything_bound = (Bcd*BcdScale(m)).^2*wA(p)*wRAA.*Runt*RuntScale(n)*...
                            P(i)*(wAR).^2*(wAP(k))^2*wRP*(wAR_prime(j))^2;
                                        

                        % Define Partition function, and P_bound states
            %             PartitionFunction_competition(:,i,j,k) = Bcd_bound + Runt_bound + Bcd_Runt_bound +...
            %                                 Pol_bound + Bcd_Pol_bound + Runt_Pol_bound + everything_bound;
                        PartitionFunction = Nothing_bound + Bcd_bound + Two_Bcd_bound + ... 
                                            Runt_bound + Bcd_Runt_bound + Two_Bcd_Runt_bound + ...
                                            Pol_bound + Bcd_Pol_bound + Two_Bcd_Pol_bound + ...
                                            Runt_Pol_bound + Bcd_Runt_Pol_bound + everything_bound;
                        RNAP_bound =        Pol_bound + Bcd_Pol_bound + Two_Bcd_Pol_bound + ...
                                            Runt_Pol_bound + Bcd_Runt_Pol_bound + everything_bound;

                        P_bound_temp = RNAP_bound./PartitionFunction;

                        % Now, let's get the features of gradients, Max, Min, Dynamic
                        % Range, Boundary position, and Slope.
                        [Max, Min, dynamicRange, boundaryPosition, slope] = getGradientFeatures (X,P_bound_temp);

            %             Maximum(i,j,k) = max;
            %             Minimum(i,j,k) = min;
            %             DynamicRange(i,j,k) = dynamicRange;
            %             BoundaryPosition(i,j,k) = boundaryPosition;
            %             Slope(i,j,k) = slope;
                        Maximum(l) = Max;
                        Minimum(l) = Min;
                        DynamicRange(l) = dynamicRange;
                        BoundaryPosition(l) = boundaryPosition;
                        Slope(l) = -1*slope; % Making the slope positive, just to compare the amplitude.
                        l=l+1;
                    end
                end
            end
        end
    end
end

%% Define Max,Boundary Position, and Slope for quenching
Maximum_quenching = Maximum;
BoundaryPosition_quenching = BoundaryPosition;
Slope_quenching = Slope;
DynamicRange_quenching = DynamicRange;

scatter3(Maximum_quenching, BoundaryPosition_quenching, Slope_quenching,1)
xlabel('Maximum')
ylabel('Boundary Position')
zlabel('Slope')
%% 3. Direct repression
P = [0.1:0.1:1,2:10];
wAR = 1; % competition, about binding
wA = exp([0 1 2 3 4 5]); % Cooperativity between Bcd molecules (binding)
wAP = exp([0 1 2 3 4 5]);%exp(5);
wRP = exp(-1*[0 1 2 3 4 5]);
wAR_prime = 1; % repression on activator's action, espeically in quenching
wRAA = 1; %repressing the Bcd's cooperativity

%PartitionFunction_competition = nan(length(X),length(P),length(wAR),length(wAP));

l = 1;
for i=1:length(P)
    for j=1:length(wRP)
        for k=1:length(wAP)
            for m=1:length(BcdScale)
                for n=1:length(RuntScale)
                    for p=1:length(wA)
                        % Let's break down each statistical weight, and think about what's
                        % dominating, to figure out a reasonable parameter range.
                        Nothing_bound = 1;
                        Bcd_bound = 2*Bcd*BcdScale(m);
                        Two_Bcd_bound = (Bcd*BcdScale(m)).^2*wA(p);
                        Runt_bound = Runt*RuntScale(n);
                        Bcd_Runt_bound = 2*Bcd*BcdScale(m).*Runt*RuntScale(n)*wAR;
                        Two_Bcd_Runt_bound = (Bcd*BcdScale(m)).^2*wA(p)*wRAA.*Runt*RuntScale(n)*(wAR^2);
                        
                        Pol_bound = P(i);
                        Bcd_Pol_bound = 2*Bcd*BcdScale(m)*P(i)*wAP(k);
                        Two_Bcd_Pol_bound = (Bcd*BcdScale(m)).^2*wA(p)*P(i)*(wAP(k)).^2;
                        Runt_Pol_bound = Runt*RuntScale(n)*P(i)*wRP(j);
                        Bcd_Runt_Pol_bound = 2*Bcd*BcdScale(m).*Runt*RuntScale(n)*...
                                            P(i)*wAR*wAP(k)*wRP(j)*wAR_prime;
                        everything_bound = (Bcd*BcdScale(m)).^2*wA(p)*wRAA.*Runt*RuntScale(n)*...
                            P(i)*(wAR).^2*(wAP(k))^2*wRP(j)*(wAR_prime)^2;
                                        

                        % Define Partition function, and P_bound states
            %             PartitionFunction_competition(:,i,j,k) = Bcd_bound + Runt_bound + Bcd_Runt_bound +...
            %                                 Pol_bound + Bcd_Pol_bound + Runt_Pol_bound + everything_bound;
                        PartitionFunction = Nothing_bound + Bcd_bound + Two_Bcd_bound + ... 
                                            Runt_bound + Bcd_Runt_bound + Two_Bcd_Runt_bound + ...
                                            Pol_bound + Bcd_Pol_bound + Two_Bcd_Pol_bound + ...
                                            Runt_Pol_bound + Bcd_Runt_Pol_bound + everything_bound;
                        RNAP_bound =        Pol_bound + Bcd_Pol_bound + Two_Bcd_Pol_bound + ...
                                            Runt_Pol_bound + Bcd_Runt_Pol_bound + everything_bound;

                        P_bound_temp = RNAP_bound./PartitionFunction;

                        % Now, let's get the features of gradients, Max, Min, Dynamic
                        % Range, Boundary position, and Slope.
                        [Max, Min, dynamicRange, boundaryPosition, slope] = getGradientFeatures (X,P_bound_temp);

            %             Maximum(i,j,k) = max;
            %             Minimum(i,j,k) = min;
            %             DynamicRange(i,j,k) = dynamicRange;
            %             BoundaryPosition(i,j,k) = boundaryPosition;
            %             Slope(i,j,k) = slope;
                        Maximum(l) = Max;
                        Minimum(l) = Min;
                        DynamicRange(l) = dynamicRange;
                        BoundaryPosition(l) = boundaryPosition;
                        Slope(l) = -1*slope; % Making the slope positive, just to compare the amplitude.
                        l=l+1;
                    end
                end
            end
        end
    end
end

%% Define Max,Boundary Position, and Slope for direct repression
Maximum_direct = Maximum;
BoundaryPosition_direct = BoundaryPosition;
Slope_direct = Slope;
DynamicRange_direct = DynamicRange;
%% Plot all cases - Max, BP, and Slope
hold on
grid on
title('Combination of P_{bound} features')
xlabel('Maximum')
ylabel('Boundary Position')
zlabel('Slope')
dotSize = 5;
scatter3(Maximum_competition, BoundaryPosition_competition, Slope_competition,dotSize,[213,108,85]/255)%red
scatter3(Maximum_anti_coop, BoundaryPosition_anti_coop, Slope_anti_coop,dotSize,[126,157,144]/255)%dark green
scatter3(Maximum_quenching, BoundaryPosition_quenching, Slope_quenching,dotSize,[234,194,100]/255)%yellow
scatter3(Maximum_direct, BoundaryPosition_direct, Slope_direct,dotSize,[108,188,233]/255)% cyan
view(-60,30)
legend('Competition','Anti-cooperativity','Quenching','Direct')
set(gca,'fontSize',20)

%% Plot all cases - Dynamic range, BP, and Slope
hold on
grid on
title('Combination of P_{bound} features')
xlabel('Dynamic range')
ylabel('Boundary Position')
zlabel('Slope')
dotSize = 5;
scatter3(DynamicRange_competition, BoundaryPosition_competition, Slope_competition,dotSize,[213,108,85]/255)%red
scatter3(DynamicRange_anti_coop, BoundaryPosition_anti_coop, Slope_anti_coop,dotSize,[126,157,144]/255)%dark green
scatter3(DynamicRange_quenching, BoundaryPosition_quenching, Slope_quenching,dotSize,[234,194,100]/255)%yellow
scatter3(DynamicRange_direct, BoundaryPosition_direct, Slope_direct,dotSize,[108,188,233]/255)% cyan
view(-60,30)
legend('Competition','Anti-cooperativity','Quenching','Direct')
set(gca,'fontSize',20)

%% Plot 2D, collased one
hold on
plot(BoundaryPosition_competition, Slope_competition./DynamicRange_competition)
plot(BoundaryPosition_anti_coop, Slope_anti_coop./DynamicRange_anti_coop)
plot(BoundaryPosition_quenching, Slope_quenching./DynamicRange_quenching)
plot(BoundaryPosition_direct, Slope_direct./DynamicRange_direct)