% Repressibility
% Explanation goes here

%% NC13
Rate_NC13 = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedInitialRate_NC13.mat');
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

%% Plot to check
hold on
plot(0:0.01:1,Rate_NC13_r0_interp)
plot(0:0.01:1,Rate_NC13_r1_interp)
plot(0:0.01:1,Rate_NC13_r2_interp)
plot(0:0.01:1,Rate_NC13_r3_interp)

title('Initial rate of RNAP loading over AP')
xlabel('AP')
ylabel('Initial rate of RNAP loading (AU/min)')
legend('r0','r1','r2','r3')

%% Logistic function fit (r0)
X = 0.2:0.01:0.47;%0:0.01:1; % X axis, EL
Y1 = Rate_NC13_r0_interp(21:48);
X1 = 0.15:0.01:0.6; % fit plotting

% Eqn : Y = L/(1 + A*exp(B*X)) + C;
% MAPPING:  L = b(1), A = b(2), B = b(3), C = b(4);
logistic_fit = @(b,X)  b(1)./(1+b(2).*exp(b(3)*X)) + b(4);
b0 = [250; 5; 1; 70];                                  % Initial Parameter Estimates
B0_13 = lsqcurvefit(logistic_fit, b0, X, Y1)

figure(1)
plot(X, Y1, 'bp')
hold on
plot(X1, logistic_fit(B0_13,X1), '-r')
hold off
grid
xlabel('AP')
ylabel('r0')
legend('Data', 'Logistic Equation Fit', 'Location','NE')

%% Logistic function fit (r1)
X = 0.2:0.01:0.57;%0:0.01:1; % X axis, EL
Y1 = Rate_NC13_r1_interp(21:58);
X1 = 0.15:0.01:0.6; % fit plotting

% Eqn : Y = L/(1 + A*exp(B*X)) + C;
% MAPPING:  L = b(1), A = b(2), B = b(3), C = b(4);
logistic_fit = @(b,X)  b(1)./(1+b(2).*exp(b(3)*X)) + b(4);
b0 = [220; 5; 1; 70];                                  % Initial Parameter Estimates
B1_13 = lsqcurvefit(logistic_fit, b0, X, Y1)

figure(1)
plot(X, Y1, 'bp')
hold on
plot(X1, logistic_fit(B1_13,X1), '-r')
hold off
grid
xlabel('AP')
ylabel('r1')
legend('Data', 'Logistic Equation Fit', 'Location','NE')

%% Logistic function fit (r2)
X = 0.27:0.01:0.47;%0:0.01:1; % X axis, EL
Y1 = Rate_NC13_r2_interp(28:48);
X1 = 0.15:0.01:0.6; % fit plotting

% Eqn : Y = L/(1 + A*exp(B*X)) + C;
% MAPPING:  L = b(1), A = b(2), B = b(3), C = b(4);
logistic_fit = @(b,X)  b(1)./(1+b(2).*exp(b(3)*(X))) + b(4);
b0 = [300; 1; 10; 10];                                  % Initial Parameter Estimates
B2_13 = lsqcurvefit(logistic_fit, b0, X, Y1)

figure(1)
plot(X, Y1, 'bp')
hold on
plot(X1, logistic_fit(B2_13,X1), '-r')
hold off
grid
xlabel('AP')
ylabel('r2')
legend('Data', 'Logistic Equation Fit', 'Location','NE')

%% Logistic function fit (r3)
X = 0.18:0.01:0.30;%0:0.01:1; % X axis, EL
Y1 = Rate_NC13_r3_interp(19:31);
X1 = 0.15:0.01:0.6; % fit plotting

% Eqn : Y = L/(1 + A*exp(B*X)) + C;
% MAPPING:  L = b(1), A = b(2), B = b(3), C = b(4);
logistic_fit = @(b,X)  b(1)./(1+b(2).*exp(b(3)*(X))) + b(4);
b0 = [240; 1; 15; 30 ];                                  % Initial Parameter Estimates
B3_13 = lsqcurvefit(logistic_fit, b0, X, Y1)

figure(1)
plot(X, Y1, 'bp')
hold on
plot(X1, logistic_fit(B3_13,X1), '-r')
plot(0:0.01:1,Rate_NC13_r3_interp,'og')
hold off
grid
xlabel('AP')
ylabel('r3')
legend('Data', 'Logistic Equation Fit','Data-full', 'Location','NE')

%% Get the logistic fitted profile
X = 0:0.01:1;
X1 = 0.15:0.01:0.6;
logistic_fitted_r0_13 = logistic_fit(B0_13,X1);
logistic_fitted_r1_13 = logistic_fit(B1_13,X1);
logistic_fitted_r2_13 = logistic_fit(B2_13,X1);
logistic_fitted_r3_13 = logistic_fit(B3_13,X1);

hold on
plot(X1,logistic_fitted_r0_13)
plot(X1,logistic_fitted_r1_13)
plot(X1,logistic_fitted_r2_13)
plot(X1,logistic_fitted_r3_13)
legend('r0','r1','r2','r3')
title('Logistic fit of Initial rate of RNAP loading')
xlabel('AP (EL)')
ylabel('Logistic fitted Initial rate (AU/min)')

%% Repressibility ( from the logistic fits)
% Let's calculate the repressibility, which is defined as
% Repressibility = 1 - ( Rate( N binding sites) / Rate ( No binding sites))
% Let's not think about the basal rate for now.
X1 = 0.15:0.01:0.6;
Repressibility_r1_13 = 1 - logistic_fitted_r1_13./logistic_fitted_r0_13;
Repressibility_r2_13 = 1 - logistic_fitted_r2_13./logistic_fitted_r0_13;
Repressibility_r3_13 = 1 - logistic_fitted_r3_13./logistic_fitted_r0_13;

Repression_r1_13 = logistic_fitted_r0_13 ./ logistic_fitted_r1_13;
Repression_r2_13 = logistic_fitted_r0_13 ./ logistic_fitted_r2_13;
Repression_r3_13 = logistic_fitted_r0_13 ./ logistic_fitted_r3_13;
% Let's plot for 20-40% of the embryo only (since these fits are reasonable
% for only this window).
AP1 = find(X1==0.2);
AP2 = find(X1==0.39);
AP2 = 26;

figure(1)
hold on 
plot(X1(AP1:AP2),Repressibility_r1_13(AP1:AP2))
plot(X1(AP1:AP2),Repressibility_r2_13(AP1:AP2))
plot(X1(AP1:AP2),Repressibility_r3_13(AP1:AP2))
title('Repressibility : 1 - Rate(N_{R}) / Rate(N_{R} = 0) @ NC 13')
xlabel('AP axis (EL)')
ylabel('Repressibility')
legend('r1','r2','r3','Location','NW')
hold off

figure(2)
hold on 
plot(X1(AP1:AP2),Repression_r1_13(AP1:AP2))
plot(X1(AP1:AP2),Repression_r2_13(AP1:AP2))
plot(X1(AP1:AP2),Repression_r3_13(AP1:AP2))
title('Repression : Rate(N_{R}=0) / Rate(N_{R}) @ NC 13')
xlabel('AP axis (EL)')
ylabel('Repressibility')
legend('r1','r2','r3','Location','NW')
hold off

%% Basal level subtraction
% There seems to be some basal level of expression, that is independent of
% Bcd. I will subtract this value from all data poins for now.
Basal_expression = min(Rate_NC13_r3_interp);

Rate_NC13_r0_interp_BGsubtracted = Rate_NC13_r0_interp - Basal_expression;
Rate_NC13_r1_interp_BGsubtracted = Rate_NC13_r1_interp - Basal_expression;
Rate_NC13_r2_interp_BGsubtracted = Rate_NC13_r2_interp - Basal_expression;
Rate_NC13_r3_interp_BGsubtracted = Rate_NC13_r3_interp - Basal_expression;

%% Repressibility
% Let's calculate the repressibility, which is defined as
% Repressibility = 1 - ( Rate( N binding sites) / Rate ( No binding sites))

Repressibility_r1 = 1 - Rate_NC13_r1_interp_BGsubtracted./Rate_NC13_r0_interp_BGsubtracted;
Repressibility_r2 = 1 - Rate_NC13_r2_interp_BGsubtracted./Rate_NC13_r0_interp_BGsubtracted;
Repressibility_r3 = 1 - Rate_NC13_r3_interp_BGsubtracted./Rate_NC13_r0_interp_BGsubtracted;

hold on 
plot(0:0.01:1,Repressibility_r1)
plot(0:0.01:1,Repressibility_r2)
plot(0:0.01:1,Repressibility_r3)
title('Repressibility : 1 - Rate(N_{R}) / Rate(N_{R} = 0) @ NC 13')
xlabel('AP axis (EL)')
ylabel('Repressibility')
legend('r1','r2','r3')



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

%% Plot to check
hold on
plot(0:0.01:1,Rate_NC14_r0_interp)
plot(0:0.01:1,Rate_NC14_r1_interp)
plot(0:0.01:1,Rate_NC14_r2_interp)
plot(0:0.01:1,Rate_NC14_r3_interp)


%% Logistic function fit (r0)
X = 0.16:0.01:0.47;%0:0.01:1; % X axis, EL
Y1 = Rate_NC14_r0_interp(17:48);
X1 = 0.15:0.01:0.6; % fit plotting

% Eqn : Y = L/(1 + A*exp(B*X)) + C;
% MAPPING:  L = b(1), A = b(2), B = b(3), C = b(4);
logistic_fit = @(b,X)  b(1)./(1+b(2).*exp(b(3)*X)) + b(4);
b0 = [1000; 5; 1; 50];                                  % Initial Parameter Estimates
B0 = lsqcurvefit(logistic_fit, b0, X, Y1)

figure(1)
plot(X, Y1, 'bp')
hold on
plot(X1, logistic_fit(B0,X1), '-r')
hold off
grid
xlabel('AP')
ylabel('r0')
legend('Data', 'Logistic Equation Fit', 'Location','NE')

%% Logistic function fit (r1)
X = 0.2:0.01:0.47;%0:0.01:1; % X axis, EL
Y1 = Rate_NC14_r1_interp(21:48);
X1 = 0.15:0.01:0.6; % fit plotting

% Eqn : Y = L/(1 + A*exp(B*X)) + C;
% MAPPING:  L = b(1), A = b(2), B = b(3), C = b(4);
logistic_fit = @(b,X)  b(1)./(1+b(2).*exp(b(3)*X)) + b(4);
b0 = [220; 5; 1; 70];                                  % Initial Parameter Estimates
B1 = lsqcurvefit(logistic_fit, b0, X, Y1)

figure(1)
plot(X, Y1, 'bp')
hold on
plot(X1, logistic_fit(B1,X1), '-r')
hold off
grid
xlabel('AP')
ylabel('r1')
legend('Data', 'Logistic Equation Fit', 'Location','NE')

%% Logistic function fit (r2)
X = 0.18:0.01:0.4;%0:0.01:1; % X axis, EL
Y1 = Rate_NC14_r2_interp(19:41);
X1 = 0.15:0.01:0.6; % fit plotting

% Eqn : Y = L/(1 + A*exp(B*X)) + C;
% MAPPING:  L = b(1), A = b(2), B = b(3), C = b(4);
logistic_fit = @(b,X)  b(1)./(1+b(2).*exp(b(3)*(X))) + b(4);
b0 = [300; 1; 10; 10];                                  % Initial Parameter Estimates
B2 = lsqcurvefit(logistic_fit, b0, X, Y1)

figure(1)
plot(X, Y1, 'bp')
hold on
plot(X1, logistic_fit(B2,X1), '-r')
hold off
grid
xlabel('AP')
ylabel('r2')
legend('Data', 'Logistic Equation Fit', 'Location','NE')

%% Logistic function fit (r3)
X = 0.16:0.01:0.55;%0:0.01:1; % X axis, EL
Y1 = Rate_NC14_r3_interp(17:56);
X1 = 0.15:0.01:0.6; % fit plotting

% Eqn : Y = L/(1 + A*exp(B*X)) + C;
% MAPPING:  L = b(1), A = b(2), B = b(3), C = b(4);
logistic_fit = @(b,X)  b(1)./(1+b(2).*exp(b(3)*(X))) + b(4);
b0 = [240; 1; 15; 30 ];                                  % Initial Parameter Estimates
B3 = lsqcurvefit(logistic_fit, b0, X, Y1)

figure(1)
plot(X, Y1, 'bp')
hold on
plot(X1, logistic_fit(B3,X1), '-r')
%plot(0:0.01:1,Rate_NC14_r3_interp,'og')
hold off
grid
xlabel('AP')
ylabel('r3')
legend('Data', 'Logistic Equation Fit', 'Location','NE')

%% Get the logistic fitted profile
X = 0:0.01:1;
X1 = 0.15:0.01:0.6;
logistic_fitted_r0 = logistic_fit(B0,X1);
logistic_fitted_r1 = logistic_fit(B1,X1);
logistic_fitted_r2 = logistic_fit(B2,X1);
logistic_fitted_r3 = logistic_fit(B3,X1);

hold on
plot(X1,logistic_fitted_r0)
plot(X1,logistic_fitted_r1)
plot(X1,logistic_fitted_r2)
plot(X1,logistic_fitted_r3)
legend('r0','r1','r2','r3')
title('Logistic fit of Initial rate of RNAP loading')
xlabel('AP (EL)')
ylabel('Logistic fitted Initial rate (AU/min)')


%% Subtract the background
% For now, I will subtract the minimal value from each fit
logistic_fitted_r0 = logistic_fitted_r0 - min(logistic_fitted_r0);
logistic_fitted_r1 = logistic_fitted_r1 - min(logistic_fitted_r1);
logistic_fitted_r2 = logistic_fitted_r2 - min(logistic_fitted_r2);
logistic_fitted_r3 = logistic_fitted_r3 - min(logistic_fitted_r3);
%% Repressibility ( from the logistic fits)
% Let's calculate the repressibility, which is defined as
% Repressibility = 1 - ( Rate( N binding sites) / Rate ( No binding sites))
% Let's not think about the basal rate for now.
X1 = 0.15:0.01:0.6;
Repressibility_r1 = 1 - logistic_fitted_r1./logistic_fitted_r0;
Repressibility_r2 = 1 - logistic_fitted_r2./logistic_fitted_r0;
Repressibility_r3 = 1 - logistic_fitted_r3./logistic_fitted_r0;

% Let's plot for 20-40% of the embryo only (since these fits are reasonable
% for only this window).
AP1 = find(X1==0.15);
AP2 = find(X1==0.6);

hold on 
plot(X1(AP1:AP2),Repressibility_r1(AP1:AP2))
plot(X1(AP1:AP2),Repressibility_r2(AP1:AP2))
plot(X1(AP1:AP2),Repressibility_r3(AP1:AP2))
title('Repressibility : 1 - Rate(N_{R}) / Rate(N_{R} = 0) @ NC 14')
xlabel('AP axis (EL)')
ylabel('Repressibility')
legend('r1','r2','r3','Location','NW')
ylim([0 1])
%% Basal level (with raw data)
% Think about somewhat basal level of Txn, we should subtract this from all
% of the rates.
% We should do a better job on estimating the basal level, for example,
% finding the basal level with BcdE1, for each construct. Ideally, they
% should be similar.
Basal_Level = min(Rate_NC14_r3);

Rate_NC14_r0 = Rate_NC14_r0 - Basal_Level;
Rate_NC14_r1 = Rate_NC14_r1 - Basal_Level;
Rate_NC14_r2 = Rate_NC14_r2 - Basal_Level;
Rate_NC14_r3 = Rate_NC14_r3 - Basal_Level;


% Smoothen the Rate_NC14_r0
Rate_NC14_r0_smooth = movmean(Rate_NC14_r0,3);
Rate_NC14_r1_smooth = movmean(Rate_NC14_r1,3);
Rate_NC14_r2_smooth = movmean(Rate_NC14_r2,3);
Rate_NC14_r3_smooth = movmean(Rate_NC14_r3,3);

% Interpolate with 1% binning
Rate_NC14_r0_interp = interp1(0:0.025:1,Rate_NC14_r0_smooth,0:0.01:1);
Rate_NC14_r1_interp = interp1(0:0.025:1,Rate_NC14_r1_smooth,0:0.01:1);
Rate_NC14_r2_interp = interp1(0:0.025:1,Rate_NC14_r2_smooth,0:0.01:1);
Rate_NC14_r3_interp = interp1(0:0.025:1,Rate_NC14_r3_smooth,0:0.01:1);

%hold on

%% Repressibility
% Side thing, but let's calculate the repressibility, which is defined as
% Repressibility = 1 - ( Rate( N binding sites) / Rate ( No binding sites))

Repressibility_r1 = 1 - Rate_NC14_r1_interp./Rate_NC14_r0_interp;
Repressibility_r2 = 1 - Rate_NC14_r2_interp./Rate_NC14_r0_interp;
Repressibility_r3 = 1 - Rate_NC14_r3_interp./Rate_NC14_r0_interp;

hold on 
plot(0:0.01:1,ones(1,101))
plot(0:0.01:1,Repressibility_r1)
plot(0:0.01:1,Repressibility_r2)
plot(0:0.01:1,Repressibility_r3)
title('Repressibility : 1 - Rate(N_{R}) / Rate(N_{R} = 0) @ NC 14')
xlabel('AP axis (EL)')
ylabel('Repressibility')
legend('Saturation','r1','r2','r3')