% Modeling for the hb P2 + r3
% Dec.2018
% The goal is to theorize the case of hb P2 + three Runt binding sites.

%% r0 case as a standard to compare with
A = [0:0.01:0.1, 0.2:0.01:1, 2:1:10, 20:10:100, 200:100:1000];
wA = exp(0.5);
%R = 1; 
P = 1;
wAP = exp(0.5);

for i = 1:length(A)
    % Only Activators bound
    A_only(i) = (1 + A(i)*wA).^6 - 1./wA;

    % Activators and RNAP
    AP(i) = P*((1 + A(i)*wA*wAP).^6 - 1./wA);
end
PartitionFunction_r0 = 1 + A_only + P + AP;

Bound_r0 = P + AP ;

P_bound_r0 = Bound_r0 ./ PartitionFunction_r0;


semilogx(A,P_bound_r0,'-')
%% Competition
R = [0:0.01:0.1, 0.2:0.01:1, 2:1:10, 20:10:100, 200:100:1000];
wA = exp(0.5);
A = 1;
P = 1;
wAP = exp(0.5);

for i = 1:length(R)
    % Either only Activators or repressors bound
    A_only(i) = (1 + A*wA).^6 - 1./wA;
    R_only(i) = (1 + R(i)).^3;

    % Activators with 1,2,3 repressors bound
    AR1(i) = R(i)*((1+A*wA).^5 - 1./wA);
    AR2(i) = R(i).^2*((1+A*wA).^4 - 1./wA);
    AR3(i) = R(i).^3*((1+A*wA).^3 - 1./wA);

    % Activators and RNAP
    AP(i) = P*((1 + A*wA*wAP).^6 - 1./wA);

    % Repressors and RNAP(including no repressor, only RNAP)
    RP(i) = P*(1 + R(i)).^3;

    % N Rep + Activators + RNAP
    APR1(i) = P.*R(i)*3*((1 + A*wA*wAP).^5 - 1./wA);
    APR2(i) = P.*R(i).^2*3*((1 + A*wA*wAP).^4 - 1./wA);
    APR3(i) = P.*R(i).^3*1*((1 + A*wA*wAP).^3 - 1./wA);
end

PartitionFunction_competition = 1 + A_only + R_only + AR1 + AR2 + AR3 + ...
            AP + RP + APR1 + APR2 + APR3;

Bound = AP + RP + APR1 + APR2 + APR3;

P_bound_competition = Bound ./ PartitionFunction_competition;


semilogx(R,P_bound_competition,'-')

%% Quenching
R = [0:0.01:0.1, 0.2:0.01:1, 2:1:10, 20:10:100, 200:100:1000];
wA = exp(0.5);
A = 1;
P = 1;
wAP = exp(0.5);
wAR = exp(-0.7);

for i = 1:length(A)
    % Either only Activators or repressors bound
    A_only(i) = (1 + A*wA).^6 - 1./wA;
    R_only(i) = (1 + R(i)).^3 -1;

    % Activators with 1,2,3 repressors bound
    AR(i) = ((1 + R(i)).^3 - 1) .* ((1 + A*wA).^6 - 1./wA);

    % Activators and RNAP
    AP(i) = P*((1 + A*wA*wAP).^6 - 1./wA);

    % Repressors and RNAP (including no repressor, only RNAP)
    RP(i) = P*(1 + R(i)).^3;

    % N Rep + Activators + RNAP
    APR(i) = P*((1 + R(i)*wAR).^3 - 1).*((1 + A*wA).^6 - 1./wA);
    
end

PartitionFunction_quenching = 1 + A_only + R_only + AR + ...
            AP + RP + APR;

Bound = AP + RP + APR;

P_bound_quenching = Bound ./ PartitionFunction_quenching;

semilogx(R,P_bound_quenching)

%% Comparison
figure(1)
hold on
semilogx(A,P_bound_r0,'-')
semilogx(A,P_bound_competition)
semilogx(A,P_bound_quenching)
set(gca, 'XScale', 'log');

title('P_{bound}')
xlabel('Activator concentration (A/K_{A})')
ylabel('P_{bound}')
legend('r0','Competition', 'Quenching')

%%
Ratio = P_bound_competition./P_bound_quenching;
figure(2)
semilogx(A,Ratio)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Change the Repressor conc. range

%% r0 case as a standard to compare with
A = [0:0.01:0.1, 0.2:0.01:1, 2:1:10, 20:10:100, 200:100:1000];
wA = exp(0.5);
%R = 1; 
P = 1;
wAP = exp(0.5);

for i = 1:length(A)
    % Only Activators bound
    A_only(i) = (1 + A(i)*wA).^6 - 1./wA;

    % Activators and RNAP
    AP(i) = P*((1 + A(i)*wA*wAP).^6 - 1./wA);
end
PartitionFunction_r0 = 1 + A_only + P + AP;

Bound_r0 = P + AP ;

P_bound_r0 = Bound_r0 ./ PartitionFunction_r0;


semilogx(A,P_bound_r0,'-')
%% Competition
A = [0:0.01:0.1, 0.2:0.01:1, 2:1:10, 20:10:100, 200:100:1000];
wA = exp(0.5);
R = 1;
P = 1;
wAP = exp(0.5);

for i = 1:length(A)
    % Either only Activators or repressors bound
    A_only(i) = (1 + A(i)*wA).^6 - 1./wA;
    R_only(i) = (1 + R).^3;

    % Activators with 1,2,3 repressors bound
    AR1(i) = R*((1+A(i)*wA).^5 - 1./wA);
    AR2(i) = R.^2*((1+A(i)*wA).^4 - 1./wA);
    AR3(i) = R.^3*((1+A(i)*wA).^3 - 1./wA);

    % Activators and RNAP
    AP(i) = P*((1 + A(i)*wA*wAP).^6 - 1./wA);

    % Repressors and RNAP(including no repressor, only RNAP)
    RP(i) = P*(1 + R).^3;

    % N Rep + Activators + RNAP
    APR1(i) = P.*R*3*((1 + A(i)*wA*wAP).^5 - 1./wA);
    APR2(i) = P.*R.^2*3*((1 + A(i)*wA*wAP).^4 - 1./wA);
    APR3(i) = P.*R.^3*1*((1 + A(i)*wA*wAP).^3 - 1./wA);
end

PartitionFunction_competition = 1 + A_only + R_only + AR1 + AR2 + AR3 + ...
            AP + RP + APR1 + APR2 + APR3;

Bound = AP + RP + APR1 + APR2 + APR3;

P_bound_competition = Bound ./ PartitionFunction_competition;


semilogx(A,P_bound_competition,'-')

%% Quenching
A = [0:0.01:0.1, 0.2:0.01:1, 2:1:10, 20:10:100, 200:100:1000];
wA = exp(0.5);
R = 1;
P = 1;
wAP = exp(0.5);
wAR = exp(-0.7);

for i = 1:length(A)
    % Either only Activators or repressors bound
    A_only(i) = (1 + A(i)*wA).^6 - 1./wA;
    R_only(i) = (1 + R).^3 -1;

    % Activators with 1,2,3 repressors bound
    AR(i) = ((1 + R).^3 - 1) .* ((1 + A(i)*wA).^6 - 1./wA);

    % Activators and RNAP
    AP(i) = P*((1 + A(i)*wA*wAP).^6 - 1./wA);

    % Repressors and RNAP (including no repressor, only RNAP)
    RP(i) = P*(1 + R).^3;

    % N Rep + Activators + RNAP
    APR(i) = P*((1 + R*wAR).^3 - 1).*((1 + A(i)*wA).^6 - 1./wA);
    
end

PartitionFunction_quenching = 1 + A_only + R_only + AR + ...
            AP + RP + APR;

Bound = AP + RP + APR;

P_bound_quenching = Bound ./ PartitionFunction_quenching;

semilogx(A,P_bound_quenching)

%% Comparison
figure(1)
hold on
semilogx(A,P_bound_r0,'-')
semilogx(A,P_bound_competition)
semilogx(A,P_bound_quenching)
set(gca, 'XScale', 'log');

title('P_{bound}')
xlabel('Activator concentration (A/K_{A})')
ylabel('P_{bound}')
legend('r0','Competition', 'Quenching')

%%
Ratio = P_bound_competition./P_bound_quenching;
figure(2)
semilogx(A,Ratio)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's think about more realistic case, like Bcd concentration gradient in
% the embryo

X = [0:0.025:1]; % Embryo AP axis
A = 10*exp(-5*X);

%% r0 case as a standard to compare with
X = [0:0.025:1]; % Embryo AP axis
A = 10*exp(-5*X);
wA = exp(0.5);
%R = 1; 
P = 1;
wAP = exp(0.5);

for i = 1:length(A)
    % Only Activators bound
    A_only(i) = (1 + A(i)*wA).^6 - 1./wA;

    % Activators and RNAP
    AP(i) = P*((1 + A(i)*wA*wAP).^6 - 1./wA);
end
PartitionFunction_r0 = 1 + A_only + P + AP;

Bound_r0 = P + AP ;

P_bound_r0 = Bound_r0 ./ PartitionFunction_r0;


plot(X,P_bound_r0,'-')
%% Competition
R = [0:0.01:0.1, 0.2:0.01:1, 2:1:10, 20:10:100, 200:100:1000];
wA = exp(0.5);
A = 1;
P = 1;
wAP = exp(0.5);

for i = 1:length(R)
    % Either only Activators or repressors bound
    A_only(i) = (1 + A*wA).^6 - 1./wA;
    R_only(i) = (1 + R(i)).^3;

    % Activators with 1,2,3 repressors bound
    AR1(i) = R(i)*((1+A*wA).^5 - 1./wA);
    AR2(i) = R(i).^2*((1+A*wA).^4 - 1./wA);
    AR3(i) = R(i).^3*((1+A*wA).^3 - 1./wA);

    % Activators and RNAP
    AP(i) = P*((1 + A*wA*wAP).^6 - 1./wA);

    % Repressors and RNAP(including no repressor, only RNAP)
    RP(i) = P*(1 + R(i)).^3;

    % N Rep + Activators + RNAP
    APR1(i) = P.*R(i)*3*((1 + A*wA*wAP).^5 - 1./wA);
    APR2(i) = P.*R(i).^2*3*((1 + A*wA*wAP).^4 - 1./wA);
    APR3(i) = P.*R(i).^3*1*((1 + A*wA*wAP).^3 - 1./wA);
end

PartitionFunction_competition = 1 + A_only + R_only + AR1 + AR2 + AR3 + ...
            AP + RP + APR1 + APR2 + APR3;

Bound = AP + RP + APR1 + APR2 + APR3;

P_bound_competition = Bound ./ PartitionFunction_competition;


semilogx(R,P_bound_competition,'-')

%% Quenching
R = [0:0.01:0.1, 0.2:0.01:1, 2:1:10, 20:10:100, 200:100:1000];
wA = exp(0.5);
A = 1;
P = 1;
wAP = exp(0.5);
wAR = exp(-0.7);

for i = 1:length(A)
    % Either only Activators or repressors bound
    A_only(i) = (1 + A*wA).^6 - 1./wA;
    R_only(i) = (1 + R(i)).^3 -1;

    % Activators with 1,2,3 repressors bound
    AR(i) = ((1 + R(i)).^3 - 1) .* ((1 + A*wA).^6 - 1./wA);

    % Activators and RNAP
    AP(i) = P*((1 + A*wA*wAP).^6 - 1./wA);

    % Repressors and RNAP (including no repressor, only RNAP)
    RP(i) = P*(1 + R(i)).^3;

    % N Rep + Activators + RNAP
    APR(i) = P*((1 + R(i)*wAR).^3 - 1).*((1 + A*wA).^6 - 1./wA);
    
end

PartitionFunction_quenching = 1 + A_only + R_only + AR + ...
            AP + RP + APR;

Bound = AP + RP + APR;

P_bound_quenching = Bound ./ PartitionFunction_quenching;

semilogx(R,P_bound_quenching)