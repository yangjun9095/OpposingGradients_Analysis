%% Modeling for one activator and one repressor binding sites (1A1R)
% December, 2018
% Yang Joon Kim

%% Partition function

%% P_{bound} for the Competition 
A = [0:0.1:1, 2:1:10,20:10:100,200:100:1000];
%A = [0:1:10^3];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [0.01 0.1 0.5 1 10 100 1000];

P = 1;
wAR = 0; % competition
wAP = 10;
wRP = 1;

PartitionFunction = nan(length(A),length(R));
RNAPBoundStates = nan(length(A),length(R));
P_bound = nan(length(A),length(R));

for i=1:length(A)
    for j=1:length(R)
        PartitionFunction(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
        RNAPBoundStates(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
    end
end
P_bound = RNAPBoundStates ./ PartitionFunction;

surf(R,A,P_bound)

title('P_{bound} for competitive repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Constraints (Overleaf eqn.12)
% There are some constraints about parameters, which is
% P*wAP >> R/K_R >> P
% Thus, I will start with set of parameters as below.

A = [0:0.1:1, 2:1:10,20:10:100];%,200:100:1000];
%A = [0:1:10^3];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [1 10 100];

P = 1;
wAR = 0; % competition
wAP = 100;
wRP = 1;

PartitionFunction = nan(length(A),length(R));
RNAPBoundStates = nan(length(A),length(R));
P_bound = nan(length(A),length(R));
for i=1:length(A)
    for j=1:length(R)
        PartitionFunction(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
        RNAPBoundStates(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
    end
end
P_bound = RNAPBoundStates ./ PartitionFunction;

surf(R,A,P_bound)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('P_{bound} for competitive repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Quenching

A = [0:0.1:1, 2:1:10,20:10:100];%,200:100:1000];
%A = [0:1:10^3];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [1 10 100];

P = 1;
wAR = 1; 
wARPrime = exp(-5);
wAP = 100;
wRP = 1;

PartitionFunction = nan(length(A),length(R));
RNAPBoundStates = nan(length(A),length(R));
P_bound = nan(length(A),length(R));
for i=1:length(A)
    for j=1:length(R)
        PartitionFunction(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP*wARPrime;
        RNAPBoundStates(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP*wARPrime;
    end
end
P_bound = RNAPBoundStates ./ PartitionFunction;

surf(R,A,P_bound)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('P_{bound} for quenching repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Direct repression

A = [0:0.1:1, 2:1:10,20:10:100];%,200:100:1000];
%A = [0:1:10^3];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [1 10 100];

P = 1;
wAR = 1; 
%wARPrime = exp(-5);
wAP = 100;
wRP = 0;

PartitionFunction = nan(length(A),length(R));
RNAPBoundStates = nan(length(A),length(R));
P_bound = nan(length(A),length(R));
for i=1:length(A)
    for j=1:length(R)
        PartitionFunction(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP*wARPrime;
        RNAPBoundStates(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
    end
end
P_bound = RNAPBoundStates ./ PartitionFunction;

surf(R,A,P_bound)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('P_{bound} for direct repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Shut down

A = [0:0.1:1, 2:1:10,20:10:100];%,200:100:1000];
%A = [0:1:10^3];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [1 10 100];

P = 1;
wAR = 0; 
%wARPrime = exp(-5);
wAP = 100;
wRP = 0;

PartitionFunction = nan(length(A),length(R));
RNAPBoundStates = nan(length(A),length(R));
P_bound = nan(length(A),length(R));
for i=1:length(A)
    for j=1:length(R)
        PartitionFunction(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
        RNAPBoundStates(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
    end
end
P_bound = RNAPBoundStates ./ PartitionFunction;

surf(R,A,P_bound)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('P_{bound} for shut down')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Let's constrain more, by fixing the repressor concentration.
A = [0:0.01:0.1,0.2:0.1:1, 2:1:10,20:10:100];%,200:100:1000];
%A = [0:1:10^3];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = 10;

% 1. Competition
P = 1;
wAR = 0; 
%wARPrime = exp(-5);
wAP = exp(5);
wRP = 1;

for i=1:length(A)
    PartitionFunction(i) = 1 + A(i) + R + A(i).*R*wAR + P + A(i)*P*wAP + R*P*wRP + A(i).*R*P*wAR*wAP*wRP;
    RNAPBoundStates(i) = P + A(i)*P*wAP + R*P*wRP + A(i).*R*P*wAR*wAP*wRP;
end
P_bound_competition = RNAPBoundStates ./ PartitionFunction;

% 2. Quenching
P = 1;
wAR = 1; 
wARPrime = exp(-5);
wAP = exp(5);
wRP = 1;

for i=1:length(A)
    PartitionFunction(i) = 1 + A(i) + R + A(i).*R*wAR + P + A(i)*P*wAP + R*P*wRP + A(i).*R*P*wAR*wAP*wRP*wARPrime;
    RNAPBoundStates(i) = P + A(i)*P*wAP + R*P*wRP + A(i).*R*P*wAR*wAP*wRP*wARPrime;
end
P_bound_quenching = RNAPBoundStates ./ PartitionFunction;

% 3. Direct repression
P = 1;
wAR = 1; 
%wARPrime = exp(-5);
wAP = exp(5);
wRP = 0;

for i=1:length(A)
    PartitionFunction(i) = 1 + A(i) + R + A(i).*R*wAR + P + A(i)*P*wAP + R*P*wRP + A(i).*R*P*wAR*wAP*wRP;
    RNAPBoundStates(i) = P + A(i)*P*wAP + R*P*wRP + A(i).*R*P*wAR*wAP*wRP;
end
P_bound_direct = RNAPBoundStates ./ PartitionFunction;

% 4. shut down 
P = 1;
wAR = 0; 
%wARPrime = exp(-5);
wAP = exp(5);
wRP = 0;

for i=1:length(A)
    PartitionFunction(i) = 1 + A(i) + R + A(i).*R*wAR + P + A(i)*P*wAP + R*P*wRP + A(i).*R*P*wAR*wAP*wRP;
    RNAPBoundStates(i) = P + A(i)*P*wAP + R*P*wRP + A(i).*R*P*wAR*wAP*wRP;
end
P_bound_shutdown = RNAPBoundStates ./ PartitionFunction;

% Plot
hold on
plot(A,P_bound_competition)
plot(A,P_bound_quenching)
plot(A,P_bound_direct)
plot(A,P_bound_shutdown)
set(gca, 'XScale', 'log');
title('P_{bound} for Activator concentrations')
xlabel('Activator (A/K_{A}) (AU)')
ylabel('P_{bound}')
legend('competition','quenching','direct','shutdown','Location','SouthEast')
standardizeFigure(gca,legend,[])

%% Let's constrain more, by fixing the activator concentration.
R = [0:0.01:0.1,0.2:0.1:1, 2:1:10,20:10:100];%,200:100:1000];
%A = [0:1:10^3];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
A = 1;

% 1. Competition
P = 1;
wAR = 0; 
%wARPrime = exp(-5);
wAP = exp(5);
wRP = 1;

for i=1:length(R)
    PartitionFunction(i) = 1 + A + R(i) + A.*R(i)*wAR + P + A*P*wAP + R(i)*P*wRP + A.*R(i)*P*wAR*wAP*wRP;
    RNAPBoundStates(i) = P + A*P*wAP + R(i)*P*wRP + A.*R(i)*P*wAR*wAP*wRP;
end
P_bound_competition = RNAPBoundStates ./ PartitionFunction;

% 2. Quenching
P = 1;
wAR = 1; 
wARPrime = exp(-5);
wAP = exp(5);
wRP = 1;

for i=1:length(R)
    PartitionFunction(i) = 1 + A + R(i) + A.*R(i)*wAR + P + A*P*wAP + R(i)*P*wRP + A*R(i)*P*wAR*wAP*wRP*wARPrime;
    RNAPBoundStates(i) = P + A*P*wAP + R(i)*P*wRP + A*R(i)*P*wAR*wAP*wRP*wARPrime;
end
P_bound_quenching = RNAPBoundStates ./ PartitionFunction;

% 3. Direct repression
P = 1;
wAR = 1; 
%wARPrime = exp(-5);
wAP = exp(5);
wRP = 0;

for i=1:length(R)
    PartitionFunction(i) = 1 + A + R(i) + A*R(i)*wAR + P + A*P*wAP + R(i)*P*wRP + A*R(i)*P*wAR*wAP*wRP;
    RNAPBoundStates(i) = P + A*P*wAP + R(i)*P*wRP + A*R(i)*P*wAR*wAP*wRP;
end
P_bound_direct = RNAPBoundStates ./ PartitionFunction;

% 4. shut down 
P = 1;
wAR = 0; 
%wARPrime = exp(-5);
wAP = exp(5);
wRP = 0;

for i=1:length(R)
    PartitionFunction(i) = 1 + A + R(i) + A*R(i)*wAR + P + A*P*wAP + R(i)*P*wRP + A*R(i)*P*wAR*wAP*wRP;
    RNAPBoundStates(i) = P + A*P*wAP + R(i)*P*wRP + A*R(i)*P*wAR*wAP*wRP;
end
P_bound_shutdown = RNAPBoundStates ./ PartitionFunction;

% Plot
hold on
plot(R,P_bound_competition)
plot(R,P_bound_quenching)
plot(R,P_bound_direct)
plot(R,P_bound_shutdown)
ylim([0 1])
set(gca, 'XScale', 'log');
title('P_{bound} for Repressor concentrations')
xlabel('Repressor (R/K_{R}) (AU)')
ylabel('P_{bound}')
legend('competition','quenching','direct','shutdown','Location','SouthWest')
standardizeFigure(gca,legend,[])