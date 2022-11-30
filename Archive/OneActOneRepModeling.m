% Theoretical Modeling, start from the simplest case of one Bcd and one
% Runt binding sites, then explicitly thinking about RNAP recruitment to
% the promoter.

% Variables : 
% A/K_A (A) = activator conc. divided by the dissociation constant of
% the activator and its binding site.
% R/K_A (R) = repressor conc. divided by the dissociation constant of
% the repressor and its binding site.
% wAR : interaction term between the activator and repressor
% wAP : interaction term between the activator and RNAP
% wRP : interaction term between the repressor and RNAP
% P : P/K_P = RNAP conc. divided by the dissociation constant of the RNAP
% and promoter.
%% Initial definition of variables (or parameters)P
clear all
X = linspace(0,1,41); 
BcdScale = 10;
RuntScale = 10;

% for Bcd, I start with this, with the length constant of 0.2, which is the value from Feng Liu's paper.
A = BcdScale * exp(-5*X); 

R = RuntScale * (-1* exp((X-0.5).^2) + 1.3)/0.3;

P = 5;

wAR = 1;
wAP = 5;
wRP = 0;

% Partition function : 
% Pbound

PartitionFunction = 1 + A + R + A.*R*wAR + P + A*P*wAP + R*P*wRP + A.*R*P*wAR*wAP*wRP;

RNAPBoundStates = P + A*P*wAP + R*P*wRP + A.*R*P*wAR*wAP*wRP;

P_bound = RNAPBoundStates ./ PartitionFunction;

% Calculating the P_bound for zero repressor site, and only one activator
% site.
%R=0;
%PartitionFunction = 1 + A + R + A.*R*wAR + P + A*P*wAP + R*P*wRP + A.*R*P*wAR*wAP*wRP;

%RNAPBoundStates = P + A*P*wAP + R*P*wRP + A.*R*P*wAR*wAP*wRP;

hold on
plot(X,PartitionFunction)
plot(X,RNAPBoundStates)
plot(X,RNAPBoundStates./PartitionFunction * 100)
legend('PartitionFunction','RNAPBound','P_Bound')
%% plot Rate_1 and Rate_0 for wR
hold on
    plot(X,A)
    plot(X,R*0.5)
    plot(X,P_bound*10)

% for i=1:length(wR)
%     plot(X,Ratio(:,:,i))
% end
title('P_bound over AP')
xlabel('AP')
ylabel('P_bound (AU)')
legend('Activator','Repressor','P_{bound}')
%,num2str(wR(1)),num2str(wR(2)),num2str(wR(3)))
standardizeFigure(gca,legend,[])



%% Explore broader parameter space (Nov.2018)
% the idea is trying to make somewhat general predictions for all possible
% activator and repressor conc. combinations. Then, try to narrow down
% realistic parameter space (sub-space)

%% P_{bound} for the Competition 
A = [1:10,20:10:100,200:100:1000];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [1 10 100 1000];

P = 1;
wAR = 0; % competition
wAP = 10;
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

title('P_{bound} for competitive repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Get the fold-change (competitive)
A = [1:10,20:10:100,200:100:1000];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [0 0 0 0];

P = 1;
wAR = 0;
wAP = 10;
wRP = 0;

PartitionFunction_R0 = nan(length(A),length(R));
RNAPBoundStates_R0 = nan(length(A),length(R));
P_bound_R0 = nan(length(A),length(R));

for i=1:length(A)
    for j=1:length(R)
        PartitionFunction_R0(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
        RNAPBoundStates_R0(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
    end
end
P_bound_R0 = RNAPBoundStates_R0 ./ PartitionFunction_R0;

Fold_change_compete_R0 = P_bound./ P_bound_R0;
%surf(R,A,Fold_change_nocompete)
%title('P_{bound} for non-competitive repression')
%xlabel('Runt(AU)')
%ylabel('Bcd(AU)')
%zlabel('Fold-change')

semilogx(A,Fold_change_compete_R0(:,1))
hold on
for i=2:length(R)
    semilogx(A,Fold_change_compete_R0(:,i))
end
title('Fold-change(competitive)')
xlabel('Bcd(AU)')
ylabel('P_{bound}')
legend('R=1','R=10','R=100','R=1000','Location','SouthEast')
standardizeFigure(gca,legend,[])
%% Explore the A, R parameter space - non-competitive repression
A = [1:10,20:10:100,200:100:1000];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [1 10 100 1000];

P = 1;
wAR = 1;
wAP = 10;
wRP = 0;

PartitionFunction_nocompete = nan(length(A),length(R));
RNAPBoundStates_nocompete = nan(length(A),length(R));
P_bound_nocompete = nan(length(A),length(R));

for i=1:length(A)
    for j=1:length(R)
        PartitionFunction_nocompete(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
        RNAPBoundStates_nocompete(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
    end
end
P_bound_nocompete = RNAPBoundStates_nocompete ./ PartitionFunction_nocompete;

surf(R,A,P_bound_nocompete)

title('P_{bound} for non-competitive repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')


%% Get the fold-change (non-competitive)
A = [1:10,20:10:100,200:100:1000];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [0 0 0 0];

P = 1;
wAR = 1;
wAP = 10;
wRP = 0;

PartitionFunction_nocompete_R0 = nan(length(A),length(R));
RNAPBoundStates_nocompete_R0 = nan(length(A),length(R));
P_bound_nocompete_R0 = nan(length(A),length(R));

for i=1:length(A)
    for j=1:length(R)
        PartitionFunction_nocompete_R0(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
        RNAPBoundStates_nocompete_R0(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
    end
end
P_bound_nocompete_R0 = RNAPBoundStates_nocompete_R0 ./ PartitionFunction_nocompete_R0;

Fold_change_nocompete = P_bound_nocompete./ P_bound_nocompete_R0;
%surf(R,A,Fold_change_nocompete)
%title('P_{bound} for non-competitive repression')
%xlabel('Runt(AU)')
%ylabel('Bcd(AU)')
%zlabel('Fold-change')

semilogx(A,Fold_change_nocompete(:,1))
hold on
for i=2:length(R)
    semilogx(A,Fold_change_nocompete(:,i))
end
title('Fold-change(non-competitive)')
xlabel('Bcd(AU)')
ylabel('P_{bound}')
legend('R=1','R=10','R=100','R=1000','Location','SouthEast')
standardizeFigure(gca,legend,[])

%% Get the ratio between competition vs non-competition

Ratio_P_bound = P_bound ./ P_bound_nocompete;

%h = axes;
%set(h,'xscale','log')
%set(h,'yscale','log')

surf(R,A,Ratio_P_bound)

title('Ratio of P_{bound} between competitive vs non-competitive')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% figure out how to constrain parameters like wAR, etc.
% A and R actual values would be on the order or 10^1 (guess...)
% Let's start with A=10, and R=10, then try to find reasonable P, wAR, wAP,
% etc.

% A=10;
% R=10;
% P = 1;
% wAR = 1;
% wAP = 10;
% wRP = 0;
% 
% PartitionFunction = 1 + A + R + A.*R*wAR + P + A*P*wAP + R*P*wRP + A.*R*P*wAR*wAP*wRP;
% 
% RNAPBoundStates = P + A*P*wAP + R*P*wRP + A.*R*P*wAR*wAP*wRP;
% 
% P_bound = RNAPBoundStates ./ PartitionFunction;

%% Plot cross-sections (competition)
h = axes;
set(h,'xscale','log')
hold on

for i=1:length(R)
    semilogx(A,P_bound(:,i))
    pause
end
ylim([0 1])
title('P_{bound} for different Repressor/Activator ratio')
xlabel('Activator concentration(AU)')
ylabel('P_{bound}')
legend('1/10','1','10','100')

standardizeFigure(gca,legend,[])

%% Plot cross-sections (non-competitive)
h = axes;
set(h,'xscale','log')
hold on

for i=1:length(R)
    semilogx(A,P_bound_nocompete(:,i))
    pause
end
ylim([0 1])
title('P_{bound} for different Repressor/Activator ratio')
xlabel('Activator concentration(AU)')
ylabel('P_{bound}')
legend('1/10','1','10','100')

standardizeFigure(gca,legend,[])
%% Realistic parameter space (of Bcd and Runt)
% Here, I will make pairs of Bcd and Runt 
% For Bcd, I start with this, with the length constant of 0.2, which is the value from Feng Liu's paper.

X =0:0.025:1; % AP axis, 2.5%
BcdScale = 100;
RuntScale = 100;

Bcd = BcdScale * exp(-5*X); 
Runt = RuntScale * (-1* exp((X-0.5).^2) + 1.3)/0.3;

% Parameter sets
P = 1;
wAR = 1;
wAP = 10;
wRP = 0;

% P_{bound}
PartitionFunction_AP = 1 + Bcd + Runt + Bcd.*Runt*wAR + P + Bcd*P*wAP + Runt*P*wRP + Bcd.*Runt*P*wAR*wAP*wRP;
RNAPBoundStates_AP = P + Bcd*P*wAP + Runt*P*wRP + Bcd.*Runt*P*wAR*wAP*wRP;

P_bound_AP = RNAPBoundStates_AP./PartitionFunction_AP;

plot(X,P_bound_AP)



%% 
hold on
for i=1:length(X)
    scatter3(Runt(i),Bcd(i),P_bound_AP(i))
end

%% 2018-11-09(Friday)
% Revisit the idea of generating polarizing plots for different scenarios
% of repression
% Plug in 5kbT for energy terms

%% Competition
A = [1:10,20:10:100,200:100:1000];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [10:10:1000];
%R = [1 10 100 1000];
eps = 5;
P = 0.1;
wAR = 0;
wAP = exp(eps);
wRP = 1;

PartitionFunction_compete = nan(length(A),length(R));
RNAPBoundStates_compete = nan(length(A),length(R));
P_bound_compete = nan(length(A),length(R));

for i=1:length(A)
    for j=1:length(R)
        PartitionFunction_compete(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
        RNAPBoundStates_compete(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP;
    end
end
P_bound_compete = RNAPBoundStates_compete ./ PartitionFunction_compete;

surf(R,A,P_bound_compete)

title('P_{bound} for competitive repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Quenching
A = [1:10,20:10:100,200:100:1000];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [10:10:1000];
eps = 3;
P = 0.1;
wAR = 1;
wAP = exp(5);
wRP = 1;
DelwAP = exp(-eps);

PartitionFunction_quench = nan(length(A),length(R));
RNAPBoundStates_quench = nan(length(A),length(R));
P_bound_quench = nan(length(A),length(R));

for i=1:length(A)
    for j=1:length(R)
        PartitionFunction_quench(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP*DelwAP;
        RNAPBoundStates_quench(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP*DelwAP;
    end
end
P_bound_quench = RNAPBoundStates_quench ./ PartitionFunction_quench;

surf(R,A,P_bound_quench)

title('P_{bound} for quenching repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Direct repression
A = [1:10,20:10:100,200:100:1000];
%A = [1 10 100 1000 10000];
%R = [1:10,20:10:100,200:100:1000];
R = [10:10:1000];

P = 0.1;
wAR = 1;
wAP = exp(5);
wRP = exp(-5);


PartitionFunction_quench = nan(length(A),length(R));
RNAPBoundStates_quench = nan(length(A),length(R));
P_bound_quench = nan(length(A),length(R));

for i=1:length(A)
    for j=1:length(R)
        PartitionFunction_direct(i,j) = 1 + A(i) + R(j) + A(i).*R(j)*wAR + P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP*DelwAP;
        RNAPBoundStates_direct(i,j) = P + A(i)*P*wAP + R(j)*P*wRP + A(i).*R(j)*P*wAR*wAP*wRP*DelwAP;
    end
end
P_bound_direct = RNAPBoundStates_direct ./ PartitionFunction_direct;

surf(R,A,P_bound_direct)

title('P_{bound} for direct repression')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('P_{bound}')

%% Ratio
Ratio_compete_quench = P_bound_compete ./ P_bound_quench;

surf(R,A,Ratio_compete_quench)

title('Ratio of P_{bound} for competition vs quenching')
xlabel('Runt(AU)')
ylabel('Bcd(AU)')
zlabel('Ratio of P_{bound}')