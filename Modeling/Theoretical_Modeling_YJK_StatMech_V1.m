% This script is for theoretical modeling for the Opposing gradient
% project. 
% Author : Yang Joon Kim (yjkim90@berkeley.edu) and Paul Jeammet
% Last Updated : Oct/2018

% I will write down some states and weights to explain Bcd's cooperative binding,
% also considering Runt's repressive action.
% In future, I can think about other models (like different modes of
% repression.
% The goal is to predict the rate of transcrition (or RNAP loading) over
% time, over all AP bins.

% For starters, I will start with several assumptions, which needs to be checked
% if they are fair.
% Assumptions :
% 1. There are 6 Bicoid binding sites on hb P2 enhancer
% 2. Runt's interaction with activator or PIC is unclear. Thus there could
% be many different scenarios, such as short-range repressor, long range
% repressor, etc. But, here, I will start with the quenching scenario,
% where Runt is not affecting the Bcd binding, but only for the Bcd's
% activation.
% 3. 


%% Data Import

Bcd=load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\Bicoid.mat');
Runt=load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\Runt.mat');

NewBcdFluo=Bcd.pchbcd;
NewRuntFluo=Runt.pchrntCut;
%New time array, with a point every 10s
NewBcdTime=(0:length(NewBcdFluo)-1)/6; 
APbin_Bcd=(1:41);

%% r0 modeling to find the set of parameters

% Define the parameters
% wB : cooperativity between Bicoid molecules
% Bcd : Bicoid concentration at a given time and space (AP bin)
% Kd : dissociation constant of Bicoid to its binding site in the enhancer,
% we assume that it's same for either A3 or X1
% P : [P]/Kp - RNAP concentration
% wBP : cooperativity of Bicoid molecule and RNAP molecule
P = 0.5;
wB = 1.5;
wBP = 1.5;
Bcd = NewBcdFluo;
Kd = 0.5;

% First, let's think about P_bound for a fixed set of parameters above.
% Then, we can think of partition function, and P_bound for each given
% APbin, and given time point. (Thus, the dimension of the matrix, Z (or
% P_bound) should be the same as NewBcdFluo)

Z_r0 = (1+P)*(1-1/wB) + (1/wB)*(1+Bcd*wB/Kd).^6 + (P/wB)*(1+Bcd*wB*wBP/Kd).^6 ;

RNAP_bound_states_weights_r0 = P*(1-1/wB) + (P/wB)*(1+Bcd*wB*wBP/Kd).^6;

P_bound_r0 = RNAP_bound_states_weights_r0./Z_r0;

plot(0:0.025:1,P_bound_r0)

%% Pick some parameter sets that looks reasonable with the data
% Let's focus on nc13, and use the initial rate of RNAL loading in nc13
% along AP, to narrow down the parameters : 


% Load the r0 datatset

%% r1 case modeling
% Define the additional parameters for the repressor (Runt)
% binding/activity first.
% wBPR : cooperativity between Bcd,Runt, and RNAP
wBPR = 0.75;
% KR : dissociation constant of Runt to its binding site
KR = 0.5;
% Runt : Runt concentration
Runt = NewRuntFluo

Z_r1 = Z_r0 + Runt/KR*(1-1/wB) + Runt/KR.*(1+Bcd/Kd*wB).^6 + P*(wBPR/wBP)*(Runt/KR).*((1+Bcd/Kd*wB*wBP).^6);

RNAP_bound_states_weights_r1 = P*(1-1/wB) + (P/wB).*(1+Bcd*wB*wBP/Kd).^6 + P*(wBPR/wBP)*(Runt/KR).*((1+Bcd/Kd*wB*wBP).^6);

P_bound_r1 = RNAP_bound_states_weights_r1./Z_r1;

plot(0:0.025:1,P_bound_r1)

%% Let's plot them over time (at one AP bin)
AP = 9;

figure(3)
hold on
plot((1:length(P_bound_r0))*0.16,P_bound_r0(:,AP))
plot((1:length(P_bound_r1))*0.16,P_bound_r1(:,AP))

title('Predicted RNAP loading rate for r0, r1')
xlabel('Time (min)')
ylabel('Predicted RNAP loading rate (AU)')
legend('r0','r1')
standardizeFigure(gca,legend,[])
%%
numStates = 2;

weight = zeros(length(NewBcdTime),length(APbin_Bcd),numStates);

weight(:,:,1) = ones(length(NewBcdTime),length(APbin_Bcd));
weight(:,:,2) = (NewBcdFluo*3).^6;


% Rates for R0
% These values are from Eric Smith paper which used MS2.V1, thus could change.
r_basal = 0; 
r_Max = 9; 

rate(1) = r_basal;
rate(2) = r_basal + r_Max;

% Assume that the Bcd and Runt has the same frame rate, as well as same
% size of the matrices. This trimming of the dataset can be easily done separately.
Rate = zeros(length(NewBcdTime),length(APbin_Bcd));
sum_weights = zeros(length(NewBcdTime),length(APbin_Bcd));

for i=1:numStates
    Rate = Rate + (rate(i)*squeeze(weight(:,:,i))); % weighted sum of rates
    sum_weights = sum_weights + squeeze(weight(:,:,i)); % sum of all weights
end
% Normalize the Rate with the partition function(sum of all weights)
Rate = Rate ./ sum_weights;



% plot 

r0= load('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient/r0');



figure(1)
t_Plot=200

X=linspace(0,100,41);


plot(X,NewBcdFluo(t_Plot,:)*8,X,NewRuntFluo(t_Plot,:)*8,X,Rate(t_Plot,:),X,r0.MeanVectorAP(round(t_Plot/4),:)/60)
legend('Bicoid','Runt','Predictive Rate','R0 data')
title('plot of input and output (simulated and mesured) data')
xlabel('AP position')
ylabel('Rate')

% figure(2)
% plot((1:length(r0.ElapsedTime)),r0.MeanVectorAP(:,11))
%%




%% Let's start with the next simplest case, r1
% Write down all states and corresponding weights
numStates = 4;

% interaction between Bcd and Runt : w (0=<w=<1, for now)
% For now, we will assume that activators and repressors can bind
% independently, and repressors act upon PIC or promoter to repress the
% transcription initiation.
w = 1; % This means there's no cooperativity

weight = zeros(length(NewBcdTime),length(APbin_Bcd),numStates);

weight(:,:,1) = ones(length(NewBcdTime),length(APbin_Bcd));
weight(:,:,2) = NewBcdFluo.^6;
weight(:,:,3) = NewRuntFluo;
weight(:,:,4) = NewBcdFluo.^6.*NewRuntFluo*w;

% Rates for R1
% These values are from Eric Smith paper which used MS2.V1, thus could change.
r_basal = 3; 
r_Max = 9; 
r_R = 3; % This is a naive guess by YJK, to make Rate = 0 in case of repressor bound.

rate(1) = r_basal;
rate(2) = r_basal + r_Max;
rate(3) = r_basal - r_R;
rate(4) = r_basal + r_Max - r_R;

% Assume that the Bcd and Runt has the same frame rate, as well as same
% size of the matrices. This trimming of the dataset can be easily done separately.
Rate = zeros(length(NewBcdTime),length(APbin_Bcd));
sum_weights = zeros(length(NewBcdTime),length(APbin_Bcd));

for i=1:numStates
    Rate = Rate + (rate(i)*squeeze(weight(:,:,i))); % weighted sum of rates
    sum_weights = sum_weights + squeeze(weight(:,:,i)); % sum of all weights
end
% Normalize the Rate with the partition function(sum of all weights)
Rate = Rate ./ sum_weights;

% Quick plot of r1

%Plot at a time point. t=200 (=33min from NC13) is good for imputs high in NC14
figure(2)
%Time point can be changed here
t_Plot=200
plot((1:41),NewBcdFluo(t_Plot,:),(1:41),NewRuntFluo(t_Plot,:),(1:41),Rate(t_Plot,:))
legend('Bicoid','Runt','Predictive Rate')
title('quick imput/input plot')

 %% For general case of NR (Number of repressor binding sites) - additive repression
% % Number of repressor binding sites on hbP2 enhancer, this value should be
% % from 0 to 3.
% NR = 3;
% % empty,activators only, all activators and repressors bound, + different
% % number of repressors in either activators are bound or not (*2)
% numStates = 3 + 2*NR;
% 
% % interaction between Bcd and Runt : w (0=<w=<1, for now)
% % For now, we will assume that activators and repressors can bind
% % independently, and repressors act upon PIC or promoter to repress the
% % transcription initiation.
% w = 1; % This means there's no cooperativity
% 
% for i=4:numStates
%     weight(1) = 1;
%     weight(2) = NewBcdFluo.^6;
%     weight(3) = NewBcdFluo.^6.*NewRuntFluo.^NR*w;
%     if i+3 <= numStates - NR
%         weight(i) = NewRuntFluo.^i;
%     else
%         weight(i) = (NewRuntFluo.^i).*(NewBcdFluo.^6);
%     end
% end
% 
% % Rates
% % These values are from Eric Smith paper which used MS2.V1, thus could change.
% r_basal = 3; 
% r_Max = 9; 
% r_R = 1; % This is a naive guess by YJK, to make Rate = 0 in case of repressor bound.
% 
% for i=1:numStates
%     rate(1) = r_basal;
%     rate(2) = r_basal + r_Max;
%     rate(3) = r_basal + r_Max - r_R;
%     
%     if i+3 <= numStates - NR
%         rate(i) = r_basal - r_R*i;
%     else
%         % Here's another assumption that the rate is additive
%         rate(i) = r_basal + r_Max - r_R*i; 
%     end
% end
% 
% % Assume that the Bcd and Runt has the same frame rate, as well as same
% % size of the matrices. This trimming of the dataset can be easily done separately.
% Rate = zeros(length(NewBcdTime),APbin_Bcd);
% sum_weights = zeros(length(NewBcdTime),APbin_Bcd);
% 
% for i=1:numStates
%     Rate = Rate + (rate(i).*weight(i)); % weighted sum of rates
%     sum_weights = sum_weights + weight(i); % sum of all weights
% end
% % Normalize the Rate with the partition function(sum of all weights)
% Rate = Rate ./ sum_weights;

% %% For general case of NR (Number of repressor binding sites) - simple repression
% % Number of repressor binding sites on hbP2 enhancer, this value should be
% % from 0 to 3.
% NR = 3;
% % empty,activators only, all activators and repressors bound, + different
% % number of repressors in either activators are bound or not (*2)
% numStates = 3 + 2*NR; 
% 
% % interaction between Bcd and Runt : w (0=<w=<1, for now)
% % For now, we will assume that activators and repressors can bind
% % independently, and repressors act upon PIC or promoter to repress the
% % transcription initiation.
% w = 1; % This means there's no cooperativity
% 
% for i=4:numStates
%     weight(1) = 1;
%     weight(2) = Bcd.^6;
%     weight(3) = Bcd.^6.*Runt.^NR*w;
%     if i+3 <= numStates - NR
%         weight(i) = Runt.^i;
%     else
%         weight(i) = (Runt.^i).*(Bcd.^6);
%     end
% end
% 
% % Rates
% % These values are from Eric Smith paper which used MS2.V1, thus could change.
% r_basal = 3; 
% r_Max = 9; 
% r_R = 3; % This is a naive guess by YJK, to make Rate = 0 in case of repressor bound.
% 
% for i=1:numStates
%     rate(1) = r_basal;
%     rate(2) = r_basal + r_Max;
%     rate(3) = r_basal + r_Max - r_R;
%     
%     if i+3 <= numStates - NR
%         rate(i) = r_basal - r_R*i;
%     else
%         % Here's another assumption that the rate is additive
%         rate(i) = r_basal + r_Max - r_R; 
%     end
% end
% 
% % Assume that the Bcd and Runt has the same frame rate, as well as same
% % size of the matrices. This trimming of the dataset can be easily done separately.
% Rate = zeros(tLength_Bcd,APbin_Bcd);
% sum_weights = zeros(tLength_Bcd,APbin_Bcd);
% 
% for i=1:numStates
%     Rate = Rate + (rate(i).*weight(i)); % weighted sum of rates
%     sum_weights = sum_weights + weight(i); % sum of all weights
% end
% % Normalize the Rate with the partition function(sum of all weights)
% Rate = Rate ./ sum_weights;

%% Generate the MS2 traces for given rate profile

%Length of gene being transcribed (P2-MS2-LacZ) is about 5.3kb
Length=5.296;
%Set an elongation rate - Garcia et al. Current Biology 2013 measured
%1.54+/-0.14 kb per minute. There is some recent infor that suggests that
%this number is only 1.2 of what it should be. 
Elongrate=1.54;
%The length of the MS2 loops
MS2=1.275;
%And the length of LacZ
LacZ=Length-MS2; 
%Set max and min initiation rates in PolII molecules/minute starting with
%17 and 5 are roughly what we see in our experiemntal data for anterior 
%and posterior regions; - Liz&Jonathan's
RateMax=17;
RateMin=3; 
% Kd=8
%Let's have DelE defined as eAccessible-eInaccessible -- a negative Del E
%would be the wildtype (with Zelda case)
%DelE=-2;
%Define the cooperativity of the bicoid molecules
%w=1;


%Kd=10+(4*d);