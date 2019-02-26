% This script is for theoretical modeling for the Opposing gradient
% project. 
% Author : Yang Joon Kim (yjkim90@berkeley.edu)
% Last Updated : June/2018

% I will write down some Hill models to explain Bcd's cooperative binding,
% also considering Runt's repressive action.
% In future, I can think about other models, 
% The goal is to predict the rate of transcrition (or RNAP loading) over
% time, over all AP bins.

% For starters, I will start with several assumptions, which needs to be checked
% if they are fair.
% Assumptions :
% 1. There are 6 Bicoid binding sites on hb P2 enhancer
% 2. Runt's interaction with activator or PIC is unclear. Thus there could
% be many different scenarios, such as short-range repressor, long range
% repressor, etc.
% 3. 

%% Load the datasets
InputData = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\InputData.mat')

Bcd = InputData.Bcd;
Runt = InputData.Runt;

NewBcdTime = Bcd.ElapsedTime;
NewBcdFluo = Bcd.Fluo;

NewRuntTime = Runt.ElapsedTime;
NewRuntFluo = Runt.Fluo;
%% r0
%% Let's start with the next simplest case, r1
% Write down all states and corresponding weights
numStates = 4;

% interaction between Bcd and Runt : w (0=<w=<1, for now)
% For now, we will assume that activators and repressors can bind
% independently, and repressors act upon PIC or promoter to repress the
% transcription initiation.
w = 1; % This means there's no cooperativity

weight = zeros(length(NewBcdTime),length(Bcd.APbinID),numStates);

weight(:,:,1) = ones(length(NewBcdTime),length(Bcd.APbinID));
weight(:,:,2) = NewBcdFluo.^6;
weight(:,:,3) = NewRuntFluo;
weight(:,:,4) = NewBcdFluo.^6.*NewRuntFluo*w;

%% Rates
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
Rate = zeros(length(NewBcdTime),APbin_Bcd);
sum_weights = zeros(length(NewBcdTime),APbin_Bcd);

for i=1:numStates
    Rate = Rate + (rate(i)*squeeze(weight(:,:,i))); % weighted sum of rates
    sum_weights = sum_weights + squeeze(weight(:,:,i)); % sum of all weights
end
% Normalize the Rate with the partition function(sum of all weights)
Rate = Rate ./ sum_weights;

%% For general case of NR (Number of repressor binding sites) - additive repression
% Number of repressor binding sites on hbP2 enhancer, this value should be
% from 0 to 3.
NR = 2;
% empty,activators only, all activators and repressors bound, + different
% number of repressors in either activators are bound or not (*2)
numStates = 3 + 2*NR;

% interaction between Bcd and Runt : w (0=<w=<1, for now)
% For now, we will assume that activators and repressors can bind
% independently, and repressors act upon PIC or promoter to repress the
% transcription initiation.
w = 1; % This means there's no cooperativity

weight = zeros(length(NewBcdTime),length(Bcd.APbinID),numStates);

weight(:,:,1) = 1;
weight(:,:,2) = NewBcdFluo.^6;
weight(:,:,3) = NewBcdFluo.^6.*NewRuntFluo.^NR*w;
    
for i=4:numStates
    if i+3 <= numStates - NR
        weight(:,:,i) = NewRuntFluo.^i;
    else
        weight(:,:,i) = (NewRuntFluo.^i).*(NewBcdFluo.^6);
    end
end

% Rates
% These values are from Eric Smith paper which used MS2.V1, thus could change.
r_basal = 3; 
r_Max = 9; 
r_R = 3; % This is a naive guess by YJK, to make Rate = 0 in case of repressor bound.

for i=1:numStates
    rate(1) = r_basal;
    rate(2) = r_basal + r_Max;
    rate(3) = r_basal + r_Max - r_R;
    
    if i+3 <= numStates - NR
        rate(i) = r_basal - r_R*i;
    else
        % Here's another assumption that the rate is additive
        rate(i) = r_basal + r_Max - r_R*i; 
    end
    
    % Make negative rates to be zero.
    if rate(i)<0
        rate(i)=0;
    end
    
end

% Assume that the Bcd and Runt has the same frame rate, as well as same
% size of the matrices. This trimming of the dataset can be easily done separately.
Rate = zeros(length(NewBcdTime),APbin_Bcd);
sum_weights = zeros(length(NewBcdTime),APbin_Bcd);

for i=1:numStates
    Rate = Rate + (rate(i).*squeeze(weight(:,:,i))); % weighted sum of rates
    sum_weights = sum_weights + squeeze(weight(:,:,i)); % sum of all weights
end
% Normalize the Rate with the partition function(sum of all weights)
Rate_two_additive = Rate ./ sum_weights;

%% For general case of NR (Number of repressor binding sites) - simple repression
% Number of repressor binding sites on hbP2 enhancer, this value should be
% from 0 to 3.
NR = 2;
% empty,activators only, all activators and repressors bound, + different
% number of repressors in either activators are bound or not (*2)
numStates = 3 + 2*NR; 

% interaction between Bcd and Runt : w (0=<w=<1, for now)
% For now, we will assume that activators and repressors can bind
% independently, and repressors act upon PIC or promoter to repress the
% transcription initiation.
w = 1; % This means there's no cooperativity

weight = zeros(length(NewBcdTime),length(Bcd.APbinID),numStates);
weight(:,:,1) = 1;
weight(:,:,2) = NewBcdFluo.^6;
weight(:,:,3) = NewBcdFluo.^6.*NewRuntFluo.^NR*w;

for i=4:numStates
    if i+3 <= numStates - NR
        weight(:,:,i) = NewRuntFluo.^i;
    else
        weight(:,:,i) = (NewRuntFluo.^i).*(NewBcdFluo.^6);
    end
end

% Rates
% These values are from Eric Smith paper which used MS2.V1, thus could change.
r_basal = 3; 
r_Max = 9; 
r_R = 3; % This is a naive guess by YJK, to make Rate = 0 in case of repressor bound.

rate(1) = r_basal;
rate(2) = r_basal + r_Max;
rate(3) = r_basal + r_Max - r_R;
    
for i=4:numStates
    if i+3 <= numStates - NR
        rate(i) = r_basal - r_R*i;
    else
        % Here's another assumption that the rate is additive
        rate(i) = r_basal + r_Max - r_R; 
    end
    
    % Make negative rates to be zero.
    if rate(i)<0
        rate(i)=0;
    end
end

% Assume that the Bcd and Runt has the same frame rate, as well as same
% size of the matrices. This trimming of the dataset can be easily done separately.
Rate = zeros(length(NewBcdTime),APbin_Bcd);
sum_weights = zeros(length(NewBcdTime),APbin_Bcd);

for i=1:numStates
    Rate = Rate + (rate(i).*squeeze(weight(:,:,i))); % weighted sum of rates
    sum_weights = sum_weights + squeeze(weight(:,:,i)); % sum of all weights
end
% Normalize the Rate with the partition function(sum of all weights)
Rate_two_simple = Rate ./ sum_weights;

%% Pick one AP bin or one time point, and compare two different scenarios (simple vs additive)
tpoint = 28; %peak of nc13
%tpoint = ; %peak of nc14

hold on
% plot(0:0.025:1,Rate(tpoint,:))
% plot(0:0.025:1,Rate_two_simple(tpoint,:))
% plot(0:0.025:1,Rate_three_simple(tpoint,:))

plot(0:0.025:1,Rate(i,:))
plot(0:0.025:1,Rate_two_additive(i,:))
plot(0:0.025:1,Rate_three_additive(i,:))

legend('','','')

%% Now, I want to compare two different scenarios in case of r1,2,3

% First, let's get the maximum RNAP loading rate in nc13 and nc14
% respectively. This is because the maximum loading rate time might not be
% synchronized. 
% For nc13, the time series is 1~100th frames
% r0
Rate_zero_nc13 = max(Rate_zero(1:100,:));
% r1 (note that there's no difference between simple vs additive in this
% case)
Rate_one_nc13 = max(Rate_one(1:100,:));
% r2
Rate_two_simple_nc13 = max(Rate_two_simple(1:100,:));
Rate_two_additive_nc13 = max(Rate_two_additive(1:100,:));
% r3
Rate_three_simple_nc13 = max(Rate_three_simple(1:100,:));
Rate_three_additive_nc13 = max(Rate_three_additive(1:100,:));

%% Plot
hold on
plot(0:0.025:1,Rate_zero_nc13)
plot(0:0.025:1,Rate_one_nc13)
plot(0:0.025:1,Rate_two_simple_nc13)
plot(0:0.025:1,Rate_two_additive_nc13)
plot(0:0.025:1,Rate_three_simple_nc13)
plot(0:0.025:1,Rate_three_additive_nc13)
legend('','','','','','')
%% Save
savePath = uigetdir; % Click the folder that you want to save the .mat file.


save([savePath,filesep,'Rate_Prediction','.mat'],...
            'Rate_zero','Rate_one','Rate_two_simple','Rate_two_additive',...
            'Rate_three_simple','Rate_three_additive')
        
%% 

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

