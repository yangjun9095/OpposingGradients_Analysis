function modeling_01_generatePredictions_hbP2
% DESCRIPTION
% Script for generating predictions for the rate of transcription as a
% function of activator and repressor concentrations.

%% Input
% activator concentration regime, repressor concentration regime,
% 

% For this, we don't know how this system works actually, meaning the way
% repressor works, also it doesn't seem like that simple (not a direct
% competition, nor a direct repression to the RNAP recruitment). But, from
% our theoretical exploration, it doesn't seem to be distinguishable by
% modeling. 
% Thus, maybe it's better to make this script to be able to compare the
% results from different models 
% : not just for adding complexity (from Hill model to individual sites), 
% but also deal with different scenarios of repression (not sure if this is necessary)

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

ColorChoice = [colorDict.purple; colorDict.green; colorDict.blue; colorDict.red; colorDict.brown]; % 4 embryos max. it could be extended easily
%% Model
% Hill model of 6 activator sites + n repressor sites
clear bound
clear P_bound_1RepSite
clear partitionFunct

A = logspace(-3,3,100);
R = logspace(-3,3,100);

p = 0.01;
wAP = 3;
%wRP = 0;
%wRP = 0.001;

for i=1:length(A)
    for j=1:length(R)
        partitionFunct(i,j) = (1+A(i).^6 + R(j) + A(i).^6*R(j)) +...
                                p*(1+A(i).^6*wAP + R(j)*wRP + A(i).^6*R(j)*wAP*wRP);
        bound(i,j) = p*(1+A(i).^6*wAP + R(j)*wRP + A(i).^6*R(j)*wAP*wRP);
    end
end

P_bound_1RepSite = bound./partitionFunct;


% plot the P_bound as a 3D plot of Activator and Repressor regime
[X,Y] = meshgrid(A,R);
interval = 4;
surf(A(1:interval:end),R(1:interval:end),P_bound_1RepSite(1:interval:end,1:interval:end))%,...
        %'FaceColor',[234,194,100]/255)

% x, y as log scale
set(gca,'XScale','log')
set(gca,'YScale','log')

% get rid of Ticks
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])

box on
grid on
alpha 0.7 % transparency
colorbar( 'off')

ylabel('Activator')
xlabel('Repressor')
zlabel('RNAP loading rate')
StandardFigure(gcf,gca)


%% Plot the actual Bicoid and Runt protein concentration gradient on top as a line.
% AP axis as from 20-60%
AP = 0.2:0.025:0.8;
% Bcd : using the length constant of 0.2
Bcd = exp(-5*AP);
Bcdscale = 25;
Bcd = Bcdscale*Bcd;

% Bcd, Runt : using the time-averaged data
filePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';
BcdData = load([filePath, filesep, 'Bcd-Averaged.mat']);
RuntData =  load([filePath, filesep, 'Runt_TimeAveraged_mixedSex_NC14.mat']);

% Average over the 10 minutes into nc14
Bcd_nc13 = BcdData.nc13;
Bcd_nc14 = BcdData.nc14;
BcdFluo = BcdData.MeanVectorAP;
BcdTime = BcdData.ElapsedTime;
BcdFluoSD = BcdData.SDVectorAP;
BcdFluoSE = BcdData.SEVectorAP;

%% generate Bcd time traces and spatial gradient plots
% Bicoid concentration time trace (20%)
APbin = 11; % 25%

hold on
errorbar(BcdTime(Bcd_nc13:end) - BcdTime(Bcd_nc13),...
            BcdFluo(Bcd_nc13:end,APbin),...
            BcdFluoSE(Bcd_nc13:end,APbin))
        
xline(18,'--')

xlim([0 60])
        
xlabel('time (min)')
ylabel('Bicoid concentration (AU)')


box on
StandardFigure(gcf,gca)

% save the plots

%% time-average over 10 min into nc14
% Bicoid concentration time trace (20%)
time_window = 10; % min
tFrame = median(diff(BcdTime));

Nstep_window = floor(time_window/tFrame);

BcdFluo_tAveraged = nanmean(BcdFluo(Bcd_nc14:Bcd_nc14+Nstep_window,:));
BcdFluo_tAveraged_SE = nanmean(BcdFluoSE(Bcd_nc14:Bcd_nc14+Nstep_window,:));

APaxis = 0:0.025:1;

hold on
errorbar(APaxis(10:end), BcdFluo_tAveraged(10:end), BcdFluo_tAveraged_SE(10:end))

xlim([0.2 0.6])
ylim([0 600])
       
xlabel('embryo length')
ylabel('Bicoid concentration (AU)')


box on
StandardFigure(gcf,gca)

%% Bcd extrapolation
BcdFluo_extrap = interp1(0.225:0.025:0.6, BcdFluo_tAveraged(10:25), 0.2:0.025:0.8, 'linear','extrap')


%%
RuntFluo_tAveraged = RuntData.AveragedFluo_tAveraged_mixed(2,:); % averaged within the 10 min into nc14
RuntFluo_tAveraged_SE = RuntData.SEFluo_tAveraged_mixed(2,:); % averaged within the 10 min into nc14

RuntFluo = RuntFluo_tAveraged;
RuntFluo_tAveraged(isnan(RuntFluo)) = [];
RuntFluo_tAveraged_SE(isnan(RuntFluo)) = [];
APaxis = 0:0.025:1;
APaxis(isnan(RuntFluo)) = [];

% extrapolate the posterior data points
RuntFluo_extrap = interp1(APaxis, RuntFluo_tAveraged, 0.2:0.025:0.8, 'linear','extrap');
RuntFluo_SE_extrap = interp1(APaxis, RuntFluo_tAveraged_SE, 0.2:0.025:0.8, 'linear','extrap');

Run = RuntFluo_extrap; % taking only 20-80% of AP axis

Run = Run/30;
Runscale = 0.01;

Run = Run*Runscale;

%% different sources of Bicoid (Liz and Jonathan)
BcdData_EE_JL = load('S:\YangJoon\Dropbox\OpposingGradient\eGFP-Bcd-From-Liz-Jonathan\BcdLevelsNC13.mat');

% extract data from 0.2-0.8 
BcdFluoTemp = nanmean(BcdData_EE_JL.BcdNC13Ant);
BcdFluoTemp_SE = nanstd(BcdData_EE_JL.BcdNC13Ant)./sqrt(42); % the number of time points are 42.
BcdFluo_NC13(1:17) = BcdFluoTemp(9:25); %20-60%
BcdFluo_NC13_SE(1:17) = BcdFluoTemp_SE(9:25);

BcdFluoTemp = nanmean(BcdData_EE_JL.BcdNC13Pos);
BcdFluoTemp_SE = nanstd(BcdData_EE_JL.BcdNC13Pos)./sqrt(42); % the number of time points are 42.
BcdFluo_NC13(18:25) = BcdFluoTemp(26:33); %62.5-80%
BcdFluo_NC13_SE(18:25) = BcdFluoTemp_SE(26:33);

%% generate plots of Bicoid and Runt spatial gradient

AP = 0.2:0.025:0.8;
% % Bcd : using the length constant of 0.2
% Bcd = exp(-5*AP);
% Bcdscale = 1500;
% Bcd = Bcdscale*Bcd;

hold on
% plot(0.2:0.025:0.8, Bcd,'color',ColorChoice(1,:),'LineWidth',2)
% plot(0.2:0.025:0.8, BcdFluo_extrap,'color',ColorChoice(1,:),'LineWidth',2)

% plot(0.2:0.025:0.8, BcdFluo_NC13*5,'color',ColorChoice(1,:),'LineWidth',2)
% plot(0.2:0.025:0.8,RuntFluo_extrap*1.5,'color',ColorChoice(2,:),'LineWidth',2)

errorbar(0.2:0.025:0.8, BcdFluo_NC13*5, BcdFluo_NC13_SE*5,'LineWidth',2)
errorbar(0.2:0.025:0.8,RuntFluo_extrap*1.5, RuntFluo_SE_extrap, 'color',ColorChoice(2,:),'LineWidth',2)

xlim([0.2 0.8])
ylim([0 600])
       
xlabel('embryo length')
ylabel('transcription factor conc.(AU)')

legend('Bicoid','Runt')
box on
StandardFigure(gcf,gca)

% save figure
figPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\Modeling';
saveas(gcf,[figPath,filesep,'inputTF_spatial_gradient','.tif']); 
saveas(gcf,[figPath,filesep,'inputTF_spatial_gradient','.pdf']); 


%% model
% Bcd : using the length constant of 0.2
Bcd = exp(-5*AP);
Bcdscale = 25;
Bcd = Bcdscale*Bcd;

Run = RuntFluo_extrap; % taking only 20-80% of AP axis
Run = Run/30;
Runscale = 0.01;
Run = Run*Runscale;


p = 0.01;
wAP = 3;
wRP = 0;
wRP = 0.001;

clear bound
clear partitionFunct
clear P_bound
for i=1:length(AP)
    
    partitionFunct(i) = (1+Bcd(i).^6 + Run(i) + Bcd(i).^6*Run(i)) +...
                            p*(1+Bcd(i).^6*wAP + Run(i)*wRP + Bcd(i).^6*Run(i)*wAP*wRP);
    bound(i) =p*(1+Bcd(i).^6*wAP + Run(i)*wRP + Bcd(i).^6*Run(i)*wAP*wRP);

end

P_bound = bound./partitionFunct;

% plot on top of the 3D surface plot
hold on
plot3(Run,Bcd,P_bound,'LineWidth',2,'Color','r')

plot(Run,Bcd,'LineWidth',2,'Color','k')

% Save the figures
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\OpposingGradientsFigures\Modeling';
%saveas(gcf,[figPath,filesep,'prediction_act_rep_3D','.tif']); 
%saveas(gcf,[figPath,filesep,'prediction_act_rep_3D','.pdf']); 

%% [Run] = 0 case (or n_R =0)
for i=1:length(AP)
    partitionFunct(i) = (1+Bcd(i).^6 ) +...
                            p*(1+Bcd(i).^6*wAP);
    bound(i) =p*(1+Bcd(i).^6*wAP );

end

P_bound_RunNull = bound./partitionFunct;

% plot on top of the 3D surface plot
hold on
plot3(Run,Bcd,P_bound_RunNull,'LineWidth',2,'Color',ColorChoice(1,:))

plot(Run,Bcd,'LineWidth',2,'Color','k')

% Save the figures
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\OpposingGradientsFigures\Modeling';
%saveas(gcf,[figPath,filesep,'prediction_act_rep_3D','.tif']); 
%saveas(gcf,[figPath,filesep,'prediction_act_rep_3D','.pdf']); 
%% Plot the predicted rate of Txn along the AP axis
figure(2)
hold on
plot(AP, P_bound,'LineWidth',2,'Color','r')
plot(AP,P_bound_RunNull,'LineWidth',2,'Color',ColorChoice(1,:))
xticks([0.2 0.3 0.4 0.5 0.6 0.7 0.8])

set(gca,'yticklabel',[])
box on

xlabel('AP axis (embryo length)')
ylabel('RNAP loading rate')
legend('WT','No Repressor','Location','SouthWest')

StandardFigure(gcf,gca)
%save the plot
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\OpposingGradientsFigures\Modeling';
saveas(gcf,[figPath,filesep,'prediction_rep_NoRep_overAP','.tif']); 
saveas(gcf,[figPath,filesep,'prediction_rep_NoRep_overAP','.pdf']); 

%% Try with 2 repressor binding sites?

clear bound
clear partitionFunct
clear P_bound
for i=1:length(AP)
    
        partitionFunct(i) = (1+Bcd(i).^6 + 2*Run(i) + Run(i)^2 + 2*Bcd(i).^6*Run(i)) +...
                                Bcd(i).^6*Run(i)^2 + ...
                                p*(1+Bcd(i).^6*wAP + 2*Run(i)*wRP + Run(i)^2*wRP^2 +...
                                2*Bcd(i).^6*Run(i)*wAP*wRP + Bcd(i).^6*Run(i)^2*wAP*wRP^2);
        bound(i) = p*(1+Bcd(i).^6*wAP + 2*Run(i)*wRP + Run(i)^2*wRP^2 +...
                                2*Bcd(i).^6*Run(i)*wAP*wRP + Bcd(i).^6*Run(i)^2*wAP*wRP^2);
end

P_bound_2sites = bound./partitionFunct;


%% Try with 3 repressor binding sites?
clear bound
clear partitionFunct
for i=1:length(AP)
    
        partitionFunct(i) = (1+Bcd(i).^6 + 3*Run(i) + 3*Run(i)^2 + Run(i)^3 +...
                                3*Bcd(i).^6*Run(i) + 3*Bcd(i).^6*Run(i)^2 + Bcd(i).^6*Run(i)^3) + ...
                                p*(1+Bcd(i).^6*wAP + 2*Run(i)*wRP + Run(i)^2*wRP^2 +...
                                2*Bcd(i).^6*Run(i)*wAP*wRP + Bcd(i).^6*Run(i)^2*wAP*wRP^2);
        bound(i) = p*(1+Bcd(i).^6*wAP + 2*Run(i)*wRP + Run(i)^2*wRP^2 +...
                                2*Bcd(i).^6*Run(i)*wAP*wRP + Bcd(i).^6*Run(i)^2*wAP*wRP^2);
end

P_bound_3sites = bound./partitionFunct;


%% plot altogether, 0,1,2,3 sites
figure(3)
hold on
plot(AP,P_bound_RunNull,'LineWidth',2,'Color',ColorChoice(1,:))
plot(AP, P_bound,'LineWidth',2,'Color','r')
plot(AP,P_bound_2sites,'LineWidth',2,'Color',ColorChoice(2,:))
plot(AP,P_bound_3sites,'LineWidth',2,'Color',ColorChoice(3,:))
xticks([0.2 0.3 0.4 0.5 0.6 0.7 0.8])

set(gca,'yticklabel',[])
box on

xlabel('AP axis (embryo length)')
ylabel('RNAP loading rate')
legend('0','1 site','2 sites','3 sites','Location','SouthWest')

StandardFigure(gcf,gca)
%save the plot
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\OpposingGradientsFigures\Modeling';
saveas(gcf,[figPath,filesep,'prediction_num_RepSites_overAP','.tif']); 
saveas(gcf,[figPath,filesep,'prediction_num_RepSites_overAP','.pdf']); 


%% 3D space prediction (input-output) : 
% Hill model of 6 activator sites + n repressor sites
% clear bound
% clear P_bound_2RepSites
% clear partitionFunct
% 
% A = logspace(-3,3,100);
% R = logspace(-3,3,100);
% 
% p = 0.01;
% wAP = 3;
% wRP = 0.001;
% 
% for i=1:length(A)
%     for j=1:length(R)
%         partitionFunct(i,j) = (1+A(i).^6 + 2*R(j) + R(j)^2 + 2*A(i).^6*R(j)) +...
%                                 A(i).^6*R(j)^2 + ...
%                                 p*(1+A(i).^6*wAP + 2*R(j)*wRP + R(j)^2*wRP^2 +...
%                                 2*A(i).^6*R(j)*wAP*wRP + A(i).^6*R(j)^2*wAP*wRP^2);
%         bound(i,j) = p*(1+A(i).^6*wAP + 2*R(j)*wRP + R(j)^2*wRP^2 +...
%                                 2*A(i).^6*R(j)*wAP*wRP + A(i).^6*R(j)^2*wAP*wRP^2);
%     end
% end
% 
% P_bound_2RepSites = bound./partitionFunct;
% 
% 
% % plot the P_bound as a 3D plot of Activator and Repressor regime
% [X,Y] = meshgrid(A,R);
% interval = 5;
% surf(A(1:interval:end),R(1:interval:end),P_bound_2RepSites(1:interval:end,1:interval:end))
% 
% % x, y as log scale
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% 
% % get rid of Ticks
% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
% set(gca,'zticklabel',[])
% 
% box on
% grid on
% 
% ylabel('Activator')
% xlabel('Repressor')
% zlabel('RNAP loading rate')
% StandardFigure(gcf,gca)

%% Try with 3 repressor binding sites?
% Hill model of 6 activator sites + n repressor sites
clear bound
clear P_bound_3RepSites
clear partitionFunct

A = logspace(-3,3,100);
R = logspace(-3,3,100);

p = 0.01;
wAP = 3;
wRP = 0.001;

for i=1:length(A)
    for j=1:length(R)
        partitionFunct(i,j) = (1+A(i).^6 + 3*R(j) + 3*R(j)^2 + R(j)^3 +...
                                3*A(i).^6*R(j) + 3*A(i).^6*R(j)^2 + A(i).^6*R(j)^3) + ...
                                p*(1+A(i).^6*wAP + 2*R(j)*wRP + R(j)^2*wRP^2 +...
                                2*A(i).^6*R(j)*wAP*wRP + A(i).^6*R(j)^2*wAP*wRP^2);
        bound(i,j) = p*(1+A(i).^6*wAP + 2*R(j)*wRP + R(j)^2*wRP^2 +...
                                2*A(i).^6*R(j)*wAP*wRP + A(i).^6*R(j)^2*wAP*wRP^2);
    end
end

P_bound_2RepSites = bound./partitionFunct;


% plot the P_bound as a 3D plot of Activator and Repressor regime
[X,Y] = meshgrid(A,R);
interval = 5;
surf(A(1:interval:end),R(1:interval:end),P_bound_2RepSites(1:interval:end,1:interval:end))

% x, y as log scale
set(gca,'XScale','log')
set(gca,'YScale','log')

% get rid of Ticks
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])

box on
grid on

ylabel('Activator')
xlabel('Repressor')
zlabel('RNAP loading rate')
StandardFigure(gcf,gca)
%% Step0. Can we do somewhat phenomenological Hill model like in Jeehae's paper?
% I don't know how to incorporate the repressor terms in here...
% How about 1/(1 + ([R]/Kr)^Nr)

%% Step1. Let's start with the simplest Hill model of Bcd activating the hb P2

end