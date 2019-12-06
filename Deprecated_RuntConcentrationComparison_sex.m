% Compare the Runt nuclear concentration between multiple embryos
% The goal of this script is to generate plots for 
% 1) Runt nuc. concentration profile for male vs female
% 2) Runt nuc. concentration embryo-to-embryo variability
% 3) Runt nuc. concentration nucleus-to-nucleus variability
function [Result] = RuntConcentrationVariability

%% Different datasets (1024x256, 200Hz, 0.85 zoom, etc., better signal to noise than above)
RuntFemale = LoadMS2Sets('Runt-1min-200Hz-Female');
RuntMale = LoadMS2Sets('Runt-1min-200Hz-Male');

%% Sanity check for AP registration
Sample = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos8\CompiledNuclei.mat')

plot(0:0.025:1, Sample.MeanVectorAP)
%% Find the peak tpoint in nc13
Data = RuntMale;
for j=1:length(Data)
    clf
    hold on
    for i=Data(j).nc14:length(Data(j).ElapsedTime)%Data(j).nc13:Data(j).nc14
        errorbar(0:0.025:1,Data(j).MeanVectorAP(i,:),Data(j).SDVectorAP(i,:))
        i
        pause
    end
end
%% Plot over AP axis - Female
Data = RuntFemale;
% % NC13
% tPoint(1) = 24; %22:25;
% tPoint(2) = 19; %17:20;
% tPoint(3) = 13; %12:14;
% tPoint(4) = 25; %24:27;

% Beginning of NC14 (kinda static)
tPoint(1) = 35;
tPoint(2) = 35;
tPoint(3) = 30;
tPoint(4) = 50;

hold on
for j=1:length(Data)
    errorbar(0:0.025:1, Data(j).MeanVectorAP(tPoint(j),:), Data(j).SDVectorAP(tPoint(j),:))
end
title(['Runt nuclear concentration over AP','@ beginning of NC14'])
xlabel('AP axis (EL)')
ylabel('Runt nuclear fluorescence (AU)')
legend('1','2','3','4','Location','northwest')
StandardFigure(gcf,gca)

%% Plot over AP axis - Male
Data = RuntMale;

% Beginning of NC14 (kinda static)
tPoint(1) = 40;
tPoint(2) = 40;
tPoint(3) = 32;
tPoint(4) = 40;
tPoint(5) = 32;

hold on
for j=1:length(Data)
    errorbar(0:0.025:1, Data(j).MeanVectorAP(tPoint(j),:), Data(j).SDVectorAP(tPoint(j),:))
end
title(['Runt nuclear concentration over AP','@ beginning of NC14'])
xlabel('AP axis (EL)')
ylabel('Runt nuclear fluorescence (AU)')
legend('1','2','3','4','5','Location','northwest')
StandardFigure(gcf,gca)
%% Part2 :  Let's plot the male/female data using the LoadDataSets.m
% Load the datasets
RuntMale = LoadMS2Sets('Runt-1min-200Hz-Male')
RuntFemale = LoadMS2Sets('Runt-1min-200Hz-Female')

%% plot the nuclear fluorescence over AP at the peak of nc 13
tDelay= 13; % frames (~min)

% Male 
hold on
for i=1:length(RuntMale)
    RuntNC13(i) = RuntMale(i).nc13;
    RuntNucFluo(i,:) = RuntMale(i).MeanVectorAP(RuntNC13(i)+tDelay,:);
    NEllipses(i,:) = RuntMale(i).NParticlesAP(RuntNC13(i)+tDelay,:);
    RuntSDFluo(i,:) = RuntMale(i).SDVectorAP(RuntNC13(i)+tDelay,:);
    RuntFluoError(i,:) = RuntNucFluo(i,:)./sqrt(NEllipses(i,:));
    errorbar(0:0.025:1,RuntNucFluo(i,:),RuntFluoError(i,:),'b')
    pause
end


% Female
hold on
for i=1:length(RuntFemale)
    RuntNC13(i) = RuntFemale(i).nc13;
    RuntNucFluo(i,:) = RuntFemale(i).MeanVectorAP(RuntNC13(i)+tDelay,:);
    NEllipses(i,:) = RuntFemale(i).NParticlesAP(RuntNC13(i)+tDelay,:);
    RuntSDFluo(i,:) = RuntFemale(i).SDVectorAP(RuntNC13(i)+tDelay,:);
    RuntFluoError(i,:) = RuntNucFluo(i,:)./sqrt(NEllipses(i,:));
    errorbar(0:0.025:1,RuntNucFluo(i,:),RuntFluoError(i,:),'r')
    pause
end
    
%FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\Runt-JB3-eGFP\Runt-MaleFemale\';
%saveas(Runt_NucFluo_overAP_nc13_fig,[FigPath 'Runt_NucFluo_overAP_nc13' , '.pdf']); 
%% Plot over Time
hold on
AP = 14;
for i=1:length(RuntMale)
    RuntNC13(i) = RuntMale(i).nc13;
    ElapsedTime = RuntMale(i).ElapsedTime;
    RuntNucFluo(i,:) = RuntMale(i).MeanVectorAP(:,AP);
    NEllipses(i,:) = RuntMale(i).NParticlesAP(:,AP);
    RuntSDFluo(i,:) = RuntMale(i).SDVectorAP(:,AP);
    RuntFluoError(i,:) = RuntNucFluo(i,:)./sqrt(NEllipses(i,:));
    errorbar(ElapsedTime(RuntNC13(i):end) - ElapsedTime(RuntNC13(i)),...
                RuntNucFluo(RuntNC13(i):end,:),RuntFluoError(RuntNC13(i):end,:),'b')
end

%% Average multiple datasets
AverageDatasets_NuclearProtein('Runt-1min-200Hz-Male','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');
AverageDatasets_NuclearProtein('Runt-1min-200Hz-Female','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');

%% Plot the Averaged traces
RuntProtein_male = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-1min-200Hz-Male-Averaged.mat')
RuntProtein_female = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-1min-200Hz-Female-Averaged.mat')

NoNBControl = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\NoNB-Averaged.mat');
%% plot
AP = 19;
APbin = (AP-1)*2.5; % percent AP bin
tWindow_female = RuntProtein_female.nc13:length(RuntProtein_female.ElapsedTime);
tWindow_male = RuntProtein_male.nc13:length(RuntProtein_male.ElapsedTime);%RuntProtein_male.nc14;
tWindow_control = NoNBControl.nc13:length(NoNBControl.ElapsedTime);%NoNBControl.nc14;

% Interpolate the female timepoint 3, since it's definitely an artefact
RuntProtein_male.MeanVectorAP(RuntProtein_male.nc13 + 3,:) = ...
             (RuntProtein_male.MeanVectorAP(RuntProtein_male.nc13 + 2,:) +...
              RuntProtein_male.MeanVectorAP(RuntProtein_male.nc13 + 4,:))/2;

hold on
errorbar(NoNBControl.ElapsedTime(tWindow_control),...
            NoNBControl.MeanVectorAP(tWindow_control,AP),...
            NoNBControl.SEVectorAP(tWindow_control,AP))

errorbar(RuntProtein_male.ElapsedTime(tWindow_male),...
            RuntProtein_male.MeanVectorAP(tWindow_male,AP),...
            RuntProtein_male.SEVectorAP(tWindow_male,AP))
        
errorbar(RuntProtein_female.ElapsedTime(tWindow_female),...
            RuntProtein_female.MeanVectorAP(tWindow_female,AP),...
            RuntProtein_female.SEVectorAP(tWindow_female,AP))


xlim([0 35])
title('Avereaged Nuclear fluorescence')
legend('NoNB','male','female')
xlabel('Time (min)')
ylabel('Nuclear fluorescence')
standardizeFigure(gca,legend,[])
saveas(gcf,['E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Runt_Nuclear_fluo_overNC1314','_',num2str(APbin),'%','.pdf'])
%% plot over AP
tPoint = 13; % at 13 min into nc13 in this case, since it's 1 min time resolution.
hold on
errorbar(0:0.025:1,NoNBControl.MeanVectorAP(tPoint,:),NoNBControl.SEVectorAP(tPoint,:))
errorbar(0:0.025:1,RuntProtein_male.MeanVectorAP(tPoint,:),RuntProtein_male.SEVectorAP(tPoint,:))
errorbar(0:0.025:1,RuntProtein_female.MeanVectorAP(tPoint,:),RuntProtein_female.SEVectorAP(tPoint,:))


title('Avereaged Nuclear fluorescence')
legend('NoNB','male','female')
xlabel('AP')
ylabel('Nuclear fluorescence')
standardizeFigure(gca,legend,[])
saveas(gcf,'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Runt_Nuclear_fluo_overAP_NC13_peak.pdf')
end