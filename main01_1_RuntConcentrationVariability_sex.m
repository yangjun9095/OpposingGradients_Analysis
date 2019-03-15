% Compare the Runt nuclear concentration between multiple embryos
% The goal of this script is to generate plots for 
% 1) Runt nuc. concentration profile for male vs female
% 2) Runt nuc. concentration embryo-to-embryo variability
% 3) Runt nuc. concentration nucleus-to-nucleus variability
function [Result] = RuntConcentrationVariability

%% Different datasets (1024x256, 200Hz, 0.85 zoom, etc., better signal to noise than above)

RuntMale1 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos7\CompiledNuclei.mat')
RuntMale2 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos22\CompiledNuclei.mat')
RuntMale3 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos23\CompiledNuclei.mat')
RuntMale4 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-02-RuntJB3-vasa-eGFP-Pos18\CompiledNuclei.mat')

RuntFemale1 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos4\CompiledNuclei.mat')
RuntFemale2 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos8\CompiledNuclei.mat')
RuntFemale3 = load('E:\YangJoon\LivemRNA\Data\Dropbox;\DynamicsResults\2018-12-04-RuntJB3-vasa-eGFP-Pos13\CompiledNuclei.mat')

%% Define varialbles to use
Time_Male1   =  RuntMale1.nc13:RuntMale1.nc14+10;
Time_Male2   =  RuntMale2.nc13:RuntMale2.nc14+10;
Time_Male3   =  RuntMale3.nc13:RuntMale3.nc14+10;
Time_Male4   =  RuntMale4.nc13:RuntMale4.nc14+10;

Time_Female1 =  RuntFemale1.nc13:RuntFemale1.nc14+10;
Time_Female2 =  RuntFemale2.nc13:RuntFemale2.nc14+10;
Time_Female3 =  RuntFemale3.nc13:RuntFemale3.nc14+2;

%% Plot Runt nuclear fluo over Time (for one AP bin)
AP = 17;
hold on
%errorbar(RuntUnknownSex1.ElapsedTime(Time_Unknown1)+5,RuntUnknownSex1.MeanVectorAP(Time_Unknown1,AP),RuntUnknownSex1.SDVectorAP(Time_Unknown1,AP))
%errorbar(RuntUnknownSex2.ElapsedTime(Time_Unknown2)-15,RuntUnknownSex2.MeanVectorAP(Time_Unknown2,AP),RuntUnknownSex2.SDVectorAP(Time_Unknown2,AP))

%RuntMale1 dataset has some weird spike in the dataset (original movie), so
%I interpolated the value.
RuntMale1.MeanVectorAP(26,:) = nanmean([RuntMale1.MeanVectorAP(25,:),RuntMale1.MeanVectorAP(27,:)]);

errorbar(RuntMale1.ElapsedTime(Time_Male1) - RuntMale1.ElapsedTime(RuntMale1.nc13),...
            RuntMale1.MeanVectorAP(Time_Male1,AP),RuntMale1.SDVectorAP(Time_Male1,AP))
errorbar(RuntMale2.ElapsedTime(Time_Male2) - RuntMale2.ElapsedTime(RuntMale2.nc13),...
            RuntMale2.MeanVectorAP(Time_Male2,AP),RuntMale2.SDVectorAP(Time_Male2,AP))
errorbar(RuntMale3.ElapsedTime(Time_Male3) - RuntMale3.ElapsedTime(RuntMale3.nc13),...
            RuntMale3.MeanVectorAP(Time_Male3,AP),RuntMale3.SDVectorAP(Time_Male3,AP))
errorbar(RuntMale4.ElapsedTime(Time_Male4) - RuntMale4.ElapsedTime(RuntMale4.nc13),...
            RuntMale4.MeanVectorAP(Time_Male4,AP),RuntMale4.SDVectorAP(Time_Male4,AP))

errorbar(RuntFemale1.ElapsedTime(Time_Female1) - RuntFemale1.ElapsedTime(RuntFemale1.nc13),...
            RuntFemale1.MeanVectorAP(Time_Female1,AP),RuntFemale1.SDVectorAP(Time_Female1,AP))
errorbar(RuntFemale2.ElapsedTime(Time_Female2) - RuntFemale2.ElapsedTime(RuntFemale2.nc13),...
            RuntFemale2.MeanVectorAP(Time_Female2,AP),RuntFemale2.SDVectorAP(Time_Female2,AP))
% errorbar(RuntFemale3.ElapsedTime(Time_Female3) - RuntFemale3.ElapsedTime(RuntFemale3.nc13),...
%             RuntFemale3.MeanVectorAP(Time_Female3,AP),RuntFemale3.SDVectorAP(Time_Female3,AP))
% errorbar(RuntFemale4.ElapsedTime(Time_Female4) - RuntFemale4.ElapsedTime(RuntFemale4.nc13),...
%             RuntFemale4.MeanVectorAP(Time_Female4,AP),RuntFemale4.SDVectorAP(Time_Female4,AP))

%ylim([0 400])

title(['Runt nuclear concentration over Time',' @ ',num2str((AP-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Runt nuclear fluorescence (AU)')
legend('Male1','Male2','Male3','Male4','Female1','Female2','Female3','Location','southeast')
standardizeFigure(gca,legend,[])
%saveas(gcf,'E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Personalcodes-YangJoon')
%% Find the peak tpoint in nc13
Data = RuntFemale3;
hold on
for i=Data.nc13:Data.nc14
    errorbar(0:0.025:1,Data.MeanVectorAP(i,:),Data.SDVectorAP(i,:))
    i
    pause
end
%% Plot over AP axis
tPoint(1) = RuntMale1.nc13+10;
tPoint(2) = RuntMale2.nc13+10;
tPoint(3) = RuntMale3.nc13+10;
tPoint(4) = RuntMale4.nc13+11;
tPoint(5) = RuntFemale1.nc13+11;
tPoint(6) = RuntFemale2.nc13+10;
%tPoint(7) = RuntFemale3.nc13+10;

hold on
errorbar(0:0.025:1, RuntMale1.MeanVectorAP(tPoint(1),:), RuntMale1.SDVectorAP(tPoint(1),:))
errorbar(0:0.025:1, RuntMale2.MeanVectorAP(tPoint(2),:), RuntMale2.SDVectorAP(tPoint(2),:))
errorbar(0:0.025:1, RuntMale3.MeanVectorAP(tPoint(3),:), RuntMale3.SDVectorAP(tPoint(3),:))
errorbar(0:0.025:1, RuntMale3.MeanVectorAP(tPoint(4),:), RuntMale3.SDVectorAP(tPoint(4),:))

errorbar(0:0.025:1, RuntFemale1.MeanVectorAP(tPoint(5),:), RuntFemale1.SDVectorAP(tPoint(5),:))
errorbar(0:0.025:1, RuntFemale2.MeanVectorAP(tPoint(6),:), RuntFemale2.SDVectorAP(tPoint(6),:))
%errorbar(0.2:0.025:1.2, RuntFemale3.MeanVectorAP(tPoint(7),:), RuntFemale3.SDVectorAP(tPoint(7),:))

title(['Runt nuclear concentration over AP','@ peak of nc13'])
xlabel('AP axis (EL)')
ylabel('Runt nuclear fluorescence (AU)')
legend('Male1','Male2','Male3','Male4','Female1','Female2','Location','northwest')
standardizeFigure(gca,legend,[])

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