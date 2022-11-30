% Compare the Runt nuclear concentration between multiple embryos
% The goal of this script is to generate plots for 
% 1) Runt nuc. concentration profile for male vs female
% 2) Runt nuc. concentration embryo-to-embryo variability
% 3) Runt nuc. concentration nucleus-to-nucleus variability
function [Result] = RuntConcentrationVariability

%% Load datasets
RuntMale1 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-14-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1-1\CompiledNuclei.mat')
RuntMale2 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-28-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1-7\CompiledNuclei.mat')

RuntFemale1 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-14-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1-2\CompiledNuclei.mat')
RuntFemale2 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-14-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1-3\CompiledNuclei.mat')
RuntFemale3 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-28-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1-5\CompiledNuclei.mat')
RuntFemale4 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-28-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1-6\CompiledNuclei.mat')

RuntUnknownSex1 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-07-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1\CompiledNuclei.mat')
%RuntUnknownSex2 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-28-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1-5\CompiledNuclei.mat')

%% Plot nuclear fluorescence (at different AP, time)
AP=13;
% Time_Male1 = RuntMale1.nc13:RuntMale1.nc14;
% Time_Female1 = RuntFemale1.nc13:RuntFemale1.nc14;
% Time_Female2 = RuntFemale2.nc13:RuntFemale2.nc14;

Time_Male1   =  RuntMale1.nc13:3:RuntMale1.nc14+50;
Time_Male2   =  RuntMale2.nc13:length(RuntMale2.ElapsedTime);
Time_Female1 =  RuntFemale1.nc13:3:RuntFemale1.nc14+50;
Time_Female2 =  RuntFemale2.nc13:RuntFemale2.nc14+50;
Time_Female3 =  RuntFemale3.nc13:RuntFemale3.nc14+20;
Time_Female4 =  RuntFemale4.nc13:RuntFemale4.nc14+20;

%Time_Unknown1 = RuntUnknownSex1.nc13:RuntUnknownSex1.nc14+50;

hold on
%errorbar(RuntUnknownSex1.ElapsedTime(Time_Unknown1)+5,RuntUnknownSex1.MeanVectorAP(Time_Unknown1,AP),RuntUnknownSex1.SDVectorAP(Time_Unknown1,AP))
%errorbar(RuntUnknownSex2.ElapsedTime(Time_Unknown2)-15,RuntUnknownSex2.MeanVectorAP(Time_Unknown2,AP),RuntUnknownSex2.SDVectorAP(Time_Unknown2,AP))

errorbar(RuntMale1.ElapsedTime(Time_Male1) - RuntMale1.ElapsedTime(RuntMale1.nc13),...
            RuntMale1.MeanVectorAP(Time_Male1,AP),RuntMale1.SDVectorAP(Time_Male1,AP))
% errorbar(RuntMale2.ElapsedTime(Time_Male2) - RuntMale2.ElapsedTime(RuntMale2.nc13),...
%             RuntMale2.MeanVectorAP(Time_Male2,AP),RuntMale2.SDVectorAP(Time_Male2,AP))

errorbar(RuntFemale1.ElapsedTime(Time_Female1) - RuntFemale1.ElapsedTime(RuntFemale1.nc13),...
            RuntFemale1.MeanVectorAP(Time_Female1,AP),RuntFemale1.SDVectorAP(Time_Female1,AP))
% errorbar(RuntFemale2.ElapsedTime(Time_Female2) - RuntFemale2.ElapsedTime(RuntFemale2.nc13),...
%             RuntFemale2.MeanVectorAP(Time_Female2,AP),RuntFemale2.SDVectorAP(Time_Female2,AP))
% errorbar(RuntFemale3.ElapsedTime(Time_Female3) - RuntFemale3.ElapsedTime(RuntFemale3.nc13),...
%             RuntFemale3.MeanVectorAP(Time_Female3,AP),RuntFemale3.SDVectorAP(Time_Female3,AP))
% errorbar(RuntFemale4.ElapsedTime(Time_Female4) - RuntFemale4.ElapsedTime(RuntFemale4.nc13),...
%             RuntFemale4.MeanVectorAP(Time_Female4,AP),RuntFemale4.SDVectorAP(Time_Female4,AP))



title(['Runt nuclear concentration over Time',' @ ',num2str((AP-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Runt nuclear fluorescence (AU)')
legend('Male1','Female1')
%legend('Male1','Male2','Female1','Female2','Female3','Female4','Location','southeast')
standardizeFigure_YJK(gca,legend,[])

%% I'll use different datasets (acquired 1 min intervals, for multiplexing the imaging)
% Now, I have roughly 10 datasets with this imaging condition, and it's
% possible to get maximum of 7 datasets per session. Also, I have some NoNB
% control datsets (~3) for this imaging condition.

% Load datasets
RuntMale1 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-RuntJB3-vasa-eGFP1-Pos4\CompiledNuclei.mat')
RuntMale2 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-RuntJB3-vasa-eGFP1-Pos15\CompiledNuclei.mat')

RuntFemale1 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-RuntJB3-vasa-eGFP1-Pos1\CompiledNuclei.mat')
RuntFemale2 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-RuntJB3-vasa-eGFP1-Pos5\CompiledNuclei.mat')
RuntFemale3 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-RuntJB3-vasa-eGFP1-Pos6\CompiledNuclei.mat')


%% Different datasets (1024x256, 200Hz, 0.85 zoom, etc., better signal to noise than above)

RuntMale1 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos7\CompiledNuclei.mat')
RuntMale2 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos22\CompiledNuclei.mat')
RuntMale3 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos23\CompiledNuclei.mat')
RuntMale4 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-12-02-RuntJB3-vasa-eGFP-Pos18\CompiledNuclei.mat')

RuntFemale1 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos4\CompiledNuclei.mat')
RuntFemale2 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-12-03-RuntJB3-vasa-eGFP-Pos8\CompiledNuclei.mat')
RuntFemale3 = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-12-04-RuntJB3-vasa-eGFP-Pos13\CompiledNuclei.mat')

%% Define varialbles to use
Time_Male1   =  RuntMale1.nc13:RuntMale1.nc14+10;
Time_Male2   =  RuntMale2.nc13:RuntMale2.nc14+10;
Time_Male3   =  RuntMale3.nc13:RuntMale3.nc14+10;
Time_Male4   =  RuntMale4.nc13:RuntMale4.nc14+10;

Time_Female1 =  RuntFemale1.nc13:RuntFemale1.nc14+10;
Time_Female2 =  RuntFemale2.nc13:RuntFemale2.nc14+10;
Time_Female3 =  RuntFemale3.nc13:RuntFemale3.nc14+2;

%% Plot Runt nuclear fluo over Time (for one AP bin)
AP = 21;
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
standardizeFigure_YJK(gca,legend,[])

%% Part2 :  Let's plot the male/female data using the LoadDataSets.m
% Load the datasets
RuntMale = LoadMS2Sets('Runt-1min-200Hz-Male')
RuntFemale = LoadMS2Sets('Runt-1min-200Hz-Female')

%% plot the nuclear fluorescence at the peak of nc 13
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
    
%% Plot the Averaged traces
RuntProtein_female = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Runt-1min-200Hz-Female-Averaged.mat');
RuntProtein_male = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Runt-1min-200Hz-Male-Averaged.mat');
NoNBControl = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\NoNB-Averaged.mat');
%% plot
AP = 17;
tWindow_female = RuntProtein_female.nc13:length(RuntProtein_female.ElapsedTime);
tWindow_male = RuntProtein_male.nc13:length(RuntProtein_male.ElapsedTime);%RuntProtein_male.nc14;
tWindow_control = NoNBControl.nc13:length(NoNBControl.ElapsedTime);%NoNBControl.nc14;

% Interpolate the female timepoint 3, since it's definitely an artefact
RuntProtein_male.MeanVectorAP(RuntProtein_male.nc13 + 3,:) = ...
             (RuntProtein_male.MeanVectorAP(RuntProtein_male.nc13 + 2,:) +...
              RuntProtein_male.MeanVectorAP(RuntProtein_male.nc13 + 4,:))/2;

hold on
errorbar(RuntProtein_female.ElapsedTime(tWindow_female),...
            RuntProtein_female.MeanVectorAP(tWindow_female,AP),...
            RuntProtein_female.SEVectorAP(tWindow_female,AP))
errorbar(RuntProtein_male.ElapsedTime(tWindow_male),...
            RuntProtein_male.MeanVectorAP(tWindow_male,AP),...
            RuntProtein_male.SEVectorAP(tWindow_male,AP))
errorbar(NoNBControl.ElapsedTime(tWindow_control),...
            NoNBControl.MeanVectorAP(tWindow_control,AP),...
            NoNBControl.SEVectorAP(tWindow_control,AP))
xlim([0 35])
title('Avereaged Nuclear fluorescence')
legend('female','male','NoNB')
xlabel('Time (min)')
ylabel('Nuclear fluorescence')
standardizeFigure_YJK(gca,legend,'purple','blue','red')
%% plot over AP
hold on
errorbar(0:0.025:1,RuntProtein_female.MeanVectorAP(13,:),RuntProtein_female.SEVectorAP(13,:))
errorbar(0:0.025:1,RuntProtein_male.MeanVectorAP(13,:),RuntProtein_male.SEVectorAP(13,:))
errorbar(0:0.025:1,NoNBControl.MeanVectorAP(13,:),NoNBControl.SEVectorAP(13,:))

title('Avereaged Nuclear fluorescence')
legend('female','male','NoNB')
xlabel('AP')
ylabel('Nuclear fluorescence')
end