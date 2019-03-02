% Bleaching Test for JB3-Runt:eGFP, and MCP-mCherry
% Yang Joon Kim (yjkim90@berkeley.edu)
% Last edited 11/19/2018

% Description
% The idea here is to use the combination of Sequantial scanning and ROI,
% to have two different exposure frequency in one imaging (one embryo).
% Thus, 

%% Check the data by plotting
% Load the dataset
% MS2Data = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-hbP2-r0-MS2V5-lacZ-2xMCP-mCherry-RuntJB3-vasa-eGFP1-20sec-BleachingTest\CompiledParticles.mat');
% RuntData = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-hbP2-r0-MS2V5-lacZ-2xMCP-mCherry-RuntJB3-vasa-eGFP1-20sec-BleachingTest\CompiledNuclei.mat');

MS2Data = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-hbP2-r0-MS2V5-lacZ-2xMCP-mCherry-RuntJB3-vasa-eGFP1-488_10uW-587_15uW\CompiledParticles.mat');
RuntData = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-11-18-hbP2-r0-MS2V5-lacZ-2xMCP-mCherry-RuntJB3-vasa-eGFP1-488_10uW-587_15uW\CompiledNuclei.mat');


Spot_ROI = MS2Data.MeanVectorAP_ROI{1,1};
Spot_nonROI = MS2Data.MeanVectorAP_nonROI;
Spot_Time = MS2Data.ElapsedTime;
Spot_Time = Spot_Time';

Spot_NROI = MS2Data.NParticlesAP_ROI{1,1};
Spot_NnonROI = MS2Data.NParticlesAP_nonROI;

Spot_SDROI = MS2Data.SDVectorAP_ROI{1,1};
Spot_SDnonROI = MS2Data.SDVectorAP_nonROI;

Spot_SEROI = Spot_SDROI./sqrt(Spot_NROI);
Spot_SEnonROI = Spot_SDnonROI./sqrt(Spot_NnonROI);

%% 
AP = 12; %out of 41 total bins

figure(1)
hold on
%plot(1:length(Spot_Time),Spot_ROI(:,AP))
%plot(1:length(Spot_Time),Spot_nonROI(:,AP))
errorbar(Spot_Time,Spot_ROI(:,AP),Spot_SEROI(:,AP))
errorbar(Spot_Time,Spot_nonROI(:,AP),Spot_SEnonROI(:,AP))
title(['MS2 Spot Fluorescence ','at AP =',num2str((AP-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Spot Fluorescence (AU)')
legend('ROI','Control')
standardizeFigure(gca,legend,[])
% figure(2)
% hold on
% %plot(1:length(Spot_Time),Spot_ROI(:,AP))
% %plot(1:length(Spot_Time),Spot_nonROI(:,AP))
% errorbar(Spot_Time,MS2Data.MeanVectorAnterior_ROI{1,1},MS2Data.SDVectorAnterior_ROI{1,1})
% errorbar(Spot_Time,MS2Data.MeanVectorAnterior_nonROI{1,1},MS2Data.SDVectorAnterior_nonROI{1,1})
% title(['MS2 Spot Fluorescence (MeanVectorAnterior(10~35% EL)'])
% xlabel('Time (min)')
% ylabel('Spot Fluorescence (AU)')
% legend('ROI','Control')

% figure(2)
% hold on
% %plot(1:length(Spot_Time),Spot_ROI(:,AP))
% %plot(1:length(Spot_Time),Spot_nonROI(:,AP))
% errorbar(Spot_Time,Spot_ROI(:,AP).*Spot_NROI(:,AP),Spot_SEROI(:,AP))
% errorbar(Spot_Time,Spot_nonROI(:,AP).*Spot_NnonROI(:,AP),Spot_SEnonROI(:,AP))
% title(['MS2 Spot Fluorescence ','at AP =',num2str((AP-1)*2.5),'%'])
% xlabel('Time (min)')
% ylabel('Number of Spots * Spot Fluorescence (AU)')
% legend('ROI','Control')

% figure(3)
% hold on
% %plot(1:length(Spot_Time),Spot_ROI(:,AP))
% %plot(1:length(Spot_Time),Spot_nonROI(:,AP))
% plot(Spot_Time,Spot_NROI(:,AP))
% plot(Spot_Time,Spot_NnonROI(:,AP))
% title(['Number of MS2 spots ','at AP =',num2str((AP-1)*2.5),'%'])
% xlabel('Time (min)')
% ylabel('Number of Spots')
% legend('ROI','Control')

%% Calculate the error bars (Nuclear fluo)
Nuclearfluo_ROI = RuntData.MeanVectorAP_ROI;
Nuclearfluo_nonROI = RuntData.MeanVectorAP_nonROI;
Nuclear_Time = RuntData.ElapsedTime;
Nuclear_Time = Nuclear_Time'

Nuclear_NROI = RuntData.NParticlesAP_ROI;
Nuclear_NnonROI = RuntData.NParticlesAP_nonROI;

Nuclear_SDROI = RuntData.SDVectorAP_ROI;
Nuclear_SDnonROI = RuntData.SDVectorAP_nonROI;

Nuclear_SEROI = Nuclear_SDROI./sqrt(Nuclear_NROI);
Nuclear_SEnonROI = Nuclear_SDnonROI./sqrt(Nuclear_NnonROI);

%% Plot the nuclear fluo
AP = 14;

figure(1) % Plotting the Nuc fluo and Spot fluo synchronously.

hold on
%plot(1:length(Spot_Time),Spot_ROI(:,AP))
%plot(1:length(Spot_Time),Spot_nonROI(:,AP))
% errorbar(Spot_Time,Spot_ROI(:,AP).*Spot_NROI(:,AP),Spot_SEROI(:,AP),'b')
% errorbar(Spot_Time,Spot_nonROI(:,AP).*Spot_NnonROI(:,AP),Spot_SEnonROI(:,AP),'r')
% title(['Nuclear and MS2 Spot Fluorescence ','at AP =',num2str((AP-1)*2.5),'%'])
 xlabel('Time (min)')
% yyaxis left
% ylabel('Spot Fluorescence (AU)')
%plot(1:length(Spot_Time),Spot_ROI(:,AP))
%plot(1:length(Spot_Time),Spot_nonROI(:,AP))
%Bcd scale
RuntScale = 10;
errorbar(Nuclear_Time,RuntScale*Nuclearfluo_ROI(:,AP),RuntScale*Nuclear_SEROI(:,AP))
errorbar(Nuclear_Time,RuntScale*Nuclearfluo_nonROI(:,AP),RuntScale*Nuclear_SEnonROI(:,AP))
title(['Nuclear Fluorescence ','at AP =',num2str((AP-1)*2.5),'%'])
%xlabel('Time (min)')
%yyaxis right
ylabel('Nuclear Fluorescence (AU)')
%legend('ROI-Spot','CTL-Spot','ROI-Nuc','CTL-Nuc')
legend('ROI-Nuc','CTL-Nuc')

%% Compare the nuclear fluorescence of ROI and non-ROI
hold on
for i=1:length(CompiledNuclei_ROI)
    plot(CompiledNuclei_ROI(i).Frames,CompiledNuclei_ROI(i).FluoMax,'r')
end

for i=1:length(CompiledNuclei_nonROI)
    plot(CompiledNuclei_nonROI(i).Frames,CompiledNuclei_nonROI(i).FluoMax,'b')
end

%% Compare the spot fluorescence - offset

%% Compare the MeanVectorAll
%% Get the Autofluorescence in nucleus and cytoplasm
AutoFluomCi = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-07-04-2xMCP-mCherry-Autofluo-mCi-ROI\CompiledNuclei.mat');

AutoFluoTime = AutoFluomCi.ElapsedTime;
AutoFluo_ROI = AutoFluomCi.MeanVectorAP_ROI;
AutoFluo_ROI_SE = AutoFluomCi.SDVectorAP_ROI./ sqrt(AutoFluomCi.NParticlesAP_ROI);
AutoFluo_nonROI = AutoFluomCi.MeanVectorAP_nonROI;
AutoFluo_nonROI_SE = AutoFluomCi.SDVectorAP_nonROI./ sqrt(AutoFluomCi.NParticlesAP_nonROI);
%% Plot
AP = 15;
hold on
errorbar(AutoFluoTime, AutoFluo_ROI(:,AP), AutoFluo_ROI_SD(:,AP))
errorbar(AutoFluoTime, AutoFluo_nonROI(:,AP), AutoFluo_nonROI_SD(:,AP))

title('Nuclear Autofluorescence over time (ROI vs nonROI)')
xlabel('Time (min)')
ylabel('Nuclear Fluorescence (AU)')
legend('ROI','non-ROI')