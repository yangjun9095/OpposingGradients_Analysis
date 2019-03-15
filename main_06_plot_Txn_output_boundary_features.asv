function main_06_plot_Txn_output_boundary_features
% This is edited from "OpposingGradientAnalysis_combined.m"
% Description : This is a collection of sub-codes to generate plots of
% analysis
% INPUTS :
% 1) Bcd, and Runt Datasets
% 2) hbP2-r0,1,2,3 MS2.V5 datasets (either averaged or single individual
% spots, etc.
% 3) 
%Combine all the info about 1)Bicoid, 2)Runt, and 3) hbP2-r0,1,2,3-MS2 in
% one plot. 

%% Load the datasets
% We will start from plotting those info along AP axis, over time (nc 14)
% First, Bicoid
Bcd = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Bcd-Averaged.mat')

% Second, Runt
RuntData = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-1min-200Hz-Female-Averaged.mat')

% Third, hbP2-r0,1,2,3-MS2.V5
% first, run the code for averaging the datasets (Note that this
% AverageDatasets.m is averaging only nc13 and nc 14 (These are saved in
% /OpposingGradient folder)
AverageDatasets('r0','savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient')
AverageDatasets('r1','savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient')
AverageDatasets('r2','savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient')
AverageDatasets('r3','savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient')
% Load the averaged rN datasets
% 1) from nc13 to nc 14 (end)
% r0Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGrad ient\Data_Processed\r0new.mat');
% r1Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\r1new.mat');
% r2Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\r2new.mat');
% r3Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\r3new.mat');

% 2) from nc12 to nc 14
% r0Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0_FromNC12.mat');
% r1Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1_FromNC12.mat');
% r2Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2_FromNC12.mat');
% r3Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3_FromNC12.mat');

r0Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0.mat');
r1Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1.mat');
r2Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2.mat');
r3Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3.mat');

%% Bcd
BcdTime = BcdData.ElapsedTime;
BcdTime = BcdData.ElapsedTime-BcdData.ElapsedTime(BcdData.nc13);

%% Runt
RuntTime = RuntData.ElapsedTime;
RuntTime = RuntData.ElapsedTime-RuntData.ElapsedTime(RuntData.nc14);

%% r0
r0Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0_FromNC12.mat')
Time_r0 = r0Data.ElapsedTime;
MeanVectorAP_r0 = r0Data.MeanVectorAP;
SDVectorAP_r0 = r0Data.SDVectorAP;
SEVectorAP_r0 = r0Data.SEVectorAP;
NParticlesAP_r0 = r0Data.NParticlesAP;
AccumulatedmRNA_r0 = r0Data.AccumulatedmRNA;
SDAccumulatedmRNA_r0 = r0Data.AccumulatedmRNA_SD;
AccumulatedmRNA_FractionON_r0 = r0Data.AccumulatedmRNA_FractionON;

% Define the time points
%nc12_r0 = r0Data.nc12;
nc13_r0 = r0Data.nc13;
nc14_r0 = r0Data.nc14;

%% r1
r1Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1_FromNC12.mat')
Time_r1 = r1Data.ElapsedTime;
MeanVectorAP_r1 = r1Data.MeanVectorAP;
SDVectorAP_r1 = r1Data.SDVectorAP;
SEVectorAP_r1 = r1Data.SEVectorAP;
NParticlesAP_r1 = r1Data.NParticlesAP;
AccumulatedmRNA_r1 = r1Data.AccumulatedmRNA;
SDAccumulatedmRNA_r1 = r1Data.AccumulatedmRNA_SD;
AccumulatedmRNA_FractionON_r1 = r1Data.AccumulatedmRNA_FractionON;

nc12_r1 = r1Data.nc12;
nc13_r1 = r1Data.nc13;
nc14_r1 = r1Data.nc14;

%% r2
r2Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2_FromNC12.mat')
Time_r2 = r2Data.ElapsedTime;
MeanVectorAP_r2 = r2Data.MeanVectorAP;
SDVectorAP_r2 = r2Data.SDVectorAP;
SEVectorAP_r2 = r2Data.SEVectorAP;
NParticlesAP_r2 = r2Data.NParticlesAP;
AccumulatedmRNA_r2 = r2Data.AccumulatedmRNA;
SDAccumulatedmRNA_r2 = r2Data.AccumulatedmRNA_SD;
AccumulatedmRNA_FractionON_r2 = r2Data.AccumulatedmRNA_FractionON;

nc12_r2 = r2Data.nc12;
nc13_r2 = r2Data.nc13;
nc14_r2 = r2Data.nc14;
%% r3
r3Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3_FromNC12.mat')
Time_r3 = r3Data.ElapsedTime;
MeanVectorAP_r3 = r3Data.MeanVectorAP;
SDVectorAP_r3 = r3Data.SDVectorAP;
SEVectorAP_r3 = r3Data.SEVectorAP;
NParticlesAP_r3 = r3Data.NParticlesAP;
AccumulatedmRNA_r3 = r3Data.AccumulatedmRNA;
SDAccumulatedmRNA_r3 = r3Data.AccumulatedmRNA_SD;
AccumulatedmRNA_FractionON_r3 = r3Data.AccumulatedmRNA_FractionON;

nc12_r3 = r3Data.nc12;
nc13_r3 = r3Data.nc13;
nc14_r3 = r3Data.nc14;

%% r3 prime (r3', scrambled Runt binding sites)
r3primeData = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3prime_FromNC12.mat')
Time_r3prime = r3primeData.ElapsedTime;
MeanVectorAP_r3prime = r3primeData.MeanVectorAP;
SDVectorAP_r3prime = r3primeData.SDVectorAP;
SEVectorAP_r3prime = r3primeData.SEVectorAP;
NParticlesAP_r3prime = r3primeData.NParticlesAP;
AccumulatedmRNA_r3prime = r3primeData.AccumulatedmRNA;
SDAccumulatedmRNA_r3prime = r3primeData.AccumulatedmRNA_SD;
AccumulatedmRNA_FractionON_r3prime = r3primeData.AccumulatedmRNA_FractionON;

nc12_r3prime = r3primeData.nc12;
nc13_r3prime = r3primeData.nc13;
nc14_r3prime = r3primeData.nc14;
%% Section 1.
% Plotting MeanVectorAP and AccumulatedmRNA over AP @ different time
% points. Just looking at how the trend looks like.
%% 1) MS2 Spot fluo (RNAP loading rate) along AP axis (over time)
% This will generate a movie of r0,1,2,3 - MeanVectorAP over AP at
% different time points.

% I want to plot the averaged MS2 spot fluorescence over AP axis (also over
% time, nc13 and nc 14). Assume that the datasets are taken with the same
% frame rate.
% Location to save the plots
Folder = 'E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\hbP2rN-SpotFluo'

% Get the longest time vector
[tLength,Index] = min([length(Time_r0),length(Time_r1),length(Time_r2),length(Time_r3),length(Time_r3prime)]);

for i=1:tLength
    clf
    hold on
    errorbar(0:0.025:1,MeanVectorAP_r0(i,:),SEVectorAP_r0(i,:))
    errorbar(0:0.025:1,MeanVectorAP_r1(i,:),SEVectorAP_r1(i,:))
    errorbar(0:0.025:1,MeanVectorAP_r2(i,:),SEVectorAP_r2(i,:))
    errorbar(0:0.025:1,MeanVectorAP_r3(i,:),SEVectorAP_r3(i,:))
    errorbar(0:0.025:1,MeanVectorAP_r3prime(i,:),SEVectorAP_r3prime(i,:))
    
    ylim([0 1200])
    title(['MS2 Spot Fluorescence @ ',num2str(floor(Time_r0(i))),'min from nc13'])
    xlabel('AP')
    ylabel('MS2 Spot Fluorescence (AU)')
    legend('r0','r1','r2','r3','r3-scr')
    standardizeFigure(gca,legend,'g','r','yellow','cyan')
    fig = gcf;
    fig.InvertHardcopy = 'off';
    pause
    %saveas(gcf, [Folder, '\', 'Time',num2str(floor(Time_r0(i))),'min', '.tif'])%[nc14Folder '\' num2str(nucleus) '.png'])

end

%% 2) mRNA Accumulation along AP axis (over time)

% Location to save the plots
Folder = 'E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\hbP2rN-AccumulatedmRNA'
% Get the longest time vector
[tLength,Index] = min([length(Time_r0),length(Time_r1),length(Time_r2),length(Time_r3),length(Time_r3prime)]);

% AccumulatedmRNA_FractionON_r0(AccumulatedmRNA_FractionON_r0==0) = nan;
% AccumulatedmRNA_FractionON_r1(AccumulatedmRNA_FractionON_r1==0) = nan;
% AccumulatedmRNA_FractionON_r2(AccumulatedmRNA_FractionON_r2==0) = nan;
% AccumulatedmRNA_FractionON_r3(AccumulatedmRNA_FractionON_r3==0) = nan;

for i=1:tLength
    clf
    hold on
    %errorbar(0:0.025:1,AccumulatedmRNA_r0(i,:),SDAccumulatedmRNA_r0(i,:))
    %errorbar(0:0.025:1,AccumulatedmRNA_r1(i,:),SDAccumulatedmRNA_r1(i,:))
    %errorbar(0:0.025:1,AccumulatedmRNA_r2(i,:),SDAccumulatedmRNA_r2(i,:))
    %errorbar(0:0.025:1,AccumulatedmRNA_r3(i,:),SDAccumulatedmRNA_r3(i,:))
    
    plot(0:0.025:1,AccumulatedmRNA_FractionON_r0(i,:))
    plot(0:0.025:1,AccumulatedmRNA_FractionON_r1(i,:))
    plot(0:0.025:1,AccumulatedmRNA_FractionON_r2(i,:))
    plot(0:0.025:1,AccumulatedmRNA_FractionON_r3(i,:))
    plot(0:0.025:1,AccumulatedmRNA_FractionON_r3prime(i,:))
    %ylim([0 30000])
    title({'Accumulated mRNA',[' @ ',num2str(floor(Time_r0(i))),'min from nc13']})
    xlabel('AP')
    ylabel('Accumulated mRNA (AU)')
    legend('r0','r1','r2','r3','r3-scr')
    standardizeFigure(gca,legend,'g','r','yellow','cyan')
    %'g','r','yellow','cyan'
    fig = gcf;
    fig.InvertHardcopy = 'off';
    pause
    %saveas(gcf, [Folder, '\', 'Time',num2str(floor(Time_r0(i))),'min', '.tif'])%[nc14Folder '\' num2str(nucleus) '.png'])
end
%% Section 2. Boundary features - Averaged Embryos
%% 1) Get the boundary profile of Accumulated mRNA
% Here, I will use FindPosition.m script that I wrote to find
% xhalf,Width,and Slope. [xhalf,Width,Slope]=FindPosition(Profile,tpoint)

[Xhalf_r0,Width_r0] = CalculateBoundaryFeatures(AccumulatedmRNA_FractionON_r0,...
                                                        Time_r0,...
                                                        nc12_r0,nc13_r0,nc14_r0);
[Xhalf_r1,Width_r1] = CalculateBoundaryFeatures(AccumulatedmRNA_FractionON_r1,...
                                                        Time_r1,...
                                                        nc12_r1,nc13_r1,nc14_r1);
[Xhalf_r2,Width_r2] = CalculateBoundaryFeatures(AccumulatedmRNA_FractionON_r2,...
                                                        Time_r2,...
                                                        nc12_r2,nc13_r2,nc14_r2);
[Xhalf_r3,Width_r3] = CalculateBoundaryFeatures(AccumulatedmRNA_FractionON_r3,...
                                                        Time_r3,...
                                                        nc12_r3,nc13_r3,nc14_r3);        
% [Xhalf_r3prime,Width_r3prime] = CalculateBoundaryFeatures(AccumulatedmRNA_FractionON_r3prime,...
%                                                         Time_r3prime,...
%                                                         nc12_r3prime,nc13_r3prime,nc14_r3prime);                                                     
                                                    
% for i=1:tLength
%     [Xhalf_r0(i),Width_r0(i),Slope_r0(i)] = FindPosition(AccumulatedmRNA_FractionON_r0,i);
%     pause(0.1)
%     [Xhalf_r1(i),Width_r1(i),Slope_r1(i)] = FindPosition(AccumulatedmRNA_FractionON_r1,i);
%     [Xhalf_r2(i),Width_r2(i),Slope_r2(i)] = FindPosition(AccumulatedmRNA_FractionON_r2,i);
%     [Xhalf_r3(i),Width_r3(i),Slope_r3(i)] = FindPosition(AccumulatedmRNA_FractionON_r3,i);
% end

%% Plot the Boundary morphology - boundary position
hold on
plot(Time_r0,Xhalf_r0)
plot(Time_r1,Xhalf_r1)
plot(Time_r2,Xhalf_r2)
plot(Time_r3,Xhalf_r3)


title('Boundary position over Time')
xlabel('Time (min)')
ylabel('Boundary position (AP)')
legend('r0','r1','r2','r3')
standardizeFigure(gca,legend,[])
fig = gcf;
fig.InvertHardcopy = 'off';

plot([Time_r0(nc13) Time_r0(nc13)],[0 0.5],'k','LineWidth',5)
plot([Time_r0(nc14) Time_r0(nc14)],[0 0.5],'k','LineWidth',5)

%saveas(gcf, ['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient', '\', 'Boundary_Position_hbP2-r0123_nc1314.tif'])%[nc14Folder '\' num2str(nucleus) '.png'])

%% Plot the Boundary width
hold on
plot(Time_r0,Width_r0)
plot(Time_r1,Width_r1)
plot(Time_r2,Width_r2)
plot(Time_r3,Width_r3)

title('Boundary width over Time')
xlabel('Time (min)')
ylabel('Boundary width')
legend('r0','r1','r2','r3')
standardizeFigure(gca,legend,[])
fig = gcf;
fig.InvertHardcopy = 'off';
plot([Time_r0(nc13) Time_r0(nc13)],[0 0.5],'k','LineWidth',5)
plot([Time_r0(nc14) Time_r0(nc14)],[0 0.5],'k','LineWidth',5)
%saveas(gcf, ['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient', '\', 'Boundary_Width_hbP2-r0123_nc1314.tif'])%[nc14Folder '\' num2str(nucleus) '.png'])

%% Section 3. Boundary characterization (total mRNA) - single, individual embryos

%% r0 dataset - individual embryos
% For now, I will just create different tabs in the DataStatus.xls.
% In future, there should be a better way to deal with multiple datasets.

%% Start with r0 datasets
% Calculate the AccumulatedmRNA, and also synchronize as the beginning of
% nc12
AverageDatasets('r0-1','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets('r0-2','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets('r0-3','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');

%% Load the datasets - r0
r0Data(1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0-1_FromNC12.mat');
r0Data(2) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0-2_FromNC12.mat');
r0Data(3) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0-3_FromNC12.mat');

[Xhalf_r0_1,Width_r0_1] = CalculateBoundary_SpotFluo(r0Data(1).AccumulatedmRNA_FractionON,...
                                                   r0Data(1).ElapsedTime,...
                                                    r0Data(1).nc12,r0Data(1).nc13,r0Data(1).nc14);
                                                
[Xhalf_r0_2,Width_r0_2] = CalculateBoundary_SpotFluo(r0Data(2).AccumulatedmRNA_FractionON,...
                                                   r0Data(2).ElapsedTime,...
                                                    r0Data(2).nc12,r0Data(2).nc13,r0Data(2).nc14); 
                                                        
[Xhalf_r0_3,Width_r0_3] = CalculateBoundary_SpotFluo(r0Data(3).AccumulatedmRNA_FractionON,...
                                                   r0Data(3).ElapsedTime,...
                                                    r0Data(3).nc12,r0Data(3).nc13,r0Data(3).nc14);      
                                                
%% Plot to quickly check - r0
hold on
plot(Xhalf_r0_1)
plot(Xhalf_r0_2)
plot(Xhalf_r0_3)

%% Plot the boundary position - r0 (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r0_1),length(Xhalf_r0_2),length(Xhalf_r0_3)]);
Xhalf_r0_1 = Xhalf_r0_1(1,1:MinLength);
Xhalf_r0_2 = Xhalf_r0_2(1,1:MinLength)
Xhalf_r0_3 = Xhalf_r0_3(1,1:MinLength)
Time_r0_Min = Time_r0(1:MinLength);

MeanXhalf_r0 = nanmean([Xhalf_r0_1;Xhalf_r0_2;Xhalf_r0_3]);
STDXhalf_r0 = nanstd([Xhalf_r0_1;Xhalf_r0_2;Xhalf_r0_3]);

hold on
errorbar(Time_r0_Min,MeanXhalf_r0,STDXhalf_r0)
plot(Time_r0,Xhalf_r0)
plot([Time_r0(nc13_r0) Time_r0(nc13_r0)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r0(nc14_r0) Time_r0(nc14_r0)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary position along Time')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')

%% Plot the boundary width - r0 (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r0_1),length(Xhalf_r0_2),length(Xhalf_r0_3)]);
Width_r0_1 = Width_r0_1(1,1:MinLength);
Width_r0_2 = Width_r0_2(1,1:MinLength)
Width_r0_3 = Width_r0_3(1,1:MinLength)
Time_r0_Min = Time_r0(1:MinLength);

MeanWidth_r0 = nanmean([Width_r0_1;Width_r0_2;Width_r0_3]);
STDWidth_r0 = nanstd([Width_r0_1;Width_r0_2;Width_r0_3]);

hold on
errorbar(Time_r0_Min,MeanWidth_r0,STDWidth_r0)
plot(Time_r0,Width_r0)
plot([Time_r0(nc13_r0) Time_r0(nc13_r0)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r0(nc14_r0) Time_r0(nc14_r0)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary width along Time')
xlabel('Time')
ylabel('Boundary width(AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')

%% Individual embryos - r1 datasets
% Calculate the AccumulatedmRNA, and also synchronize as the beginning of
% nc12
AverageDatasets_FineTuning('r1-1','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets_FineTuning('r1-2','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets_FineTuning('r1-3','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');

%% Load the datasets - r1
r1Data(1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1-1_FromNC12.mat');
r1Data(2) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1-2_FromNC12.mat');
r1Data(3) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1-3_FromNC12.mat');

[Xhalf_r1_1,Width_r1_1] = CalculateBoundary_SpotFluo(r1Data(1).AccumulatedmRNA_FractionON,...
                                                   r1Data(1).ElapsedTime,...
                                                    r1Data(1).nc12,r1Data(1).nc13,r1Data(1).nc14);
                                                
[Xhalf_r1_2,Width_r1_2] = CalculateBoundary_SpotFluo(r1Data(2).AccumulatedmRNA_FractionON,...
                                                   r1Data(2).ElapsedTime,...
                                                    r1Data(2).nc12,r1Data(2).nc13,r1Data(2).nc14); 
                                                        
[Xhalf_r1_3,Width_r1_3] = CalculateBoundary_SpotFluo(r1Data(3).AccumulatedmRNA_FractionON,...
                                                   r1Data(3).ElapsedTime,...
                                                    r1Data(3).nc12,r1Data(3).nc13,r1Data(3).nc14);
                                                
%% Plot the boundary position (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r1_1),length(Xhalf_r1_2),length(Xhalf_r1_3)]);
Xhalf_r1_1 = Xhalf_r1_1(1,1:MinLength);
Xhalf_r1_2 = Xhalf_r1_2(1,1:MinLength)
Xhalf_r1_3 = Xhalf_r1_3(1,1:MinLength)
Time_r1_Min = Time_r1(1:MinLength);

MeanXhalf_r1 = nanmean([Xhalf_r1_1;Xhalf_r1_2;Xhalf_r1_3]);
STDXhalf_r1 = nanstd([Xhalf_r1_1;Xhalf_r1_2;Xhalf_r1_3]);

hold on
errorbar(Time_r1_Min,MeanXhalf_r1,STDXhalf_r1)
plot(Time_r1,Xhalf_r1)
plot([Time_r1(nc13_r1) Time_r1(nc13_r1)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r1(nc14_r1) Time_r1(nc14_r1)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary position along Time')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Plot the boundary width (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r1_1),length(Xhalf_r1_2),length(Xhalf_r1_3)]);
Width_r1_1 = Width_r1_1(1,1:MinLength);
Width_r1_2 = Width_r1_2(1,1:MinLength)
Width_r1_3 = Width_r1_3(1,1:MinLength)
Time_r1_Min = Time_r1(1:MinLength);

MeanWidth_r1 = nanmean([Width_r1_1;Width_r1_2;Width_r1_3]);
STDWidth_r1 = nanstd([Width_r1_1;Width_r1_2;Width_r1_3]);

hold on
errorbar(Time_r1_Min,MeanWidth_r1,STDWidth_r1)
plot(Time_r1,Width_r1)
plot([Time_r1(nc13_r1) Time_r1(nc13_r1)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r1(nc14_r1) Time_r1(nc14_r1)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary width along Time')
xlabel('Time')
ylabel('Boundary width(AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Individual embryos vs averaged embryos : r2 datasets
% Calculate the AccumulatedmRNA, and also synchronize as the beginning of
% nc12
AverageDatasets_FineTuning('r2-1','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets_FineTuning('r2-2','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets_FineTuning('r2-3','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');

%% Load the datasets - r2
r2Data(1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2-1_FromNC12.mat');
r2Data(2) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2-2_FromNC12.mat');
r2Data(3) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2-3_FromNC12.mat');

[Xhalf_r2_1,Width_r2_1] = CalculateBoundary_SpotFluo(r2Data(1).AccumulatedmRNA_FractionON,...
                                                   r2Data(1).ElapsedTime,...
                                                    r2Data(1).nc12,r2Data(1).nc13,r2Data(1).nc14);
                                                
[Xhalf_r2_2,Width_r2_2] = CalculateBoundary_SpotFluo(r2Data(2).AccumulatedmRNA_FractionON,...
                                                   r2Data(2).ElapsedTime,...
                                                    r2Data(2).nc12,r2Data(2).nc13,r2Data(2).nc14); 
                                                        
[Xhalf_r2_3,Width_r2_3] = CalculateBoundary_SpotFluo(r2Data(3).AccumulatedmRNA_FractionON,...
                                                   r2Data(3).ElapsedTime,...
                                                    r2Data(3).nc12,r2Data(3).nc13,r2Data(3).nc14);
                                                
%% Plot the boundary position (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r2_1),length(Xhalf_r2_2),length(Xhalf_r2_3)]);
Xhalf_r2_1 = Xhalf_r2_1(1,1:MinLength);
Xhalf_r2_2 = Xhalf_r2_2(1,1:MinLength)
Xhalf_r2_3 = Xhalf_r2_3(1,1:MinLength)
Time_r2_Min = Time_r2(1:MinLength);

MeanXhalf_r2 = nanmean([Xhalf_r2_1;Xhalf_r2_2;Xhalf_r2_3]);
STDXhalf_r2 = nanstd([Xhalf_r2_1;Xhalf_r2_2;Xhalf_r2_3]);

hold on
errorbar(Time_r2_Min,MeanXhalf_r2,STDXhalf_r2)
plot(Time_r2,Xhalf_r2)
plot([Time_r2(nc13_r2) Time_r2(nc13_r2)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r2(nc14_r2) Time_r2(nc14_r2)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary position along Time')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Plot the boundary width (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r2_1),length(Xhalf_r2_2),length(Xhalf_r2_3)]);
Width_r2_1 = Width_r2_1(1,1:MinLength);
Width_r2_2 = Width_r2_2(1,1:MinLength)
Width_r2_3 = Width_r2_3(1,1:MinLength)
Time_r2_Min = Time_r2(1:MinLength);

MeanWidth_r2 = nanmean([Width_r2_1;Width_r2_2;Width_r2_3]);
STDWidth_r2 = nanstd([Width_r2_1;Width_r2_2;Width_r2_3]);

hold on
errorbar(Time_r2_Min,MeanWidth_r2,STDWidth_r2)
plot(Time_r2,Width_r2)
plot([Time_r2(nc13_r2) Time_r2(nc13_r2)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r2(nc14_r2) Time_r2(nc14_r2)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary width along Time')
xlabel('Time')
ylabel('Boundary width(AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')

%% Individual embryos vs averaged embryos : r3 datasets
% Calculate the AccumulatedmRNA, and also synchronize as the beginning of
% nc12
AverageDatasets_FineTuning('r3-1','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets_FineTuning('r3-2','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets_FineTuning('r3-3','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');

%% Load the datasets - r3
r3Data(1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3-1_FromNC12.mat');
r3Data(2) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3-2_FromNC12.mat');
r3Data(3) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3-3_FromNC12.mat');

[Xhalf_r3_1,Width_r3_1] = CalculateBoundary_SpotFluo(r3Data(1).AccumulatedmRNA_FractionON,...
                                                   r3Data(1).ElapsedTime,...
                                                    r3Data(1).nc12,r3Data(1).nc13,r3Data(1).nc14);
                                                
[Xhalf_r3_2,Width_r3_2] = CalculateBoundary_SpotFluo(r3Data(2).AccumulatedmRNA_FractionON,...
                                                   r3Data(2).ElapsedTime,...
                                                    r3Data(2).nc12,r3Data(2).nc13,r3Data(2).nc14); 
                                                        
[Xhalf_r3_3,Width_r3_3] = CalculateBoundary_SpotFluo(r3Data(3).AccumulatedmRNA_FractionON,...
                                                   r3Data(3).ElapsedTime,...
                                                    r3Data(3).nc12,r3Data(3).nc13,r3Data(3).nc14);
                                                
%% Plot the boundary position (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r3_1),length(Xhalf_r3_2),length(Xhalf_r3_3)]);
Xhalf_r3_1 = Xhalf_r3_1(1,1:MinLength);
Xhalf_r3_2 = Xhalf_r3_2(1,1:MinLength)
Xhalf_r3_3 = Xhalf_r3_3(1,1:MinLength)
Time_r3_Min = Time_r3(1:MinLength);

MeanXhalf_r3 = nanmean([Xhalf_r3_1;Xhalf_r3_2;Xhalf_r3_3]);
STDXhalf_r3 = nanstd([Xhalf_r3_1;Xhalf_r3_2;Xhalf_r3_3]);

hold on
errorbar(Time_r3_Min,MeanXhalf_r3,STDXhalf_r3)
plot(Time_r3,Xhalf_r3)
plot([Time_r3(nc13_r3) Time_r3(nc13_r3)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r3(nc14_r3) Time_r3(nc14_r3)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary position along Time')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Plot the boundary width (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r3_1),length(Xhalf_r3_2),length(Xhalf_r3_3)]);
Width_r3_1 = Width_r3_1(1,1:MinLength);
Width_r3_2 = Width_r3_2(1,1:MinLength)
Width_r3_3 = Width_r3_3(1,1:MinLength)
Time_r3_Min = Time_r3(1:MinLength);

MeanWidth_r3 = nanmean([Width_r3_1;Width_r3_2;Width_r3_3]);
STDWidth_r3 = nanstd([Width_r3_1;Width_r3_2;Width_r3_3]);

hold on
errorbar(Time_r3_Min,MeanWidth_r3,STDWidth_r3)
plot(Time_r3,Width_r3)
plot([Time_r3(nc13_r3) Time_r3(nc13_r3)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r3(nc14_r3) Time_r3(nc14_r3)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary width along Time')
xlabel('Time')
ylabel('Boundary width(AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Plot all boundary positions (r0,1,2,3)
hold on
plot([Time_r0(nc13_r0) Time_r0(nc13_r0)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r0(nc14_r0) Time_r0(nc14_r0)],[0 0.6],'k')%,'LineWidth',5)

Nstep = 3;
errorbar(Time_r0_Min(1:Nstep:end),MeanXhalf_r0(1:Nstep:end),STDXhalf_r0(1:Nstep:end),'capsize',0)%,'b','LineWidth',2)
errorbar(Time_r1_Min(1:Nstep:end),MeanXhalf_r1(1:Nstep:end),STDXhalf_r1(1:Nstep:end),'capsize',0)%,'r','LineWidth',2)
errorbar(Time_r2_Min(1:Nstep:end),MeanXhalf_r2(1:Nstep:end),STDXhalf_r2(1:Nstep:end),'capsize',0)%,'y','LineWidth',2)
errorbar(Time_r3_Min(1:Nstep:end),MeanXhalf_r3(1:Nstep:end),STDXhalf_r3(1:Nstep:end),'capsize',0)%,'g','LineWidth',2)

%plot(Time_r0,Xhalf_r0)%,'b','LineWidth',3)
%plot(Time_r1,Xhalf_r1)%,'r','LineWidth',3)
%plot(Time_r2,Xhalf_r2)%,'y','LineWidth',3)
%plot(Time_r3,Xhalf_r3)%,'g','LineWidth',3)


ylim([0 0.6])
title('Boundary position along Time-Accumulated mRNA')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('nc13','nc14','r0','r1','r2','r3')%,...
        %'r0-pooled','r1-pooled','r2-pooled','r3-pooled')
standardizeFigure_YJK(gca,legend,'red','yellow','cyan','magenta','black','black')%'red','yellow','cyan','magenta',
%% Save values
% 1) Mean Xhalf(mRNA)
MeanXhalf_r0_mRNA = MeanXhalf_r0;
MeanXhalf_r1_mRNA = MeanXhalf_r1;
MeanXhalf_r2_mRNA = MeanXhalf_r2;
MeanXhalf_r3_mRNA = MeanXhalf_r3;

% 2) STD Xhalf (mRNA)
STDXhalf_r0_mRNA = STDXhalf_r0;
STDXhalf_r1_mRNA = STDXhalf_r1;
STDXhalf_r2_mRNA = STDXhalf_r2;
STDXhalf_r3_mRNA = STDXhalf_r3;
%% Plot all boundary widths (r0,1,2,3)
hold on
plot([Time_r0(nc13_r0) Time_r0(nc13_r0)],[0 0.5],'k')
plot([Time_r0(nc14_r0) Time_r0(nc14_r0)],[0 0.5],'k')

Nstep = 3;
errorbar(Time_r0_Min(1:Nstep:end),MeanWidth_r0(1:Nstep:end),STDWidth_r0(1:Nstep:end),'capsize',0)
errorbar(Time_r1_Min(1:Nstep:end),MeanWidth_r1(1:Nstep:end),STDWidth_r1(1:Nstep:end),'capsize',0)
errorbar(Time_r2_Min(1:Nstep:end),MeanWidth_r2(1:Nstep:end),STDWidth_r2(1:Nstep:end),'capsize',0)
errorbar(Time_r3_Min(1:Nstep:end),MeanWidth_r3(1:Nstep:end),STDWidth_r3(1:Nstep:end),'capsize',0)

%plot(Time_r0,Width_r0)
%plot(Time_r1,Width_r1)
%plot(Time_r2,Width_r2)
%plot(Time_r3,Width_r3)


ylim([0 0.5])
title('Boundary width along Time')
xlabel('Time')
ylabel('Boundary width (AP)')
legend('nc13','nc14','r0','r1','r2','r3')%,...
        %'r0-pooled','r1-pooled','r2-pooled','r3-pooled')
standardizeFigure_YJK(gca,legend,'red','yellow','cyan','magenta','black','black')%'red','yellow','cyan','magenta',

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 4. Boundary features for MS2 Spot fluorescence (RNAP loading rate) instead of total mRNA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Boundary features of RNAP loading rates (Averaged embryos)
% Here, I will use FindPosition.m scrit that I wrote to find
% xhalf,Width,and Slope. [xhalf,Width,Slope]=FindPosition(Profile,tpoint)

[Xhalf_r0,Width_r0] = CalculateBoundary_SpotFluo(MeanVectorAP_r0,...
                                                        Time_r0,...
                                                        nc12_r0,nc13_r0,nc14_r0);
[Xhalf_r1,Width_r1] = CalculateBoundary_SpotFluo(MeanVectorAP_r1,...
                                                        Time_r1,...
                                                        nc12_r1,nc13_r1,nc14_r1);
[Xhalf_r2,Width_r2] = CalculateBoundary_SpotFluo(MeanVectorAP_r2,...
                                                        Time_r2,...
                                                        nc12_r2,nc13_r2,nc14_r2);
[Xhalf_r3,Width_r3] = CalculateBoundary_SpotFluo(MeanVectorAP_r3,...
                                                        Time_r3,...
                                                        nc12_r3,nc13_r3,nc14_r3);                                                 
%% Plot the Boundary morphology - boundary position
hold on
plot(Time_r0,Xhalf_r0)
plot(Time_r1,Xhalf_r1)
plot(Time_r2,Xhalf_r2)
plot(Time_r3,Xhalf_r3)


title('Boundary position over Time')
xlabel('Time (min)')
ylabel('Boundary position (AP)')
legend('r0','r1','r2','r3')
standardizeFigure(gca,legend,[])
fig = gcf;
fig.InvertHardcopy = 'off';

% plot([Time_r0(r0Data.nc13) Time_r0(r0Data.nc13)],[0 0.6],'k','LineWidth',5)
% plot([Time_r0(r0Data.nc14) Time_r0(r0Data.nc14)],[0 0.6],'k','LineWidth',5)
standardizeFigure(gca,legend,[])
%saveas(gcf, ['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient', '\', 'Boundary_Position_hbP2-r0123_nc1314.tif'])%[nc14Folder '\' num2str(nucleus) '.png'])

%% Plot the Boundary width
hold on
plot(Time_r0,Width_r0)
plot(Time_r1,Width_r1)
plot(Time_r2,Width_r2)
plot(Time_r3,Width_r3)

title('Boundary width over Time')
xlabel('Time (min)')
ylabel('Boundary width')
legend('r0','r1','r2','r3')
standardizeFigure(gca,legend,[])
fig = gcf;
fig.InvertHardcopy = 'off';
plot([Time_r0(r0Data.nc13) Time_r0(r0Data.nc13)],[0 0.5],'k','LineWidth',5)
plot([Time_r0(r0Data.nc14) Time_r0(r0Data.nc14)],[0 0.5],'k','LineWidth',5)
%saveas(gcf, ['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient', '\', 'Boundary_Width_hbP2-r0123_nc1314.tif'])%[nc14Folder '\' num2str(nucleus) '.png'])

%% Section 5. Boundary features from individual embryos 
%% Boundary features from individual embryos - MeanVectorAP (r0)
r0Data(1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0-1_FromNC12.mat');
r0Data(2) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0-2_FromNC12.mat');
r0Data(3) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0-3_FromNC12.mat');

[Xhalf_r0_1,Width_r0_1] = CalculateBoundary_SpotFluo(r0Data(1).MeanVectorAP,...
                                                   r0Data(1).ElapsedTime,...
                                                    r0Data(1).nc12,r0Data(1).nc13,r0Data(1).nc14);
                                                
[Xhalf_r0_2,Width_r0_2] = CalculateBoundary_SpotFluo(r0Data(2).MeanVectorAP,...
                                                   r0Data(2).ElapsedTime,...
                                                    r0Data(2).nc12,r0Data(2).nc13,r0Data(2).nc14); 
                                                        
[Xhalf_r0_3,Width_r0_3] = CalculateBoundary_SpotFluo(r0Data(3).MeanVectorAP,...
                                                   r0Data(3).ElapsedTime,...
                                                    r0Data(3).nc12,r0Data(3).nc13,r0Data(3).nc14);
                                                
%% Plot the boundary position (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r0_1),length(Xhalf_r0_2),length(Xhalf_r0_3)]);
Xhalf_r0_1 = Xhalf_r0_1(1,1:MinLength);
Xhalf_r0_2 = Xhalf_r0_2(1,1:MinLength)
Xhalf_r0_3 = Xhalf_r0_3(1,1:MinLength)
Time_r0_Min = Time_r0(1:MinLength);

MeanXhalf_r0 = nanmean([Xhalf_r0_1;Xhalf_r0_2;Xhalf_r0_3]);
STDXhalf_r0 = nanstd([Xhalf_r0_1;Xhalf_r0_2;Xhalf_r0_3]);

hold on
errorbar(Time_r0_Min,MeanXhalf_r0,STDXhalf_r0)
plot(Time_r0,Xhalf_r0)
plot([Time_r0(nc13_r0) Time_r0(nc13_r0)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r0(nc14_r0) Time_r0(nc14_r0)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary position along Time - r0')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('Individual','Pooled')

standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Plot the boundary width (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r0_1),length(Xhalf_r0_2),length(Xhalf_r0_3)]);
Width_r0_1 = Width_r0_1(1,1:MinLength);
Width_r0_2 = Width_r0_2(1,1:MinLength)
Width_r0_3 = Width_r0_3(1,1:MinLength)
Time_r0_Min = Time_r0(1:MinLength);

MeanWidth_r0 = nanmean([Width_r0_1;Width_r0_2;Width_r0_3]);
STDWidth_r0 = nanstd([Width_r0_1;Width_r0_2;Width_r0_3]);

hold on
errorbar(Time_r0_Min,MeanWidth_r0,STDWidth_r0)
plot(Time_r0,Width_r0)
plot([Time_r0(nc13_r0) Time_r0(nc13_r0)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r0(nc14_r0) Time_r0(nc14_r0)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary width along Time - r0')
xlabel('Time')
ylabel('Boundary width(AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% r1
%% Boundary features from individual embryos - MeanVectorAP (r1)
r1Data(1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1-1_FromNC12.mat');
r1Data(2) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1-2_FromNC12.mat');
r1Data(3) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1-3_FromNC12.mat');

[Xhalf_r1_1,Width_r1_1] = CalculateBoundary_SpotFluo(r1Data(1).MeanVectorAP,...
                                                   r1Data(1).ElapsedTime,...
                                                    r1Data(1).nc12,r1Data(1).nc13,r1Data(1).nc14);
                                                
[Xhalf_r1_2,Width_r1_2] = CalculateBoundary_SpotFluo(r1Data(2).MeanVectorAP,...
                                                   r1Data(2).ElapsedTime,...
                                                    r1Data(2).nc12,r1Data(2).nc13,r1Data(2).nc14); 
                                                        
[Xhalf_r1_3,Width_r1_3] = CalculateBoundary_SpotFluo(r1Data(3).MeanVectorAP,...
                                                   r1Data(3).ElapsedTime,...
                                                    r1Data(3).nc12,r1Data(3).nc13,r1Data(3).nc14);
                                                
%% Plot the boundary position (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r1_1),length(Xhalf_r1_2),length(Xhalf_r1_3)]);
Xhalf_r1_1 = Xhalf_r1_1(1,1:MinLength);
Xhalf_r1_2 = Xhalf_r1_2(1,1:MinLength)
Xhalf_r1_3 = Xhalf_r1_3(1,1:MinLength)
Time_r1_Min = Time_r1(1:MinLength);

MeanXhalf_r1 = nanmean([Xhalf_r1_1;Xhalf_r1_2;Xhalf_r1_3]);
STDXhalf_r1 = nanstd([Xhalf_r1_1;Xhalf_r1_2;Xhalf_r1_3]);

hold on
errorbar(Time_r1_Min,MeanXhalf_r1,STDXhalf_r1)
plot(Time_r1,Xhalf_r1)
plot([Time_r1(nc13_r1) Time_r1(nc13_r1)],[0 0.5],'k')%,'LineWidth',5)
plot([Time_r1(nc14_r1) Time_r1(nc14_r1)],[0 0.5],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary position along Time - r1')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('Individual','Pooled')

standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Plot the boundary width (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r1_1),length(Xhalf_r1_2),length(Xhalf_r1_3)]);
Width_r1_1 = Width_r1_1(1,1:MinLength);
Width_r1_2 = Width_r1_2(1,1:MinLength)
Width_r1_3 = Width_r1_3(1,1:MinLength)
Time_r1_Min = Time_r1(1:MinLength);

MeanWidth_r1 = nanmean([Width_r1_1;Width_r1_2;Width_r1_3]);
STDWidth_r1 = nanstd([Width_r1_1;Width_r1_2;Width_r1_3]);

hold on
errorbar(Time_r1_Min,MeanWidth_r1,STDWidth_r1)
plot(Time_r1,Width_r1)
plot([Time_r1(nc13_r1) Time_r1(nc13_r1)],[0 0.5],'k')%,'LineWidth',5)
plot([Time_r1(nc14_r1) Time_r1(nc14_r1)],[0 0.5],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary width along Time - r1')
xlabel('Time')
ylabel('Boundary width(AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Individual embryos vs averaged embryos : r2 datasets
% % Calculate the AccumulatedmRNA, and also synchronize as the beginning of
% % nc12
% AverageDatasets_FineTuning('r2-1','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
% AverageDatasets_FineTuning('r2-2','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
% AverageDatasets_FineTuning('r2-3','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
%% Boundary features from individual embryos - MeanVectorAP (r2)
r2Data(1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2-1_FromNC12.mat');
r2Data(2) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2-2_FromNC12.mat');
r2Data(3) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2-3_FromNC12.mat');

[Xhalf_r2_1,Width_r2_1] = CalculateBoundary_SpotFluo(r2Data(1).MeanVectorAP,...
                                                   r2Data(1).ElapsedTime,...
                                                    r2Data(1).nc12,r2Data(1).nc13,r2Data(1).nc14);
                                                
[Xhalf_r2_2,Width_r2_2] = CalculateBoundary_SpotFluo(r2Data(2).MeanVectorAP,...
                                                   r2Data(2).ElapsedTime,...
                                                    r2Data(2).nc12,r2Data(2).nc13,r2Data(2).nc14); 
                                                        
[Xhalf_r2_3,Width_r2_3] = CalculateBoundary_SpotFluo(r2Data(3).MeanVectorAP,...
                                                   r2Data(3).ElapsedTime,...
                                                    r2Data(3).nc12,r2Data(3).nc13,r2Data(3).nc14);
                                                
%% Plot the boundary position (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r2_1),length(Xhalf_r2_2),length(Xhalf_r2_3)]);
Xhalf_r2_1 = Xhalf_r2_1(1,1:MinLength);
Xhalf_r2_2 = Xhalf_r2_2(1,1:MinLength)
Xhalf_r2_3 = Xhalf_r2_3(1,1:MinLength)
Time_r2_Min = Time_r2(1:MinLength);

MeanXhalf_r2 = nanmean([Xhalf_r2_1;Xhalf_r2_2;Xhalf_r2_3]);
STDXhalf_r2 = nanstd([Xhalf_r2_1;Xhalf_r2_2;Xhalf_r2_3]);

hold on
errorbar(Time_r2_Min,MeanXhalf_r2,STDXhalf_r2)
plot(Time_r2,Xhalf_r2)
plot([Time_r2(nc13_r2) Time_r2(nc13_r2)],[0 0.5],'k')%,'LineWidth',5)
plot([Time_r2(nc14_r2) Time_r2(nc14_r2)],[0 0.5],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary position along Time - r2')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('Individual','Pooled')

standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Plot the boundary width (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r2_1),length(Xhalf_r2_2),length(Xhalf_r2_3)]);
Width_r2_1 = Width_r2_1(1,1:MinLength);
Width_r2_2 = Width_r2_2(1,1:MinLength)
Width_r2_3 = Width_r2_3(1,1:MinLength)
Time_r2_Min = Time_r2(1:MinLength);

MeanWidth_r2 = nanmean([Width_r2_1;Width_r2_2;Width_r2_3]);
STDWidth_r2 = nanstd([Width_r2_1;Width_r2_2;Width_r2_3]);

hold on
errorbar(Time_r2_Min,MeanWidth_r2,STDWidth_r2)
plot(Time_r2,Width_r2)
plot([Time_r2(nc13_r2) Time_r2(nc13_r2)],[0 0.5],'k')%,'LineWidth',5)
plot([Time_r2(nc14_r2) Time_r2(nc14_r2)],[0 0.5],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary width along Time - r2')
xlabel('Time')
ylabel('Boundary width(AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')

%% Individual embryos vs averaged embryos : r3 datasets
% Calculate the AccumulatedmRNA, and also synchronize as the beginning of
% nc12
AverageDatasets_FineTuning('r3-1','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets_FineTuning('r3-2','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets_FineTuning('r3-3','NC',12,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');

%% Boundary features from individual embryos - MeanVectorAP (r3)
r3Data(1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3-1_FromNC12.mat');
r3Data(2) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3-2_FromNC12.mat');
r3Data(3) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3-3_FromNC12.mat');

[Xhalf_r3_1,Width_r3_1] = CalculateBoundary_SpotFluo(r3Data(1).MeanVectorAP,...
                                                   r3Data(1).ElapsedTime,...
                                                    r3Data(1).nc12,r3Data(1).nc13,r3Data(1).nc14);
                                                
[Xhalf_r3_2,Width_r3_2] = CalculateBoundary_SpotFluo(r3Data(2).MeanVectorAP,...
                                                   r3Data(2).ElapsedTime,...
                                                    r3Data(2).nc12,r3Data(2).nc13,r3Data(2).nc14); 
                                                        
[Xhalf_r3_3,Width_r3_3] = CalculateBoundary_SpotFluo(r3Data(3).MeanVectorAP,...
                                                   r3Data(3).ElapsedTime,...
                                                    r3Data(3).nc12,r3Data(3).nc13,r3Data(3).nc14);
                                                
%% Plot the boundary position (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r3_1),length(Xhalf_r3_2),length(Xhalf_r3_3)]);
Xhalf_r3_1 = Xhalf_r3_1(1,1:MinLength);
Xhalf_r3_2 = Xhalf_r3_2(1,1:MinLength)
Xhalf_r3_3 = Xhalf_r3_3(1,1:MinLength)
Time_r3_Min = Time_r3(1:MinLength);

MeanXhalf_r3 = nanmean([Xhalf_r3_1;Xhalf_r3_2;Xhalf_r3_3]);
STDXhalf_r3 = nanstd([Xhalf_r3_1;Xhalf_r3_2;Xhalf_r3_3]);

hold on
errorbar(Time_r3_Min,MeanXhalf_r3,STDXhalf_r3)
plot(Time_r3,Xhalf_r3)
plot([Time_r3(nc13_r3) Time_r3(nc13_r3)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r3(nc14_r3) Time_r3(nc14_r3)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary position along Time - r3')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('Individual','Pooled')

standardizeFigure_YJK(gca,legend,'black','black','red','cyan')
%% Plot the boundary width (average over individual embryo, vs from averaged traces)
MinLength = min([length(Xhalf_r3_1),length(Xhalf_r3_2),length(Xhalf_r3_3)]);
Width_r3_1 = Width_r3_1(1,1:MinLength);
Width_r3_2 = Width_r3_2(1,1:MinLength)
Width_r3_3 = Width_r3_3(1,1:MinLength)
Time_r3_Min = Time_r3(1:MinLength);

MeanWidth_r3 = nanmean([Width_r3_1;Width_r3_2;Width_r3_3]);
STDWidth_r3 = nanstd([Width_r3_1;Width_r3_2;Width_r3_3]);

hold on
errorbar(Time_r3_Min,MeanWidth_r3,STDWidth_r3)
plot(Time_r3,Width_r3)
plot([Time_r3(nc13_r3) Time_r3(nc13_r3)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r3(nc14_r3) Time_r3(nc14_r3)],[0 0.6],'k')%,'LineWidth',5)
ylim([0 0.6])
title('Boundary width along Time - r3')
xlabel('Time')
ylabel('Boundary width(AP)')
legend('Individual','Pooled')
standardizeFigure_YJK(gca,legend,'black','black','red','cyan')

%% Plot all boundary positions (r0,1,2,3)
hold on
plot([Time_r0(nc13_r0) Time_r0(nc13_r0)],[0 0.6],'k')%,'LineWidth',5)
plot([Time_r0(nc14_r0) Time_r0(nc14_r0)],[0 0.6],'k')%,'LineWidth',5)

Nstep = 3;
errorbar(Time_r0_Min(1:Nstep:end),MeanXhalf_r0(1:Nstep:end),STDXhalf_r0(1:Nstep:end),'capsize',0)%,'b','LineWidth',2)
errorbar(Time_r1_Min(1:Nstep:end),MeanXhalf_r1(1:Nstep:end),STDXhalf_r1(1:Nstep:end),'capsize',0)%,'r','LineWidth',2)
errorbar(Time_r2_Min(1:Nstep:end),MeanXhalf_r2(1:Nstep:end),STDXhalf_r2(1:Nstep:end),'capsize',0)%,'y','LineWidth',2)
errorbar(Time_r3_Min(1:Nstep:end),MeanXhalf_r3(1:Nstep:end),STDXhalf_r3(1:Nstep:end),'capsize',0)%,'g','LineWidth',2)

%plot(Time_r0,Xhalf_r0)%,'b','LineWidth',3)
%plot(Time_r1,Xhalf_r1)%,'r','LineWidth',3)
%plot(Time_r2,Xhalf_r2)%,'y','LineWidth',3)
%plot(Time_r3,Xhalf_r3)%,'g','LineWidth',3)


ylim([0 0.6])
title('Boundary position along Time - Mean Spot Fluorescence')
xlabel('Time')
ylabel('Boundary Position (AP)')
legend('nc13','nc14','r0','r1','r2','r3',...
        'r0-pooled','r1-pooled','r2-pooled','r3-pooled')
standardizeFigure_YJK(gca,legend,'red','yellow','cyan','magenta','black','black')%,'red','yellow','cyan','magenta','black','black')
%% Plot all boundary widths (r0,1,2,3)
hold on
plot([Time_r0(nc13) Time_r0(nc13)],[0 0.5],'k')
plot([Time_r0(nc14) Time_r0(nc14)],[0 0.5],'k')

errorbar(Time_r0_Min,MeanWidth_r0,STDWidth_r0)
errorbar(Time_r1_Min,MeanWidth_r1,STDWidth_r1)
errorbar(Time_r2_Min,MeanWidth_r2,STDWidth_r2)
errorbar(Time_r3_Min,MeanWidth_r3,STDWidth_r3)

plot(Time_r0,Width_r0)
plot(Time_r1,Width_r1)
plot(Time_r2,Width_r2)
plot(Time_r3,Width_r3)


ylim([0 0.5])
title('Boundary width along Time')
xlabel('Time')
ylabel('Boundary width (AP)')
legend('nc13','nc14','r0','r1','r2','r3',...
        'r0-pooled','r1-pooled','r2-pooled','r3-pooled')
standardizeFigure_YJK(gca,legend,'red','yellow','cyan','magenta','red','yellow','cyan','magenta','black','black')


%% Save values
% 1) Mean Xhalf(mRNA)
MeanXhalf_r0_Spots = MeanXhalf_r0;
MeanXhalf_r1_Spots = MeanXhalf_r1;
MeanXhalf_r2_Spots = MeanXhalf_r2;
MeanXhalf_r3_Spots = MeanXhalf_r3;

% 2) STD Xhalf (mRNA)
STDXhalf_r0_Spots = STDXhalf_r0;
STDXhalf_r1_Spots = STDXhalf_r1;
STDXhalf_r2_Spots = STDXhalf_r2;
STDXhalf_r3_Spots = STDXhalf_r3;
end
 %% Old script for making a movie of AP profile (Bcd, Runt, and hb P@-N
% % First, make a time series of nc 14 with 10sec interval
% TimeSeries = 0:0.5:90; % min
% 
% Folder = 'E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\hbP2-rN-AccumulatedmRNA-withBcd'
% for i=1:length(TimeSeries)
%     time = TimeSeries(i);
%     clf
%     hold on
%     % Bicoid conc. plot over AP
%     % Find the time point close to the TimeSeries
%     Bcdscale = 8;
%     [Dummy,BcdIndex]=min((BcdTime-time).^2);
%     
%     errorbar(0:0.025:1,Bcdscale*BcdData.MeanVectorAP(BcdIndex,:),...
%         Bcdscale*BcdData.SDVectorAP(BcdIndex,:)./sqrt(BcdData.NParticlesAP(BcdIndex,:)))
%     
%     % Runt conc. plot over AP
%     % Find the time point close to the TimeSeries
% %      [Dummy,RuntIndex]=min((RuntTime-time).^2);
% %      RuntScale = 1;
% % %     1) Nuclear Fluo
% %     errorbar(0:0.01:1,RuntScale*(RuntData.MeanVectorAP(RuntIndex,:)-1700),...
% %         RuntScale*RuntData.SDVectorAP(RuntIndex,:)./sqrt(RuntData.NParticlesAP(RuntIndex,:)))
% %     2) Nuc Fluo - CytoFluo / KmCh
% %       errorbar(0:0.01:1,TF_Nuclei(RuntIndex,:),TFError(RuntIndex,:))
% %     3) Just subtracting the nuclear fluo of No NB data.
% %      errorbar(0:0.01:1,TF_mCh(RuntIndex,:),TF_mCh_Error(RuntIndex,:))
% 
%     % r0 (MS2)
%     % Find the time point close to the TimeSeries
%     [Dummy,r0Index]=min((Time_r0-time).^2);
%     Scale_mRNA = 1/20;
%     % 1) Instantaneous RNAP loading rate (Spot fluo)
% %     errorbar(0:0.025:1,r0.MeanVectorAP(r0Index,:),...
% %         r0.SDVectorAP(r0Index,:)./sqrt(r0.NParticlesAP(r0Index,:)))
%     % 2) Accumulated mRNA (accumulated spot fluo)
%      errorbar(0:0.025:1,Scale_mRNA*(AccumulatedmRNA_r0(r0Index,:)),...
%                         Scale_mRNA*SDAccumulatedmRNA_r0(r0Index,:))
%      
%      % r1 (MS2)
%         [Dummy,r1Index]=min((Time_r1-time).^2);
%     Scale_mRNA = 1/20;
%     % 1) Instantaneous RNAP loading rate (Spot fluo)
% %     errorbar(0:0.025:1,r1.MeanVectorAP(r1Index,:),...
% %         r1.SDVectorAP(r1Index,:)./sqrt(r1.NParticlesAP(r1Index,:)))
%     % 2) Accumulated mRNA (accumulated spot fluo)
%      errorbar(0:0.025:1,Scale_mRNA*AccumulatedmRNA_r1(r1Index,:),...
%                         Scale_mRNA*SDAccumulatedmRNA_r1(r1Index,:))
%      
%     % r2 (MS2)
%         [Dummy,r2Index]=min((Time_r2-time).^2);
%     Scale_mRNA = 1/20;
%     % 1) Instantaneous RNAP loading rate (Spot fluo)
% %     errorbar(0:0.025:1,r2.MeanVectorAP(r2Index,:),...
% %         r2.SDVectorAP(r2Index,:)./sqrt(r2.NParticlesAP(r2Index,:)))
%     % 2) Accumulated mRNA (accumulated spot fluo)
%      errorbar(0:0.025:1,Scale_mRNA*AccumulatedmRNA_r2(r2Index,:),...
%                         Scale_mRNA*SDAccumulatedmRNA_r2(r2Index,:))
% 
%     % r3 (MS2)
%     % Find the time point close to the TimeSeries
%     [Dummy,r3Index]=min((Time_r3-time).^2);
%     % 1) Spot Fluorescence (instantaneous RNAP loading rate)
% %     errorbar(0:0.025:1,r3Data.MeanVectorAP(r3Index,:),...
% %         r3Data.SDVectorAP(r3Index,:)./sqrt(r3Data.NParticlesAP(r3Index,:)))
%     % 2) Accumulated mRNA (accumulated spot fluo)
%      errorbar(0:0.025:1,Scale_mRNA*AccumulatedmRNA_r3(r3Index,:),...
%                         Scale_mRNA*SDAccumulatedmRNA_r3(r3Index,:))
%     %xlim([0.15 0.7])
%     ylim([0 1500])
%     title(['Bcd,and Accumulated mRNA @ ',num2str(TimeSeries(i)),'min'])
%     xlabel('AP')
%     ylabel('Fluorescence (AU)')
%     legend('Bicoid','r0','r1','r2','r3')
%     standardizeFigure(gca,legend,'g','r','yellow','cyan')
%     %pause
%     fig = gcf;
%     fig.InvertHardcopy = 'off';
%     saveas(gcf, [Folder, '\', 'Time',num2str(i),'min', '.tif'])%[nc14Folder '\' num2str(nucleus) '.png'])
%     
% end




