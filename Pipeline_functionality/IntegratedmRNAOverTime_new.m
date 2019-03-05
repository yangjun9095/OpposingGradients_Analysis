function IntegratedmRNAOverTime
% Analysis on Integrated mRNA over Time 
% Here, the slope is Txn rate
% Also, note that I used MeanVectorAP for accumulating, so I did not take
% into account the FractionON

%% Load the datasets
r0Data = load('r0_FromNC12.mat');
r1Data = load('r1_FromNC12.mat');
r2Data = load('r2_FromNC12.mat');
r3Data = load('r3_FromNC12.mat');
%% plotting
AP = 19;
hold on
errorbar(r0Data.ElapsedTime,r0Data.AccumulatedmRNA(:,AP),r0Data.AccumulatedmRNA_SD(:,AP))
errorbar(r1Data.ElapsedTime,r1Data.AccumulatedmRNA(:,AP),r1Data.AccumulatedmRNA_SD(:,AP))
errorbar(r2Data.ElapsedTime,r2Data.AccumulatedmRNA(:,AP),r2Data.AccumulatedmRNA_SD(:,AP))
errorbar(r3Data.ElapsedTime,r3Data.AccumulatedmRNA(:,AP),r3Data.AccumulatedmRNA_SD(:,AP))

%xlim([0 100])
title(['Integrated mRNA over Time @ AP =',num2str((AP-1)*2.5),' %'])
xlabel('Time (min)')
ylabel('Integrated mRNA (AU)')
legend('r0','r1','r2','r3')
standardizeFigure_YJK(gca,legend,'red','cyan','yellow','lightblue')
% for i=10:21
%     plot(ElapsedTime,AccumulatedmRNA(:,i))
%     pause
% end







end