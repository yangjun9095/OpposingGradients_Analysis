function [Time,Fold_change] = CalculateFoldChange_Rate(data1,data2)
% This is a function to calculate the Fold-change of Mean spot fluorescence
% (MeanVectorAP) over time.
% INPUT : MeanVectorAP1 and MeanVectorAP2 also with ElapsedTime for
% corresponding ones.
% For this, I will use the Averaged datasets (by AverageDatasets), which
% has synchronized the datasets as the beginning of nc12 (or nc13)

%% Load the Datasets
data1 = 'r0'
data2 = 'r2'
Dataset1 = load([data1,'_FromNC12.mat']);
Dataset2 = load([data2,'_FromNC12.mat']);

%% Define the useful fields
Time1 = Dataset1.ElapsedTime;
Time2 = Dataset2.ElapsedTime;

MeanVectorAP1 = Dataset1.MeanVectorAP;
MeanVectorAP2 = Dataset2.MeanVectorAP;

nc13 = Dataset1.nc13;
nc14 = Dataset1.nc14;

% Find the shorter time point
tLength = min([length(Time1),length(Time2)]);
ShortIndex = find([length(Time1),length(Time2)]==tLength);

% Truncate the longer series to have the same t length.
if ShortIndex == 1
    MeanVectorAP2 = MeanVectorAP2(1:tLength,:);
    Time2 = Time2(1:tLength);
elseif ShortIndex ==2
    MeanVectorAP1 = MeanVectorAP1(1:tLength,:);
    Time1 = Time1(1:tLength);
end

% Calculate the Fold-change
Fold_change = MeanVectorAP2./MeanVectorAP1;
%Fold_change_13 = MeanVectorAP2_13./MeanVectorAP1_13;
%% Plot to check
APbin = 10;

hold on
% plot(Time1,MeanVector1(:,APbin))
% plot(Time2,MeanVector2(:,APbin))
plot(1:27,Fold_change_13(:,APbin))
title(['Fold-change over Time (MeanVectorAP)-',data2,' / ',data1,' @ AP bin = ',num2str((APbin-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Fold-change')

%% Save the fields
Time = Time1;
Fold_change = Fold_change;
MeanVectorAP1 = MeanVectorAP1;
MeanVectorAP2 = MeanVectorAP2;
%save('Time','Fold_change','MeanVectorAP1','MeanVectorAP2',)

%% 
%% Plot to check
% APbin = 11;
% hold on
% plot(1:length(Fold_change_r0r1),Fold_change_r0r1(:,APbin))
% plot(1:length(Fold_change_r0r2),Fold_change_r0r2(:,APbin))
% plot(1:length(Fold_change_r0r3),Fold_change_r0r3(:,APbin))
end