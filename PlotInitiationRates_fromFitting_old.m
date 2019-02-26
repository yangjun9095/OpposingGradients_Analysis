function PlotInitiationRates_fromFitting
% YJK : This is my old script used for plotting MeanFits.mat for nc13, or
% nc14...
%% MeanFit
%% Load the data set
clear all

mRNAr0=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r0/MeanFits.mat')

mRNAr1=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r1/MeanFits.mat')

mRNAr2=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r2/MeanFits.mat');

mRNAr3=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r3/MeanFits.mat');

%% FitResults
r0=mRNAr0.FitResults;
r1=mRNAr1.FitResults;
r2=mRNAr2.FitResults;
r3=mRNAr3.FitResults;
%% Plotting Mean Rate of transcription at nc 13 along the AP aixs

for i=1:length(r0)
    if r0(i,2).Approved==1
        r00(i)=r0(i,2).RateFit;
    elseif i<22
        r00(i)=nan;
    else
        r00(i)=0;
    end
end

for i=1:length(r1)
    if r1(i,2).Approved==1
        r11(i)=r1(i,2).RateFit;
    elseif i<10
        r11(i)=nan;
    else
        r11(i)=0;
    end
end

for i=1:length(r2)
    if r2(i,2).Approved==1
        r22(i)=r2(i,2).RateFit;
    elseif i<10
        r22(i)=nan;
    else
        r22(i)=0;
    end
end

for i=1:length(r3)
    if r3(i,2).Approved==1
        r33(i)=r3(i,2).RateFit;
    elseif i<10
        r33(i)=nan;
    else
        r33(i)=0;
    end
end

%% Averaging over 3 AP bins
for j=2:length(r0)-1
    r000(j)=nanmean(r00(j-1:j+1));
end

for j=2:length(r1)-1
    r111(j)=nanmean(r11(j-1:j+1));
end

for j=2:length(r2)-1
    r222(j)=nanmean(r22(j-1:j+1));
end

for j=2:length(r3)-1
    r333(j)=nanmean(r33(j-1:j+1));
end

r000(41)=0;
r111(41)=0;
r222(41)=0;
r333(41)=0;
%% 
hold on
plot(0:0.025:1,r000,'k')
plot(0:0.025:1,r111,'b')
plot(0:0.025:1,r222,'g')
plot(0:0.025:1,r333,'r')

xlim([0 1])
ylim([0 250])
title('Rate of Transcription for different number of repressor binding sites')
xlabel('AP')
ylabel('Rate of Transcription (AU)')
legend('0','1','2','3')

set(gca,'Fontsize',30)
hold off
end