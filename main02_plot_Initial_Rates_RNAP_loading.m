function main02_plot_Initial_Rates_RNAP_loading

%Edited from the script written by Paul, 05/30/2018.

%This code is using the "mean fits" files generated by my version of FitMeanAPSymetric 
%to plot the initial rates along APaxis. 
%The T_on (time of activation of the transcription,extrapolated from the initial rate)
%can also be plotted.

%% Load the datasets (using LoadMS2Sets.m)
Data_r0 = LoadMS2Sets('r0')
Data_r1 = LoadMS2Sets('r1')
Data_r2 = LoadMS2Sets('r2')
Data_r3 = LoadMS2Sets('r3')
Data_r3_prime = LoadMS2Sets('r3prime')

%% Extract useful fields
function [fittedRate,fittedRateSD,fittedTon] = Extract_Fields_MeanFits(Data)
    % Data is the compiled datasets by LoadMS2Sets.m
    
    % From this, let's construct 3D matrices for fitted Rate and fitted T on,
    % which have dimensions like (AP,NC,index of embryo)

    % Initialize the matrices, fill with Nans.
    fittedRate = nan(41,3,length(Data));
    fittedRateSD = nan(41,3,length(Data));
    fittedTon = nan(41,3,length(Data));

    for i=1:length(Data)
        MeanFits = Data(i).MeanFits;
        for j=1:length(MeanFits)
            for NC = 1:3 % NC12 to NC14 by default
                if ~isempty(MeanFits(j,NC).RateFit) & (MeanFits(j,NC).Approved==1)
                    fittedRate(j,NC,i) = MeanFits(j,NC).RateFit;
                    fittedRateSD(j,NC,i) = MeanFits(j,NC).SDRateFit;
                    fittedTon(j,NC,i) = MeanFits(j,NC).TimeStart;
                end
            end
        end
    end
end

%% Extract the fitted values from all of my datasets
[fittedRate_r0,fittedRateSD_r0,fittedTon_r0] = Extract_Fields_MeanFits(Data_r0);
[fittedRate_r1,fittedRateSD_r1,fittedTon_r1] = Extract_Fields_MeanFits(Data_r1);
[fittedRate_r2,fittedRateSD_r2,fittedTon_r2] = Extract_Fields_MeanFits(Data_r2);
[fittedRate_r3,fittedRateSD_r3,fittedTon_r3] = Extract_Fields_MeanFits(Data_r3);
[fittedRate_r3_prime,fittedRateSD_r3_prime,fittedTon_r3_prime] = Extract_Fields_MeanFits(Data_r3_prime);

%% plot individual fitted curves to check
NC = 2; %NC13
hold on
for i=1:length(Data_r0)
    errorbar(0:0.025:1,fittedRate_r0(:,NC,i),fittedRateSD_r0(:,NC,i))
end

for i=1:length(Data_r1)
    errorbar(0:0.025:1,fittedRate_r1(:,NC,i),fittedRateSD_r1(:,NC,i))
end

hold on
for i=1:length(Data_r2)
    errorbar(0:0.025:1,fittedRate_r2(:,NC,i),fittedRateSD_r2(:,NC,i))
    pause
end

hold on
for i=1:length(Data_r3)
    errorbar(0:0.025:1,fittedRate_r3(:,NC,i),fittedRateSD_r3(:,NC,i))
    pause
end

hold on
for i=1:length(Data_r3_prime)
    errorbar(0:0.025:1,fittedRate_r3_prime(:,NC,i),fittedRateSD_r3_prime(:,NC,i))
    pause
end

%% Calculate the average using nanmean, nanstd
average_fittedRate_r0 = nanmean(fittedRate_r0,3);
SEM_fittedRate_r0 = nanstd(fittedRate_r0,0,3)/sqrt(length(Data_r0));

average_fittedRate_r1 = nanmean(fittedRate_r1,3);
SEM_fittedRate_r1 = nanstd(fittedRate_r1,0,3)/sqrt(length(Data_r1));

average_fittedRate_r2 = nanmean(fittedRate_r2,3);
SEM_fittedRate_r2 = nanstd(fittedRate_r2,0,3)/sqrt(length(Data_r2));

average_fittedRate_r3 = nanmean(fittedRate_r3,3);
SEM_fittedRate_r3 = nanstd(fittedRate_r3,0,3)/sqrt(length(Data_r3));

average_fittedRate_r3_prime = nanmean(fittedRate_r3_prime,3);
SEM_fittedRate_r3_prime = nanstd(fittedRate_r3_prime,0,3)/sqrt(length(Data_r3_prime));

%% Plot the averaged fittedRate (initial rate of RNAP loading), and SEM

% NC12
InitialRate_NC12_figure = figure;
nc = 1; % NC12
hold on
errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc))
errorbar(0:0.025:1,average_fittedRate_r1(:,nc),SEM_fittedRate_r1(:,nc))
errorbar(0:0.025:1,average_fittedRate_r2(:,nc),SEM_fittedRate_r2(:,nc))
errorbar(0:0.025:1,average_fittedRate_r3(:,nc),SEM_fittedRate_r3(:,nc))
errorbar(0:0.025:1,average_fittedRate_r3_prime(:,nc),SEM_fittedRate_r3_prime(:,nc))

legend('r0','r1','r2','r3','r3-mutated')
xlabel('AP Position')
ylabel('Initial rate (AU/min)')
title('Initial rate of RNAP loading along AP axis, at NC 14')
standardizeFigure(gca,legend,[])

FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\InitialLoadingRates_r0123_fitted\';
saveas(InitialRate_NC12_figure,[FigPath 'InitialRate' , '_NC12' , '.pdf']); 

% NC13
InitialRate_NC13_figure = figure;
nc = 2; % NC13
hold on
errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc))
errorbar(0:0.025:1,average_fittedRate_r1(:,nc),SEM_fittedRate_r1(:,nc))
errorbar(0:0.025:1,average_fittedRate_r2(:,nc),SEM_fittedRate_r2(:,nc))
errorbar(0:0.025:1,average_fittedRate_r3(:,nc),SEM_fittedRate_r3(:,nc))
errorbar(0:0.025:1,average_fittedRate_r3_prime(:,nc),SEM_fittedRate_r3_prime(:,nc))

legend('r0','r1','r2','r3','r3-mutated')
xlabel('AP Position')
ylabel('Initial rate (AU/min)')
title('Initial rate of RNAP loading along AP axis, at NC 13')
standardizeFigure(gca,legend,[])

FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\InitialLoadingRates_r0123_fitted\';
saveas(InitialRate_NC13_figure,[FigPath 'InitialRate' , '_NC13' , '.pdf']); 

% NC14
InitialRate_NC14_figure = figure;
nc = 2; % NC13
hold on
errorbar(0:0.025:1,average_fittedRate_r0(:,nc),SEM_fittedRate_r0(:,nc))
errorbar(0:0.025:1,average_fittedRate_r1(:,nc),SEM_fittedRate_r1(:,nc))
errorbar(0:0.025:1,average_fittedRate_r2(:,nc),SEM_fittedRate_r2(:,nc))
errorbar(0:0.025:1,average_fittedRate_r3(:,nc),SEM_fittedRate_r3(:,nc))
errorbar(0:0.025:1,average_fittedRate_r3_prime(:,nc),SEM_fittedRate_r3_prime(:,nc))

legend('r0','r1','r2','r3','r3-mutated')
xlabel('AP Position')
ylabel('Initial rate (AU/min)')
title('Initial rate of RNAP loading along AP axis, at NC 14')
standardizeFigure(gca,legend,[])

FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\InitialLoadingRates_r0123_fitted\';
saveas(InitialRate_NC14_figure,[FigPath 'InitialRate' , '_NC14' , '.pdf']); 


%% Save the fitted initial rate and SEM 
save('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedInitialRate_trapezoidalfit.mat',...
        'average_fittedRate_r0','SEM_fittedRate_r0',...
        'average_fittedRate_r1','SEM_fittedRate_r1',...
        'average_fittedRate_r2','SEM_fittedRate_r2',...
        'average_fittedRate_r0','SEM_fittedRate_r0')







%%%%%%%%%%%%%%%%%%%%%Old script by Paul%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 %% Load the datasets

%Data for r0 on Data(1,:). Change the N value, number of datasets
Data(1,1)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-03-23-hbP2-r0-MS2V5-lacZ\MeanFits.mat');
Data(1,2)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-04-23-hbP2-r0-MS2V5-lacZ\MeanFits.mat');
Data(1,3)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-04-24-hbP2-r0-MS2V5-lacZ\MeanFits.mat');
Data(1,4)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-04-25-hbP2-r0-MS2V5-lacZ\MeanFits.mat');
N(1)=4;

%Data for r1
Data(2,1)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-04-30-hbP2-r1-MS2V5-lacZ\MeanFits.mat');
Data(2,2)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-05-01-hbP2-r1-MS2V5-lacZ\MeanFits.mat');
Data(2,3)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-05-04-hbP2-r1-MS2V5-lacZ\MeanFits.mat');
N(2)=3;

%Data for r2
Data(3,1)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-04-08-hbP2-r2-MS2V5-lacZ\MeanFits.mat');
Data(3,2)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-05-09-hbP2-r2-MS2V5-lacZ\MeanFits.mat');
Data(3,3)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-05-09-hbP2-r2-MS2V5-lacZ-2\MeanFits.mat');
N(3)=3;

%Data for r3
Data(4,1)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-03-23-hbP2-r3-MS2V5-lacZ\MeanFits.mat');
Data(4,2)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-03-25-hbP2-r3-MS2V5-lacZ\MeanFits.mat');
Data(4,3)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-04-13-hbP2-r3-MS2V5-lacZ\MeanFits.mat');
N(4)=3;

%Data for r3_prime
Data(5,1)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-07-18-hbP2-r3prime-MS2V5-lacZ\MeanFits.mat');
Data(5,2)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-07-20-hbP2-r3prime-MS2V5-lacZ\MeanFits.mat');
Data(5,3)=load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-07-20-hbP2-r3prime-MS2V5-lacZ2\MeanFits.mat');
N(5)=3;



%% Fill the variables to plot

% T0(Construct, DataSet, AP, NC), same for the others
T0=nan(5,4,41,3);
Rate=nan(5,4,41,3);
SD_rate=nan(5,4,41,3);


Nb_Approved_Rate=zeros(5,41,3);


for i = 1:5 %Construct r0,1,2,3,3' loop
    for j = 1:N(i) % Different number of dataSets for each construct
        for NC= 12:14
            for AP = 1:41
                if Data(i,j).FitResults(AP,NC-11).Approved == 1
                    
                    T0(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).TimeStart;
                    Rate(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).RateFit; % RateFit1 for manually fitted values
                    SD_rate(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).SDRateFit;
                    
                    Nb_Approved_Rate(i,AP,NC-11)=Nb_Approved_Rate(i,AP,NC-11)+1;
                end
            end
        end
    end
end

%% Plotting ! - for Rate of Txn initiation

%NC 12


X=0:0.025:1;

NC=1
figure(1)

% for i=1:N(1)
%     Y=squeeze(Rate(1,i,:,NC));
%     Z=squeeze(SD_rate(1,i,:,NC));
%     p(i)=errorbar(X,Y,Z);
%     p(i).Marker='*';
%     p(i).LineStyle='none';
%     p(i).Color='red';%[1,i/(1.5*N0),0]
%     %p(i).legend ='r0';
%     hold on
% end
% for i=1:N(2)
%     Y=squeeze(Rate(2,i,:,NC));
%     Z=squeeze(SD_rate(2,i,:,NC));
%     q(i)=errorbar(X,Y,Z);
%     q(i).Marker='+';
%     q(i).LineStyle='none';
%     q(i).Color='magenta';%[1,i/(1.5*N0),0]
%     %q(i).legend ='r0';
%     hold on
% end
% for i=1:N(3)
%     Y=squeeze(Rate(3,i,:,NC));
%     Z=squeeze(SD_rate(3,i,:,NC));
%     r(i)=errorbar(X,Y,Z);
%     r(i).Marker='.';
%     r(i).LineStyle='none';
%     r(i).Color='blue';%[1,i/(1.5*N0),0]
%     %r(i).legend ='r0';
%     %hold on
% end
% for i=1:N(4)
%     Y=squeeze(Rate(4,i,:,NC));
%     Z=squeeze(SD_rate(4,i,:,NC));
%     s(i)=errorbar(X,Y,Z);
%     s(i).Marker='s';
%     s(i).LineStyle='none';
%     s(i).Color='cyan';%[0,i/N3,1]
%     s(i).legend ='r3';
% end

for i=1:5
    Y=(nanmean(squeeze((Rate(i,:,:,NC)))));
    Z=(nanstd(squeeze((Rate(i,:,:,NC)))));
    Mean12(i,1:41)=(Y);
    STD12(i,1:41)=(Z); 
    SEM12(i,1:41) = STD12(i,1:41) / sqrt(N(i));
end

% for i=1:4
%     for j=1:41
%         Mean12(i,j)=Mean(i,j)./Nb_Approved_Rate(i,j,NC);
%     end
% end

%t=plot(X,Mean12(1,:),'-r',X,Mean12(2,:),'-m',X,Mean12(3,:),'-b',X,Mean12(4,:),'-c','LineWidth',2);
InitialRate_NC12_figure = figure;
hold on
errorbar(X,Mean12(1,:),SEM12(1,:))
errorbar(X,Mean12(2,:),SEM12(2,:))
errorbar(X,Mean12(3,:),SEM12(3,:))
errorbar(X,Mean12(4,:),SEM12(4,:))
errorbar(X,Mean12(5,:),SEM12(5,:))
ylim([0 400])
%legend(t,{'r0, mean','r1, mean','r2, mean','r3, mean'})
legend('r0','r1','r2','r3','r3-mutated')
xlabel('AP Position (%)')
ylabel('Initial rate')
title('Initial rate of r0, r1, r2, r3, and r3(mutated) along AP axis, at NC 12')
standardizeFigure(gca,legend,[])

%FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\InitialLoadingRates_r0123_fitted\';
%saveas(InitialRate_NC12_figure,[FigPath 'InitialRate' , '_NC12' , '.pdf']); 
%save('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\AveragedInitialRate_NC12.mat','Mean12','STD12')

%% NC13


NC=2
%figure(2)

% for i=1:N(1)
%     Y=squeeze(Rate(1,i,:,NC));
%     Z=squeeze(SD_rate(1,i,:,NC));
%     p(i)=errorbar(X,Y,Z);
%     p(i).Marker='*';
%     p(i).LineStyle='none';
%     p(i).Color='red';%[1,i/(1.5*N0),0]
%     %p(i).legend ='r0';
%     hold on
% end
% for i=1:N(2)
%     Y=squeeze(Rate(2,i,:,NC));
%     Z=squeeze(SD_rate(2,i,:,NC));
%     q(i)=errorbar(X,Y,Z);
%     q(i).Marker='+';
%     q(i).LineStyle='none';
%     q(i).Color='magenta';%[1,i/(1.5*N0),0]
%     %q(i).legend ='r0';
%     hold on
% end
% for i=1:N(3)
%     Y=squeeze(Rate(3,i,:,NC));
%     Z=squeeze(SD_rate(3,i,:,NC));
%     r(i)=errorbar(X,Y,Z);
%     r(i).Marker='.';
%     r(i).LineStyle='none';
%     r(i).Color='blue';%[1,i/(1.5*N0),0]
%     %r(i).legend ='r0';
%     hold on
% end
% for i=1:N(4)
%     Y=squeeze(Rate(4,i,:,NC));
%     Z=squeeze(SD_rate(4,i,:,NC));
%     s(i)=errorbar(X,Y,Z);
%     s(i).Marker='s';
%     s(i).LineStyle='none';
%     s(i).Color='cyan';%[0,i/N3,1]
%     %s(i).legend ='r3';
% end
% 
% for i=1:4
%     Y=(nansum(squeeze((Rate(i,:,:,NC)))));
%     Mean(i,1:41)=(Y);
% end
% 
% for i=1:4
%     for j=1:41
%         Mean13(i,j)=Mean(i,j)./Nb_Approved_Rate(i,j,NC);
%     end
% end
% 
% 
% %save('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Paul_codes\Fitting initial rate/Mean13.mat','Mean')
% t=plot(X,Mean13(1,:),'-r',X,Mean13(2,:),'-m',X,Mean13(3,:),'-b',X,Mean13(4,:),'-c','LineWidth',2);
% ylim([0 400])
% legend(t,{'r0, mean','r1, mean','r2, mean','r3, mean'})
% xlabel('AP Position (%)')
% ylabel('Initial rate')
% title('Initial rate of r0, r1, r2 and r3 along AP axis, at NC 13')


%standardizeFigure_YJK(gca,legend,'blue')
%saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\Initial_rate/NC13.tif')
for i=1:5
    Y=(nanmean(squeeze((Rate(i,:,:,NC)))));
    Z=(nanstd(squeeze((Rate(i,:,:,NC)))));
    Mean13(i,1:41)=(Y);
    STD13(i,1:41)=(Z); 
    SEM13(i,1:41) = STD13(i,1:41) / sqrt(N(i));
end

% fix the artefact?
Mean13(4,19) = 30;

InitialRate_NC13_figure = figure;
hold on
errorbar(X,Mean13(1,:),SEM13(1,:))
errorbar(X,Mean13(2,:),SEM13(2,:))
errorbar(X,Mean13(3,:),SEM13(3,:))
errorbar(X,Mean13(4,:),SEM13(4,:))
errorbar(X,Mean13(5,:),SEM13(5,:))
ylim([0 400])

legend('r0','r1','r2','r3','r3(mutated)')
xlabel('AP Position (%)')
ylabel('Initial rate')
title('Initial rate of r0, r1, r2,r3 and r3(mutated) along AP axis, at NC 13')
standardizeFigure(gca,legend,[])

FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\InitialLoadingRates_r0123_fitted\';
saveas(InitialRate_NC13_figure,[FigPath 'InitialRate' , '_NC13' , '.pdf']); 
save('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\AveragedInitialRate_NC13.mat','Mean13','STD13')

%% NC14
NC=3
%figure(3)

% for i=1:N(1)
%     Y=squeeze(Rate(1,i,:,NC));
%     Z=squeeze(SD_rate(1,i,:,NC));
%     p(i)=errorbar(X,Y,Z);
%     p(i).Marker='*';
%     p(i).LineStyle='none';
%     p(i).Color='red';%[1,i/(1.5*N0),0]
%     %p(i).legend ='r0';
%     hold on
% end
% for i=1:N(2)
%     Y=squeeze(Rate(2,i,:,NC));
%     Z=squeeze(SD_rate(2,i,:,NC));
%     q(i)=errorbar(X,Y,Z);
%     q(i).Marker='+';
%     q(i).LineStyle='none';
%     q(i).Color='magenta';%[1,i/(1.5*N0),0]
%     %q(i).legend ='r0';
%     hold on
% end
% for i=1:N(3)
%     Y=squeeze(Rate(3,i,:,NC));
%     Z=squeeze(SD_rate(3,i,:,NC));
%     r(i)=errorbar(X,Y,Z);
%     r(i).Marker='.';
%     r(i).LineStyle='none';
%     r(i).Color='blue';%[1,i/(1.5*N0),0]
%     %r(i).legend ='r0';
%     hold on
% end
% 
% for i=1:N(4)
%     Y=squeeze(Rate(4,i,:,NC));
%     Z=squeeze(SD_rate(4,i,:,NC));
%     s(i)=errorbar(X,Y,Z);
%     s(i).Marker='s';
%     s(i).LineStyle='none';
%     s(i).Color='cyan';%[0,i/N3,1]
%     %s(i).legend ='r3';
% end
% 
% for i=1:4
%     Y=(nansum(squeeze((Rate(i,:,:,NC)))));
%     Mean(i,1:41)=(Y);
% end
% 
% for i=1:4
%     for j=1:41
%         Mean14(i,j)=Mean(i,j)./Nb_Approved_Rate(i,j,NC);
%     end
% end
% 
% t=plot(X,Mean14(1,:),'-r',X,Mean14(2,:),'-m',X,Mean14(3,:),'-b',X,Mean14(4,:),'-c','LineWidth',2);
% ylim([0 400])
% legend(t,{'r0, mean','r1, mean','r2, mean','r3, mean'})
% xlabel('AP Position (%)')
% ylabel('Initial rate')
% title('Initial rate of r0, r1, r2 and r3 along AP axis, at NC 14')


%standardizeFigure_YJK(gca,legend,'blue')
%saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\Initial_rate/NC14.tif')


for i=1:5
    Y=(nanmean(squeeze((Rate(i,:,:,NC)))));
    Z=(nanstd(squeeze((Rate(i,:,:,NC)))));
    Mean14(i,1:41)=(Y);
    STD14(i,1:41)=(Z); 
    SEM14(i,1:41) = STD14(i,1:41) / sqrt(N(i));
end

InitialRate_NC14_figure = figure;
hold on
errorbar(X,Mean14(1,:),SEM14(1,:))
errorbar(X,Mean14(2,:),SEM14(2,:))
errorbar(X,Mean14(3,:),SEM14(3,:))
errorbar(X,Mean14(4,:),SEM14(4,:))
errorbar(X,Mean14(5,:),SEM14(5,:))
ylim([0 400])

legend('r0','r1','r2','r3','r3(mutated)')
xlabel('AP Position (%)')
ylabel('Initial rate')
title('Initial rate of r0, r1, r2,r3, and r3(mutated) along AP axis, at NC 14')
standardizeFigure(gca,legend,[])

 FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\InitialLoadingRates_r0123_fitted\';
 saveas(InitialRate_NC14_figure,[FigPath 'InitialRate' , '_NC14' , '.pdf']); 
 save('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\AveragedInitialRate_NC14.mat','Mean14','STD14')
%% save
%% Part2. Fitting the initial loading rates (mean) with theoretical models
% Here, I will use Mean12, Mean13, and Mean14 for fitting with some
% expression for different constructs with different numbers of binding
% sites. Note that the Mean12,13, and 14 means that initial rate of RNAP
% loading fit from individual embryos, then averaged.

% Load the Bcd, Runt data
Bicoid = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed/Bcd-Averaged');
Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-05-24-Runt-JB3-MCP-mCherry-vasa-eGFP1\CompiledNuclei');

% First, we should build a prediction matrix using handful of parameters.
% Use Rate_r0.m for optimizing the parameters, r_max, r_basal, and K_d
% using lsqnonlin.

% Second, we need to decide which nc we're going to model.
% One good starting place would be nc13, although we still don't have a
% good predictive power on the whole timecourse. Or, the beginning of nc14
% might be good as well.

% Now, I will start from 
BcdFluo = Bicoid.MeanVectorAP(Bicoid.nc14+20,:);
RuntFluo = Runt.MeanVectorAP(Runt.nc14+20,:);

% plot to check
hold on
plot(BcdFluo)
plot(RuntFluo)
title('Bcd and Runt profile over AP')
xlabel('AP')
ylabel('Protein concentration (AU)')
legend('Bcd','Runt')

%% limit the AP bins for fitting
APbinstart = 11;
APbinend = 26;

% Initial condition
R_max=500; % maximum
R_bas=100; % basal, probably due to ectopic expression...
Kd = 0.1;
% fitting for r0 first,
p0=[0.1,R_max,R_bas];
%lb=[0,R_max,R_bas];
lb = [0, 0, 0];
ub=[1,R_max,R_bas];
options = optimset('Display','iter');

f= @(p)Rate_r0(p,BcdFluo(APbinstart:APbinend))-Mean14(1,APbinstart:APbinend);

P = lsqnonlin(f, p0, lb, ub, options);
%% Part3. Rates of RNAP loading for trapezium
% added by Yang Joon Kim, 9/1/2018
% The goal is to plot two different RNAP loading rates for trapezium.
% The fitting is done by FitTiltedTrapezoid_Mean.m script

%% Load the datasets
Data(1,1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\r0_InitialFit.mat');
Data(2,1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\r1_InitialFit.mat');
Data(3,1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\r2_InitialFit.mat');
Data(4,1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\r3_InitialFit.mat');

%% Fill the variables to plot

% T0(Construct, DataSet, AP, NC), same for the others
T0=nan(4,1,41,3);
Rate1=nan(4,1,41,3);
Rate2=nan(4,1,41,3);
Rate3=nan(4,1,41,3);
%SD_rate=nan(4,4,41,3);


%Nb_Approved_Rate=zeros(4,41,3);


for i = 1:4 %Construct r0 to r3 loop
    for j = 1:length(Data(i,:)) % Different number of dataSets for each construct
        for NC= 12:14
            for AP = 1:41
                if Data(i).FitResults(AP,NC-11).Approved == 1
                    
                    T0(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).TimeStart;
                    Rate1(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).RateFit1;
                    Rate2(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).RateFit2;
                    Rate3(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).RateFit3;
                    %SD_rate(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).SDRateFit;
                    
                    %Nb_Approved_Rate(i,AP,NC-11)=Nb_Approved_Rate(i,AP,NC-11)+1;
                end
            end
        end
    end
end




%% Plotting for Rate of Txn initiation (Rate1)

%NC 12

X=0:0.025:1;

NC=2
figure(1)
hold on
for i=1:length(Data(1,:))
    Y=squeeze(Rate1(1,i,:,NC));
    %Z=squeeze(SD_rate(1,i,:,NC));
    %p(i)=errorbar(X,Y,Z);
    p(i)=plot(X,Y,'-o');
    %p(i).Marker='*';
    %p(i).LineStyle='none';
    p(i).Color='red';%[1,i/(1.5*N0),0]
    %p(i).legend ='r0';
    hold on
end
for i=1:length(Data(2,:))
    Y=squeeze(Rate1(2,i,:,NC));
    %Z=squeeze(SD_rate(2,i,:,NC));
    %q(i)=errorbar(X,Y,Z);
    q(i)=plot(X,Y,'-o');
    %q(i).Marker='+';
    %q(i).LineStyle='none';
    q(i).Color='magenta';%[1,i/(1.5*N0),0]
    %q(i).legend ='r1';
    hold on
end
for i=1:length(Data(3,:))
    Y=squeeze(Rate1(3,i,:,NC));
    %Z=squeeze(SD_rate(3,i,:,NC));
    %r(i)=errorbar(X,Y,Z);
    r(i)=plot(X,Y,'-o');
    %r(i).Marker='.';
    %r(i).LineStyle='none';
    r(i).Color='blue';%[1,i/(1.5*N0),0]
    %r(i).legend ='r2';
    hold on
end
for i=1:length(Data(4,:))
    Y=squeeze(Rate1(4,i,:,NC));
    %Z=squeeze(SD_rate(4,i,:,NC));
    %s(i)=errorbar(X,Y,Z);
    s(i)=plot(X,Y,'-o');
    %s(i).Marker='s';
    %s(i).LineStyle='none';
    s(i).Color='cyan';%[0,i/N3,1]
    %s(i).legend ='r3';
end

title(['Initial RNAP loading rate in NC',num2str(NC+11)])
xlabel('AP')
ylabel('RNAP loading rate (AU)')
legend('r0','r1','r2','r3')
standardizeFigure_YJK(gca,legend,'red','yellow','cyan','magenta')

% for i=1:4
%     Y=(nansum(squeeze((Rate(i,:,:,NC)))));
%     Mean(i,1:41)=(Y);
% end

% for i=1:4
%     for j=1:41
%         Mean(i,j)=Mean(i,j)./Nb_Approved_Rate(i,j,NC);
%     end
% end

% t=plot(X,Mean(1,:),'-r',X,Mean(2,:),'-m',X,Mean(3,:),'-b',X,Mean(4,:),'-c','LineWidth',2);
% ylim([0 400])
% legend(t,{'r0, mean','r1, mean','r2, mean','r3, mean'})
% xlabel('AP Position (%)')
% ylabel('Initial rate')
% title('Initial rate of r0, r1, r2 and r3 along AP axis, at NC 12')


%standardizeFigure_YJK(gca,legend,'blue')
%saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\Initial_rate/NC12.tif')


%% NC13


NC=2
figure(2)

for i=1:N(1)
    Y=squeeze(Rate(1,i,:,NC));
    Z=squeeze(SD_rate(1,i,:,NC));
    p(i)=errorbar(X,Y,Z);
    p(i).Marker='*';
    p(i).LineStyle='none';
    p(i).Color='red';%[1,i/(1.5*N0),0]
    %p(i).legend ='r0';
    hold on
end
for i=1:N(2)
    Y=squeeze(Rate(2,i,:,NC));
    Z=squeeze(SD_rate(2,i,:,NC));
    q(i)=errorbar(X,Y,Z);
    q(i).Marker='+';
    q(i).LineStyle='none';
    q(i).Color='magenta';%[1,i/(1.5*N0),0]
    %q(i).legend ='r0';
    hold on
end
for i=1:N(3)
    Y=squeeze(Rate(3,i,:,NC));
    Z=squeeze(SD_rate(3,i,:,NC));
    r(i)=errorbar(X,Y,Z);
    r(i).Marker='.';
    r(i).LineStyle='none';
    r(i).Color='blue';%[1,i/(1.5*N0),0]
    %r(i).legend ='r0';
    hold on
end
for i=1:N(4)
    Y=squeeze(Rate(4,i,:,NC));
    Z=squeeze(SD_rate(4,i,:,NC));
    s(i)=errorbar(X,Y,Z);
    s(i).Marker='s';
    s(i).LineStyle='none';
    s(i).Color='cyan';%[0,i/N3,1]
    %s(i).legend ='r3';
end

for i=1:4
    Y=(nansum(squeeze((Rate(i,:,:,NC)))));
    Mean(i,1:41)=(Y);
end

for i=1:4
    for j=1:41
        Mean(i,j)=Mean(i,j)./Nb_Approved_Rate(i,j,NC);
    end
end


save('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Paul_codes\Fitting initial rate/Mean13.mat','Mean')
t=plot(X,Mean(1,:),'-r',X,Mean(2,:),'-m',X,Mean(3,:),'-b',X,Mean(4,:),'-c','LineWidth',2);
ylim([0 400])
legend(t,{'r0, mean','r1, mean','r2, mean','r3, mean'})
xlabel('AP Position (%)')
ylabel('Initial rate')
title('Initial rate of r0, r1, r2 and r3 along AP axis, at NC 13')


standardizeFigure_YJK(gca,legend,'blue')
saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\Initial_rate/NC13.tif')


%NC14
NC=3
figure(3)

for i=1:N(1)
    Y=squeeze(Rate(1,i,:,NC));
    Z=squeeze(SD_rate(1,i,:,NC));
    p(i)=errorbar(X,Y,Z);
    p(i).Marker='*';
    p(i).LineStyle='none';
    p(i).Color='red';%[1,i/(1.5*N0),0]
    %p(i).legend ='r0';
    hold on
end
for i=1:N(2)
    Y=squeeze(Rate(2,i,:,NC));
    Z=squeeze(SD_rate(2,i,:,NC));
    q(i)=errorbar(X,Y,Z);
    q(i).Marker='+';
    q(i).LineStyle='none';
    q(i).Color='magenta';%[1,i/(1.5*N0),0]
    %q(i).legend ='r0';
    hold on
end
for i=1:N(3)
    Y=squeeze(Rate(3,i,:,NC));
    Z=squeeze(SD_rate(3,i,:,NC));
    r(i)=errorbar(X,Y,Z);
    r(i).Marker='.';
    r(i).LineStyle='none';
    r(i).Color='blue';%[1,i/(1.5*N0),0]
    %r(i).legend ='r0';
    hold on
end

for i=1:N(4)
    Y=squeeze(Rate(4,i,:,NC));
    Z=squeeze(SD_rate(4,i,:,NC));
    s(i)=errorbar(X,Y,Z);
    s(i).Marker='s';
    s(i).LineStyle='none';
    s(i).Color='cyan';%[0,i/N3,1]
    %s(i).legend ='r3';
end

for i=1:4
    Y=(nansum(squeeze((Rate(i,:,:,NC)))));
    Mean(i,1:41)=(Y);
end

for i=1:4
    for j=1:41
        Mean(i,j)=Mean(i,j)./Nb_Approved_Rate(i,j,NC);
    end
end

t=plot(X,Mean(1,:),'-r',X,Mean(2,:),'-m',X,Mean(3,:),'-b',X,Mean(4,:),'-c','LineWidth',2);
ylim([0 400])
legend(t,{'r0, mean','r1, mean','r2, mean','r3, mean'})
xlabel('AP Position (%)')
ylabel('Initial rate')
title('Initial rate of r0, r1, r2 and r3 along AP axis, at NC 14')


standardizeFigure_YJK(gca,legend,'blue')
saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\Initial_rate/NC14.tif')

save('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Paul_codes\Fitting initial rate/Mean14.mat','Mean')

%% NC 13
% NC=2
% figure(2)
% for i=1:N0
%     p(i)=errorbar(X,Rate_r0(i,:,NC),SD_r0(i,:,NC))
%     p(i).Marker='*';
%     p(i).LineStyle='none'
%     p(i).Color='red';%[1,i/(1.5*N0),0]
%     hold on
% end
% for j=1:N3
%     q(j)=errorbar(X,Rate_r3(j,:,NC),SD_r3(j,:,NC))
%     q(j).Marker='o';
%     q(j).LineStyle='none'
%     q(j).Color='blue'; %[0,j/N3,1]
% end
% Y0=(nansum(Rate_r0(:,:,NC)));
% Y3=(nansum(Rate_r3(:,:,NC)));
% 
% 
% for i=1:41
% Y0(i)=Y0(i)./N_r0(i,NC);
% Y3(i)=Y3(i)./N_r3(i,NC);
% end
% r=plot(X,Y0,'-r',X,Y3,'-b','LineWidth',2);
% ylim([0 400])
% legend(r,{'r0, mean','r3, mean'})
% xlabel('AP Position (%)')
% ylabel('Initial rate')
% title('Initial rate of r0 and r3 along AP axis, at NC 13')
% 
% standardizeFigure_YJK(gca,legend,'blue')
% saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\Initial_rate/NC13.tif')
% 
% 
% 
% %% NC 14
% NC=3
% figure(3)
% 
% for i=1:N0
%     p(i)=errorbar(X,Rate_r0(i,:,NC),SD_r0(i,:,NC))
%     p(i).Marker='*';
%     p(i).LineStyle='none'
%     p(i).Color='red'; %[1,i/(1.5*N0),0]
%     hold on
% end
% for j=1:N3
%     q(j)=errorbar(X,Rate_r3(j,:,NC),SD_r3(j,:,NC))
%     q(j).Marker='o';
%     q(j).LineStyle='none'
%     q(j).Color='blue';%[0,j/N3,1]
% end
% Y0=(nansum(Rate_r0(:,:,NC)));
% Y3=(nansum(Rate_r3(:,:,NC)));
% 
% 
% for i=1:41
% Y0(i)=Y0(i)./N_r0(i,NC);
% Y3(i)=Y3(i)./N_r3(i,NC);
% end
% r=plot(X,Y0,'-r',X,Y3,'-b','LineWidth',2);
% ylim([0 400])
% legend(r,{'r0, mean','r3, mean'})
% xlabel('AP Position (%)')
% ylabel('Initial rate')
% title('Initial rate of r0 and r3 along AP axis, at NC 14')
% 
% standardizeFigure_YJK(gca,legend,'b','r')
% saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\Initial_rate/NC14.tif')

% %% Plot for T_ON
% % YJK : Here, we need to get average, and SD of T_ON over embryos.
% 
% X=0:0.025:1;
% 
% NC=1
% figure(1)
% for i=1:N0
%     p(i)=plot(X,T0_r0(i,:,NC),'*');
%     p(i).Marker='*';
%     %p(i).LineStyle='none'
%     p(i).Color='red';%[1,i/(1.5*N0),0]
%     %p(i).legend ='r0';
%     hold on
% end
% for j=1:N3
%     q(j)=plot(X,T0_r3(j,:,NC),'o');
%     q(j).Marker='o';
%     %q(j).LineStyle='none'
%     q(j).Color='blue';%[0,j/N3,1]
%     %q(i).legend ='r3';
% end
% Y0=(nansum(T0_r0(:,:,NC)));
% stdY0 = std(T0_r0(:,:,NC));
% Y3=(nansum(T0_r3(:,:,NC)));
% stdY3 = std(T0_r3(:,:,NC));
% 
% 
% for i=1:41
% Y0(i)=Y0(i)./N_r0(i,NC);
% Y3(i)=Y3(i)./N_r3(i,NC);
% end
% 
% errorbar(X,Y0,stdY0,'-r','LineWidth',2);
% errorbar(X,Y3,stdY3,'-b','LineWidth',2);
% %ylim([0 400])
% legend({'r0, mean','r3, mean'})
% xlabel('AP Position (%)')
% ylabel('T_{ON} (min)')
% title('T_{ON} of r0 and r3 along AP axis, at NC 12')
% 
% standardizeFigure_YJK(gca,legend,'blue')
% saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\T_ON/NC12.tif')
% 
% %% NC 13
% NC=2
% figure(2)
% for i=1:N0
%     p(i)=plot(X,T0_r0(i,:,NC),'*');
%     p(i).Marker='*';
%     %p(i).LineStyle='none'
%     p(i).Color='red';%[1,i/(1.5*N0),0]
%     hold on
% end
% for j=1:N3
%     q(j)=plot(X,T0_r3(j,:,NC),'o');
%     q(j).Marker='o';
%     %q(j).LineStyle='none'
%     q(j).Color='blue'; %[0,j/N3,1]
% end
% 
% Y0=(nansum(T0_r0(:,:,NC)));
% stdY0 = std(T0_r0(:,:,NC));
% Y3=(nansum(T0_r3(:,:,NC)));
% stdY3 = std(T0_r3(:,:,NC));
% 
% 
% for i=1:41
% Y0(i)=Y0(i)./N_r0(i,NC);
% Y3(i)=Y3(i)./N_r3(i,NC);
% end
% 
% errorbar(X,Y0,stdY0,'-r','LineWidth',2);
% errorbar(X,Y3,stdY3,'-b','LineWidth',2);
% %ylim([0 400])
% legend({'r0, mean','r3, mean'})
% xlabel('AP Position (%)')
% ylabel('T_{ON} (min)')
% title('T_{ON} of r0 and r3 along AP axis, at NC 13')
% 
% standardizeFigure_YJK(gca,legend,'blue')
% saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\T_ON/NC13.tif')
% 
% 
% 
% %% NC 14
% NC=3
% figure(3)
% 
% for i=1:N0
%     p(i)=plot(X,T0_r0(i,:,NC),'*');
%     p(i).Marker='*';
%     %p(i).LineStyle='none'
%     p(i).Color='red'; %[1,i/(1.5*N0),0]
%     hold on
% end
% for j=1:N3
%     q(j)=plot(X,T0_r3(j,:,NC),'o');
%     q(j).Marker='o';
%     %q(j).LineStyle='none'
%     q(j).Color='blue';%[0,j/N3,1]
% end
% 
% Y0=(nansum(T0_r0(:,:,NC)));
% stdY0 = std(T0_r0(:,:,NC));
% Y3=(nansum(T0_r3(:,:,NC)));
% stdY3 = std(T0_r3(:,:,NC));
% 
% 
% for i=1:41
% Y0(i)=Y0(i)./N_r0(i,NC);
% Y3(i)=Y3(i)./N_r3(i,NC);
% end
% 
% errorbar(X,Y0,stdY0,'-r','LineWidth',2);
% errorbar(X,Y3,stdY3,'-b','LineWidth',2);
% %ylim([0 400])
% legend({'r0, mean','r3, mean'})
% xlabel('AP Position (%)')
% ylabel('T_{ON} (min)')
% title('T_{ON} of r0 and r3 along AP axis, at NC 14')
% 
% standardizeFigure_YJK(gca,legend,'b','r')
% saveas(gcf,'E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Plots-PaulJ\T_ON/NC14.tif')
% 
end