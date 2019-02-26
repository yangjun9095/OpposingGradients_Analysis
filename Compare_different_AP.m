%% plot at differents AP to see time cycles correlations

r0=load('r0.mat');
r1=load('r1.mat');
r2=load('r2.mat');
r3=load('r3.mat');

%% Build Matrix with right dimension

%Longest Elapsed Time is r0. To automatized the process, compare the
%datasets. 

L = length(r0.ElapsedTime);
ElapsedTimeTotal = r0.ElapsedTime;


r0MeanVectorAP = zeros(L,41);
r0MeanVectorAP(1:length(r0.ElapsedTime),:) = r0.MeanVectorAP;
r0SDVectorAP = zeros(L,41);
r0SDVectorAP(1:length(r0.ElapsedTime),:) = r0.SDVectorAP;

r1MeanVectorAP = zeros(L,41);
r1MeanVectorAP(1:length(r1.ElapsedTime),:) = r1.MeanVectorAP;
r1SDVectorAP = zeros(L,41);
r1SDVectorAP(1:length(r1.ElapsedTime),:) = r1.SDVectorAP;

r2MeanVectorAP = zeros(L,41);
r2MeanVectorAP(1:length(r2.ElapsedTime),:) = r2.MeanVectorAP;
r2SDVectorAP = zeros(L,41);
r2SDVectorAP(1:length(r2.ElapsedTime),:) = r2.SDVectorAP;

r3MeanVectorAP = zeros(L,41);
r3MeanVectorAP(1:length(r3.ElapsedTime),:) = r3.MeanVectorAP;
r3SDVectorAP = zeros(L,41);
r3SDVectorAP(1:length(r3.ElapsedTime),:) = r3.SDVectorAP;


%% Ploting
%Chose an AP position

AP_Pos = 20;

t=ElapsedTimeTotal;
p=plot(t,r0MeanVectorAP(:,AP_Pos),t,r1MeanVectorAP(:,AP_Pos),t,r2MeanVectorAP(:,AP_Pos),t,r3MeanVectorAP(:,AP_Pos),'LineWidth',1);

title('Mean fluorescence over time in different constructs at AP Position 20')
xlabel('Time')
ylabel('Mean fluorescence')
legend('r0','r1','r2','r3')

%saveas(gcf,'AP20.tif')
