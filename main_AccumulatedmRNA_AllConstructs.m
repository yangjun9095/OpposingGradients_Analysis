function main02_21_AccumulatedmRNA_AllConstructs
% This is a main script for 
% 1)calculating the accumulated mRNA for different constructs
% : r0,1,2,3, r1 variants, r2 variants
% 2) comparing the ratio between NC13 and NC14. Check whether we can ignore
% the NC13 reasonably.

%% Step1. Calculate the Accumulated mRNA at each time point
% Caveats : from NC13 to the end time points
% Script : AccumulatedmRNA.m which utilizes the result "DataType.mat" file
% produced from the AverageDatasets.m
dataTypes = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3'};
MinParticles = 2;

% initialize the structure to save the results
AccumulatedmRNA = cell(8,1);
AccumulatedmRNA_SD = cell(8,1);
AccumulatedmRNA_SE= cell(8,1);
AccumulatedmRNA_individual = cell(8,1);
AccumulatedmRNA_SD_individual = cell(8,1);
% Extract NC frames
NC12 = nan(8,1);
NC13 = nan(8,1);
NC14 = nan(8,1);

% for loop for all data types
for i=1:length(dataTypes)
    [AccumulatedmRNA{i,1}, AccumulatedmRNA_SD{i,1}, AccumulatedmRNA_SE{i,1},...
     NC12(i), NC13(i), NC14(i),...
     AccumulatedmRNA_individual{i,1}, AccumulatedmRNA_SD_individual{i,1}] =...
     accumulatemRNA(dataTypes{i}, MinParticles);
end

%% Step2. Compare the accumulated mRNA at different time points
% Here, the starting point and end point are both important in comparing
% the accumulated mRNA profiles.

% Notes. All of the datasets processed by AverageDatasets.m are
% synchronized with the beginning of the NC13.
% We will start by thinking about several time windows

% 1) [NC14, NC14+5min]
tFrame = 0.6888; %(min) ~41sec/frame
tInterval = round(5/tFrame); % converting 10 minutes into frames
tWindow1 = [NC14 NC14+tInterval];

% 2) [NC14, NC14+20]
tInterval = round(20/tFrame); % converting 10 minutes into frames
tWindow2 = [NC14 NC14+tInterval];

% 3) [NC14, NC14+30] : during NC14 early-mid
tInterval = round(30/tFrame); % converting 10 minutes into frames
tWindow3 = [NC14 NC14+tInterval];

% 4) [NC13, NC14] : during NC13
tWindow4 = [NC13+5 NC14];

% 5) [NC13, NC14+30]
tInterval = round(30/tFrame); % converting 10 minutes into frames
tWindow5 = [NC13+5 NC14+tInterval];

tWindows{1,1} = tWindow1;
tWindows{2,1} = tWindow2;
tWindows{3,1} = tWindow3;
tWindows{4,1} = tWindow4;
tWindows{5,1} = tWindow5;

%% Calculate the accumulated mRNA within a given time window [t1 t2]
% Basically, subtract the accumulated mRNA at t2 with mRNA at t1
% AccumulatedmRNA{i}(t2,:) - AccumulatedmRNA{i}(t1,:)
% In detail, I'll use AccumulatedmRNA_individual (from individual embryos),
% then calculate this net mRNA during a time window, then get average/SEM
% over multiple embryos.

% Create a structure whose rows are constructs (r0-r3), and whose columns
% are different time windows
for i=1:length(dataTypes) % for different constrcuts (rows)
    
    for t=1:length(tWindows) % for different time windows (column)
        % initialize the variables
        clear accumulatedmRNA_temp
        clear accumulatedmRNA_mean_temp
        clear accumulatedmRNA_SEM_temp
        % Subtract the accumulated mRNA (t2) - (t1), for each AP bin, for each
        % embryo.
        
        % Extract the tWindow from tWindows : tWindows{time_window_index}(construct,:)
        tWindow = tWindows{t}(i,:); %t column, i th embryo
        accumulatedmRNA_temp = AccumulatedmRNA_individual{i}(tWindow(2),:,:) - ...
            AccumulatedmRNA_individual{i}(tWindow(1),:,:);

        % Average/SEM over embryos
        [~,~,numEmbryos] = size(accumulatedmRNA_temp);
        accumulatedmRNA_mean_temp = nanmean(accumulatedmRNA_temp,3);
        accumulatedmRNA_SEM_temp = nanstd(accumulatedmRNA_temp,0,3)./sqrt(numEmbryos);

        % Save into a structure
        accumulatedmRNA_mean{i,t} = accumulatedmRNA_mean_temp;
        accumulatedmRNA_SEM{i,t} = accumulatedmRNA_SEM_temp;
    end
end

%% Generate plots for different time windows
% 1) NC14 : NC14+5 mins
APaxis = 0:0.025:1;
figure(1)
hold on
for i=1:4%length(dataTypes)
    errorbar(APaxis, accumulatedmRNA_mean{i,1}, accumulatedmRNA_SEM{i,1},'Color',ColorChoice(i,:))
end

ylim([0 4500])
StandardFigure(gcf,gca)
%FigPath = ''
%% 2) NC14 : NC14+30mins
% figure(2)
APaxis = 0:0.025:1;
tIndex = 3; % 3rd time window (NC14:NC14+30min)
hold on
for i=1:4%length(dataTypes)
    errorbar(APaxis, accumulatedmRNA_mean{i,tIndex}, accumulatedmRNA_SEM{i,tIndex},'Color',ColorChoice(i,:))
end
StandardFigure(gcf,gca)
%% NC13 only for all constructs
APaxis = 0:0.025:1;
tIndex = 4; % 
hold on
for i=[1,2,3,4]%1:length(dataTypes)
    errorbar(APaxis, accumulatedmRNA_mean{i,tIndex}, accumulatedmRNA_SEM{i,tIndex},'Color',ColorChoice(i,:))
end
StandardFigure(gcf,gca)
%% Plot time-evolution of accumulated mRNA
% Start with r0
hold on
for t=1:length(tWindows)
    errorbar(APaxis, accumulatedmRNA_mean{1,t}, accumulatedmRNA_SEM{1,t})
    pause
end

%% Calculate the accumulated mRNA (t2-t1) in a different way
% Subtract using already averaged profiles, then use the SD to calculate
% the new error bar
% AccumulatedmRNA, SD, and numEmbryos needed
tWindow = tWindows{5};
for i=1:length(dataTypes)
    time = tWindow(i,:);
    Acc_mRNA_tWindow{i} = AccumulatedmRNA{i}(time(2),:) -...
        AccumulatedmRNA{i}(time(1),:);
    Acc_mRNA_tWindow_SD{i} = sqrt((AccumulatedmRNA_SD{i}(time(2),:)).^2 +...
               (AccumulatedmRNA_SD{i}(time(1),:)).^2);
end

%% Plot
APaxis = 0:0.025:1;
hold on
for i=1:length(dataTypes)
    errorbar(APaxis, Acc_mRNA_tWindow{i}, Acc_mRNA_tWindow_SD{i},'Color',ColorChoice(i,:))
    pause
end

%% Part2. Accumulated mRNA in a different way : Use IntegratemRNA script.
%% Load the datasets
Data_r0 = LoadMS2Sets('r0-new','dontCompare')
Data_r1 = LoadMS2Sets('r1-new','dontCompare')
Data_r2 = LoadMS2Sets('r2-new','dontCompare')
Data_r3 = LoadMS2Sets('r3-new','dontCompare')

% Variants
% r1 variants
Data_r1_close = LoadMS2Sets('r1-close','dontCompare')
Data_r1_mid = LoadMS2Sets('r1-mid','dontCompare')
% r2 variants
Data_r2_close = LoadMS2Sets('r2_1+2','dontCompare')
Data_r2_far = LoadMS2Sets('r2_1+3','dontCompare')


%% Integrate mRNA
[TotalProd_r0,TotalProdError_r0,TotalProdN_r0,...
    MeanTotalProd_r0,SDTotalProd_r0,SETotalProd_r0]=IntegratemRNA(Data_r0,2,2)

[TotalProd_r1,TotalProdError_r1,TotalProdN_r1,...
    MeanTotalProd_r1,SDTotalProd_r1,SETotalProd_r1]=IntegratemRNA(Data_r1,2,2)

[TotalProd_r2,TotalProdError_r2,TotalProdN_r2,...
    MeanTotalProd_r2,SDTotalProd_r2,SETotalProd_r2]=IntegratemRNA(Data_r2,2,2)

[TotalProd_r3,TotalProdError_r3,TotalProdN_r3,...
    MeanTotalProd_r3,SDTotalProd_r3,SETotalProd_r3]=IntegratemRNA(Data_r3,2,2)

% r1-variants
[TotalProd_r1_close,TotalProdError_r1_close,TotalProdN_r1_close,...
    MeanTotalProd_r1_close,SDTotalProd_r1_close,SETotalProd_r1_close]=IntegratemRNA(Data_r1_close,2,2)

[TotalProd_r1_mid,TotalProdError_r1_mid,TotalProdN_r1_mid,...
    MeanTotalProd_r1_mid,SDTotalProd_r1_mid,SETotalProd_r1_mid]=IntegratemRNA(Data_r1_mid,2,2)

% r2-variants
[TotalProd_r2_close,TotalProdError_r2_close,TotalProdN_r2_close,...
    MeanTotalProd_r2_close,SDTotalProd_r2_close,SETotalProd_r2_close]=IntegratemRNA(Data_r2_close,2,2)

[TotalProd_r2_far,TotalProdError_r2_far,TotalProdN_r2_far,...
    MeanTotalProd_r2_far,SDTotalProd_r2_far,SETotalProd_r2_far]=IntegratemRNA(Data_r2_far,2,2)

%% Plot
APaxis = 0:0.025:1;
NC = 14;
hold on
errorbar(APaxis, MeanTotalProd_r0(:,NC),SETotalProd_r0(:,NC),'Color',ColorChoice(1,:))
errorbar(APaxis, MeanTotalProd_r1(:,NC),SETotalProd_r1(:,NC),'Color',ColorChoice(2,:))
errorbar(APaxis, MeanTotalProd_r2(:,NC),SETotalProd_r2(:,NC),'Color',ColorChoice(3,:))
errorbar(APaxis, MeanTotalProd_r3(:,NC),SETotalProd_r3(:,NC),'Color',ColorChoice(4,:))

% errorbar(APaxis, MeanTotalProd_r1_close(:,NC),SETotalProd_r1_close(:,NC),'Color',ColorChoice(5,:))
% errorbar(APaxis, MeanTotalProd_r1_mid(:,NC),SETotalProd_r1_mid(:,NC),'Color',ColorChoice(6,:))
% 
% errorbar(APaxis, MeanTotalProd_r2_close(:,NC),SETotalProd_r2_close(:,NC),'Color',ColorChoice(7,:))
% errorbar(APaxis, MeanTotalProd_r2_far(:,NC),SETotalProd_r2_far(:,NC),'Color',ColorChoice(8,:))
%% Calculate the amount during NC13 to NC14
Scale_NC13 = 0.5;
Scale_NC14 = 1;

integratedmRNA_r0 = TotalProd_r0(:,:,13)*Scale_NC13 + ...
                    TotalProd_r0(:,:,14)*Scale_NC14;
                
integratedmRNA_r1 = TotalProd_r1(:,:,13)*Scale_NC13 + ...
                    TotalProd_r1(:,:,14)*Scale_NC14;
                
integratedmRNA_r2 = TotalProd_r2(:,:,13)*Scale_NC13 + ...
                    TotalProd_r2(:,:,14)*Scale_NC14;
                
integratedmRNA_r3 = TotalProd_r3(:,:,13)*Scale_NC13 + ...
                    TotalProd_r3(:,:,14)*Scale_NC14;
                
integratedmRNA_r3prime = TotalProd_r3prime(:,:,13)*Scale_NC13 + ...
                    TotalProd_r3prime(:,:,14)*Scale_NC14;
                
% Error estimation
integratedmRNA_error_r0 = sqrt(TotalProd_r0(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r0(:,:,14).^2*Scale_NC14);
                        
integratedmRNA_error_r1 = sqrt(TotalProd_r1(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r1(:,:,14).^2*Scale_NC14);
                        
integratedmRNA_error_r2 = sqrt(TotalProd_r2(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r2(:,:,14).^2*Scale_NC14);
                        
integratedmRNA_error_r3 = sqrt(TotalProd_r3(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r3(:,:,14).^2*Scale_NC14);
                        
integratedmRNA_error_r3prime = sqrt(TotalProd_r3prime(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r3prime(:,:,14).^2*Scale_NC14);

%% Plot the IntegratedmRNA from all different constructs (NC13 + NC14)


% r0
figure_integratedmRNA_r0 = figure;
hold on
for i=1:length(integratedmRNA_r0(:,1))
    errorbar(0:0.025:1, integratedmRNA_r0(i,:), integratedmRNA_error_r0(i,:))
end

% r1
figure_integratedmRNA_r1 = figure;
hold on
for i=1:length(integratedmRNA_r1(:,1))
    errorbar(0:0.025:1, integratedmRNA_r1(i,:), integratedmRNA_error_r1(i,:))
end       

% r2
figure_integratedmRNA_r2 = figure;
hold on
for i=1:length(integratedmRNA_r2(:,1))
    errorbar(0:0.025:1, integratedmRNA_r2(i,:), integratedmRNA_error_r2(i,:))
end

% r3
figure_integratedmRNA_r3 = figure;
hold on
for i=1:length(integratedmRNA_r3(:,1))
    errorbar(0:0.025:1, integratedmRNA_r3(i,:), integratedmRNA_error_r3(i,:))
end

% r3 prime
figure_integratedmRNA_r3prime = figure;
hold on
for i=1:length(integratedmRNA_r3prime(:,1))
    errorbar(0:0.025:1, integratedmRNA_r3prime(i,:), integratedmRNA_error_r3prime(i,:))
end

%% Optional 
%% Averaging over multiple embryos
Averaged_integratedmRNA_r0 = nanmean(integratedmRNA_r0);
Averaged_integratedmRNA_r1 = nanmean(integratedmRNA_r1);
Averaged_integratedmRNA_r2 = nanmean(integratedmRNA_r2);
Averaged_integratedmRNA_r3 = nanmean(integratedmRNA_r3);
Averaged_integratedmRNA_r3prime = nanmean(integratedmRNA_r3prime);

SEM_integratedmRNA_r0 = nanstd(integratedmRNA_r0,[],1)./sqrt(length(integratedmRNA_r0(:,1)));
SEM_integratedmRNA_r1 = nanstd(integratedmRNA_r1,[],1)./sqrt(length(integratedmRNA_r1(:,1)));
SEM_integratedmRNA_r2 = nanstd(integratedmRNA_r2,[],1)./sqrt(length(integratedmRNA_r2(:,1)));
SEM_integratedmRNA_r3 = nanstd(integratedmRNA_r3,[],1)./sqrt(length(integratedmRNA_r3(:,1)));
SEM_integratedmRNA_r3prime = nanstd(integratedmRNA_r3prime,[],1)./sqrt(length(integratedmRNA_r3prime(:,1)));

%% Plot to check (averaged accumulated mRNA profile)
APaxis = 0:0.025:1;

hold on
errorbar(APaxis, Averaged_integratedmRNA_r0, SEM_integratedmRNA_r0)
errorbar(APaxis, Averaged_integratedmRNA_r1, SEM_integratedmRNA_r1)
errorbar(APaxis, Averaged_integratedmRNA_r2, SEM_integratedmRNA_r2)
errorbar(APaxis, Averaged_integratedmRNA_r3, SEM_integratedmRNA_r3)
errorbar(APaxis, Averaged_integratedmRNA_r3prime, SEM_integratedmRNA_r3prime)

title('Accumulated mRNA')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('0','1','2','3','3(mutated)')

StandardFigure(gcf,gca)

% Save the plots
% File path
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA'
saveas(gcf,[FigPath,filesep,'AccumulatedmRNA_AllConstructs(SmallLab)-MS2MCP-Averaged.tif'])
saveas(gcf,[FigPath,filesep,'AccumulatedmRNA_AllConstructs(SmallLab)-MS2MCP-Averaged.pdf'])


end