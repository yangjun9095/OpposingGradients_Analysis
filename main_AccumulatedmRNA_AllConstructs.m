function main_AccumulatedmRNA_AllConstructs
% This is a main script for 
% 1)calculating the accumulated mRNA for different constructs
% : r0,1,2,3, r1 variants, r2 variants
% 2) comparing the ratio between NC13 and NC14. Check whether we can ignore
% the NC13 reasonably.

%% Step1. Calculate the Accumulated mRNA at each time point
% Caveats : from NC13 to the end time points
% Script : AccumulatedmRNA.m which utilizes the result "DataType.mat" file
% produced from the AverageDatasets.m
dataTypes = {'r0-new','r1-new','r2-new','r3-new',...
                'r1-close','r1-mid','r2_1+2','r2_1+3',...
                'r0_RuntNull','r3_RuntNull'};
% Note.
% Nulls only have NC14, thus having some NaNs at the
% beginning. It's easier to calculate this using means, rather than
% individuals.

MinParticles = 1;

% initialize the structure to save the results
AccumulatedmRNA = cell(10,1);
AccumulatedmRNA_SD = cell(10,1);
AccumulatedmRNA_SE= cell(10,1);
AccumulatedmRNA_individual = cell(10,1);
AccumulatedmRNA_SD_individual = cell(10,1);
% Extract NC frames
NC12 = nan(10,1);
NC13 = nan(10,1);
NC14 = nan(10,1);
% Length of the time frames
Length = nan(10,1);

% for loop for all data types
for i=1:length(dataTypes)
    [AccumulatedmRNA{i,1}, AccumulatedmRNA_SD{i,1}, AccumulatedmRNA_SE{i,1},...
     NC12(i), NC13(i), NC14(i),...
     AccumulatedmRNA_individual{i,1}, AccumulatedmRNA_SD_individual{i,1}] =...
     accumulatemRNA(dataTypes{i}, MinParticles);
    % Define the length of the frames for each dataset
    [Length(i),~] = size(AccumulatedmRNA{i,1});
end


%% Option
% In case I already have run the step above for calculating the accumulated
% mRNA, I'll just load the results from the "OpposingGradients_ProcessedData\AveragedDatasets_Feb2020"
% folder path.
resultPath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020'

% Load the saved variables
% AccumulatedmRNA, SD, SE,  AccumualtedmRNA_individual, SD_individual

%% Step2. Compare the accumulated mRNA at different time points
% Here, the starting point and end point are both important in comparing
% the accumulated mRNA profiles.

% Notes. All of the datasets processed by AverageDatasets.m are
% synchronized with the beginning of the NC13.
% We will start by thinking about several time windows

% Caveat : The null datasets are shorter than the others, as it's taken up
% to 30 minutes or so. We need to take this into account when we compare
% with different datasets.

% 1) [NC14, NC14+10min]
tFrame = 0.6888; %(min) ~41sec/frame
tInterval = round(10/tFrame); % converting 10 minutes into frames
tWindow1 = [NC14 NC14+tInterval];

% 2) [NC14, NC14+20]
tInterval = round(20/tFrame); % converting 20 minutes into frames
tWindow2 = [NC14 NC14+tInterval];


% 3) [NC13, NC14] : during NC13
tWindow3 = [NC13+5 NC14];

% % optional) [NC14, end of measurement] : 
% This is specifically for the nulls, as there's NaNs at the beginning
% often times not exactly until the same frames.
% Length = Length' % transpose
tWindow_end = [NC14 Length];
tWindow_end(9,1) = 9; % [000], Runt null
tWindow_end(10,1) = 7; % [111], Runt null
% For nulls, because the measurement was done beginning of nc14 mostly, so
% we can also grab the end points.

% 4) [NC14, NC14+40] : during NC14 
tInterval = round(40/tFrame); % converting 40 minutes into frames
tWindow4 = [NC14+5 NC14+tInterval];

% 5) [NC13, NC14+40]
% tInterval = round(40/tFrame); % converting 30 minutes into frames
% tWindow5 = [NC13+5 NC14+tInterval];

tWindows{1,1} = tWindow1;
tWindows{2,1} = tWindow2;
tWindows{3,1} = tWindow3;
tWindows{4,1} = tWindow_end;
%tWindows{4,1} = tWindow4;
%tWindows{5,1} = tWindow5;

%% Calculate the accumulated mRNA within a given time window [t1 t2]
% Basically, subtract the accumulated mRNA at t2 with mRNA at t1
% AccumulatedmRNA{i}(t2,:) - AccumulatedmRNA{i}(t1,:)
% In detail, I'll use AccumulatedmRNA_individual (from individual embryos),
% then calculate this net mRNA during a time window, then get average/SEM
% over multiple embryos.

% Create a structure whose rows are constructs (r0-r3), and whose columns
% are different time windows
for i=1:length(dataTypes) % for different constrcuts (rows)
    
    for t=1:length(tWindows)-2 %"-2" is for taking exceptions for Nulls. 
        % for different time windows (column)
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

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
purple = [171,133,172]/255;
colorDict.purple = (4*purple - [1,1,1])/3;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
brown = [179,155,142]/255;
colorDict.brown = (2*brown - [1,1,1])/1;

colorDict.darkgreen = [126,157,144]/255;
colorDict.lightgreen = [205,214,209]/255;
thickpink = [132,27,69]/255;
colorDict.thickpink = (3*thickpink + [1,1,1]) / 4; % adding white

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.magenta; colorDict.thickpink]; 
%% Generate plots for different time windows
% 1) NC14 : NC14+10 mins
APaxis = 0:0.025:1;
figure(1)
hold on
for i=1:4%length(dataTypes)
    errorbar(APaxis, accumulatedmRNA_mean{i,1}, accumulatedmRNA_SEM{i,1},'Color',ColorChoice(i,:),'LineWidth',2)
end

legend('[0,0,0]','[1,0,0]','[0,1,1]','[1,1,1]')
ylim([0 max(accumulatedmRNA_mean{2,1})+10000])
xlim([0.2 0.6])
StandardFigure(gcf,gca)

DropboxPath = 'S:/YangJoon/Dropbox';
% Save the plot
figPath = [DropboxPath,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA\NewConstructs'];
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r0123_10minIntoNC14_20%_60%' , '.tif']); 
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r0123_10minIntoNC14_20%_60%' , '.pdf']); 
%% 2) NC14 : NC14+40mins
% figure(2)
APaxis = 0:0.025:1;
tIndex = 3; % 3rd time window (NC14:NC14+30min)
hold on
for i=1:4%length(dataTypes)
    errorbar(APaxis, accumulatedmRNA_mean{i,tIndex}, accumulatedmRNA_SEM{i,tIndex},'Color',ColorChoice(i,:),'LineWidth',2)
end

legend('[0,0,0]','[1,0,0]','[0,1,1]','[1,1,1]')
ylim([0 max(accumulatedmRNA_mean{2,tIndex})+40000]) % since usually r1([1,0,0]) seems to be higher than r0,2,3.
xlim([0.2 0.6])
StandardFigure(gcf,gca)

DropboxPath = 'S:/YangJoon/Dropbox';
% Save the plot
figPath = [DropboxPath,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA\NewConstructs'];
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r0123_40minIntoNC14_20%_60%' , '.tif']); 
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r0123_40minIntoNC14_20%_60%' , '.pdf']); 
%% NC13 only for all constructs
APaxis = 0:0.025:1;
tIndex = 4; % 
hold on
for i=[1,2,3,4]%1:length(dataTypes)
    errorbar(APaxis, accumulatedmRNA_mean{i,tIndex}, accumulatedmRNA_SEM{i,tIndex},'Color',ColorChoice(i,:),'LineWidth',2)
end
legend('[0,0,0]','[1,0,0]','[0,1,1]','[1,1,1]')
ylim([0 max(accumulatedmRNA_mean{2,tIndex})+20000]) % since usually r1([1,0,0]) seems to be higher than r0,2,3.
xlim([0.2 0.6])
StandardFigure(gcf,gca)

DropboxPath = 'S:/YangJoon/Dropbox';
% Save the plot
figPath = [DropboxPath,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA\NewConstructs'];
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r0123_NC13_20%_60%' , '.tif']); 
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r0123_NC13_20%_60%' , '.pdf']); 
 %% Plot time-evolution of accumulated mRNA
% % Start with r0
% hold on
% for t=1:length(tWindows)
%     errorbar(APaxis, accumulatedmRNA_mean{1,t}, accumulatedmRNA_SEM{1,t})
%     pause
% end

%% Plot for nulls [000] 
APaxis = 0:0.025:1;
tIndex = 4; % 
hold on
for i=[1,9] % [000], index from the AccumulatedmRNA defined on the top
    errorbar(APaxis, accumulatedmRNA_mean{i,tIndex}, accumulatedmRNA_SEM{i,tIndex},'Color',ColorChoice(i,:),'LineWidth',2)
end

% [000], Runt null
% j=9;
% errorbar(APaxis, AccumulatedmRNA{j,1}(end,:),AccumulatedmRNA_SE{j,1}(end,:),'Color',ColorChoice(j,:),'LineWidth',2)

legend('[000]','[000], Runt null')
ylim([0 max(AccumulatedmRNA{i,1}(end,:))]) 
xlim([0.2 0.6])

xlabel('AP axis (EL)')
ylabel('accumulated mRNA (AU)')

StandardFigure(gcf,gca)

DropboxPath = 'S:/YangJoon/Dropbox';
% Save the plot
figPath = [DropboxPath,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA\NewConstructs'];
saveas(gcf,[figPath,filesep,'accumulatedmRNA_[000]_Runtnull_NC14_30min_20%_60%' , '.tif']); 
saveas(gcf,[figPath,filesep,'accumulatedmRNA_[000]_Runtnull_NC14_30min_20%_60%' , '.pdf']); 

%% Plot for nulls [111] 
APaxis = 0:0.025:1;
tIndex = 4; % Check the indexing above.
hold on
for i=4 % [111], index from the AccumulatedmRNA defined on the top
    errorbar(APaxis, accumulatedmRNA_mean{i,tIndex}, accumulatedmRNA_SEM{i,tIndex},'Color',ColorChoice(i,:),'LineWidth',2)
end

%[111], Runt null
% Grab the last time point, as all the nulls are only calculated for NC14.
j=10;
errorbar(APaxis, AccumulatedmRNA{j,1}(end,:),AccumulatedmRNA_SE{j,1}(end,:),'Color',ColorChoice(j,:),'LineWidth',2)

legend('[111]','[111], Runt null')
%ylim([0 max(AccumulatedmRNA{10,1}(end,:))+10000]) 
xlim([0.2 0.6])

xlabel('AP axis (EL)')
ylabel('accumulated mRNA (AU)')

StandardFigure(gcf,gca)

DropboxPath = 'S:/YangJoon/Dropbox';
% Save the plot
%figPath = [DropboxPath,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA\NewConstructs'];
%saveas(gcf,[figPath,filesep,'accumulatedmRNA_[111]_Runtnull_NC14_30min_20%_60%' , '.tif']); 
%saveas(gcf,[figPath,filesep,'accumulatedmRNA_[111]_Runtnull_NC14_30min_20%_60%' , '.pdf']); 

%% Plot for r1 variants
APaxis = 0:0.025:1;
tIndex = 3; % 3rd time window (NC14:NC14+30min)
hold on
for i=[1,2,5,6]
    errorbar(APaxis, accumulatedmRNA_mean{i,tIndex}, accumulatedmRNA_SEM{i,tIndex},'Color',ColorChoice(i,:),'LineWidth',2)
end

legend('[0,0,0]','[1,0,0]','[0,0,1]','[0,1,0]')
ylim([0 max(accumulatedmRNA_mean{5,tIndex})+50000]) % since usually r1-close([0,0,1]) seems to be higher than r0,2,3.
xlim([0.2 0.6])
StandardFigure(gcf,gca)

DropboxPath = 'S:/YangJoon/Dropbox';
% Save the plot
figPath = [DropboxPath,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA\NewConstructs'];
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r1_variants_40minIntoNC14_20%_60%' , '.tif']); 
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r1_variants_40minIntoNC14_20%_60%' , '.pdf']); 
%% Plot for r2 variants
APaxis = 0:0.025:1;
tIndex = 3; % 3rd time window (NC14:NC14+30min)
hold on
for i=[1,3,7,8]
    errorbar(APaxis, accumulatedmRNA_mean{i,tIndex}, accumulatedmRNA_SEM{i,tIndex},'Color',ColorChoice(i,:),'LineWidth',2)
end

legend('[0,0,0]','[0,1,1]','[1,1,0]','[1,0,1]')
ylim([0 max(accumulatedmRNA_mean{5,tIndex})+50000]) % since usually r1-close([0,0,1]) seems to be higher than r0,2,3.
xlim([0.2 0.6])
StandardFigure(gcf,gca)

DropboxPath = 'S:/YangJoon/Dropbox';
% Save the plot
figPath = [DropboxPath,filesep,'Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA\NewConstructs'];
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r2_variants_40minIntoNC14_20%_60%' , '.tif']); 
saveas(gcf,[figPath,filesep,'accumulatedmRNA_r2_variants_40minIntoNC14_20%_60%' , '.pdf']); 





%%










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

%% Plot (sanity check)
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