function main_compile_AccumulatedmRNA_AllConstructs
% This is a main script for 
% 1)calculating the accumulated mRNA for different constructs
% : r0,1,2,3, r1 variants, r2 variants
% 2) comparing the ratio between NC13 and NC14. Check whether we can ignore
% the NC13 reasonably.

%% Load datasets into a structure (a master one to save)

% This assumes that the DataType is the name of the constructs in the
% DataStatus.xlsx tab, for example, r0-new

DropboxPath = 'S:/YangJoon/Dropbox/OpposingGradient';
% If already processed these steps, I'll just load the AccumulatedData
load([DropboxPath, filesep, 'OpposingGradients_ProcessedData', filesep, 'AccumulatedData.mat']);

filePath = [DropboxPath,filesep,'OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'];

AccumulatedData{1,1} = 'DataType';

DataTypes = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};

% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};
                

% Define a master structure
for i=1:length(DataTypes)
    AccumulatedData{i+1,1} = DataTypes{i}; % DataType
end

%% Generate a master structure with fields for all DataTypes
% First, initialize the structure with the field names
AccumulatedData{1,2} = 'AccumulatedmRNA';
AccumulatedData{1,3} = 'AccumulatedmRNA_SD';
AccumulatedData{1,4} = 'AccumulatedmRNA_SE';
AccumulatedData{1,5} = 'NC14';
AccumulatedData{1,6} = 'AccumulatedmRNA_individual';
AccumulatedData{1,7} = 'AccumulatedmRNA_SD_individual';

AccumulatedData{1,8} = 'AccumulatedmRNA_tWindow_mean';
AccumulatedData{1,9} = 'AccumulatedmRNA_tWindow_SEM';

AccumulatedData{1,10} = 'AccumulatedmRNA_tWindow_mean_avg';
AccumulatedData{1,11} = 'AccumulatedmRNA_tWindow_SEM_avg';

%% Step0. Run AverageDatasets.m for all DataTypes
% the default folder is AveragedDatasets_Feb2020

%% Step1. Calculate the Accumulated mRNA at each time point
% Caveats : from NC13 to the end time points
% Script : AccumulatedmRNA.m which utilizes the result "DataType.mat" file
% produced from the AverageDatasets.m

% Note.
% Nulls only have NC14, thus having some NaNs at the
% beginning. It's easier to calculate this using means, rather than
% individuals.

MinParticles = 1;

% for loop for all data types
for i=[1,9]%1:length(DataTypes)
    
    % initialize the varialbes
%     vars = 
    
    % calculate the accumulated mRNA using a custom function,
    % accumulatemRNA.m
    [AccumulatedmRNA, AccumulatedmRNA_SD, AccumulatedmRNA_SE,...
     ~, ~, NC14,...
     AccumulatedmRNA_individual, AccumulatedmRNA_SD_individual] =...
     accumulatemRNA(DataTypes{i}, MinParticles);
 
    % Save the result into the structure, into corresponding fields
    AccumulatedData{i+1,2} = AccumulatedmRNA;
    AccumulatedData{i+1,3} = AccumulatedmRNA_SD;
    AccumulatedData{i+1,4} = AccumulatedmRNA_SE;
    AccumulatedData{i+1,5} = NC14;
    AccumulatedData{i+1,6} = AccumulatedmRNA_individual;
    AccumulatedData{i+1,7} = AccumulatedmRNA_SD_individual;

    % Calculate the accumulated mRNA over a time window (NC14-NC14+30 min)
    tStart = NC14;
    
    % end time point (Note that for [000], Runt WT/nulls, we will only
    % consider up to 25min
    if i==1
        tEnd = NC14+round(25/0.6833);
    elseif i==9
        tEnd = round(13/0.6833); % [000], Runt null is 2 min faster
        
    elseif i==10 || i==2
        tEnd = NC14+round(20/0.6833);
    else
        tEnd = NC14+round(20/0.6833);
    end
    
    tWindow = [tStart tEnd];
    
    % Method1. Using individual embryos, then average
    % All the Runt null data starts from NC14, thus no need to subtract at
    % all.
    if i>=9
         accumulatedmRNA_temp = AccumulatedmRNA_individual(tWindow(2),:,:);
    else
    accumulatedmRNA_temp = AccumulatedmRNA_individual(tWindow(2),:,:) - ...
        AccumulatedmRNA_individual(tWindow(1),:,:);
    end

    % Average/SEM over embryos
    [~,~,numEmbryos] = size(accumulatedmRNA_temp);
    accumulatedmRNA_mean_temp = nanmean(accumulatedmRNA_temp,3);
    accumulatedmRNA_SEM_temp = nanstd(accumulatedmRNA_temp,0,3)./sqrt(numEmbryos);

    % Save into a structure
    AccumulatedData{i+1,8} = accumulatedmRNA_mean_temp;
    AccumulatedData{i+1,9} = accumulatedmRNA_SEM_temp;

    % Method2. Using the averaged accumulated mRNA, subtract from one time
    % point to the other.
    
    % All the Runt null data starts from NC14, thus no need to subtract at
    % all.
    if i>=9
         accumulatedmRNA_avg_temp = AccumulatedmRNA(tWindow(2),:);
         accumulatedmRNA_SEM_temp_avg = AccumulatedmRNA_SE(tWindow(2),:);
    else
        accumulatedmRNA_avg_temp = AccumulatedmRNA(tWindow(2),:) - ...
                                        AccumulatedmRNA(tWindow(1),:);
        accumulatedmRNA_SEM_temp_avg = sqrt(AccumulatedmRNA_SE(tWindow(2),:).^2 +...
                                            AccumulatedmRNA_SE(tWindow(1),:).^2);
    end
        % Save into a structure
    AccumulatedData{i+1,10} = accumulatedmRNA_avg_temp;
    AccumulatedData{i+1,11} = accumulatedmRNA_SEM_temp_avg;
end

%% 


%% Save the master structure into the Dropbox

save('S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AccumulatedData.mat',...
        'AccumulatedData')

    
    
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD PARTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end