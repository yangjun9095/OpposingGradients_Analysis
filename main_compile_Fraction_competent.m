function main_compile_Fraction_competent

%% Description

%% Load the ID variables
% Load the compiledData
DropboxFolder = 'S:/YangJoon/Dropbox/OpposingGradient';
% load data structure
load([DropboxFolder,filesep,'OpposingGradients_ProcessedData',filesep,'compiledData.mat'])

% compiledData is full up to 15th column.


%% Make a structure that compiles all information
DataTypesForFit = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};

% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};
                
% SAve the DataType and LoadMS2Sets into a structure, named compiledData
% The first row will contain the field names
compiledData{1,16} = 'Fraction_Competent_mean';
compiledData{1,17} = 'Fraction_Competent_SEM';
compiledData{1,18} = 'Fraction_Competent_individual';

% Loop through all DataTypes to calculate the fraction competent for each dataType
% Note that 
for DataType = [1:8,10,11,13:1:16]%1:length(DataTypesForFit)
    vars = {'FractionON_individual','FractionON_Average',...
            'FractionON_Average_SEM','numEmbryos'};
    clear vars
    
    % Calculate the averaged/individual fraction competent using the
    % Average_FractionON function
    [FractionON_individual, FractionON_Average,...
            FractionON_Average_SEM, numEmbryos] = Average_FractionON(compiledData{DataType+1,1})%,'skipFigure');
    compiledData{DataType+1,16} = FractionON_Average;
    compiledData{DataType+1,17} = FractionON_Average_SEM;
    compiledData{DataType+1,18} = FractionON_individual;
end

%% Runt null datasets without Histone markers

DropboxPath = 'S:/YangJoon/Dropbox/OpposingGradient';
filePath = [DropboxPath,filesep,'OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'];

% Runt null datasets ([0,0,0] and [1,1,1])
r0NullData = load([filePath, filesep, 'r0_RuntNull.mat']);
r3NullData = load([filePath, filesep, 'r3_RuntNull.mat']);

%% Extract useful info from the [000]/[111] Runt null datasets
% [000], Runt null
r0Null_Time = r0NullData.ElapsedTime;
r0Null_nc14 = r0NullData.nc14;
r0Null_FractionON_mean = nanmean(r0NullData.FractionON_individual,3);
r0Null_FractionON_STD = nanstd(r0NullData.FractionON_individual,0,3);
[~,~,r0Null_numEmbryos] = size(r0NullData.FractionON_individual);
r0Null_FractionON_SEM = r0Null_FractionON_STD./sqrt(r0Null_numEmbryos);

% Since the fraction on is #spots/APbinArea, we want to rescale this using
% the 95% percentile value of FractionON
[row,col] = size(r0Null_FractionON_mean);
r0Null_scale = prctile(reshape(r0Null_FractionON_mean,[1,row*col]),99);
r0Null_FractionON_mean_rescaled = r0Null_FractionON_mean./r0Null_scale;
r0Null_FractionON_SEM_rescaled = r0Null_FractionON_SEM./r0Null_scale;

% We only have up until ~30 minutes for the [000], Runt nullfor the plotting purpose, let's add some zero values after the time point
% in here, let's hardcode for the [000] datasets, ~40 min into nc14
r0Null_tMax = 45; % min
tRes = median(diff(r0Null_Time));
nLength = r0Null_tMax/tRes;
nLength_add = ceil(nLength - length(r0Null_Time));
[~,nAPbins] = size(r0NullData.FractionON);
vec_add = nan(nLength_add,nAPbins);
tvec = tRes*[1:nLength_add] + r0Null_Time(end);

r0Null_Time_elongated = [r0Null_Time, tvec];
r0Null_FractionON_mean_rescaled_elongated = [r0Null_FractionON_mean_rescaled; vec_add];


% [111], Runt null
r3Null_Time = r3NullData.ElapsedTime;
r3Null_nc14 = r3NullData.nc14;
r3Null_FractionON_mean = nanmean(r3NullData.FractionON_individual,3);
r3Null_FractionON_STD = nanstd(r3NullData.FractionON_individual,0,3);
[~,~,r3Null_numEmbryos] = size(r3NullData.FractionON_individual);
r3Null_FractionON_SEM = r3Null_FractionON_STD./sqrt(r3Null_numEmbryos);

% Since the fraction on is #spots/APbinArea, we want to rescale this using
% the 95% percentile value of FractionON
[row,col] = size(r3Null_FractionON_mean);
r3Null_scale = prctile(reshape(r3Null_FractionON_mean,[1,row*col]),99);
r3Null_FractionON_mean_rescaled = r3Null_FractionON_mean./r3Null_scale;
r3Null_FractionON_SEM_rescaled = r3Null_FractionON_SEM./r3Null_scale;

% Pick the maximum time point when the fracion competent is maximum
% r0 Null
[~, maxIndex] = max(r0Null_FractionON_mean_rescaled);
tMax_r0Null = 10; % 10th

FractionON_mean_r0Null = r0Null_FractionON_mean_rescaled(tMax_r0Null,:);
FractionON_SEM_r0Null = r0Null_FractionON_SEM_rescaled(tMax_r0Null,:);

% Convert zeros to NaNs
FractionON_mean_r0Null(FractionON_mean_r0Null==0) = nan;
FractionON_SEM_r0Null(FractionON_SEM_r0Null==0) = nan;

% r3 Null
[~, maxIndex] = max(r3Null_FractionON_mean_rescaled);
tMax_r3Null = 13; % 10th

FractionON_mean_r3Null = r3Null_FractionON_mean_rescaled(tMax_r3Null,:);
FractionON_SEM_r3Null = r3Null_FractionON_SEM_rescaled(tMax_r3Null,:);

% Convert zeros to NaNs
FractionON_mean_r3Null(FractionON_mean_r3Null==0) = nan;
FractionON_SEM_r3Null(FractionON_SEM_r3Null==0) = nan;

% save into the structure
compiledData{10,16} = FractionON_mean_r0Null';
compiledData{10,17} = FractionON_SEM_r0Null';

compiledData{13,16} = FractionON_mean_r3Null';
compiledData{13,17} = FractionON_SEM_r3Null';

%% Color definition
% This is defining the line color
% We have 8 distinct datasets, with or without Runt protein.
% I think selecting 8 distinguishable color sets, then changing the
% brightness by either adding/subtracting white would be a better idea than
% selecting 16 different color sets.

colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.purple = [171,133,172]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

%colorDict.magenta = [208,109,171]/255;
%colorDict.lightBlue = [115,142,193]/255;
colorDict.lightgreen = [205,214,209]/255;
colorDict.pink = [232,177,157]/255;
colorDict.thickpink = [132,27,69]/255;

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.darkgreen; colorDict.thickpink]; 

% For now, I'll add white (color+[111])/2 to make thinner color (for the
% Runt nulls)























% DESCRIPTION
% Compile the Fraction ON (Fraction of Active Nuclei) in different
% definitions, then generate/save the plots.
% This function uses sub-function called Average_FractionON(DataType), and
% also the default number of nuclei for thresholding is 2, since there
% aren't many nuclei in nc12. This can be improved in future for different
% numbers for different NCs.

% Note. Fraction ON is defined as in Garcia, 2013, which means that the
% ratio of active nuclei (whether nuclei was ever turned on during one
% nuclear cycle) to the total number of nuclei. 
% Thus,  accumrate tracking of the nuclei is important in this calculation
%
% ARGUMENTS
% DataType: Name of the tab in DataStatus that will be compiled. Make sure
% to put "READY" for the datasets that you want to compile. 
% [Any other parameters you have]
%
% OPTIONS
%
% OUTPUT
% [Any output that your script provides with a description of what it is. If 
% applicable, please note where this output is saved.]
%
% Author (contact): Yang Joon Kim (yjkim90@berkeley.edu)
% Created: 03/11/2019
% Last Updated: 04/02/2019
%
% Documented by: Yang Joon Kim (yjkim90@berkeley.edu)

%[Leave at least one blank line (without a comment symbol) before beginning
%any additonal commenting or beginning the code. Everything contained in
%the inital comment block will be included; everything after the first 
%blank, uncommented line will be excluded. Thus, do not add any blank,
%uncommented lines within the documentation header.]

%[CODE]
%% Calculate the Fraction ON using Average_FractionON function
[FractionON_individual_r0,FractionON_Average_r0,...
            FractionON_Average_Error_r0,...
                numEmbryos_r0] = Average_FractionON('r0-new','skipFigure');
        
[FractionON_individual_r1,FractionON_Average_r1,...
            FractionON_Average_Error_r1,...
                numEmbryos_r1] = Average_FractionON('r1-new','skipFigure');
        
[FractionON_individual_r2,FractionON_Average_r2,...
            FractionON_Average_Error_r2,...
                numEmbryos_r2] = Average_FractionON('r2-new','skipFigure');
        
[FractionON_individual_r3,FractionON_Average_r3,...
            FractionON_Average_Error_r3,...
                numEmbryos_r3] = Average_FractionON('r3-new','skipFigure');
        
% [FractionON_individual_r3prime,FractionON_Average_r3prime,...
%             FractionON_Average_Error_r3prime,...
%                 numEmbryos_r3prime] = Average_FractionON('r3prime','skipFigure');

% r1-close
[FractionON_individual_r1_close,FractionON_Average_r1_close,...
            FractionON_Average_Error_r1_close,...
                numEmbryos_r1_close] = Average_FractionON('r1-close','skipFigure');
% r1-mid            
[FractionON_individual_r1_mid,FractionON_Average_r1_mid,...
            FractionON_Average_Error_r1_mid,...
                numEmbryos_r1_mid] = Average_FractionON('r1-mid','skipFigure');

% r2-close
[FractionON_individual_r2_close,FractionON_Average_r2_close,...
            FractionON_Average_Error_r2_close,...
                numEmbryos_r2_close] = Average_FractionON('r2_1+2','skipFigure');
% r2-far
[FractionON_individual_r2_far,FractionON_Average_r2_far,...
            FractionON_Average_Error_r2_far,...
                numEmbryos_r2_far] = Average_FractionON('r2_1+3','skipFigure');

%% Runt nulls

%% Put all processed Fraction ON data into a structure for easier plotting

% 1) individual Fraction ON
FractionON_individual{1} = FractionON_individual_r0;
FractionON_individual{2} = FractionON_individual_r1;
FractionON_individual{3} = FractionON_individual_r2;
FractionON_individual{4} = FractionON_individual_r3;
% FractionON_individual{5} = FractionON_individual_r3prime;
FractionON_individual{5} = FractionON_individual_r1_close;
FractionON_individual{6} = FractionON_individual_r1_mid;
FractionON_individual{7} = FractionON_individual_r2_close;
FractionON_individual{8} = FractionON_individual_r2_far;

% 2) Fraction ON averaged over multiple embryos
FractionON_Average{1} = FractionON_Average_r0;
FractionON_Average{2} = FractionON_Average_r1;
FractionON_Average{3} = FractionON_Average_r2;
FractionON_Average{4} = FractionON_Average_r3;
% FractionON_Average{5} = FractionON_Average_r3prime;
FractionON_Average{5} = FractionON_Average_r1_close;
FractionON_Average{6} = FractionON_Average_r1_mid;
FractionON_Average{7} = FractionON_Average_r2_close;
FractionON_Average{8} = FractionON_Average_r2_far;

% 3) Fraction ON SEM (over multiple embryos, SD / sqrt(number of embryos))
FractionON_SEM{1} = FractionON_Average_Error_r0;
FractionON_SEM{2} = FractionON_Average_Error_r1;
FractionON_SEM{3} = FractionON_Average_Error_r2;
FractionON_SEM{4} = FractionON_Average_Error_r3;
% FractionON_SEM{5} = FractionON_Average_Error_r3prime;
FractionON_SEM{5} = FractionON_Average_Error_r1_close;
FractionON_SEM{6} = FractionON_Average_Error_r1_mid;
FractionON_SEM{7} = FractionON_Average_Error_r2_close;
FractionON_SEM{8} = FractionON_Average_Error_r2_far;

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;
colorDict.lightgreen = [205,214,209]/255;
% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.magenta; colorDict.lightgreen]; 
%% Plot the Fraction ON from all constructs for NC12-NC14
% Define the path to save plots
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\FractionON';

for NC=12:14 % NC
    FractionON_figure(NC-11) = figure(NC-11);
    hold on
    % Loop over different constructs (r0,1,2,3,3')
    for i=1:length(FractionON_individual)
        % Extract the Fraction ON from the structure
        FractionON_individual_temp = FractionON_individual{i};
        FractionON_Average_temp = FractionON_Average{i};
        FractionON_SEM_temp = FractionON_SEM{i};
        
        [~,~,numEmbryos] = size(FractionON_individual_temp);
        
        % individual embryo (Fraction ON)
%         for j=1:numEmbryos
%             plot(0:0.025:1,FractionON_individual_temp(:,NC-11,j),'-o','Color',ColorChoice(i,:))
%         end
        %H(i) = shadedErrorBar(0:0.025:1,FractionON_Average_temp(:,NC-11),FractionON_SEM_temp(:,NC-11),'lineprops',{lineColor(i),'markerfacecolor',ColorChoice(i,:)})
        errorbar(0:0.025:1,FractionON_Average_temp(:,NC-11),FractionON_SEM_temp(:,NC-11),'Color',ColorChoice(i,:))
    end
    %h(i) = H(i).mainLine;

    hold off
    title(['Fraction of Active Nuclei',' @ NC ',num2str(NC)])
    xlabel('AP (EL)')
    ylabel('Fraction ON')
    xlim([0.15 0.6])
    ylim([0 1.2])
    %legend(h,'r0','r1','r2','r3','r3prime')
    legend('r0','r1','r2','r3')%,'r3prime')
    StandardFigure(FractionON_figure(NC-11),FractionON_figure(NC-11).CurrentAxes)
    %standardizeFigure_YJK(gca,legend,[])
    % Save the plots
    saveas(FractionON_figure(NC-11),[FigPath,filesep,'FractionON_NC',num2str(NC),'.tif']); 
    saveas(FractionON_figure(NC-11),[FigPath,filesep,'FractionON_NC',num2str(NC),'.pdf']); 
end

%% r1 variants with r0 (r1, r2, r3) - NC14
hold on
% Loop over different constructs
for i=[1,2,5,6]
    % Extract the Fraction ON from the structure
    FractionON_individual_temp = FractionON_individual{i};
    FractionON_Average_temp = FractionON_Average{i};
    FractionON_SEM_temp = FractionON_SEM{i};

    [~,~,numEmbryos] = size(FractionON_individual_temp);

    % individual embryo (Fraction ON)
%         for j=1:numEmbryos
%             plot(0:0.025:1,FractionON_individual_temp(:,NC-11,j),'-o','Color',ColorChoice(i,:))
%         end
    %H(i) = shadedErrorBar(0:0.025:1,FractionON_Average_temp(:,NC-11),FractionON_SEM_temp(:,NC-11),'lineprops',{lineColor(i),'markerfacecolor',ColorChoice(i,:)})
    errorbar(0:0.025:1,FractionON_Average_temp(:,3),FractionON_SEM_temp(:,3),'Color',ColorChoice(i,:))
end
%h(i) = H(i).mainLine;

hold off
title(['Fraction of Active Nuclei',' @ NC14'])
xlabel('AP (EL)')
ylabel('Fraction ON')
xlim([0.15 0.6])
ylim([0 1.2])
%legend(h,'r0','r1','r2','r3','r3prime')
legend('r0','r1[1,0,0]','r1[0,0,1]','r1[0,1,0]','Location','SouthWest')%,'r3prime')
StandardFigure(gcf,gca)
%standardizeFigure_YJK(gca,legend,[])
% Save the plots
% saveas(gcf,[FigPath,filesep,'FractionON_NC14_AllConstructs','.tif']); 
% saveas(gcf,[FigPath,filesep,'FractionON_NC14_AllConstructs','.pdf']); 

%% r1-mid sanity check
FractionON_r1_mid_NC14 = squeeze(FractionON_individual_r1_mid(:,3,:));
hold on
for i=1:6
    plot(0:0.025:1, FractionON_r1_mid_NC14(:,i))
    pause
    
end
%% Double check individual embryos
% r0
hold on
for i=1:5
    plot(0:0.025:1, squeeze(FractionON_individual_r0(:,3,i)))
    pause
end
%% Save the calculated fields

% Define the folder structure
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
DataPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

save([DataPath,filesep,'FractionON_r0123.mat'],...
    'FractionON_Average','FractionON_SEM','FractionON_individual');

end
