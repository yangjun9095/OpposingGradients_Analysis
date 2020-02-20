function main02_02_plot_FractionON
%
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
% % Old r0, for NC14
% [FractionON_individual_r0,FractionON_Average_r0,...
%             FractionON_Average_Error_r0,...
%                 numEmbryos_r0] = Average_FractionON('r0','skipFigure');
        
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
% r2-close (r2-[1,1,0]
[FractionON_individual_r2_close,FractionON_Average_r2_close,...
            FractionON_Average_Error_r2_close,...
                numEmbryos_r2_close] = Average_FractionON('r2_1+2','skipFigure');
% r2-far (r2-[1,0,1])
[FractionON_individual_r2_far,FractionON_Average_r2_far,...
            FractionON_Average_Error_r2_far,...
                numEmbryos_r2_far] = Average_FractionON('r2_1+3','skipFigure');
%% Put all processed Fraction ON data into a structure for easier plotting
% filePath = '/Users/yangjoonkim/Dropbox/OpposingGradient/OpposingGradients_ProcessedData/FractionON_Feb2020'
% load([filePath, filesep, 'r0-new_FractionON.mat'])

% 1) individual Fraction ON
FractionON_individual{1} = FractionON_individual_r0;
FractionON_individual{2} = FractionON_individual_r1;
FractionON_individual{3} = FractionON_individual_r2;
FractionON_individual{4} = FractionON_individual_r3;
% FractionON_individual{5} = FractionON_individual_r3prime;
FractionON_individual{5} = FractionON_individual_r1_close(:,:,1:5);
FractionON_individual{6} = FractionON_individual_r1_mid(:,:,4:6);
FractionON_individual{7} = FractionON_individual_r2_close(:,:,1:3);
FractionON_individual{8} = FractionON_individual_r2_far(:,:,[3,4,6]);

% 2) Fraction ON averaged over multiple embryos
% FractionON_Average{1} = FractionON_Average_r0;
% FractionON_Average{2} = FractionON_Average_r1;
% FractionON_Average{3} = FractionON_Average_r2;
% FractionON_Average{4} = FractionON_Average_r3;
% FractionON_Average{5} = FractionON_Average_r3prime;

% 3) Fraction ON SEM (over multiple embryos, SD / sqrt(number of embryos))
% FractionON_SEM{1} = FractionON_Average_Error_r0;
% FractionON_SEM{2} = FractionON_Average_Error_r1;
% FractionON_SEM{3} = FractionON_Average_Error_r2;
% FractionON_SEM{4} = FractionON_Average_Error_r3;
% FractionON_SEM{5} = FractionON_Average_Error_r3prime;

%% Re-calculate the mean/SEM for selected embryos (the ones that have normal Fraction ON s, not all 1s)
% Loop through all different constructs, calculate the mean/SEM
for i=1:length(FractionON_individual)
    % initialize temporary variables
    clear FractionON_ind_temp
    clear FractionON_mean_temp
    clear FractionON_SEM_temp

    FractionON_ind_temp = FractionON_individual{i}; % APbins(41) x NCs(3) x numEmbryos
    [~,~,numEmbryos] = size(FractionON_ind_temp); % number of embryos
    FractionON_mean_temp = nanmean(FractionON_ind_temp(:,3,:),3); % only extract NC14
    FractionON_SEM_temp = nanstd(FractionON_ind_temp(:,3,:),[],3) ./ sqrt(numEmbryos);
    
    FractionON_Average{i} =  FractionON_mean_temp;
    FractionON_SEM{i} = FractionON_SEM_temp;
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
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;
colorDict.lightgreen = [205,214,209]/255;
colorDict.thickpink = [132,27,69]/255;
% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.magenta; colorDict.thickpink]; 

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

%% Plot different combinations
FigPath = '/Users/yangjoonkim/Dropbox/Garcia Lab/Figures/Opposing Gradients/Data/FractionON';

%% 1) 1 binding site at different positions
APaxis = 0:0.025:1;
hold on
% r0
errorbar(0:0.025:1,FractionON_Average{1},FractionON_SEM{1},'Color',ColorChoice(1,:))
% r1 : [1,0,0]->2, [0,0,1]->5, [0,1,0]->6
errorbar(0:0.025:1,FractionON_Average{2},FractionON_SEM{2},'Color',ColorChoice(2,:))
errorbar(0:0.025:1,FractionON_Average{5},FractionON_SEM{5},'Color',ColorChoice(5,:))
errorbar(0:0.025:1,FractionON_Average{6},FractionON_SEM{6},'Color',ColorChoice(6,:))

errorbar(0:0.025:1,FractionON_Average{4},FractionON_SEM{4},'Color',ColorChoice(4,:))

xlim([0.15 0.6])
xticks([0.15 0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])

legend('r0','r1','r1-close','r1-mid', 'r3')
xlabel('AP Position')
ylabel('fraction of active nuclei')

StandardFigure(gcf, gca)

saveas(gcf,[FigPath, filesep, 'FractionON_r1_variants' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath, filesep, 'FractionON_r1_variants' , '_NC14' , '.pdf']); 

%%  2) 2 binding sites at different positionsAPaxis = 0:0.025:1;
hold on
% r0
errorbar(0:0.025:1,FractionON_Average{1},FractionON_SEM{1},'Color',ColorChoice(1,:))
% r2 : [0,1,1]->3, [1,1,0]->7, [1,0,1]->8
errorbar(0:0.025:1,FractionON_Average{3},FractionON_SEM{3},'Color',ColorChoice(3,:))
errorbar(0:0.025:1,FractionON_Average{7},FractionON_SEM{7},'Color',ColorChoice(7,:))
errorbar(0:0.025:1,FractionON_Average{8},FractionON_SEM{8},'Color',ColorChoice(8,:))

errorbar(0:0.025:1,FractionON_Average{4},FractionON_SEM{4},'Color',ColorChoice(4,:))

xlim([0.15 0.6])
xticks([0.15 0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])

legend('r0','r2','r2-close','r2-far', 'r3')
xlabel('AP Position')
ylabel('fraction of active nuclei')

StandardFigure(gcf, gca)

saveas(gcf,[FigPath, filesep, 'FractionON_r2_variants' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath, filesep, 'FractionON_r2_variants' , '_NC14' , '.pdf']); 

%% Different combinations of 1, 2 sites
%% 1st, [1,0,0] (r1-original) + [0,1,0](r1-mid) -> [1,1,0] (r2_close)
hold on
% r0
errorbar(0:0.025:1,FractionON_Average{1},FractionON_SEM{1},'Color',ColorChoice(1,:))

errorbar(0:0.025:1,FractionON_Average{2},FractionON_SEM{2},'Color',ColorChoice(2,:))
errorbar(0:0.025:1,FractionON_Average{6},FractionON_SEM{6},'Color',ColorChoice(6,:))
errorbar(0:0.025:1,FractionON_Average{7},FractionON_SEM{7},'Color',ColorChoice(7,:))

errorbar(0:0.025:1,FractionON_Average{4},FractionON_SEM{4},'Color',ColorChoice(4,:))

xlim([0.15 0.6])
xticks([0.15 0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])

legend('r0','r1[1,0,0]','r1[0,1,0]','r2[1,1,0]', 'r3')
xlabel('AP Position')
ylabel('fraction of active nuclei')

StandardFigure(gcf, gca)

saveas(gcf,[FigPath, filesep, 'FractionON_1+1=2_[1,0,0]+[0,1,0]=[1,1,0]' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath, filesep, 'FractionON_1+1=2_[1,0,0]+[0,1,0]=[1,1,0]' , '_NC14' , '.pdf']); 
%% 2nd, [1,0,0] (r1-original) + [0,0,1](r1-close) -> [1,0,1] (r2_far)
hold on
% r0
errorbar(0:0.025:1,FractionON_Average{1},FractionON_SEM{1},'Color',ColorChoice(1,:))

errorbar(0:0.025:1,FractionON_Average{2},FractionON_SEM{2},'Color',ColorChoice(2,:))
errorbar(0:0.025:1,FractionON_Average{5},FractionON_SEM{5},'Color',ColorChoice(5,:))
errorbar(0:0.025:1,FractionON_Average{8},FractionON_SEM{8},'Color',ColorChoice(8,:))

errorbar(0:0.025:1,FractionON_Average{4},FractionON_SEM{4},'Color',ColorChoice(4,:))

xlim([0.15 0.6])
xticks([0.15 0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])

legend('r0','r1[1,0,0]','r1[0,0,1]','r2[1,0,1]', 'r3')
xlabel('AP Position')
ylabel('fraction of active nuclei')

StandardFigure(gcf, gca)

saveas(gcf,[FigPath, filesep, 'FractionON_1+1=2_[1,0,0]+[0,0,1]=[1,0,1]' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath, filesep, 'FractionON_1+1=2_[1,0,0]+[0,0,1]=[1,0,1]' , '_NC14' , '.pdf']); 
%% 3rd, [0,1,0] (r1-mid) + [0,0,1](r1-close) -> [0,1,1] (r2-original)
hold on
% r0
errorbar(0:0.025:1,FractionON_Average{1},FractionON_SEM{1},'Color',ColorChoice(1,:))

errorbar(0:0.025:1,FractionON_Average{6},FractionON_SEM{6},'Color',ColorChoice(6,:))
errorbar(0:0.025:1,FractionON_Average{5},FractionON_SEM{5},'Color',ColorChoice(5,:))
errorbar(0:0.025:1,FractionON_Average{3},FractionON_SEM{3},'Color',ColorChoice(3,:))

errorbar(0:0.025:1,FractionON_Average{4},FractionON_SEM{4},'Color',ColorChoice(4,:))

xlim([0.15 0.6])
xticks([0.15 0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])

legend('r0','r1[0,1,0]','r1[0,0,1]','r2[0,1,1]', 'r3')
xlabel('AP Position')
ylabel('fraction of active nuclei')

StandardFigure(gcf, gca)

saveas(gcf,[FigPath, filesep, 'FractionON_1+1=2_[0,1,0]+[0,0,1]=[0,1,1]' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath, filesep, 'FractionON_1+1=2_[0,1,0]+[0,0,1]=[0,1,1]' , '_NC14' , '.pdf']); 
%% 1+2 = 3?
%% [0,1,0](r1-mid) + [1,0,1](r2_1+3) -> [1,1,1] (r3)
hold on
% r0
errorbar(0:0.025:1,FractionON_Average{1},FractionON_SEM{1},'Color',ColorChoice(1,:))

errorbar(0:0.025:1,FractionON_Average{6},FractionON_SEM{6},'Color',ColorChoice(6,:))
errorbar(0:0.025:1,FractionON_Average{8},FractionON_SEM{8},'Color',ColorChoice(8,:))

errorbar(0:0.025:1,FractionON_Average{4},FractionON_SEM{4},'Color',ColorChoice(4,:))

xlim([0.15 0.6])
xticks([0.15 0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])

legend('r0','r1[0,1,0]','r2[1,0,1]','r3')
xlabel('AP Position')
ylabel('fraction of active nuclei')

StandardFigure(gcf, gca)

saveas(gcf,[FigPath, filesep, 'FractionON_1+2=3_[0,1,0]+[1,0,1]=[1,1,1]' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath, filesep, 'FractionON_1+2=3_[0,1,0]+[1,0,1]=[1,1,1]' , '_NC14' , '.pdf']); 
%% [0,0,1](r1-close) + [1,1,0](r2_1+2) -> [1,1,1] (r3)
hold on
% r0
errorbar(0:0.025:1,FractionON_Average{1},FractionON_SEM{1},'Color',ColorChoice(1,:))

errorbar(0:0.025:1,FractionON_Average{5},FractionON_SEM{5},'Color',ColorChoice(5,:))
errorbar(0:0.025:1,FractionON_Average{7},FractionON_SEM{7},'Color',ColorChoice(7,:))

errorbar(0:0.025:1,FractionON_Average{4},FractionON_SEM{4},'Color',ColorChoice(4,:))

xlim([0.15 0.6])
xticks([0.15 0.2 0.3 0.4 0.5 0.6])
ylim([0 1.2])

legend('r0','r1[0,0,1]','r2[1,1,0]','r3')
xlabel('AP Position')
ylabel('fraction of active nuclei')

StandardFigure(gcf, gca)

saveas(gcf,[FigPath, filesep, 'FractionON_1+2=3_[0,0,1]+[1,1,0]=[1,1,1]' , '_NC14' , '.tif']); 
saveas(gcf,[FigPath, filesep, 'FractionON_1+2=3_[0,0,1]+[1,1,0]=[1,1,1]' , '_NC14' , '.pdf']); 

%% Double check individual embryos
% % r0
% hold on
% for i=1:5
%     plot(0:0.025:1, squeeze(FractionON_individual_r0(:,3,i)))
%     pause
% end
%% Save the calculated fields

% Define the folder structure
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
%     DetermineLocalFolders;
% DataPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
DataPath = '/Users/yangjoonkim/Dropbox/OpposingGradient/OpposingGradients_ProcessedData/FractionON_Feb2020';

save([DataPath,filesep,'FractionON_AllConstructs_Feb_2020.mat'],...
    'FractionON_Average','FractionON_SEM','FractionON_individual');

%% Load the calculated fraction on

end
