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
                numEmbryos_r0] = Average_FractionON('r0-new-female','skipFigure');
        
[FractionON_individual_r1,FractionON_Average_r1,...
            FractionON_Average_Error_r1,...
                numEmbryos_r1] = Average_FractionON('r1-new-female','skipFigure');
        
[FractionON_individual_r2,FractionON_Average_r2,...
            FractionON_Average_Error_r2,...
                numEmbryos_r2] = Average_FractionON('r2-new-female','skipFigure');
        
[FractionON_individual_r3,FractionON_Average_r3,...
            FractionON_Average_Error_r3,...
                numEmbryos_r3] = Average_FractionON('r3-new-female','skipFigure');
        
% [FractionON_individual_r3prime,FractionON_Average_r3prime,...
%             FractionON_Average_Error_r3prime,...
%                 numEmbryos_r3prime] = Average_FractionON('r3prime','skipFigure');

%% Put all processed Fraction ON data into a structure for easier plotting

% 1) individual Fraction ON
FractionON_individual{1} = FractionON_individual_r0;
FractionON_individual{2} = FractionON_individual_r1;
FractionON_individual{3} = FractionON_individual_r2;
FractionON_individual{4} = FractionON_individual_r3;
% FractionON_individual{5} = FractionON_individual_r3prime;

% 2) Fraction ON averaged over multiple embryos
FractionON_Average{1} = FractionON_Average_r0;
FractionON_Average{2} = FractionON_Average_r1;
FractionON_Average{3} = FractionON_Average_r2;
FractionON_Average{4} = FractionON_Average_r3;
% FractionON_Average{5} = FractionON_Average_r3prime;

% 3) Fraction ON SEM (over multiple embryos, SD / sqrt(number of embryos))
FractionON_SEM{1} = FractionON_Average_Error_r0;
FractionON_SEM{2} = FractionON_Average_Error_r1;
FractionON_SEM{3} = FractionON_Average_Error_r2;
FractionON_SEM{4} = FractionON_Average_Error_r3;
% FractionON_SEM{5} = FractionON_Average_Error_r3prime;

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

ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow; colorDict.red; colorDict.brown]; % 4 embryos max. it could be extended easily
%lineColor = ['m', 'b', 'y', 'r', 'g'];
%% Plot the Fraction ON from all constructs for NC12-NC14
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
    xlim([0.2 0.8])
    ylim([0 1.2])
    %legend(h,'r0','r1','r2','r3','r3prime')
    legend('r0','r1','r2','r3')%,'r3prime')
    standardizeFigure_YJK(gca,legend,[])

end
%% Save figures
for i=1:3 % NC
    standardizeFigure_YJK(FractionON_figure(i).CurrentAxes,FractionON_figure(i).legend,[])

end

%% Save the calculated fields

% Define the folder structure
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
DataPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

save([DataPath,filesep,'FractionON_r01233prime_20190406_Preliminary.mat'],...
    'FractionON_Average','FractionON_SEM','FractionON_individual');