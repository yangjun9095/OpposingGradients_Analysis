function  [FractionON_individual,FractionON_Average,...
            FractionON_Average_Error, numEmbryos] = Average_FractionON(DataType,varargin)
%
% DESCRIPTION
% Compile the Fraction ON (Fraction of Active Nuclei) in different
% definitions, then generate/save the plots.

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
% N_nuclei_thresh
% skipFigure
% OUTPUT
% [Any output that your script provides with a description of what it is. If 
% applicable, please note where this output is saved.]
%
% Author (contact): Yang Joon Kim (yjkim90@berkeley.edu)
% Created: 03/11/2019
% Last Updated: 03/11/2019
%
% Documented by: Yang Joon Kim (yjkim90@berkeley.edu)

%[Leave at least one blank line (without a comment symbol) before beginning
%any additonal commenting or beginning the code. Everything contained in
%the inital comment block will be included; everything after the first 
%blank, uncommented line will be excluded. Thus, do not add any blank,
%uncommented lines within the documentation header.]

%[CODE]
%% Checking Options
% Checking Varargin 
if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'skipFigure')
            skipFigure=1;
        elseif strcmpi(varargin{i},'N_nuclei_thresh')
            N_nuclei_thresh=varargin{i+1};
        end
    end
end
%% Load the datasets
Data = LoadMS2Sets(DataType,'dontCompare');

%% Fraction ON 

numEmbryos = length(Data);

% Here, I'll make an assumption that if the total number of nuclei in one
% AP bin is too small, then the Fraction can be biased. Thus, I will set
% this threshold as N_nuclei_thresh = 3 or 4;
N_nuclei_thresh = 2;
N_filter = zeros(41,3);
%(1) Averaging the Fraction ON from each embryo.
FractionON_individual = zeros(41,3,numEmbryos);
for i=1:numEmbryos
    clear N_TotalNuclei
    clear N_filter
    N_filter(:,1) = Data(i).TotalEllipsesAP(:,1) > 2;%N_nuclei_thresh;
    N_filter(:,2) = Data(i).TotalEllipsesAP(:,2) > 4;%N_nuclei_thresh;
    N_filter(:,3) = Data(i).TotalEllipsesAP(:,3) > 6;%N_nuclei_thresh;
    
    N_TotalNuclei = Data(i).TotalEllipsesAP.*N_filter;
    FractionON_individual(:,:,i) = Data(i).EllipsesOnAP{1,1}./N_TotalNuclei;
    FractionON_individual(FractionON_individual==inf) = nan;
end
FractionON_Average = nanmean(FractionON_individual,3);
FractionON_Average_Error = nanstd(FractionON_individual,[],3) / sqrt(numEmbryos);
%% Check by plotting (FractionON)
if ~skipFigure
    for NC=12:14 % NC

        figure(NC-11)
        hold on
        for j=1:numEmbryos
            plot(0:0.025:1,FractionON_individual(:,NC-11,j),'-o')
        end
            shadedErrorBar(0:0.025:1,FractionON_Average(:,NC-11),FractionON_Average_Error(:,NC-11))
        hold off
        title(['Fraction of Active Nuclei',' @ NC ',num2str(NC)])
        xlabel('AP (EL)')
        ylabel('Fraction ON')
        xlim([0.2 0.8])
        ylim([0 1.2])
        %legend('embryo1','embryo2','Average')
        StandardFigure(gcf,gca)
    end
end
%% Save plots

%% Save useful fields
savePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
save([savePath,filesep,DataType,'_FractionON.mat'],...
    'FractionON_individual','FractionON_Average','FractionON_Average_Error','numEmbryos','N_nuclei_thresh')
end