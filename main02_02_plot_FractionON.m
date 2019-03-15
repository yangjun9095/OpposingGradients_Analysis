function  main02_02_plot_FractionON(DataType,varargin)
%  main02_02_plot_FractionON(DataType,varargin)
%
% DESCRIPTION
% Compile the Fraction ON (Fraction of Active Nuclei) in different
% definitions, then generate/save the plots.
% Note. Fraction ON is defined as in Garcia, 2013, which means that the
% ratio of active nuclei (whether nuclei was ever turned on during one
% nuclear cycle) to the total number of nuclei. 
% Thus, in some sense, accumrate counting of the number of nuclei is
% important.
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
% Last Updated: 03/11/2019
%
% Documented by: Yang Joon Kim (yjkim90@berkeley.edu)

%[Leave at least one blank line (without a comment symbol) before beginning
%any additonal commenting or beginning the code. Everything contained in
%the inital comment block will be included; everything after the first 
%blank, uncommented line will be excluded. Thus, do not add any blank,
%uncommented lines within the documentation header.]

%[CODE]
%% Load the datasets
Data = LoadMS2Sets(DataType);

%% Fraction ON 

numEmbryos = length(Data);
% Here, I'll make an assumption that if the total number of nuclei in one
% AP bin is too small, then the Fraction can be biased. Thus, I will set
% this threshold as N_nuclei_thresh = 3 or 4;
N_nuclei_thresh = 3;

%(1) Averaging the Fraction ON from each embryo.
FractionON_individual = zeros(41,3,numEmbryos);
for i=1:numEmbryos
    clear N_TotalNuclei
    clear N_filter
    N_filter = Data(i).TotalEllipsesAP > N_nuclei_thresh;
    N_TotalNuclei = Data(i).TotalEllipsesAP > N_nuclei_thresh;
    FractionON_individual(:,:,i) = Data(i).EllipsesOnAP./N_TotalNuclei;
    FractionON_individual(FractionON_individual==inf) = nan;
end
FractionON_Average = nanmean(FractionON_individual,3);
FractionON_Average_Error = nanstd(FractionON_individual,[],3) / sqrt(numEmbryos);
%  %% Check by plotting (FractionON)
% NC = 12;
% hold on
% plot(0:0.025:1,FractionON_individual(:,NC-11,1))
% plot(0:0.025:1,FractionON_individual(:,NC-11,2))
% 
% errorbar(0:0.025:1,FractionON_Average(:,NC-11),FractionON_Average_Error(:,NC-11))
% 
% legend('embryo1','embryo2','Average')
end