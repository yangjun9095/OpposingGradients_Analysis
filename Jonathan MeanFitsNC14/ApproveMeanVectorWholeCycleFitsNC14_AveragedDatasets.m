function ApproveMeanVectorWholeCycleFitsNC14_AveragedDatasets(DataType,varargin)
%Last updated: 4/16/2020 by Jonathan Liu

%Analyzes saved whole cycle NC14 fit results from FitMeanVectorWholeCycleNC14. The user has the option of
%approving or rejecting the results of each AP position's fit. The
%approved/rejected results are saved in the same .mat file as the whole cycle fit results.

%The fit results should be saved as a .mat file with 2 stuctures,
%ResultsPlot and FitResults. Each should contain a
%structure of size 1xN, where N is the number of AP bins in the dataset.
%Additionally, each AP bin possesses a structure field called ApprovedFits.
%The values of ApprovedFits correspond to as follows:
%   1: approved
%   0: uncurated (can be approved if ApproveAll is used in post-analysis)
%   -1: rejected

%Variable input arguments:
% DataType

%% Input arguments
LoadPrefix = false; %By default, user selects which Prefix to load.
UseLocalMovieDatabase = false; %By default, don't use a local MovieDatabase

for i=1:length(varargin)
    if strcmpi(varargin{i},'filepath')
        filepath = varargin{i+1};
    end
end
%% Load results
%Get the default folders for Prefix loading
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders;
        
%Load fit results                                
m = load([filepath,'MeanWholeCycleFitsNC14_',DataType,'.mat']);
N_bins = length(m.ResultsPlot); %Number of AP positions in this dataset

%Extract the data into temporary variables.
ResultsPlot = m.ResultsPlot;
FitResults = m.FitResults;

%Exit out if no fits
if isempty(FitResults)
    error('Fit results are empty.');
end

%% Approve/reject fits

%User keypress options
InputKey = {'a';'r';',';'.';'k';'l';'j';'f';'g';'x'};
Function = {'Approve this fit';'Reject this fit';'Previous particle';...
    'Next particle';'Previous approved particle';'Next approved particle';...
    'Jump to a specific particle';'Flag as bursty';'Unflag as bursty';'Exit and save'};
t = table(InputKey,Function);

%Figure color map
colormap = {'red',[0.94 0.94 0.94],'green'}; %Approved/Uncurated/Rejected colors

running = true; %Keep running program until stopped

%Find first nc and particle that has data
i = 1; %Start with first indexed particle that has data

firstflag = false;
while ~firstflag
    if ~isempty(ResultsPlot(i).t_plot)
        firstflag = true;
    else
        if i < length(ResultsPlot)
            i = i+1;
        else
            warning('These fit results appear to be empty.');
        end
    end
end


close all
f = figure('Name','Linear Fit Results','Position',[200 100 800 600]);

while running
    
clf(f); %Clear figures

%Extract plotting variables from results
t_plot = ResultsPlot(i).t_plot;
MS2_plot = ResultsPlot(i).MS2_plot;
SD_MS2_plot = ResultsPlot(i).SD_MS2_plot;
NParticles_plot = ResultsPlot(i).NParticles_plot;
SE_MS2_plot = SD_MS2_plot./sqrt(NParticles_plot); %standard derror
t_fit = ResultsPlot(i).t_fit;
WholeCycleFit = ResultsPlot(i).WholeCycleFit;

%Check to see if this AP bin/nuclear cycle has any results
emptyflag = false;
if isempty(t_plot)
    emptyflag = true;
end

%Extract inference results
basal = FitResults(i).basal;
t_on = FitResults(i).t_on;
t_dwell = FitResults(i).t_dwell;
t_active = FitResults(i).t_active;
tau = FitResults(i).tau;
rate = FitResults(i).rate;
resnorm = FitResults(i).resnorm;
APPos = FitResults(i).APPos;

%% Plot fit results
figure(f);

if ~emptyflag
    hold on
    errorbar(t_plot,MS2_plot,SE_MS2_plot,'ko','CapSize',0,'DisplayName','Fluorescence data');
    plot(t_fit,WholeCycleFit,'r-','DisplayName','Whole cycle fit');
    hold off

    xlim([t_fit(1), t_fit(end)]);
    line(xlim,basal*[1,1],'Color','magenta','LineStyle','--',...
        'DisplayName','Inferred basal fluorescence');
    line(t_on*[1,1],ylim,'Color','blue','LineStyle','--',...
        'DisplayName','Inferred time on');
    line((t_on+t_dwell)*[1,1],ylim,'Color','green','LineStyle','--',...
        'DisplayName','Inferred time of steady state');
    line((t_on+t_active)*[1,1],ylim,'Color','red','LineStyle','--',...
        'DisplayName','Inferred time of promoter shutoff');
end
xlabel('Time since nuclear cycle start (min)');
ylabel('Fluorescence (AU)');
legend('Location','northeast');
title({['Whole cycle fit: Particle ',num2str(i),' of ',num2str(N_bins),...
    ', AP Position ',num2str(APPos)],...
    ['Onset time = ',num2str(t_on), ' min, Rate = ',num2str(rate),' AU/min'],...
    ['Prefix = ',Prefix],['Dwell time = ',num2str(t_dwell),' min, Active time = ',...
    num2str(t_active),' min, tau = ',num2str(tau),' min']});
f.Color = colormap{FitResults(i).ApprovedFits+2}; %Set color depending on approval status

%% User options (approve/reject, change nucleus)
disp(t); %Display options
exitflag = false; %Loop keypress query until valid exit keypress

while ~exitflag
figure(f);
waitforbuttonpress; %User input to press a key
key = f.CurrentCharacter; %Last pressed key

if strcmp(key,'a')
    FitResults(i).ApprovedFits = 1;
    disp('Approved');
elseif strcmp(key,'r')
    FitResults(i).ApprovedFits = -1;
    disp('Rejected');
elseif strcmp(key,',')
    if i > 1
        i = i - 1;
        disp('Switching to previous AP position');
        exitflag = true;
    elseif i == 1
        disp('Already at first inferred AP position!');
    end
elseif strcmp(key,'.')
    if i < N_bins
        i = i + 1;
        disp('Switching to next AP position');
        exitflag = true;
    elseif i == N_bins
        disp('Already at last AP position!');
    end
elseif strcmp(key,'k')
    if i > 1
        Approved = [FitResults.ApprovedFits];
        Approved = Approved(1:(i-1));
        indToSwitchTo = find(Approved == 1,1,'last');
        if ~isempty(indToSwitchTo)
            i = indToSwitchTo;
            disp('Switching to previous approved particle');
            exitflag = true;
        else
            disp('No approved particles before this one!');
        end
    elseif i == 1
        disp('Already at first inferred particle!');
    end
elseif strcmp(key,'l')
    if i < N_bins
        Approved = [FitResults.ApprovedFits];
        Approved = Approved((i+1):N_bins);
        indToSwitchTo = find(Approved == 1,1,'first') + i;
        if ~isempty(indToSwitchTo)
            i = indToSwitchTo;
            disp('Switching to next approved particle');
            exitflag = true;
        else
            disp('No approved particles after this one!');
        end
    elseif i == N_bins
        disp('Already at first inferred particle!');
    end
elseif strcmp(key,'j')
    j = input('Enter in particle index to jump to:');
    if j >= 1 && j <= N_bins
        i = j;
        disp(['Switching to particle ',num2str(i)]);
        exitflag = true;
    else
        disp(['Error: please enter an integer between 1 and ',num2str(N_bins)]);
    end
elseif strcmp(key,'x')
    disp('Exiting and saving results')
    exitflag = true;
    running = 0;
end
if ~isempty(FitResults(i).ApprovedFits)
    f.Color = colormap{FitResults(i).ApprovedFits+2}; %Set color depending on approval status
end
end

end

%Save approved fits results
m.FitResults = FitResults;
save([DropboxFolder,filesep,Prefix,'\MeanWholeCycleFitsNC14.mat'],'FitResults','ResultsPlot','Prefix');

disp('Results saved.');

close all;

end