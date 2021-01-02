function ApproveSingleParticleLinearFits(varargin)
%Last updated: 3/3/2020 by Jonathan Liu

%Analyzes saved linear fit results from FitSingleParticleLinearRise. The user has the option of
%approving or rejecting the results of each single nucleus fit. The
%approved/rejected results are saved in the same .mat file as the linear fit results.

%The linear fit results should be saved as a .mat file with 2 stuctures,
%ResultsPlot and FitResults. Each should contain a
%structure of size 1xN, where N is the number of AP bins in the dataset.
%Additionally, each AP bin possesses a structure field called ApprovedFits.
%The values of ApprovedFits correspond to as follows:
%   1: approved
%   0: uncurated (can be approved if ApproveAll is used in post-analysis)
%   -1: rejected

%Variable input arguments:
%   Prefix: Prefix string. If none chosen, user has the option to select
%           using a dialog menu.

%% Input arguments
LoadPrefix = true; %By default, user selects which Prefix to load.
for i=1:length(varargin)
    if strcmpi(varargin{i},'Prefix')
        Prefix = varargin{i+1};
        LoadPrefix = false;
    end
end

%% Load results
%Get the default folders for Prefix loading
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders;

if LoadPrefix
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end

%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);
        
%Load MCMC results                                
m = load([DropboxFolder,filesep,Prefix,'\LinearRiseFits.mat']);
N_particles = length(m.ResultsPlot); %Number of particles in this dataset

%Extract the data into temporary variables.
ResultsPlot = m.ResultsPlot;
FitResults = m.FitResults;

%% Approve/reject fits

%User keypress options
InputKey = {'a';'r';',';'.';'j';'x'};
Function = {'Approve this fit';'Reject this fit';'Previous particle';...
    'Next particle';'Jump to a specific particle';'Exit and save'};
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
t_fit = ResultsPlot(i).t_fit;
LinearFit = ResultsPlot(i).LinearFit;

%Check to see if this AP bin/nuclear cycle has any results
emptyflag = false;
if isempty(t_plot)
    emptyflag = true;
end

%Extract inference results
basal = FitResults(i).basal;
t_on = FitResults(i).t_on;
t_steady = FitResults(i).t_steady;
rate = FitResults(i).rate;
resnorm = FitResults(i).resnorm;
MeanAP = FitResults(i).MeanAP;

%% Plot fit results
figure(f);

if ~emptyflag
    hold on
    plot(t_plot,MS2_plot,'ko','DisplayName','Fluorescence data');
    plot(t_fit,LinearFit,'r-','DisplayName','Linear fit');
    hold off

    xlim([t_fit(1), t_fit(end)]);
    line(xlim,basal*[1,1],'Color','magenta','LineStyle','--',...
        'DisplayName','Inferred basal fluorescence');
    line(t_on*[1,1],ylim,'Color','blue','LineStyle','--',...
        'DisplayName','Inferred time on');
    line(t_steady*[1,1],ylim,'Color','green','LineStyle','--',...
        'DisplayName','Inferred time of steady state');
end
xlabel('Time since nuclear cycle start (min)');
ylabel('Fluorescence (AU)');
legend('Location','southeast');
title({['Linear fit: Particle ',num2str(i),' of ',num2str(N_particles),...
    ', AP Position ',num2str(MeanAP)],...
    ['Onset time = ',num2str(t_on), ' min, Rate = ',num2str(rate),' AU/min'],...
    ['Prefix = ',Prefix]});
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
        disp('Switching to previous particle');
        exitflag = true;
    elseif i == 1
        disp('Already at first inferred particle!');
    end
elseif strcmp(key,'.')
    if i < N_particles
        i = i + 1;
        disp('Switching to next particle');
        exitflag = true;
    elseif i == N_particles
        disp('Already at last particle!');
    end
elseif strcmp(key,'j')
    j = input('Enter in particle index to jump to:');
    if j >= 1 && j <= N_particles
        i = j;
        disp(['Switching to particle ',num2str(i)]);
        exitflag = true;
    else
        disp(['Error: please enter an integer between 1 and ',num2str(N_particles)]);
    end
elseif strcmp(key,'x')
    disp('Exiting and saving results')
    exitflag = true;
    running = 0;
end
f.Color = colormap{FitResults(i).ApprovedFits+2}; %Set color depending on approval status
end

end

%Save approved fits results
m.FitResults = FitResults;
save([DropboxFolder,filesep,Prefix,'\LinearRiseFits.mat'],'FitResults','ResultsPlot','Prefix');

disp('Results saved.');

close all;

end