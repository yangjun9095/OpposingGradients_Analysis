function main_08_compare_MS2Traces_single_vs_mean_MultipleEmbryos(DataType, varargin)
% DESCRIPTION
% This script is for comparing the individual MS2 traces with the mean, to
% see if the average represents the individual trace features well.

% Plus, on top of main08_02 function, this function plots Mean vs single traces
% for multiple embryos (for specific NC, and AP bin) altogether.
% This will be helpful for looking at embryo-to-embryo variability easily.

% INPUT :
% 1) DataType : I might need to edit this so that I can compile all single
% MS2 traces from the same Data Type. 
% OUTPUT :

% OPTIONS : 

% Checking Varargin 
if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'NC')
            NC=varargin{i+1};
%         elseif strcmpi(varargin{i},'Index')
%             Index=varargin{i+1};
        end
    end
end

%% Load the datasets
% Define the folder structure
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

Data = LoadMS2Sets(DataType);
%Index = 1; % Default, the index of the dataset in the DataStatus.xlsx tab

% Generate folder to save the data
DataFolder=['E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\MS2traces_single_average\',DataType,'_40sec'];
mkdir(DataFolder)
FigPath = DataFolder;

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.brown = [179.155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

ColorChoice = [colorDict.blue; colorDict.red; colorDict.green; colorDict.yellow;...
                colorDict.magenta; colorDict.cyan; colorDict.blue; colorDict.red]; % 4 embryos max. it could be extended easily
lineColor = ['b', 'r', 'g', 'y','m','c','b','r'];
%% Sort particles into corresponding AP bin, and NC, for all embryos,
% thus I need to mark which embryo this is from. Embryo index
% Make a cell structure to save the sorted out particles, this has rows of
% NC (13 and 14), and columns of AP bins (2.5% binned)

%Particles_pooled = cell(2,41); 


for m=1:length(Data)
    % Initialize some of the variables
    clear compiledParticles
    clear CompiledParticles
    %clf
    
    Index = m; % the index of the dataset in the DataStatus.xlsx tab
    Data_Prefix = Data(Index).Prefix; % This should be expanded to all indices of the Datasets.

    compiledParticles = load([DropboxFolder,filesep,Data_Prefix,filesep,'CompiledParticles.mat']);
    CompiledParticles = compiledParticles.CompiledParticles{1,1};

    %% Extract useful fields
    NC = 13; % This should be an option in future.

    % Define nuclear cycles (beginning frame)
    try 
        NC12 = compiledParticles.nc12;
    catch
        NC12 = nan;
    end
    NC13 = compiledParticles.nc13;
    NC14 = compiledParticles.nc14;
    ElapsedTime = compiledParticles.ElapsedTime; % min
    TimeStart13 = ElapsedTime(NC13); % min
    TimeStart14 = ElapsedTime(NC14); % min
    
    TimeRange_NC13 = NC13:NC14; % define the frames for NC13
    TimeRange_NC14 = NC14:length(ElapsedTime); % define the frames for NC13
    
    % Define the time range depending on NC
    if NC ==13
        TimeRange = TimeRange_NC13;
        TimeStart = TimeStart13;
    elseif NC == 14
        TimeRange = TimeRange_NC14;
        TimeStart = TimeStart14;
    end

    %% Loop through all CompiledParticles, then sort out for specific NC, and AP, 
    % then plot their individual fluo traces, also with MeanVectorAP for that NC, and AP


    %l=1;% index of figures
    % Loop through all AP bins, then generate plots for each AP bin.
    for i= 9:17 % 20-40% of AP axis % 1:length(compiledParticles.APbinID)-1  %[9,13,17]

        AP = compiledParticles.APbinID(i) * 100; % percent of EL
        
        % Create figure window that could be added in each loop for embryos

        Particle_trace_figure(i-8) = figure(i-8)
        hold on
        
        % Only take the AP bins where there are actual values, not Nans.
        if sum(~isnan(compiledParticles.MeanVectorAP{1}(:,i))) > 0 

            % Loop through all particles, then sort out the particles corresponding
            % to that AP bin.

            % Indexing for the number of sorted out particles.
            % This was for the gradation coloring of single traces
            k=1;
            
            for j=1:length(CompiledParticles)
            % check if the AP bin matches && NC
             % NC 13 only (Note that I'm using the frame index to cut
             % out the ones before NC13, and also after NC13.)
                if CompiledParticles(j).MeanAP>compiledParticles.APbinID(i) && ...
                         CompiledParticles(j).MeanAP<compiledParticles.APbinID(i+1) && ...
                         CompiledParticles(j).Frame(1) > NC13 && ...
                             CompiledParticles(j).Frame(1) < NC14 + 1

                     % Define the fields that are needed for plotting the single
                     % MS2 traces. (Note. This is assuming that the particle
                     % tracking is really great, thus if I see some crappy traces,
                     % then I might need to go back to the particle tracking.)
                     clear Frames
                     clear SpotFluo
                     clear SpotTime
                     Frames = CompiledParticles(j).Frame;
                     SpotTime = ElapsedTime(Frames) - TimeStart; % To start T=0 from the previous mitosis.
                     SpotFluo = CompiledParticles(j).Fluo;

                     plot(SpotTime,SpotFluo,'-','color',ColorChoice(m,:))
                     k=k+1;  % count the index of particles   
                end
            end
            % Now, plot the Mean spot fluorescence (MeanVectorAP, and
            % SDVectorAP (SEVectorAP) in that AP bin, using errorbar

%             shadedErrorBar(ElapsedTime(TimeRange) - TimeStart,...
%                         compiledParticles.MeanVectorAP{1}(TimeRange,i),...
%                         compiledParticles.SDVectorAP{1}(TimeRange,i)./...
%                         sqrt(compiledParticles.NParticlesAP{1}(TimeRange,i)),'lineprops',{lineColor(m),'markerfacecolor',ColorChoice(m,:)})

            errorbar(ElapsedTime(TimeRange) - TimeStart,...
                         compiledParticles.MeanVectorAP{1}(TimeRange,i),...
                         compiledParticles.SDVectorAP{1}(TimeRange,i)./...
                         sqrt(compiledParticles.NParticlesAP{1}(TimeRange,i)),'color',ColorChoice(m,:),'LineWidth',3)
                     
            xlim([ElapsedTime(NC13)-TimeStart ElapsedTime(NC14)-TimeStart])
            %title({'MS2 spot fluorescence : single vs Mean ';[' AP = ',num2str(AP),'%']})
            title({'Mean MS2 Spot fluorescence over time';[' AP = ',num2str(AP),'%']})
            xlabel('Time into NC13(min)')
            ylabel('Spot fluorescence (AU)')

            %l = l+1;

        else
            display(['APbin ',num2str(AP),'% does not have non-Nan values'])

        end
    end
end
hold off

%% Make the figure look better & save
for i=1:length(Particle_trace_figure)
    StandardFigure(Particle_trace_figure(i),Particle_trace_figure(i).CurrentAxes)

    % saveas(Particle_trace_figure(i),[FigPath,filesep, 'Embryos_together','_single_vs_averaged_MS2_traces at ' ,num2str((i+8-1)*2.5), '%_NC',num2str(NC) , '.tif']); 
    % saveas(Particle_trace_figure(i),[FigPath,filesep, 'Embryos_together','_single_vs_averaged_MS2_traces at ' ,num2str((i+8-1)*2.5), '%_NC',num2str(NC) , '.pdf']); 
    saveas(Particle_trace_figure(i),[FigPath,filesep, 'Embryos_together','Mean_MS2_traces at ' ,num2str((i+8-1)*2.5), '%_NC',num2str(NC) , '.tif']); 
    saveas(Particle_trace_figure(i),[FigPath,filesep, 'Embryos_together','Mean_MS2_traces at ' ,num2str((i+8-1)*2.5), '%_NC',num2str(NC) , '.pdf']); 
    %pause
end

%% Make a montage of all images from a single embryo 
% ( for all APbins in NC13, for example)
FigurePath = [FigPath,filesep];
fileFolder = FigurePath;
dirOutput = dir(fullfile(fileFolder,'Embryos_together*.tif'));
fileNames = string({dirOutput.name});

Montage_single_vs_average = montage(fileNames)%'Size',[3 3]) % ,'Size',[1 9]) % option for m x n matrix
% r1 : [3 5 7  11 13 15 20 22 24] for 25, 30, 35% AP bins

% save the figure
saveas(Montage_single_vs_average,[FigurePath 'Montage_',DataType,'_embryos_together ', '_NC13' , '.pdf']); 
end