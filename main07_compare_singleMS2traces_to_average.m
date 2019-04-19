function main_07_compare_singleMS2traces_to_average(DataType, varargin)
% DESCRIPTION
% This script is for comparing the individual MS2 traces with the mean, to
% see if the average represent the individual trace features well.

% INPUT :
% 1) Prefix : Prefix = 'Dataset', as defined for our image analysis pipeline.
% 2) DataType : I might need to edit this so that I can compile all single
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

for m=1:length(Data)
    % Initialize some of the variables
    clear compiledParticles
    clear CompiledParticles
    clf
    
    Index = m; % the index of the dataset in the DataStatus.xlsx tab
    Data_Prefix = Data(Index).Prefix; % This should be expanded to all indices of the Datasets.

    compiledParticles = load([DropboxFolder,filesep,Data_Prefix,filesep,'CompiledParticles.mat']);
    CompiledParticles = compiledParticles.CompiledParticles{1,1};

    %% Extract useful fields
    NC = 13; % This should be an option in future.

    % Define nuclear cycles (beginning frame)
    NC12 = compiledParticles.nc12;
    NC13 = compiledParticles.nc13;
    NC14 = compiledParticles.nc14;
    ElapsedTime = compiledParticles.ElapsedTime; % min
    TimeStart = ElapsedTime(NC13); % min

    TimeRange_NC13 = NC13:NC14; % define the frames for NC13



 
    %% Loop through all CompiledParticles, then sort out for specific NC, and AP, 
    % then plot their individual fluo traces, also with MeanVectorAP for that NC, and AP

    l=1;% index of figures
    
    % Loop through all AP bins, then generate plots for each AP bin.
    for i= 9:17 % 20-40% of AP axis % 1:length(compiledParticles.APbinID)-1  

        AP = compiledParticles.APbinID(i) * 100; % percent of EL

        % Only take the AP bins where there are actual values, not Nans.
        if sum(~isnan(compiledParticles.MeanVectorAP{1}(:,i))) > 0 

            % Loop through all particles, then sort out the particles corresponding
            % to that AP bin.
            %clear Particle_trace_figure
            %clear h
            Particle_trace_figure = figure;
            %figure(i)
            hold on

            % Indexing for the number of sorted out particles.
            k=1;
            % color map
            colormap(jet(256)); % or viridis, inferno, magma, or plasma
            cmap=colormap;
            jStart = 1;
            jEnd = 20; %length(CompiledParticles);
            Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

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

                     plot(SpotTime,SpotFluo,'-','color',Color(k-jStart+1,:))
                     k=k+1;  % count the index of particles   
                end
            end
            % Now, plot the Mean spot fluorescence (MeanVectorAP, and
            % SDVectorAP (SEVectorAP) in that AP bin, using errorbar
            shadedErrorBar(ElapsedTime(TimeRange_NC13) - TimeStart,...
                        compiledParticles.MeanVectorAP{1}(TimeRange_NC13,i),...
                        compiledParticles.SDVectorAP{1}(TimeRange_NC13,i)./...
                        sqrt(compiledParticles.NParticlesAP{1}(TimeRange_NC13,i)),'lineprops',{'r','markerfacecolor','k'})
            hold off
            %h(k+1) = H.mainLine;

            xlim([ElapsedTime(NC13)-TimeStart ElapsedTime(NC14)-TimeStart])

            title({'MS2 spot fluorescence : single vs Mean ';[' AP = ',num2str(AP),'%']})
            xlabel('Time into NC13(min)')
            ylabel('Spot fluorescence (AU)')
            %legend([h(k) H.mainLine],'single','Mean')

            % pause
            % Make the figure look better
            StandardFigure(Particle_trace_figure,gca)

            % count the index of figure (for montage)
            l=l+1;
            ImageList(l) = gcf; % Save the current figure in ImageList with Index of l

            % save the figure
            saveas(Particle_trace_figure,[FigPath,filesep, 'Embryo',num2str(Index),'_single_vs_averaged_MS2_traces at ' ,num2str(AP), '%_NC',num2str(NC) , '.pdf']); 

        else
            display(['APbin ',num2str(AP),'% does not have non-Nan values'])
        end
    end
end
%% Make a montage of all images from a single embryo 
% ( for all APbins in NC13, for example)
FigPath = [FigPath,filesep];
fileFolder = FigPath;
dirOutput = dir(fullfile(fileFolder,'Embryo*.pdf'));
fileNames = string({dirOutput.name});

Montage_single_vs_average = montage(fileNames([1 5 9 10 14 18 19 23 27]))%'Size',[3 3]) % ,'Size',[1 9]) % option for m x n matrix
% r1 : [3 5 7  11 13 15 20 22 24] for 25, 30, 35% AP bins

% save the figure
%saveas(Montage_single_vs_average,[FigPath 'Montage_',DataType,'_embryos ', '_NC13' , '.pdf']); 

