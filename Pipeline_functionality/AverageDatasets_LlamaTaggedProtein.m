%% Averaging multiple datasets
function AverageDatasets_NuclearProtein(DataType,varargin)
% Author : Yang Joon Kim (yjkim90@berkeley.edu)
% This is edited from Meghan's CombineMultipleEmbryos.m script

% DESCRIPTION
% This function has input of datatype in DataStatus.xls, grabs all datasets
% in that tab,and calculates 
% 1) Averaged  Nuclear fluorescence (weighted sum), Standard
% Deviation, and the total number of nuclei from multiple embryos in
% (nc12), nc13 and nc14.


% Note. This is assuming that you're interested in nc12, nc13 and nc14.
% If you don't have the whole nc12, you might need to comment that dataset
% out in the DataStatus.xlsx
%
% OPTIONS
% 'NC', N : designate the nuclear cycle to begin avearging, this can be 12,13 or
% 14.

% PARAMETERS
% DataType: This is a string that is identical to the name of the tab in
% dataStatus.xlsx that you wish to analyze.
%
% OUTPUT
% Variables for plotting, or more detailed analysis with 
% 1) the Averaged nuclear fluorescence over time. Save as 'Name of the DataType'.mat file
% (nc12, nc13, nc14, NParticlesAP_BGsubtracted,MeanVectorsAP, SDVectorAP_BGsubtracted, ElapsedTime,
%  AccumulatedmRNA, FractionON, AccumulatedmRNA_FractionON  ) 


[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

Data = LoadMS2Sets(DataType);

% Save path option
savePath = '';

for i=1:length(varargin)
    if strcmpi(varargin{i}, 'savePath')
        savePath = varargin{i+1};
    end
end

NC = 13; % Default
for i = 1:length(varargin)
   if strcmpi(varargin{i}, 'NC')
        NC = varargin{i+1};  
   end
end

numEmbryos=length(Data);

%Find the total number of frames for each embryo
numFrames = zeros(1, numEmbryos);
maxAPIndex = zeros(1, numEmbryos);
maxTime = zeros(1, numEmbryos);
for i = 1:numEmbryos
    if NC~=12
        numFrames(i) = size(Data(i).ElapsedTime, 2);
        nc13(i) = Data(i).nc13;
        length_total_NC(i) = numFrames(i) - nc13(i)+1; % length of frames from nc13 to the last frame
        maxAPIndex(i) = length(Data(i).APbinID);  %Data(i).MaxAPIndex; % Usually 41, in 2.5% binning
        maxTime(i) = Data(i).ElapsedTime(numFrames(i));
    elseif NC==12
        numFrames(i) = size(Data(i).ElapsedTime, 2);
        nc12(i) = Data(i).nc12;
        length_total_NC(i) = numFrames(i) - nc12(i)+1; % length of frames from nc13 to the last frame
        maxAPIndex(i) = length(Data(i).APbinID); %Data(i).MaxAPIndex; % Usually 41, in 2.5% binning
        maxTime(i) = Data(i).ElapsedTime(numFrames(i));
    end
end

%Store the number of AP bins (this should always be 41).
numAPBins = maxAPIndex(1);

%% Synchronize the vectors as the beginning of the (nc 12), 13, and 14
% This could be edited to include even earlier nuclear cycles in the future.
% For now, nc12,nc13 and nc14 might be good enough.

% There are cases that the lengths of nuclear cycles were different (June,2018)
% To average multiple embryos properly, I will separate each nuclear cycle
% separately, since there's no transcription during the mitosis, it's fine
% to synchronize the datasets based on the beginning of the interphase
% (Use the APDivision.mat)

% Define the new ElapsedTime vector for the combined embryo. 
% The new ElapsedTime should start with the beginning of nc12, and also has
% the length of the frames of nc12 ~ nc14 (of the longest dataset for each APbin), all the
% empty values can be plugged with Nans.

% This ElapsedTime variable has evenly spaced time points estimated from diff(ElapsedTime)
% (This is assumption that we took the data with negligible time between serieses, which is pretty fair)

%% For all AP bins, define the new time vector
if NC==12
    Length12 = zeros(numEmbryos,numAPBins);
end
Length13 = zeros(numEmbryos,numAPBins);
Length14 = zeros(numEmbryos,numAPBins);

for j=1:numAPBins
    % For all embryos, go through to find the longest nuclear cycle (number
    % of frames)
    for i=1:numEmbryos
        % Define the nuclear cycle (In case we start from nc 12)
        if NC==12
            nc12(i,j) = Data(i).APDivision(12,j);
        end
            nc13(i,j) = Data(i).APDivision(13,j);
            nc14(i,j) = Data(i).APDivision(14,j);
            
        % Calculate the number of frames (during nc)
        if NC==12 && nc12(i,j)~=0
            Length12 (i,j) = nc13(i,j) - nc12(i,j);
            Length13 (i,j) = nc14(i,j) - nc13(i,j);
            Length14 (i,j) = numFrames(i) - nc14(i,j);
        elseif NC==13 && nc13(i,j)~=0
            Length13 (i,j) = nc14(i,j) - nc13(i,j);
            Length14 (i,j) = numFrames(i) - nc14(i,j);
        end
    end
    % Find the maximum length for each cycle
        numFrames13(j) = max(Length13(:,j));
        numFrames14(j) = max(Length14(:,j));
        TotalFrames(j) = numFrames13(j) + numFrames14(j);
    if NC==12
        numFrames12(j) = max(Length12(:,j));
        TotalFrames(j) = numFrames12(j) + numFrames13(j) + numFrames14(j);
    end
   
end
 
% NewFrameLength = max(TotalFrames);

%% Define empty matrices (filled with Nans)

if NC==12
    L12 = max(max(Length12));
else
    L12 = 0;
end
L13 = max(max(Length13));
L14 = max(max(Length14)); % Get the minimum for now, we can fix this better later.

% NC12 (Optional)
if NC==12
    MeanVectorAP_BGsubtracted_12 = NaN(L12,numAPBins,numEmbryos);
    SDVectorAP_BGsubtracted_12 = NaN(L12,numAPBins,numEmbryos);
    NParticlesAP_BGsubtracted_12 = NaN(L12,numAPBins,numEmbryos);
end

% NC13
MeanVectorAP_BGsubtracted_13 = NaN(L13,numAPBins,numEmbryos);
SDVectorAP_BGsubtracted_13 = NaN(L13,numAPBins,numEmbryos);
NParticlesAP_BGsubtracted_13 = NaN(L13,numAPBins,numEmbryos);

% NC14
MeanVectorAP_BGsubtracted_14 = NaN(L14,numAPBins,numEmbryos);
SDVectorAP_BGsubtracted_14 = NaN(L14,numAPBins,numEmbryos);
NParticlesAP_BGsubtracted_14 = NaN(L14,numAPBins,numEmbryos);


% Total matrices
NewFrameLength = L12+L13+L14;

MeanVectorAP_BGsubtracted = NaN(NewFrameLength,numAPBins,numEmbryos);
SDVectorAP_BGsubtracted = NaN(NewFrameLength,numAPBins,numEmbryos);
NParticlesAP_BGsubtracted = NaN(NewFrameLength,numAPBins,numEmbryos);


% Synchornize all fields as all of them starts from nc 13 (or 12)
% This synchronization should be done for each AP bin, since they might
% have different anaphase time point.

for j=1:numAPBins       
    for i=1:numEmbryos
        if NC==12 && nc12(i,j)~=0
            % Sync the fields for each nc
            % NC12
            MeanVectorAP_BGsubtracted_12(1:L12,j,i) = Data(i).MeanVectorAP_BGsubtracted(nc12(i,j):nc12(i,j)+L12-1,j);
            SDVectorAP_BGsubtracted_12(1:L12,j,i) = Data(i).SDVectorAP_BGsubtracted(nc12(i,j):nc12(i,j)+L12-1,j);
            NParticlesAP_BGsubtracted_12(1:L12,j,i) = Data(i).NParticlesAP_BGsubtracted(nc12(i,j):nc12(i,j)+L12-1,j);
           
            % NC13
            MeanVectorAP_BGsubtracted_13(1:L13,j,i) = Data(i).MeanVectorAP_BGsubtracted(nc13(i,j):nc13(i,j)+L13-1,j);
            SDVectorAP_BGsubtracted_13(1:L13,j,i) = Data(i).SDVectorAP_BGsubtracted(nc13(i,j):nc13(i,j)+L13-1,j);
            NParticlesAP_BGsubtracted_13(1:L13,j,i) = Data(i).NParticlesAP_BGsubtracted(nc13(i,j):nc13(i,j)+L13-1,j);
       
            % NC14                          
            MeanVectorAP_BGsubtracted_14(1:numFrames(i)-nc14(i,j),j,i) = Data(i).MeanVectorAP_BGsubtracted(nc14(i,j):numFrames(i)-1,j);
            SDVectorAP_BGsubtracted_14(1:numFrames(i)-nc14(i,j),j,i) = Data(i).SDVectorAP_BGsubtracted(nc14(i,j):numFrames(i)-1,j);
            NParticlesAP_BGsubtracted_14(1:numFrames(i)-nc14(i,j),j,i) = Data(i).NParticlesAP_BGsubtracted(nc14(i,j):numFrames(i)-1,j);

        elseif nc13(i,j)~=0 && nc14(i,j)~=0

            % NC13
            MeanVectorAP_BGsubtracted_13(1:L13,j,i) = Data(i).MeanVectorAP_BGsubtracted(nc13(i,j):nc13(i,j)+L13-1,j);
            SDVectorAP_BGsubtracted_13(1:L13,j,i) = Data(i).SDVectorAP_BGsubtracted(nc13(i,j):nc13(i,j)+L13-1,j);
            NParticlesAP_BGsubtracted_13(1:L13,j,i) = Data(i).NParticlesAP_BGsubtracted(nc13(i,j):nc13(i,j)+L13-1,j);

            % NC14                          
            MeanVectorAP_BGsubtracted_14(1:numFrames(i)-nc14(i,j),j,i) = Data(i).MeanVectorAP_BGsubtracted(nc14(i,j):numFrames(i)-1,j);
            SDVectorAP_BGsubtracted_14(1:numFrames(i)-nc14(i,j),j,i) = Data(i).SDVectorAP_BGsubtracted(nc14(i,j):numFrames(i)-1,j);
            NParticlesAP_BGsubtracted_14(1:numFrames(i)-nc14(i,j),j,i) = Data(i).NParticlesAP_BGsubtracted(nc14(i,j):numFrames(i)-1,j);
        end
    end
end

% % Take the most frequent value of dT from the ElapsedTime. It's because the
% % dT can be different in case we stopped and restarted the movie.
deltaT = mode(diff(Data(1).ElapsedTime)); 
ElapsedTime = deltaT*(0:NewFrameLength-1);

%% Average all fields at each time point
%% Make Nans as zeros
% MeanVectorAP_BGsubtracted(isnan(MeanVectorAP_BGsubtracted)) = 0;
% SDVectorAP_BGsubtracted(isnan(SDVectorAP_BGsubtracted)) = 0;
% NParticlesAP_BGsubtracted(isnan(NParticlesAP_BGsubtracted)) = 0;

if NC==12
    MeanVectorAP_BGsubtracted_12(isnan(MeanVectorAP_BGsubtracted_12)) = 0;
    SDVectorAP_BGsubtracted_12(isnan(SDVectorAP_BGsubtracted_12)) = 0;
    NParticlesAP_BGsubtracted_12(isnan(NParticlesAP_BGsubtracted_12)) = 0;

end

MeanVectorAP_BGsubtracted_13(isnan(MeanVectorAP_BGsubtracted_13)) = 0;
SDVectorAP_BGsubtracted_13(isnan(SDVectorAP_BGsubtracted_13)) = 0;
NParticlesAP_BGsubtracted_13(isnan(NParticlesAP_BGsubtracted_13)) = 0;

MeanVectorAP_BGsubtracted_14(isnan(MeanVectorAP_BGsubtracted_14)) = 0;
SDVectorAP_BGsubtracted_14(isnan(SDVectorAP_BGsubtracted_14)) = 0;
NParticlesAP_BGsubtracted_14(isnan(NParticlesAP_BGsubtracted_14)) = 0;

%% Concatenate the MeanVectorAP_BGsubtracted_12,13,and 14 into MeanVectorAP_BGsubtracted 
% (same for SD, NParticles, and FractionON_Instant
if NC==12
    MeanVectorAP_BGsubtracted = cat(1,MeanVectorAP_BGsubtracted_12,MeanVectorAP_BGsubtracted_13,MeanVectorAP_BGsubtracted_14);
    SDVectorAP_BGsubtracted = cat(1,SDVectorAP_BGsubtracted_12,SDVectorAP_BGsubtracted_13,SDVectorAP_BGsubtracted_14);
    NParticlesAP_BGsubtracted = cat(1,NParticlesAP_BGsubtracted_12,NParticlesAP_BGsubtracted_13,NParticlesAP_BGsubtracted_14);
elseif NC==13
    MeanVectorAP_BGsubtracted = cat(1,MeanVectorAP_BGsubtracted_13,MeanVectorAP_BGsubtracted_14);
    SDVectorAP_BGsubtracted = cat(1,SDVectorAP_BGsubtracted_13,SDVectorAP_BGsubtracted_14);
    NParticlesAP_BGsubtracted = cat(1,NParticlesAP_BGsubtracted_13,NParticlesAP_BGsubtracted_14);
else 
    warning('This part is left as an option. You can designate earlier cycles by editing this code.')

end

%% Plot to check before the averaging (Save the MeanVectorAP_BGsubtracted, etc. from individual embryos for future plots)
%AP = 10; % You can change this
%hold on
% for i=1:numEmbryos
%     errorbar(ElapsedTime,MeanVectorAP_BGsubtracted(:,AP,i),SDVectorAP_BGsubtracted(:,AP,i))
% end
MeanVectorAP_BGsubtracted_individual= MeanVectorAP_BGsubtracted;
SDVectorAP_BGsubtracted_individual = SDVectorAP_BGsubtracted;
NParticlesAP_BGsubtracted_individual = NParticlesAP_BGsubtracted;

%% Plot to check individual embryos - over time (optional)
%  APbin = 17;
% % 
% hold on
% for i=1:length(Data)
%     errorbar(ElapsedTime, MeanVectorAP_BGsubtracted(:,APbin,i), SDVectorAP_BGsubtracted(:,APbin,i))
%     pause
% end

%% Plot to check individual embryos - over AP (optional)
% APaxis = 0:0.025:1;
% Tpoint = 7;
% 
% hold on
% for i=1:length(Data)
%     errorbar(APaxis, MeanVectorAP_BGsubtracted(Tpoint,:,i), SDVectorAP_BGsubtracted(Tpoint,:,i))
%     pause
% end
%% Averaging
sumMean = zeros(NewFrameLength,numAPBins);
sumSD = zeros(NewFrameLength,numAPBins);
sumNParticles = zeros(NewFrameLength,numAPBins);

for i=1:numEmbryos
    sumMean = sumMean + squeeze(MeanVectorAP_BGsubtracted(:,:,i).*NParticlesAP_BGsubtracted(:,:,i));
    sumSD = sumSD + squeeze(SDVectorAP_BGsubtracted(:,:,i).^2.*NParticlesAP_BGsubtracted(:,:,i));
    sumNParticles = sumNParticles + squeeze(NParticlesAP_BGsubtracted(:,:,i));
end
    
MeanVectorAP_BGsubtractedTemp = sumMean./sumNParticles;
SDVectorAP_BGsubtractedTemp = sqrt(sumSD./sumNParticles);
NParticlesAP_BGsubtractedTemp = sumNParticles;

MeanVectorAP_BGsubtracted = MeanVectorAP_BGsubtractedTemp;
SDVectorAP_BGsubtracted = SDVectorAP_BGsubtractedTemp;
NParticlesAP_BGsubtracted = NParticlesAP_BGsubtractedTemp;
SEVectorAP = SDVectorAP_BGsubtracted./sqrt(NParticlesAP_BGsubtracted); % Standard error of mean (SD / sqrt(number of observation(nuclei))
ElapsedTime = ElapsedTime;

%% Plot averaged Nuclear fluo traces along with individual embryos
% EmbryosLegend = {'Embryo1','Embryo2','Embryo3','Embryo4','Embryo5'};
% EmbryosLegend = EmbryosLegend(1:numEmbryos)
% for AP = 1:41
%     APbin = (AP-1)*2.5;
%     clf
%     hold on
%     for i=1:numEmbryos
%         errorbar(ElapsedTime,MeanVectorAP_BGsubtracted_forPlot(:,AP,i),SDVectorAP_BGsubtracted_forPlot(:,AP,i))
%     end
% 
%     errorbar(ElapsedTime,MeanVectorAP_BGsubtracted(:,AP),SEVectorAP(:,AP))
% 
%     title({'Averaged MS2 spot fluorescence';['@ APbin = ',num2str(APbin),'%']})
%     xlabel('Time (min)')
%     ylabel('MS2 Spot Fluorescence (AU)')
%     legend(EmbryosLegend,'Mean')
%     yMax = max(max(MeanVectorAP_BGsubtracted_forPlot(:,AP,:))) + max(max(SDVectorAP_BGsubtracted_forPlot(:,AP,:)));
%     if yMax==0
%         yMax=1;
%     end
%     ylim([0 yMax])
%     
%     standardizeFigure_YJK(gca,legend,'yellow','cyan','magenta','lightblue','red')
%     saveas(gcf,['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Figures-OpposingGradients\hbP2-r0123-AveragedMS2Traces\',...
%             DataType,'AveragedTraces_IndividualEmbryos-AP=',num2str(APbin),'%'],'tif');
% end

%% Convert Zeros to NaNs
MeanVectorAP_BGsubtracted_individual(MeanVectorAP_BGsubtracted_individual==0) = nan;
SDVectorAP_BGsubtracted_individual(SDVectorAP_BGsubtracted_individual==0) = nan;
NParticlesAP_BGsubtracted_individual(NParticlesAP_BGsubtracted_individual==0) = nan;

MeanVectorAP_BGsubtracted(MeanVectorAP_BGsubtracted==0) = nan;
SDVectorAP_BGsubtracted(SDVectorAP_BGsubtracted==0) = nan;
SEVectorAP(SEVectorAP==0) = nan;
NParticlesAP_BGsubtracted(NParticlesAP_BGsubtracted==0) = nan;
%% Define the fields that needs to be saved
% 1) Nuclear cycle (all AP bins are synchronized)
if NC==12
    nc12 = 1;
    nc13 = nc12 + L12;
elseif NC==13
    nc12 = 0;
    nc13 = 1;
end
    %nc13 = nc12 + L12;
    nc14 = nc13 + L13;
    
% 2) 
%% Save the fields in .mat file
    if ~isempty(savePath)
        save([savePath,filesep,DataType,'_BGsubtracted-Averaged.mat'],...
            'MeanVectorAP_BGsubtracted','SDVectorAP_BGsubtracted','SEVectorAP','NParticlesAP_BGsubtracted','ElapsedTime',...
            'MeanVectorAP_BGsubtracted_individual','SDVectorAP_BGsubtracted_individual','NParticlesAP_BGsubtracted_individual',...
            'nc12', 'nc13', 'nc14')
    else
        warning('Define the File Path in the argument above')
    end
end