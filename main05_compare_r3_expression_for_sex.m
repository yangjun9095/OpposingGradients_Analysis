% Compare the transcription of r3 (hbP2 + 3 Runt binding sites) from the male/female 

%% Load the Datasets
% Eventually, I need to use the LoadMS2Sets, and spreadsheet tabs to pool
% multiple datasets. But for now, I will just manually load the datasets.
r3_Male1 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-11-hbP2-r3-MS2V5-20sec-8\CompiledParticles.mat')
r3_Male2 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-11-hbP2-r3-MS2V5-20sec-9\CompiledParticles.mat')

r3_Female1 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-11-hbP2-r3-MS2V5-20sec-10\CompiledParticles.mat')
r3_Female2 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-11-hbP2-r3-MS2V5-20sec-12\CompiledParticles.mat')
r3_Female3 = load('E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\2018-12-12-hbP2-r3-MS2V5-20sec-15\CompiledParticles.mat')

%% Extract the fields that are needed for plotting

% ElapsedTime
Time_Male1 = r3_Male1.ElapsedTime;
Time_Male2 = r3_Male2.ElapsedTime;
Time_Female1 = r3_Female1.ElapsedTime;
Time_Female2 = r3_Female2.ElapsedTime;
Time_Female3 = r3_Female3.ElapsedTime;

% MeanVectorAP
MeanFluo_Male1 = cell2mat(r3_Male1.MeanVectorAP);
MeanFluo_Male2 = cell2mat(r3_Male2.MeanVectorAP);
MeanFluo_Female1 = cell2mat(r3_Female1.MeanVectorAP);
MeanFluo_Female2 = cell2mat(r3_Female2.MeanVectorAP);
MeanFluo_Female3 = cell2mat(r3_Female3.MeanVectorAP);

% SDVectorAP
SDFluo_Male1 = cell2mat(r3_Male1.SDVectorAP);
SDFluo_Male2 = cell2mat(r3_Male2.SDVectorAP);
SDFluo_Female1 = cell2mat(r3_Female1.SDVectorAP);
SDFluo_Female2 = cell2mat(r3_Female2.SDVectorAP);
SDFluo_Female3 = cell2mat(r3_Female3.SDVectorAP);

% nc
nc13_Male1 = r3_Male1.nc13;
nc13_Male2 = r3_Male2.nc13;
nc13_Female1 = r3_Female1.nc13;
nc13_Female2 = r3_Female2.nc13;
nc13_Female3 = r3_Female3.nc13;

%%
AP = 15;
hold on
errorbar(Time_Male1(nc13_Male1:end) - Time_Male1(nc13_Male1), MeanFluo_Male1(nc13_Male1:end, AP),...
            SDFluo_Male1(nc13_Male1:end, AP));
errorbar(Time_Male2(nc13_Male2:end) - Time_Male2(nc13_Male2), MeanFluo_Male2(nc13_Male2:end, AP),...
            SDFluo_Male2(nc13_Male2:end, AP));
errorbar(Time_Female1(nc13_Female1:end) - Time_Female1(nc13_Female1), MeanFluo_Female1(nc13_Female1:end, AP),...
            SDFluo_Female1(nc13_Female1:end, AP));
errorbar(Time_Female2(nc13_Female2:end) - Time_Female2(nc13_Female2), MeanFluo_Female2(nc13_Female2:end, AP),...
            SDFluo_Female2(nc13_Female2:end, AP));
errorbar(Time_Female3(nc13_Female3:end) - Time_Female3(nc13_Female3), MeanFluo_Female3(nc13_Female3:end, AP),...
            SDFluo_Female3(nc13_Female3:end, AP));
        
title('r3 MS2 spot fluorescence over time (nc13)')
xlabel('Time (min)')
ylabel('MS2 spot fluorescence (AU)')
legend('Male1','Male2','Female1','Female2','Female3')

%% Let's use the LoadMS2Sets
r3Male = LoadMS2Sets('r3-male');
r3Female = LoadMS2Sets('r3-female');

%% (optional) Load the Runt protein dataset to show the gradient simultaneously
RuntProtein_female = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-1min-200Hz-Female-Averaged.mat');
RuntProtein_male = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-1min-200Hz-Male-Averaged.mat');
%% 1) Mean Spot fluorescence
% Use MeanVectorAP
AP =13;
hold on
for i=1:length(r3Male)
    nc13 = r3Male(i).nc13;
    nc14 = r3Male(i).nc14;
    Time = r3Male(i).ElapsedTime;
    MeanSpotFluo = cell2mat(r3Male(i).MeanVectorAP);
    NSpots = cell2mat(r3Male(i).NParticlesAP);
    SDSpotFluo = cell2mat(r3Male(i).SDVectorAP);
    errorbar(Time(nc13:nc14) - Time(nc13), MeanSpotFluo(nc13:nc14,AP),...
                SDSpotFluo(nc13:nc14,AP),'b')
end


for i=1:length(r3Female)
    nc13 = r3Female(i).nc13;
    nc14 = r3Female(i).nc14;
    Time = r3Female(i).ElapsedTime;
    MeanSpotFluo = cell2mat(r3Female(i).MeanVectorAP);
    NSpots = cell2mat(r3Female(i).NParticlesAP);
    SDSpotFluo = cell2mat(r3Female(i).SDVectorAP);
    errorbar(Time(nc13:nc14) - Time(nc13), MeanSpotFluo(nc13:nc14,AP),...
                SDSpotFluo(nc13:nc14,AP),'r')
end
            
    
    
%% 2) Initial rate of RNAP loading (initial slope)
% Do this in two ways. The first way is using the FitMeanAPSymmetric to get the slope
% for the Mean spot fluo, and the second is using the single trace fitting,
% to extract one initial slope per particle, then get the average of those.

%% 3) Fraction ON
% Male
for i=1:length(r3Male)
    FractionON_instant = cell2mat(r3Male(i).NParticlesAP)./8; %number of nuclei in one AP bin
    FractionON_13_Male(i,:) = max(FractionON_instant(r3Male(i).nc13:r3Male(i).nc14,:));
end
FractionON_13_Male(FractionON_13_Male==0) = nan;
% Female
for i=1:length(r3Female)
    FractionON_instant = cell2mat(r3Female(i).NParticlesAP)./8; %number of nuclei in one AP bin
    FractionON_13_Female(i,:) = max(FractionON_instant(r3Female(i).nc13:r3Female(i).nc14,:));
end
FractionON_13_Female(FractionON_13_Female==0) = nan;
FractionON_13_Female(:,8) = nan;
% Plot the Fraction ON for male and female (r3)
%yyaxis left
hold on
plothandle(1) = errorbar(0:0.025:1,nanmean(FractionON_13_Male),nanstd(FractionON_13_Male)./sqrt(length(r3Male)),'-')
plothandle(2) = errorbar(0:0.025:1,nanmean(FractionON_13_Female),nanstd(FractionON_13_Female)./sqrt(length(r3Female)),'-')
xlim([0.15 0.6])
title('Fraction ON')
xlabel('AP')
ylabel('FractionON')
legend('Fraction On_{Male}','Fraction On_{Female}')
standardizeFigure(gca,legend,[])
hold off

% % Plot the Runt concentration profile over AP
yyaxis right
hold on
plothandle(3) = errorbar(0:0.025:1,RuntProtein_female.MeanVectorAP(13,:),RuntProtein_female.SEVectorAP(13,:),'-')
plothandle(4) = errorbar(0:0.025:1,RuntProtein_male.MeanVectorAP(13,:),RuntProtein_male.SEVectorAP(13,:),'-')

ylabel('Runt nuclear fluorescence (AU)')
legend('Fraction On_{Male}','Fraction On_{Female}','Runt_{female}','Runt_{male}')
standardizeFigure(gca,legend,[])
hold off
%% 4) Single traces
% Define AP bin (and it's in nc13)
AP = 15;
% Data = 
% Particles = 
%% 4) Single traces from multiple datasets
% Define AP bin (and it's in nc13)
AP = 15;
k=1; % counting the number of particles
hold on
for i=1:length(r3Female)
    clear Particles
    Particles = cell2mat(r3Female(i).CompiledParticles);
    NC13 = r3Female(i).nc13;
    for j=1:length(Particles)

        if Particles(j).nc==13 && length(Particles(j).Frame)>3 &&...
                Particles(j).MeanAP >(AP-1)*0.025 && Particles(j).MeanAP < (AP*0.025)
            k=k+1
            clear Frame
            clear SpotFluo
            Frame = Particles(i).Frame;
            SpotFluo = Particles(i).Fluo;
            %SpotFluoError = Particles(i).FluoError;
            %errorbar(Frame,SpotFluo,SpotFluoError)
            plot(Frame-NC13+1,SpotFluo)
            pause
        end

    %     NucFrame = Nuclei(find(Nuclei.schnitz ==i)).Frames;
    %     NucFluo = Nuclei(Nuclei.schnitz == i).MaxFluo;
    
        %hold on
        %errorbar(Frame,SpotFluo,SpotFluoError)
        %errorbar(Frame,SpotFluo,SpotFluoError)
        %plot(NucFrame,NucFluo)
        %pause
    end
end

%% Plot single MS2 spot trace at specific AP bin along with Mean
% We will start with looking at one AP bin, and one nc

%AP = 10;
Data.MeanVectorAP = cell2mat(Data.MeanVectorAP);
Data.SDVectorAP = cell2mat(Data.SDVectorAP);
Data.NParticlesAP = cell2mat(Data.NParticlesAP);

%%
nc = 13;
for AP = 1:41
    APbin = (AP-1)*0.025; % Percentage. For particles, we filter by (APbin - 0.025) <APPosition < APbin
    if ~isnan(Data.APbinArea(AP))
%       nc = 13;
%         if nc==13
            clf
            hold on
            for i=1:length(Particles)
                if Particles(i).MeanAP < APbin && Particles(i).MeanAP > APbin-0.025
                    if Particles(i).nc == nc && length(Particles(i).Frame)>2
                        plot(Data.ElapsedTime(Particles(i).Frame), Particles(i).Fluo,'-o')
                        %pause
                    end
                end
            end

            % Plot the MeanVectorAP at that AP bin for that nc
            if nc==12
                NC = Data.nc12:Data.nc13;
            elseif nc==13
                NC = Data.nc13:Data.nc14;
            elseif nc==14
                NC = Data.nc14:length(Data.ElapsedTime);
            elseif nc==11
                NC = Data.nc11:Data.nc12;
            else
                warning('Expand this code for the earlier nc')
            end

            errorbar(Data.ElapsedTime(NC),Data.MeanVectorAP(NC,AP),Data.SDVectorAP(NC,AP),'r')
            % Label the figure properly and save
            title ({['Individual MS2 Spot Fluo traces'];[' @ APbin =',num2str(APbin*100),'%',' NC=',num2str(nc)]})
            xlabel('Time (min)')
            ylabel('MS2 Spot Fluorescence (AU)')
            %standardizeFigure_YJK(gca,[],[])
            pause
            %saveas(gcf,['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Figures-OpposingGradients\hbP2-r0123-AveragedMS2Traces\SingleMS2TracesWithAverage-20180509-hbP2-r2\',...
            %            'AP=',num2str(AP),' NC=',num2str(nc)],'tif')
        %end
    end
end
%% 5) Averaged MS2 traces
AverageDatasets('r3-male','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');
AverageDatasets('r3-female','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient');

%% 
Male_Averaged = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3-male.mat');
Female_Averaged = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3-female.mat');

%% plot
AP = 17;
hold on
errorbar(Male_Averaged.ElapsedTime(Male_Averaged.nc13:Male_Averaged.nc14),...
            Male_Averaged.MeanVectorAP(Male_Averaged.nc13:Male_Averaged.nc14,AP),...
            Male_Averaged.SEVectorAP(Male_Averaged.nc13:Male_Averaged.nc14,AP))
errorbar(Female_Averaged.ElapsedTime(Female_Averaged.nc13:Female_Averaged.nc14),...
            Female_Averaged.MeanVectorAP(Female_Averaged.nc13:Female_Averaged.nc14,AP),...
            Female_Averaged.SEVectorAP(Female_Averaged.nc13:Female_Averaged.nc14,AP))
xlim([0 14])
title(['Averaged MS2 spot fluorescence over NC 13',' @ AP = ',num2str((AP-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Averaged MS2 spot fluorescence (AU)')
legend('Male','Female')
standardizeFigure(gca,legend,'red','blue')

%% 6) Fitted to MeanVectorAP (FitMeanAPSymmetric)

% male
%% Initial rate of RNAP loading (Fitted by FitMeanAPSymmetric) (written by Yang Joon Kim)
% Last Updated : 12/29/2018

% Load the fitted fields
dataset = 'r3-male'
data = LoadMS2Sets(dataset);
nSets = length(data);
ap = data(1).APbinID;
numAPBins = length(ap);
Prefix = cell(1, nSets);

%% First, let's only load nonzero RateFits, and only nc13 
% this NC should be made as optional later.
RateFit = nan(numAPBins,nSets);
SDRateFit = nan(numAPBins,nSets);
for i=1:nSets
    %MeanFits = load(['E:\YangJoon\LivemRNA\Data\DynamicsResults\',data(i).Prefix,'\MeanFitsV2.mat']);
    MeanFits = data(i).MeanFits; % manual fit (tiltedTrapezoid)
    for j=1:numAPBins % 41
        %MeanFitsTemp = MeanFits(j,2);
        if ~isempty(MeanFits(j,2).RateFit) && ~MeanFits(j,2).RateFit == 0
            RateFit(j,i) = MeanFits(j,2).RateFit;
            SDRateFit(j,i) = MeanFits(j,2).SDRateFit;
%             if isempty(MeanFits(j,2).SDRateFit)
%                 SDRateFit(j,i) = nan;
%             else
%                 SDRateFit(j,i) = MeanFits(j,2).SDRateFit;
%             end
        end
    end
end

% Average the fitted rates over multiple embryos
AvgRateFit = nanmean(RateFit,2);
STDRateFit = nanstd(RateFit,0,2);
%% plot
hold on
for i=1:nSets
    errorbar(ap,RateFit(:,i),SDRateFit(:,i))
end
% Averaged rate, with std
errorbar(ap,AvgRateFit(:),STDRateFit(:))

standardizeFigure(gca,legend,[])

title('Initial rate of RNAP loading over AP')
xlabel('AP')
ylabel('Inital rate of RNAP loading')
legend('embryo1','embryo2','embryo3','embryo4','Average')

%% male data saved temporarily
RateFit_male = AvgRateFit;
STDRateFit_male = STDRateFit;
NDatasets_male = nSets;
% female
%% Initial rate of RNAP loading (Fitted by FitMeanAPSymmetric) (written by Yang Joon Kim)
% Last Updated : 12/29/2018

% Load the fitted fields
dataset = 'r3-female'
data = LoadMS2Sets(dataset);
nSets = length(data);
ap = data(1).APbinID;
numAPBins = length(ap);
Prefix = cell(1, nSets);

%% First, let's only load nonzero RateFits, and only nc13 
% this NC should be made as optional later.
RateFit = nan(numAPBins,nSets);
SDRateFit = nan(numAPBins,nSets);
for i=1:nSets
    %MeanFits = load(['E:\YangJoon\LivemRNA\Data\DynamicsResults\',data(i).Prefix,'\MeanFitsV2.mat']);
    %MeanFits = MeanFits.FitResults;
    MeanFits = data(i).MeanFits; % manual fit (tiltedTrapezoid)
    for j=1:numAPBins % 41
        %MeanFitsTemp = MeanFits(j,2);
        if ~isempty(MeanFits(j,2).RateFit) && ~MeanFits(j,2).RateFit == 0
            RateFit(j,i) = MeanFits(j,2).RateFit;
            SDRateFit(j,i) = MeanFits(j,2).SDRateFit;
%             if isempty(MeanFits(j,2).SDRateFit1)
%                 SDRateFit(j,i) = nan;
%             else
%                 SDRateFit(j,i) = MeanFits(j,2).SDRateFit1;
%             end
        end
    end
end

% Average the fitted rates over multiple embryos
AvgRateFit = nanmean(RateFit,2);
STDRateFit = nanstd(RateFit,0,2);
%% plot
hold on
for i=1:nSets
    errorbar(ap,RateFit(:,i),SDRateFit(:,i))
end
% Averaged rate, with std
errorbar(ap,AvgRateFit(:),STDRateFit(:))

standardizeFigure(gca,legend,[])

title('Initial rate of RNAP loading over AP')
xlabel('AP')
ylabel('Inital rate of RNAP loading')
%legend('embryo1','embryo2','embryo3','embryo4','Average')
%% female data saved temporarily
RateFit_female = AvgRateFit;
STDRateFit_female = STDRateFit;
NDatasets_female = nSets;
%% Plot
yyaxis left
hold on
errorbar(ap,RateFit_male,STDRateFit_male./sqrt(NDatasets_male),'-')
errorbar(ap,RateFit_female,STDRateFit_female./sqrt(NDatasets_female),'-')
xlim([0.15 0.6])
title('Initial rate of RNAP loading over AP')
xlabel('AP')
ylabel('Inital rate of RNAP loading')
%legend('male','female')
standardizeFigure(gca,legend,[])

% % Plot the Runt concentration profile over AP
yyaxis right
hold on
plothandle(3) = errorbar(0:0.025:1,RuntProtein_female.MeanVectorAP(13,:),RuntProtein_female.SEVectorAP(13,:),'-')
plothandle(4) = errorbar(0:0.025:1,RuntProtein_male.MeanVectorAP(13,:),RuntProtein_male.SEVectorAP(13,:),'-')

ylabel('Runt nuclear fluorescence (AU)')
legend('Loading Rate_{Male}','Loading Rate_{Female}','Runt_{female}','Runt_{male}')
standardizeFigure_YJK(gca,legend,[])
hold off

%% Step2. Compare the Fraction ON after the nuclear tracking curation with Tr2D.
%% Load datasets
% Only load the datasets that are after the Tr2D tracking / compiling
% Note that these datasets are valid in nc13, but not all of them are
% valide before/after then, since I mostly took only nc13.
r3_male = LoadMS2Sets('r3-male');
r3_female = LoadMS2Sets('r3-female');

%% Calculate the Fraction ON
% Here, I'll make an assumption that if the total number of nuclei in one
% AP bin is too small, then the Fraction can be biased. Thus, I will set
% this threshold as N_nuclei_thresh = 4;
N_nuclei_thresh = 3;

% Go through each dataset, then calculate the Fraction ON in nc13
for i=1:length(r3_male)
    % filter out the AP bins that has less than N_nuclei_thresh nuclei
    clear Nfilter
    clear N_total_nuclei
    Nfilter = r3_male(i).TotalEllipsesAP(:,2) > N_nuclei_thresh;
    N_total_nuclei = r3_male(i).TotalEllipsesAP(:,2).*Nfilter;
    FractionON_individual_male(:,i) = r3_male(i).EllipsesOnAP{1,1}(:,2)./N_total_nuclei ;
    FractionON_individual_male(FractionON_individual_male==inf) = nan;
end

FractionON_male = nanmean(FractionON_individual_male,2);
FractionON_SEM_male = nanstd(FractionON_individual_male,0,2)./sqrt(length(r3_male));

for i=1:length(r3_female)
    % filter out the AP bins that has less than N_nuclei_thresh nuclei
    clear Nfilter
    Nfilter = r3_female(i).TotalEllipsesAP(:,2) > N_nuclei_thresh;
    N_total_nuclei = r3_female(i).TotalEllipsesAP(:,2).*Nfilter;
    FractionON_individual_female(:,i) = r3_female(i).EllipsesOnAP{1,1}(:,2) ./N_total_nuclei ;
    FractionON_individual_female(FractionON_individual_female==inf) = nan;
end

FractionON_female = nanmean(FractionON_individual_female,2);
FractionON_SEM_female = nanstd(FractionON_individual_female,0,2)./sqrt(length(r3_female));

hold on
errorbar(0:0.025:1,FractionON_male,FractionON_SEM_male)
errorbar(0:0.025:1,FractionON_female,FractionON_SEM_female)
% pause
% for i=1:length(r3_male)
%     plot(0:0.025:1, FractionON_individual_male(:,i),'b')
%     pause
% end
% for i=1:length(r3_female)
%     plot(0:0.025:1, FractionON_individual_female(:,i),'r')
%     pause
% end

title('Fraction ON in NC13 for different sex (r3)')
xlabel('AP axis (Embryo Length)')
ylabel('Fraction ON')
legend('male','female')
% save the plots

%% save the processed results
save('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\FractionON_sex_afterTr2D.mat',...
        'FractionON_female','FractionON_male',...
        'FractionON_individual_female','FractionON_individual_male',...
        'FractionON_SEM_female','FractionON_SEM_male')

