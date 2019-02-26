function hbP2_r0123_Variability (varargin)
% This code is for the analysis of variability (noise) between nuclei in
% the same AP bin, for hbP2-r0,1,2,3 constructs.

% The idea is to see the variability of mRNA (Integrated fluorescence)
% 1) For different AP bins (different Bcd, Runt concentration)
% 2) For different numbers of Runt binding sites
% 3) For different cycles (nc12, nc13, and nc14)
clear all
% Load the dataset
%Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-03-23-hbP2-r3-MS2V5-lacZ\CompiledParticles.mat')
%Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-04-25-hbP2-r0-MS2V5-lacZ\CompiledParticles.mat')
r0Data = LoadMS2Sets('r0');
r3Data = LoadMS2Sets('r3');
%% Plot Accumulated Fluo (mRNA) of single particles, over AP (at certain NC)
NC=13;
%Particles = Data.CompiledParticles;
%hold on
% for i=1:length(Particles)
%     if Particles(i).nc==NC && ~isempty(Particles(i).TotalmRNA)
%         errorbar(Particles(i).MeanAP,Particles(i).TotalmRNA,Particles(i).TotalmRNAError)
%     end
% end
k=1;
hold on
for i=1:length(r0Data)
    Particles=r0Data(i).CompiledParticles;
    for j=1:length(Particles)
        if Particles(j).nc==NC && ~isempty(Particles(j).TotalmRNA) &&...
                Particles(j).TotalmRNA<=14000 && Particles(j).TotalmRNA>=1000
            errorbar(Particles(j).MeanAP,Particles(j).TotalmRNA,Particles(j).TotalmRNAError)
            % Define a new Particle field for the ones that passed through
            % the thresholding.
            NewParticles(k).MeanAP = Particles(j).MeanAP;
            NewParticles(k).TotalmRNA = Particles(j).TotalmRNA;
            NewParticles(k).TotalmRNAError = Particles(j).TotalmRNAError;
            k = k+1;
        end
    end
end
    
title(['Total mRNA (single particle) over AP during ','NC=',num2str(NC)])
xlabel('AP')
ylabel('Total mRNA(AU)')
%% variability of RNAP Loading rate

%% Variability of mRNA (accumulated) 
% Note that this result could be compared with Shawn Little's paper(Cell, 2013)
% To calculate the variability of mRNA (produced in single nucleus), I need
% to track single traces pretty accurately. mRNA1(t), mRNA2(t),...,etc.

% For now, I will pick some AP bins, with well-tracked MS2 spot traces
% Question : How should I think about the lineages? Do I need His-iRFP?
% Even then, do I need to assume that the mRNA is separated equally during
% mitosis, or exported from the nucleus during mitosis?

% Let's start from considering nc13 only.
%% (Test) Pick one AP bin, and check if the particle tracking is done well
% Note. Particles.mat and CompiledParticles.mat are different.
% I will start with CompiledParticles.

% Definitely, I need to be careful for CheckParticleTracking.
% I should go back, and manually curate the missing spots.
%Particles = Data.CompiledParticles;

% NewParticles
for AP = 1:41
    Index =[]; % Get the indices of particles that satisfies conditions below.
    OriginalIndex = [];
    for i=1:length(NewParticles)
        %clf
        % First, let's start with nc 13
        %if Particles(i).nc==13
            % Pick some AP bins, let's start with 0.2~0.3 AP bin
            if (Particles(i).MeanAP <(AP*0.025)) && Particles(i).MeanAP > (AP-1)*0.025
                %if sum(Particles(i).FrameApproved)>5 % At least 5 frames should exist.
                    %plot(Particles(i).Frame,Particles(i).Fluo,'-o')
                    Index = [Index i];
                    %OriginalIndex = [OriginalIndex Particles(i).OriginalParticle;]
                    %pause
                %end
            end
        %end
    end
    ParticleIndex(AP).Index = Index;
    %ParticlesIndex(AP).OriginalIndex = OriginalIndex;
end
%% Check the Particle's position
% xpos =[];
% ypos= [];
% for i=1:length(Index)
%     j=Index(i);
%     xpos(i) = nanmean(Particles(j).xPos);
%     ypos(i) = nanmean(Particles(j).yPos);
%     
%     plot(xpos,ypos,'o','MarkerSize',11)
%     pause
% end

%% TotalmRNA and TotalmRNAError @ end of nc
% Be careful for the particles that survived during the mitosis, since they
% can be tracked for two nuclear cycles.
for AP=1:41
    Index = ParticleIndex(AP).Index;
    SortedParticles = NewParticles(Index);
    TotalmRNA = [];
    TotalmRNAError = [];

    % Use the TotalmRNA calculated by CompileParticles
    for i=1:length(SortedParticles)
        TotalmRNA = [TotalmRNA SortedParticles(i).TotalmRNA];
        TotalmRNAError = [TotalmRNAError SortedParticles(i).TotalmRNAError];
    end
    mRNA(AP).TotalmRNA = TotalmRNA;
    mRNA(AP).TotalmRNAError = TotalmRNAError;
    mRNA(AP).STD_mRNA = std(TotalmRNA);
    mRNA(AP).Var_mRNA = var(TotalmRNA);
    mRNA(AP).MeanmRNA = nanmean(TotalmRNA);
    mRNA(AP).FanoFactor = var(TotalmRNA)./nanmean(TotalmRNA);
    % Error of the STD should be estimated later.
end

%% Plot the Variability (Standard Deviation)
for AP=1:41
    STD_mRNA(AP) = mRNA(AP).STD_mRNA;
    FanoFactor(AP) = mRNA(AP).FanoFactor;
end

%% Plot the Fano factor
plot(0:0.025:1,FanoFactor)
title('Fano factor (in nc13) over AP')
xlabel('AP')
ylabel('Fano Factor (mRNA)')

%% For each Particle, calculate the accumulated mRNA at each time point.

% SortedParticles = Particles(Index);
% 
% % Define a cell to save all variables.
% mRNA ={};
% for i=1:length(SortedParticles)
%     mRNA(i).Frame = SortedParticles(i).Frame;
%     for j=2:length(SortedParticles(i).Frame)
%         mRNA(i).TotalmRNA(j) = trapz(SortedParticles(i).Frame(1:j),SortedParticles(i).Fluo(1:j));
%     end
% end

%% plot to check the mRNA accumulation
% for i=1:length(SortedParticles)
%     plot(mRNA(i).Frame,mRNA(i).TotalmRNA,'-o')
%     pause
% end
%% Match the Accumulated mRNA with time points
% Generate the time frame series

%% Calculate the mean and SEM of TotalmRNA
% MeanTotalmRNA = mean(TotalmRNA);
% STDTotalmRNA = std(TotalmRNA);
% SEMTotalmRNA = STDTotalmRNA / sqrt(length(TotalmRNA));

%% variability of protein (inferred)

%% variability of protein (measured by NB signal)
% For this, I can use the SDVectorAP and NParticles in CompiledNuclei.mat
%% Compare the level of variabilities at different steps of Central Dogma
end