function CheckSingleTraces
% Author : Yang Joon Kim (yjkim90@berkeley.edu)
% This code is plotting individual MS2 spot traces along with the Averaged
% MS2 trace, to see the trend.
% So, we will start from importing one embryo, then expand this code for
% averaged over multiple embryos later.

%% Load the data set - this needs to be edited for conveniences later.
clear all
Data = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-12-11-hbP2-r3-MS2V5-20sec-8\CompiledParticles.mat');
Particles = cell2mat(Data.CompiledParticles);
%Protein = load('CompiledNuclei.mat');
%Ellipses = load('Ellipses.mat');
%Particles = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-04-30-hbP2-r1-MS2V5-lacZ\Particles.mat');
%Spots = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-04-30-hbP2-r1-MS2V5-lacZ\Spots.mat');

%Load the schnitzcells
%Schnitzcells = load('*lin.mat');
%load('2017-03-06-Hb-P2P-MS2V5-NB-MCP-NLS-mCherry_lin.mat')

%% Plot all MS2 spot (particle) traces
% Particles - individual MS2 spots
% CompiledNuclei - individual nuclei
%Particles = mRNA.CompiledParticles;
%Nuclei = Protein.CompiledNuclei;

for i=1:length(Particles)
    Frame = Particles(i).Frame;
    SpotFluo = Particles(i).Fluo;
    %SpotFluoError = Particles(i).FluoError;
%     NucFrame = Nuclei(find(Nuclei.schnitz ==i)).Frames;
%     NucFluo = Nuclei(Nuclei.schnitz == i).MaxFluo;
    
    hold on
    plot(Frame,SpotFluo)
    %errorbar(Frame,SpotFluo,SpotFluoError) 
    %plot(NucFrame,NucFluo)
    pause
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
%% Additional part
% compare single traces between male and females
% First, load the datasets

end