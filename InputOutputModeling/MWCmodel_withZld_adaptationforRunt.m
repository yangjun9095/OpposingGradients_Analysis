%function[MeanVectorAP,DelT]=MWCmodel_withRunt(DelE,Kd,Kdz,w,BcdData,RuntData,trigger,varargin)
% For this script we are going to simulate the fluorescent traces that we
% expect to see based on the experimentally measured Bicoid and Runt concentrations.
%  We are assuming here that all of the nuclei are "on." Later we could add
%  in some probabilities if we want. We will also model the effect of
%  dynamic Zld concentrations by using a single trace of Zld and assuming
%  each AP position has that profile (this is okay since Zld seems to be
%  homogeneous throughout the embryo).

% The inputs to the model are:
% DelE - energy of inaccessible state - energy of accessible state.
% Kdz - Dissociation constant of Zld binding to the gene.
% Kd - Dissociation constant of Bcd binding to the gene.
% w - higher-order cooperativity constant.
% trigger - trigger time (min) - we assume transcription to be forbidden
% before this time due to e.g. mitotic repression
% model - specified type of binding model:
%   'one' = one or more
%   'all' = all or nothing
%   'add' = additive
%         Additional options:
%   'static' - use static inputs as mean over time

% Here are some assumptions we will make:
    %DeltaT is the time for a PolII molecule to tranverse the length of
    %gene - this is the time during which it will be localized to our site
    %of transcription. 
     
    %The minimum number of PolII required to be able to detect the
    %fluorescence above the background is 3 molecules. 
    
%So given the experimental traces of [Bcd] we are going to generate
%predicted fluorescence traces and then fit with FitMeanAP to determine
%an experimental t_on.
Kd=13;
%Hill coefficient-starting with 6
n=6;
%Disocciation constat - starting with 2600 
% Kd=12;
%Length of gene being transcribed (P2-MS2-LacZ) is about 5.3kb
Length=5.296;
%Set an elongation rate - Garcia et al. Current Biology 2013 measured
%1.54+/-0.14 kb per minute. There is some recent infor that suggests that
%this number is only 1.2 of what it should be. 
Elongrate=1.54;
%The length of the MS2 loops
MS2=1.275;
%And the length of LacZ
LacZ=Length-MS2; 
%Set max and min initiation rates in PolII molecules/minute starting with
%17 and 5 are roughly what we see in our experiemntal data for anterior 
%and posterior regions;
RateMax=17;
RateMin=3; 
% Kd=8
%Let's have DelE defined as eAccessible-eInaccessible -- a negative Del E
%would be the wildtype (with Zelda case)
%DelE=-2;
%Define the cooperativity of the bicoid molecules
%w=1;


%Kd=10+(4*d);

%% 
%Load the bicoid concentration data

% Load in nuclear cycle 13 bicoid data
load('BcdLevelsNC13.mat');
BcdNC13Ant = BcdData.BcdNC13Ant;
BcdNC13Mid = BcdData.BcdNC13Mid;
BcdNC13Pos = BcdData.BcdNC13Pos;

TimeNC13Ant = BcdData.TimeNC13Ant;
TimeNC13Mid = BcdData.TimeNC13Mid;
TimeNC13Pos = BcdData.TimeNC13Pos;

% This data is broken into three compartments:
% Anterior:     columns 9-25
% Middle:       columns 13-29
% Posterior:    columns 22-37

%How many times steps do we want?
TimeSteps=200;

%created vector from 0 to 1 broken into 41 bins
Pos=linspace(0,1,41);  

%first column of usable bicoid data
First=9;

%last column of usable bicoid data
Last=37;

%range of usable bicoid data
NumCol=Last-First+1;


% For the anterior compartment, we'll use columns 9 to 13
% For the middle area, we'll use columns 14 to 28
% For the posterior area, we'll use columns 29 to 37

% Create 

DelT = nan(1,41); %Create timestep array per AP bin

for i=First:13
    TimeEnd=max(TimeNC13Ant);
    % Create time matrix
    Time(:,i)=linspace(0,TimeEnd,TimeSteps)';
    % Interpolate between timepoints to approximate temporally continuous bicoid data
    pchbcd(:,i)=pchip(TimeNC13Ant(1,:),BcdNC13Ant(:,i),Time(:,i));
    % Create matrix that records the size of timesteps from the Time matrix
    DelT(1,i)=Time(2,i)-Time(1,i);
end

% Operate same as loop above, but for Middle range of Bicoid data
for i=14:28
    TimeEnd=max(TimeNC13Mid);
    Time(:,i)=linspace(0,TimeEnd,TimeSteps)';
    pchbcd(:,i)=pchip(TimeNC13Mid(1,:),BcdNC13Mid(:,i),Time(:,i));
    DelT(1,i)=Time(2,i)-Time(1,i);
end

% Operates same as loop above, but for Posterior range of Bicoid data
for i=29:Last
    TimeEnd=max(TimeNC13Pos);
    Time(:,i)=linspace(0,TimeEnd,TimeSteps)';
    pchbcd(:,i)=pchip(TimeNC13Pos(1,:),BcdNC13Pos(:,i),Time(:,i));
    DelT(1,i)=Time(2,i)-Time(1,i);
end

for i=First:Last;
    for j=51:55;
        Time(j,i)=Time(j-1,i)+DelT(i);
    end
end

%populate columns 1-8 with NaNs so they don't affect our Rate matrix
% in the all columns for rows 1:First-1, fill the verticle columns with NaNs
pchbcd(:,1:First-1)=nan(size(pchbcd,1),First-1);
pchbcd(:,Last+1:41)=nan(size(pchbcd,1),41-Last);

pchbcd=0.17*pchbcd; %Rescale to get the right arbitrary units.

% If static option selected, use mean Bcd over time as input
% if any(strcmp(varargin,'static'))
%     meanbcd = nanmean(pchbcd,1);
%     for i = 1:size(pchbcd,1)
%         pchbcd(i,:) = meanbcd;
%     end
% end

%Assume 10 Zelda binding sites.
%nz = 10;

% Load the Runt data. Right now the data is the same for all AP bins, using
%an interpolation of traces from the Runt female (averaged) datasets. 
% For the Runt null simulations, we pass in a "fake" dataset with all zero values, to make the
%code usable for both data types.
RuntData = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-1min-200Hz-Female-Averaged.mat')
RuntNC13 = RuntData.MeanVectorAP(RuntData.nc13:RuntData.nc14,:);

TimeNC13Runt = RuntData.ElapsedTime(RuntData.nc13:RuntData.nc14) - RuntData.ElapsedTime(RuntData.nc13);


% We will be using column 15 as the proxy Zld profile for each AP position.
pchRunt = nan(size(Time,1),41);
% Interpolate the Zelda profile to match the Bcd timing.
for i=First:Last
    pchRunt(:,i)=pchip(TimeNC13Runt(1,:),RuntNC13(:,15),Time(:,i));
end

% Runt background fluo (free eGFP) should be subtracted properly.
% For now, I'll just subtract the minimal value, assuming that the
% effective Runt concentration in mitosis is zero...

% % If static option selected, use mean Runt over time as input
% if any(strcmp(varargin,'static'))
%     meanRunt = nanmean(pchRunt,1);
%     for i = 1:size(pchRunt,1)
%         pchRunt(i,:) = meanRunt;
%     end
% end
% plot(Time(10:end,:),pchbcd(10:end,:))
% 
% plot(Pos(1:37),max(pchbcd(10:end,:)))
%%
%Now we need to make a rate matrix corresponding to the dynamic bicoid and
%Runt concentration data (using the pchbcd pchRunt matrix). 

%Starting with the basic equation for the rate.  Calling pRate in case we
%add other Rate equations later.

%Enumerate the different states and their weights.
% 
% pAcc=((pchbcd/Kd).^n)./(exp(DelE)+(pchbcd/Kd).^n);
% 
% Rate=((pchbcd/Kd).^n)./(1+(pchbcd/Kd).^n)*(RateMax-RateMin)+RateMin;
% 
%Weights of Bicoid bound states.
P_0=1;
P_1=nchoosek(6,1)*(pchbcd/Kd);
P_2=nchoosek(6,2)*(w*(pchbcd/Kd).^2);
P_3=nchoosek(6,3)*((w^2)*(pchbcd/Kd).^3);
P_4=nchoosek(6,4)*((w^3)*(pchbcd/Kd).^4);
P_5=nchoosek(6,5)*((w^4)*(pchbcd/Kd).^5);
P_6=nchoosek(6,6)*((w^5)*(pchbcd/Kd).^6);

%Rates of transcribing state, depending on model

if any(strcmp(varargin,'one'))
    r = [RateMin, RateMax, RateMax, RateMax, RateMax, RateMax, RateMax];
elseif any(strcmp(varargin,'all'))
    r = [RateMin, RateMin, RateMin, RateMin, RateMin, RateMin, RateMax];
elseif any(strcmp(varargin,'add'))
    r = zeros(1,7);
    for x = 0:6
        r(x+1) = RateMin + (RateMax - RateMin)*(x/6);
    end
end
    
R_0 = r(1);
R_1 = r(2);
R_2 = r(3);
R_3 = r(4);
R_4 = r(5);
R_5 = r(6);
R_6 = r(7);
%Calculate the partition function.

PartFun=P_Inacc+((1+pchzld/Kdz).^nz).*(P_0+P_1+P_2+P_3+P_4+P_5+P_6);

%Get the total rate matrix.
pRate = (P_0*R_0 + ((1+pchzld/Kdz).^nz) .* (P_1*R_1+P_2*R_2+P_3*R_3+P_4*R_4+P_5*R_5+P_6*R_6))./PartFun;


%Now plot the rates vs time if you want
% % 
% figure('Name','Initiation Rates');
% for i=First:Last;
% plot(Time,pRate(:,i)); 
% hold on;
% end

%Calculate the number of PolII on the DNA at any given time

%For each time point, lets calculate the number of PolIIs on the DNA.  
%% 


%Clearing stuff
clear PolOcc;
clear PolPos;
clear InstRate;
clear TotOcc;
clear Polocctemp;



    %PolOcc is a matrix of Pol II occupation.  This matrix will have an entry for 
%each position along the embryo and each time point.  This entry represents
%the PolII loaded at that time point. This will stay constant after the Pol 
%II is loaded until the time at which the Pol II falls off the gene.  
%At this time, and after, the entry will be zero.
%Establish the firt time point with a nan vector
NoCol=NumCol;
PolOcc=nan(1,Last);

%PolPos is a companion matrix to PolOcc that record the position
%of the pol II loaded during each time point. The position of the Pol II
%starts at zero for the time point at which it's loaded that updates at
%each time point. 
PolPos=zeros(1,Last);


%Start with a loop over time steps.
for i=2:TimeSteps;
%     For each time step, loop over each postion along the embryo(for which
%     we have Bcd data)
    for k=1:Last;
      
         %Find the Rate for the ith time point and kth position in the embryo
           %Programing in a shutoff of new Pol loading at 18-3=15 minutes
           %(t(1,76)). This is the approx length of nc13. If simulating other cycles change this
           %shutoff time
    if Time(i,k) <= 15
        InstRate=pRate(i,k);  
    else
        InstRate=NaN;
    end
    
    %Trigger initiation in beginning (e.g. due to mitotic repression)
    if Time(i,k)<=trigger
        InstRate=NaN;
    end
    
 
      
   
    %The time difference between each time point is always DelT minutes.
    %Update the pol II occupation matrix
    %First determine the number of Pol II that load during this time point
    PolOcc(i,k)=InstRate*DelT(1,k);
    %Set the newly loaded Pol II's position as zero. 
    PolPos(i,k)=0;
    
        
    %End the position loop
    end
    
    %Now that we have assigned occupations and postions for every position
    %in the ith time point, update all of the postions and occupations for
    %every timepoint up to i-1.
    %Establish a temporary occupation matrix to save the info in the 
    %occupation matrixt before we update it.
    Polocctemp=PolOcc;
    %Loop over all time up to i-1 (we want to update occ and pos of the
    %entries established at earlier time points. 
   for j=1:i-1;
   
        %Dont forget to loop over every position along the AP axis for the
        %updates
        for l=1:Last;
        % We move foward polymerases by (elongation rate x deltaT)
        if max(PolOcc(:,l))>0
            % move forward already loaded polymerases by their elongation
            % rate * time-step
            PolPos(j,l)=PolPos(j,l)+DelT(1,l)*Elongrate;

        end
            %Update the positions but moving the previously loaded Pol II
            %along the reporter gene
            %If there's nothign there, position will remain 0.
%             if PolOcc(j,l)==0;
%                 PolPos(j,l)=0;
%                 %If it's not a number, set the position to zero.
%             else
%                 if isnan(PolOcc(j,1));
%                     PolPos(j,1)=0;
                    %Now if it is a number and is not zero, update the
                    %position by multiplying the elongation rate by DelT.
        %If the updated position of the Pol II is further than the length
        %of the reporter gene, then remove it from the matrix and set the
        %occupation to zero.
        
        if PolPos(j,l)>=Length;
        PolOcc(j,l)=0;
        %Or else only the position is updated and the occupation stays the
        %same.
        end
        %End the postion loop
        end
   %End the time (up to i-1) loop
   end

         %We need to temporarily calibrate the occupation (for fluoresence calc later)
         %for the PolII that are not past the MS2 DNA
   %Loop over all time
   for j=1:i;
       %Loop over all positions
       for l=1:Last;
           %Determine if the Pol II are within the MS2 sequence.
           if PolPos(j,l)<=MS2;
               %Use the fraction of MS2 traversed to calibrate the signal
               %(determined by the number of Pol II for now. 
                Polocctemp(j,l)=Polocctemp(j,l)*(PolPos(j,l)/MS2);
           else
           %End the if statement    
           end
       %End loop over postion    
       end
   %End loop over time    
   end

   
   %Now compute a matrix which stores the total number of Pol II on the
   %gene at each time point and at each position by using nansum (ignores
   %nans)   
 
    TotOcc(i,:)=nansum(Polocctemp);
   %Now end the time loop.
end

%Plot the total Pol II occupation 

% 
% figure('Name','Pol II Occupancies');
% plot(Time,TotOcc);


%Now create a matrix to simulate the fluorescence at each AP position
%Start with the TotOcc matrix
Fluor=TotOcc;

%Take out anything with less than 4 Pol II (again looping over all time and
%all positions

for i=1:TimeSteps;
    for p=1:Last;
    if Fluor(i,p)<=4;
        Fluor(i,p)=NaN;
    else
        Fluor(i,p)=Fluor(i,p)-1;
    end
    end
end

% %Plot the fluorescence
% figure('Name','Fluorescence Signals');
% plot(Time,Fluor);

%Fill out the Fluor matrix with NaN vectors so there are 41 position 
%columns again.

%Create a nan vector to add to the r.h. side of the Fluor matrix.
ZeroVec=zeros(TimeSteps,1);

 %Now loop over the positions in which there is no data.
for i=1:First-1;
    PosFluor(:,i)=ZeroVec;
end

%Set the poisitions with data to equal the Fluor matrix
PosFluor(:,First:Last)=Fluor(:,First:Last);

%Again adding the nan vector to postions in which there is no data.
for i=Last+1:41;
     PosFluor(:,i)=ZeroVec;
end

%plot(Pos,max(PosFluor));

%Now let normalize the PosFluor matrix so the fluorescence intensity is
%something similar to our fluorescence data which maxes out around 1400ish
%I've commented most of this out because we only need the scaling factor,
%but in the future, with new Bcd datasets we'll want to use these
%calculations. 

%First set this max 
% MaxFluor=1400;
% 
% % %Check the max value of the Fluor matrix
% % MaxSim=max(max(Fluor));
% % 
% % %Divide to find the scaling factor
% % ScaFac=MaxFluor/MaxSim;
% 
% %The scaling factor below was determined by some tests independent of the
% %calculations above.  The max value of the fluor matrix in this case is
% %55.8007
ScaFac=25;
% 
% %Now scale the entire matrix
SimuFluor(:,:)=ScaFac*PosFluor;
% ScaTotOcc(:,:)=ScaFac*TotOcc;
% 
% %Clearing just in case.
% clear TotScaPolOcc;
% %Filling in positions around the scaled occupation matrix with the nan
% %vector.
% for i=1:First-1;
%     TotScaPolOcc(:,i)=NanVec;
% end
% 
% TotScaPolOcc(:,First:Last)=ScaTotOcc;
% 
% for i=Last+1:41;
%      TotScaPolOcc(:,i)=NanVec;
% end
% 
% figure('Name','Scaled Fluorescence Signals');
% plot(Time,SimuFluor(:,9));



%Make an NParticlesAP matrix by
%%%We don't need these NParticles.  We actually need at least 2 or 3
%%%particles for each AP position we're fitting to have them show up in the
%%%fitting script
Threevector=ones(TimeSteps,1)*3;
NanVec=nan(TimeSteps,1);

for i=1:First-1
    NParticlesAP(:,1)=NanVec;
end

for i=First:Last
    NParticlesAP(:,i)=Threevector;
end

for i=Last+1:41
    NParticlesAP(:,i)=NanVec;
end


% more stuff to make curve fitting work
Factor=3;
SDVectorAP=SimuFluor/3;
MeanVectorAP=SimuFluor;
nc12=0;
nc13=1;
nc14=TimeSteps;
APbinID=linspace(0,1,41);
ElapsedTime=Time(:,10)';

%Old things that saved information
%{
%save
%save('kinetic-simulation-first-test','APbinID','ElapsedTime','SDVectorAP','MeanVectorAP','NParticlesAP','nc12','nc13','nc14');

save(['ThermoSimpleAct_Rates_-2_w1_KDB=' num2str(Kd)], 'APbinID','ElapsedTime','SDVectorAP','MeanVectorAP','NParticlesAP','nc12','nc13','nc14');

%create an APDivision file
Zvec=zeros(14,41);

for j=First:Last
    APDivision(13,j)=1;
    APDivision(14,j)=TimeSteps;
end

%Add zeros to columns 37-41 to fix error with fit software
APDivision(:,Last+1:41)=zeros(14,41-Last);

%NParticlesAP(45,:)=zeros(1,41);

save('APDivision','APDivision')

%}
% end








