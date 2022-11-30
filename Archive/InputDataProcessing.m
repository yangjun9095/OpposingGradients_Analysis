% Data processing for Bcd and Runt inputs

% Bcd and Runt profiles are matrices, which have the dimension of Time x AP
% bins. This means that we should have the same frame rate, to make it
% easier in prediction. (or interpolation, if it's not that bursty)
% Bcd : [Bcd]/Kd, Bicoid concentration divided by Kd of Bicoid
% Runt : [Runt]/K_R, Runt concentration divided by K_R
BcdData = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\eGFP-Bcd-From-Liz-Jonathan\BcdGFPMid.mat')
Bcd = BcdData.DataBcd; % This has all fields like nc, MeanVectorAP, etc.
Kd = 1; % This value sholud be changed later.
BcdFluo = Bcd.MeanVectorAP / Kd;

Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-04-29-Runt-JB3-MCP-mCherry-vasa-eGFP1-edited-2.5%APbins\CompiledNuclei.mat')
KR = 1; % This value should be changed later.
RuntFluo = Runt.MeanVectorAP / KR;

% NR : number of repressor binding sites
% nR : number of repressor bound
% r_basal : basal rate of RNAP loading, in the empty enhancer.
% r_Max : maximum rate of RNAP loading, in case of 6 Bcd binding
% r_R : decrease in the rate of RNAP loading upon Runt binding.

% Define the size of the Bcd and Runt, assuming that they have the same
% size for now.
[tLength_Bcd,APbin_Bcd] = size(BcdFluo);
[tLength_Runt,APbin_Runt] = size(RuntFluo);

%% Trim the Bcd and Runt data for the same size
% First, we ultimately need a nice interpolation, in case Bcd and Runt data
% were taken with different time resolution, which is actually the case.
% Second, we need to trim down the late nc 14 for the shorter one.

% Third, we will start from only thinking about nc13 and nc14 for now.
% We can include nc12 or earlier later.

% Bcd
BcdTime = Bcd.ElapsedTime(Bcd.nc13:end)-Bcd.ElapsedTime(Bcd.nc13);
BcdTime = BcdTime';
BcdInterval = mode(diff(BcdTime))*60; % [sec]

BcdFluo = BcdFluo(Bcd.nc13:end,:);
BcdFluo(isnan(BcdFluo)) = 0;
% Runt
RuntTime = Runt.ElapsedTime(Runt.nc13:end)-Runt.ElapsedTime(Runt.nc13);
RuntTime = RuntTime';
RuntInterval = mode(diff(RuntTime))*60; % [sec]

RuntFluo = RuntFluo(Runt.nc13:end,:);
RuntFluo(isnan(RuntFluo)) = 0;
%%
% For now, we have 30 seconds for Bcd, and 41(40) seconds for Runt.
% We will interpolate these datasets into 10 seconds interval.

% Time interpolation (Define new time vectors)
NewBcdTime = (0:0.1667:BcdTime(end));
NewRuntTime = (0:0.1667:RuntTime(end));

NewBcdFluo = zeros(length(NewBcdTime),length(Bcd.APbinID));
NewRuntFluo = zeros(length(NewRuntTime),length(Runt.APbinID));

for i=1:length(Bcd.APbinID)
    NewBcdFluo(:,i) = pchip(BcdTime',BcdFluo(:,i),NewBcdTime);
end

for i=1:length(Runt.APbinID)
    NewRuntFluo(:,i) = pchip(RuntTime,RuntFluo(:,i),NewRuntTime);
end

%% Check the interpolation (Bcd)
hold on
plot(BcdTime,BcdFluo,'-or')
plot(NewBcdTime,NewBcdFluo,'-ob')

%% Note that the Runt Background is not subtracted here.
% For now, I will just subtract a constant (that seems to be background),
% which is 95 for all time points.
% This background subtraction should be revisited later, in a more serious
% way. For example, we can subtract the nuclear fluorescence in the very
% anterior (or posterior) at each time point. Or, we can also use the
% formula in the LlamaTag paper using the cytoplasmic fluorescence to infer
% the free fluorophore in the nucleus.
RuntFluo = RuntFluo - 95;
NewRuntFluo = NewRuntFluo - 95;

%% Check the interpolation (Runt)
hold on
plot(RuntTime,RuntFluo,'-or')
plot(NewRuntTime,NewRuntFluo,'-ob')

%% Normalization of two inputs
% Since we don't know how to scale those two inputs, I will just normalize
% them with their maximal intensity for now. We can also do this above with
% Kd and KR
BcdMax = max(max(NewBcdFluo));
RuntMax = max(max(NewRuntFluo));

NewBcdFluo = NewBcdFluo/BcdMax;
NewRuntFluo = NewRuntFluo/RuntMax;
%% Now, plug in those two inputs (Bcd and Runt) into the equation
% This is a subtle point, but we need to set the dimension of both inputs
% same for easier calculation.
% Find out which one has shorter time length.
% For now, I will truncate the Bcd data to fit to the Runt data.
NewBcdTime = NewBcdTime(1:length(NewRuntTime));
NewBcdFluo = NewBcdFluo(1:length(NewRuntTime),:);


%% Saving the processed data

% Define the path
savePath = uigetdir; % Click the folder that you want to save the .mat file.

Bcd.Fluo = NewBcdFluo;
Bcd.ElapsedTime = NewBcdTime;

Runt.Fluo = NewRuntFluo;
Runt.ElapsedTime = NewRuntTime;

save([savePath,filesep,'InputData-MidBcd','.mat'],...
            'Bcd','Runt')


