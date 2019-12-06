function Supple_05_LlamaTag_TF_BGsubtraction
%% DESCRIPTION

%% Load the datasets
RuntMale = LoadMS2Sets('Runt-1min-200Hz-Male')
RuntFemale = LoadMS2Sets('Runt-1min-200Hz-Female')


%% Process background subtraction using processLlamaTagData_subtractBGFluo
for i=1:length(RuntMale)
    DatasetName = RuntMale(i).SetName;
    Prefix = DatasetName(11:end-1); % extracting the Prefix = ''
    processLlamaTagData_subtractBGFluo(Prefix);
end

for i=1:length(RuntFemale)
    DatasetName = RuntFemale(i).SetName;
    Prefix = DatasetName(11:end-1); % extracting the Prefix = ''
    processLlamaTagData_subtractBGFluo(Prefix);
end

%% Take background subtracted nuclear fluo to get mean,std over APbins.
% Make a function for this task, so that I can make things modular.
% CompileNuclearLlamaTaggedProtein.m
%CompileNuclearLlamaTaggedProtein(Prefix)
for i=1:length(RuntMale)
    DatasetName = RuntMale(i).SetName;
    Prefix = DatasetName(11:end-1); % extracting the Prefix = ''
    CompileNuclearLlamaTaggedProtein(Prefix);
end

for i=1:length(RuntFemale)
    DatasetName = RuntFemale(i).SetName;
    Prefix = DatasetName(11:end-1); % extracting the Prefix = ''
    CompileNuclearLlamaTaggedProtein(Prefix);
end


end


