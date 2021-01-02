function FitP2PLinearRise(DataStatusTab)
%Fits single particle MS2 data from P2P datasets to linear rise model.
%Input in the name of the data status tab to fit the entire dataset.

P2P = LoadMS2Sets(DataStatusTab);

for i = 1:length(P2P)
    if ~isempty(P2P(i).Prefix)
        FitSingleParticleLinearRise('Prefix',P2P(i).Prefix);
    else
        Prefix = P2P(i).SetName(10:end-1);
        FitSingleParticleLinearRise('Prefix',Prefix);
    end
end
end