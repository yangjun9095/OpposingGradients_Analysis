% Analyze Single trace fits
function Analyze_single_fits

% This script is for analysis of initial fits for single MS2 traces.
% Here, I will use this to compare the expression of r3 for males vs
% females
 
% Use LoadMS2Sets in the future, for now, I will just use Particles.mat
Folder = 'E:\YangJoon\LivemRNA\Data\DynamicsResults';
Male = LoadMS2Sets('r3-male');
Female = LoadMS2Sets('r3-female');

for i=1:length(Male)
    Particles{i} = load([Folder,filesep,Male(i).Prefix,filesep,'Particles.mat']);
    MaleParticles{i} = Particles{i}.Particles;
end

for i=1:length(Female)
    Particles{i} = load([Folder,filesep,Female(i).Prefix,filesep,'Particles.mat']);
    FemaleParticles{i} = Particles{i}.Particles;
end
%% Put all particles into one cell structure, with their AP position, fitted slope and T_ON
k=0;
for i=1:length(Male)
    for j=1:length(MaleParticles{i})
        if ~isempty(MaleParticles{i}(j).fitApproved)
            k=k+1;
            maleParticles(k).APpos = nanmean(MaleParticles{i}(j).APpos);
            maleParticles(k).fittedSlope = MaleParticles{i}(j).fittedSlope;
            maleParticles(k).tON = MaleParticles{i}(j).fittedTON;
        end
    end
end

k=0;
for i=1:length(Female)
    for j=1:length(FemaleParticles{i})
        if ~isempty(FemaleParticles{i}(j).fittedSlope)
            k=k+1;
            femaleParticles(k).APpos = nanmean(FemaleParticles{i}(j).APpos);
            femaleParticles(k).fittedSlope = FemaleParticles{i}(j).fittedSlope;
            femaleParticles(k).tON = FemaleParticles{i}(j).fittedTON;
        end
    end
end

%% Plot single fits over AP
hold on
for i=1:length(maleParticles)
    plot(maleParticles(i).APpos,maleParticles(i).fittedSlope,'ob')
end
    
for i=1:length(femaleParticles)
    plot(femaleParticles(i).APpos,femaleParticles(i).fittedSlope,'or')
end

title('Initial RNAP loading rates (from single traces)')
xlabel('AP axis')
ylabel('Initial RNAP loading rate (AU/min)')


%% Taking average of these single fits
% Here, the binning might change the results, I will start with 1%
X = 0:0.01:1;
%X = 0:0.025:1;


APfilter_male = zeros(length(maleParticles),length(X));
APfilter_female = zeros(length(femaleParticles),length(X));

% Make matrices of fittedSlope
for i=1:length(maleParticles)
    singlefits_male(i) = maleParticles(i).fittedSlope;
end

for i=1:length(femaleParticles)
    singlefits_female(i) = femaleParticles(i).fittedSlope;
end

% Make AP filter
for i=2:length(X)
    for j=1:length(maleParticles)
        if maleParticles(j).APpos > X(i-1) && maleParticles(j).APpos < X(i)
            APfilter_male(j,i-1) = 1;
        else
            APfilter_male(j,i-1) = 0;
        end
    end
    
    for k=1:length(femaleParticles)
        if femaleParticles(k).APpos > X(i-1) && femaleParticles(k).APpos < X(i)
            APfilter_female(k,i-1) = 1;
        else
            APfilter_female(k,i-1) = 0;
        end
    end
end

%% Calculate the mean and std
clear N_male
clear N_female
clear sum_male
clear sum_female

mean_singlefits_male = nan(1,length(X));
mean_singlefits_female = nan(1,length(X));

sem_singlefits_male = nan(1,length(X));
sem_singlefits_female = nan(1,length(X));

N_male = nan(1,length(X));
N_female = nan(1,length(X));

sum_male = nan(1,length(X));
sum_female = nan(1,length(X));


for i=1:length(X)-1
    filterd_singlefits_male = transpose(singlefits_male).*APfilter_male(:,i);
    filterd_singlefits_male(filterd_singlefits_male==0) = nan;
    sum_male(i) = nansum(filterd_singlefits_male);
    N_male(i) = sum(APfilter_male(:,i));
    if N_male(i) == 0
        N_male(i) = 1;
    end
    mean_singlefits_male(i) = sum_male(i) / N_male(i);
    sem_singlefits_male(i) = nanstd(filterd_singlefits_male)./sqrt(N_male(i));
end
mean_singlefits_male(mean_singlefits_male==0) = nan;
sem_singlefits_male(std_singlefits_male==0) = nan;


for i=1:length(X)-1
    filterd_singlefits_female = transpose(singlefits_female).*APfilter_female(:,i);
    filterd_singlefits_female(filterd_singlefits_female==0) = nan;
    sum_female(i) = nansum(filterd_singlefits_female);
    N_female(i) = sum(APfilter_female(:,i));
    if N_female(i) == 0
        N_female(i) = 1;
    end
    mean_singlefits_female(i) = sum_female(i) / N_female(i);
    sem_singlefits_female(i) = nanstd(filterd_singlefits_female)./sqrt(N_female(i));
end
mean_singlefits_female(mean_singlefits_female==0) = nan;
sem_singlefits_female(std_singlefits_female==0) = nan;

%% 
hold on
errorbar(X, mean_singlefits_male, sem_singlefits_male)
errorbar(X, mean_singlefits_female, sem_singlefits_female)
xlim([0.15 0.6])
title('Averaged initial loading rates (2.5% binning)')
xlabel('AP axis')
ylabel('Initial rate of RNAP loading (AU/min)')
legend('male','female')
standardizeFigure_YJK(gca,legend,[])
end