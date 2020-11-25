function [results,chain,s2chain]  = fitstuff_mcmc2glob(varargin)

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')

expmnt = "affinities"; %affinities, phases or phaff
md = "simpleweak"; %simpleweak, simpleweakdimer, repression, tfdriven, artifact, fourstate
metric = "fraction"; %fraction, fluo
lsq = false;
noOff = true;
nSimu = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats
minKD = 200;
maxKD = 1E4;
minw = 1E-2; %1E-2
maxw = 1E1; %1E2
minR = 10;
maxR = 1E3;
displayFigures = true;
wb = true;
fixedKD = NaN; %if this value isn't nan, this value will determine the fixed KD parameter and KD won't be fitted
fixedOffset = NaN;
fixedR = NaN;
fixedw = NaN;
enhancerSubset = {};
scoreSubset = [];
positionSubset = [];
useBatches = true; %fit all the data across embryos and not just means 
noAverage = false;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if expmnt == "phaff"
    lsq = false;
    noOff = true;
end

if noAverage
    useBatches = false;
end

enhancers_1dg = {'1Dg11'};
enhancers_aff =  {'1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
enhancers_ph = {'1Dg-5', '1Dg-8D'};
scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
positions = [0, -5, -8]';

if strcmpi(expmnt, 'affinities')
    enhancers = [enhancers_1dg, enhancers_aff];
elseif strcmpi(expmnt, 'phases')
    enhancers = [enhancers_1dg, enhancers_ph];
elseif expmnt=="phaff"
    enhancers = [enhancers_1dg, enhancers_aff, enhancers_ph];
end

if ~isempty(enhancerSubset)
    enhancers = enhancerSubset;
    scores= scoreSubset;
    positions = positionSubset;
end


%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = getXRange(enhancers, expmnt);

nSets = length(enhancers);
xo = {};
yo = {};
xs = {};
ys = {};
dsid = [];
T = [];
Y = [];
T_batch = [];
Y_batch = [];
dsid_batch = [];
for k = 1:nSets
    cond = strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k});
    xo{k} = dorsalResultsDatabase.dorsalFluoBins(cond);
    if metric == "fraction"
        yo{k} = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
        yo_batch{k} = dorsalResultsDatabase.fracFluoEmbryo(cond, :);
    elseif metric == "fluo"
        yo{k} = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond);
        yo_batch{k} = dorsalResultsDatabase.allMaxFluoEmbryo(cond, :);
    end
    
     xo_batch{k} = repmat(xo{k}, [1, size(yo_batch{k}, 2)]);
    [xs_batch{k}, ys_batch{k}]= processVecs(xo_batch{k}, yo_batch{k}, xrange(k, :));
    
    [xs{k}, ys{k}]= processVecs(xo{k}, yo{k}, xrange(k, :));
    dsid = [dsid; k*ones(size(xs{k}))];
    dsid_batch = [dsid_batch; k*ones(size(xo{k}))];

    if isnan(xrange(k, 1))
        xrange(k, 1) = min(xo{k});
    end
    if isnan(xrange(k, 2))
        xrange(k, 2) = max(xo{k});
    end
    
    assert(~isempty(xs{k}));
    
    T = [T; xs{k}];
    Y = [Y; ys{k}];
    T_batch = [T_batch; xo{k}];
    Y_batch = [Y_batch; ys_batch{k}];
end

if ~noAverage
    data.ydata = [T, Y];
    data.dsid = dsid;
    data.X =  [T dsid];
end

data_batch = {};
for k = 1:size(Y_batch, 2)
    data_batch{k}.ydata = [T_batch Y_batch(:, k)];
    data_batch{k}.X =  [T_batch dsid_batch];
    data_batch{k}.dsid = dsid_batch;
end

if useBatches
    disp('Doing batched fits');
    data = data_batch;
end

%%


%rate, kd, hill, y offset
y_max = nanmax(Y(:));
x_max = max(T);

[p0, lb, ub, names] = getInits(expmnt, md, metric ,x_max, y_max, nSets,...
    'minKD', minKD, 'maxKD', maxKD, 'minw', minw, 'maxw', maxw, 'minR', minR, 'maxR', maxR, 'subset', ~isempty(enhancerSubset));


if lsq
    % Refine the first guess for the parameters with fminseacrh and calculate residual variance as an estimate of the model error variance.
    [k0, mse] = globfit2('expmnt',expmnt, 'metric', metric, 'md', md, 'maxKD', maxKD, 'displayFigures', displayFigures);
else
    k0 = p0;
end


if noOff
    %ignore the offset in the initial parameters and bounds
    p0(names=="offset") = [];
    lb(names=="offset") = [];
    ub(names=="offset") = [];
    k0(names=="offset") = [];
    names(names=="offset") = [];
end

% put the initial parameters and bounds in a form that the mcmc function
% accepts
params = cell(1, length(k0));
pri_mu = NaN; %default prior gaussian mean
pri_sig = Inf; %default prior gaussian variance
localflag = 0; %is this local to this dataset or shared amongst batches?


for i = 1:length(k0)
    
    targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.
    
    if ~isnan(fixedKD) && contains(names(i), "KD")
        k0(i) = fixedKD;
        targetflag = 0;
    end
    
    if ~isnan(fixedOffset) && contains(names(i), "off")
        k0(i) = fixedOffset;
        targetflag = 0;
    end
    
    if ~isnan(fixedR) && contains(names(i), "R")
        k0(i) = fixedR;
        targetflag = 0;
    end
    
    if ~isnan(fixedw) && contains(names(i), "w")
        k0(i) = fixedw;
        targetflag = 0;
    end
    
    params{1, i} = {names(i),k0(i), lb(i), ub(i), pri_mu, pri_sig, targetflag, localflag};
    
end


model = struct;

simpleWeakOptions = struct('noOff', noOff, 'fraction', metric=="fraction",...
    'dimer', contains(md, "dimer"), 'expmnt', expmnt);

if noAverage
    load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'combinedCompiledProjects_allEnhancers')
    ccp = combinedCompiledProjects_allEnhancers(strcmpi({combinedCompiledProjects_allEnhancers.dataSet}, '1Dg11_2xDl' )  &...
        [combinedCompiledProjects_allEnhancers.cycle]==12);
    ccp = ccp(cellfun(@any, {ccp.particleFrames}));
    data.ydata = [ [ccp.dorsalFluoFeature]' , cellfun(@max, {ccp.particleFluo3Slice})'];
    simpleWeakOptions.onedsid = true;
end

mdl = @(x, p) simpleweak(x, p, simpleWeakOptions);

%leaving this here in case it'll be useful in the future
%     model.ssfun = @(params, data) sum( (data.ydata(:,2)-mdl(data.X(:,1), params)).^2 );

model.modelfun   = mdl;  %use mcmcrun generated ssfun 

if lsq
    model.sigma2 = mse;
end

options.drscale = 5; % a high value (5) is important for multimodal parameter spaces
options.waitbar = wb; %the waitbar is rate limiting sometimes
options.nsimu = nSimu; %should be between 1E3 and 1E6
options.updatesigma = 1; %honestly don't know what this does

rng(1,'twister'); %set the rng seed so we get the same results every run of this function

if lsq
    [results,chain,s2chain] = mcmcrun(model,data,params,options);
else
    %we're gonna run this three times and use the initial results of one
    %run as conditions for the next. this is an alternative when common least
    %squares gives results too poor to initialize with
    results = [];
    [results,~,~,~]=mcmcrun(model,data,params,options,results);
    [results,~,~,~]=mcmcrun(model,data,params,options,results);
    [results,chain,s2chain,~]=mcmcrun(model,data,params,options,results);
end

if displayFigures
    
    burnInTime = .25; %let's burn the first 25% of the chain just in case
    chain = chain(round(burnInTime*nSimu):nSimu, :);
    if ~isempty(s2chain)
        s2chain = s2chain(round(.25*nSimu):nSimu, :);
    end
    
    chainfig = figure(); clf
    mcmcplot(chain,[],results,'chainpanel')
    
    % Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.
    %geweke is a measure of whether things converged between 0 and 1.
    chainstats(chain,results)
    
    
    
    %ideally, these guys look like ellipses. if certain parameters give weird
    %shapes, it might mean those parameters should be removed from the model if
    %possible
    pairFig = figure; clf
    mcmcplot(chain,[],results,'pairs', .5);
    %
    
    
    figure(4); clf
    til = tiledlayout(1, nSets);
    dsid2 = [];
    
    xx = (0:10:max(T))';
    
    for k = 1:nSets
        dsid2 = [dsid2; k*ones(length(xx), 1)];
    end
    X2 = [repmat(xx, nSets, 1), dsid2];
    
     if results.nbatch ~= 1
        for k = 1:results.nbatch
          data_mcmcpred{k} =   repmat(xx, nSets, 1);
        end
    else
        data_mcmcpred = repmat(xx, nSets, 1);
    end
    
    
    %get the prediction intervals for the parameters and function vals
    out = mcmcpred(results,chain,[],data_mcmcpred, mdl);
    
    % mcmcpredplot(out);
    nn = (size(out.predlims{1}{1},1) + 1) / 2;
    plimi = out.predlims{1};
    yl = plimi{1}(1,:);
    yu = plimi{1}(2*nn-1,:);
    
    km = mean(chain);
    ks = std(chain);
    
    %add back the parameters not included in the fits so we can add the
    %values to titles, etc.
    missingInd = find(~ismember(names, string(results.names)));
    km(missingInd) = p0(missingInd);
    ks(missingInd) = NaN;
    
    
    
    yf = plimi{1}(nn,:);
    
    for i = 1:nSets
        
        yy = yf(X2(:, 2)==i);
        yyl = yl(X2(:, 2)==i);
        yyu = yu(X2(:, 2)==i);
        nexttile;
        fillyy(xx,yyl,yyu,[0.9 0.9 0.9]);
        hold on
        plot(xo{i}, yo{i}, 'o-', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'r');
        plot( xo{i}(xo{i} >= xrange(i, 1) & xo{i} <=xrange(i, 2) ), yo{i}(xo{i} >=  xrange(i, 1) & xo{i} <= xrange(i, 2)),'o-','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Color', 'r');
        plot(xx,yy,'-k')
        xlim([0,3500])
        
        if expmnt ~= "phaff"
            titleCell = enhancers{i};
            for t = 1:length(names)
                hasNoDigits = cellfun(@isempty, regexp(names', '\d'));
                if hasNoDigits(t) || contains(names(t), num2str(i))
                    titleCell = [ titleCell; join([names(t), ' = ', num2str(round2(km(t))), ' \pm ', num2str(round2(ks(t)))]) ];
                end
            end
        elseif expmnt == "phaff"
            naff = 7;
            nph = 3;
            if i <= naff
                titleCell = {enhancers{i}, [' \omega'' = ', num2str(round2(km(1))), ' \pm ', num2str(round2(ks(1)))],...
                    [ 'KD = ', num2str(round2(km(nph+i))), ' \pm ', num2str(round2(ks(nph+i)))]};
            elseif i > naff
                titleCell = {enhancers{i}, [' \omega'' = ', num2str(round2(km(i - (naff-1) ))), ' \pm ', num2str(round2(ks(i - (naff-1))))],...
                    [ 'KD = ', num2str(round2(km(nph+1))), ' \pm ', num2str(round2(ks(nph+1)))]};
            end
        end
        
        if metric=="fraction"
            ylim([0, 1])
        elseif metric=="fluo"
            ylim([0, y_max])
        end
        title(titleCell);
    end
    
        figure;
        if expmnt == "phaff"
            
            tilo = tiledlayout('flow');
            nexttile;
            errorbar(scores, km(nph+1:nph+naff), ks(nph+1:nph+naff));
            ylabel('KD (au)')
            xlabel('affinity (Patser score)')
            
            nexttile;
            errorbar(positions, km(1:nph), ks(1:nph));
            ylabel('w'' (au)')
            xlabel('position (bp)')
            xlim([-10, 2])
        else
            
            if expmnt == "affinities"
                hasKD = contains(names, 'kd', 'IgnoreCase', true);
                errorbar(scores, km(hasKD), ks(hasKD));
                ylabel('KD (au)')
                xlabel('affinity (Patser score)')
            elseif expmnt == "phases"
                scores = positions;
                hasw = contains(names, 'w', 'IgnoreCase', true);
                errorbar(scores, km(hasw), ks(hasw));
                ylabel('w'' (au)')
                xlabel('position (bp)')
                xlim([-10, 2])
            end
            
    end
    
    
    covfig = figure;
    tiledlayout('flow')
    nexttile;
    cv = @(x, y) sqrt(abs(x)) ./ sqrt((y'*y)); %sketch normalized covariance i made up
    imagesc(cv(results.cov, results.mean));
    colorbar;
    ylabel('parameter 1')
    xlabel('parameter 2')
    title('Norm. covariance matrix of the parameters');
    colormap(viridis);
    nexttile;
    rho = @(x, y) x ./ (y'*y); %pearson's correlation coefficient
    imagesc(rho(results.cov, sqrt(diag(results.cov))));
    colorbar;
    ylabel('parameter 1')
    xlabel('parameter 2')
    title('Correlation coefficient');
    colormap(viridis);
    
end

end
%%