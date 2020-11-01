function [b, MSE] = globfit2(varargin)
% Set up data so that Y is a function of T with a specific functional form,
% but there are multiple groups and one parameter varies across groups.

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')

expmnt = "affinities";
md = "simpleweak";
metric = "fraction";
displayFigures = true;
minKD = 0;
maxKD = 2E4;
%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if strcmpi(expmnt, 'affinities')
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
elseif strcmpi(expmnt, 'phases')
    enhancers = {'1Dg11', '1Dg-5', '1Dg-8D'};
    scores = [0, -5, -8]';
end



%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
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
for k = 1:nSets
    cond = strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k});
    xo{k} = dorsalResultsDatabase.dorsalFluoBins(cond);
    if metric == "fraction"
        yo{k} = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
    elseif metric == "fluo"
        %         yo{k} = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond);
        yo{k} = dorsalResultsDatabase.meanallmrnasEmbryo(cond);
    end
    [xs{k}, ys{k}]= processVecs(xo{k}, yo{k}, xrange(k, :));
    dsid = [dsid; k*ones(size(xs{k}))];
    
    if isnan(xrange(k, 1))
        xrange(k, 1) = min(xo{k});
    end
    if isnan(xrange(k, 2))
        xrange(k, 2) = max(xo{k});
    end
    
    assert(~isempty(xs{k}));

    T = [T; xs{k}];
    Y = [Y; ys{k}];
end

% if displayFigures
%     scatterFig = figure;
%     gscatter(T,Y,dsid)
% end



y_max = nanmax(Y(:));
x_max = max(T);

[p0, lb, ub] = getInits(expmnt, md, metric ,x_max, y_max, nSets, 'minKD', minKD, 'maxKD', maxKD);





%%
% Pack up the time and dataset id variables into X for later unpacking
X = [T dsid];

optimoptions = optimset('TolFun',1E-6, 'MaxIter', 1E6, 'MaxFunEvals', 1E5);

if expmnt == "affinities" && md=="simpleweak" && metric=="fraction"
    mdlo = @subfun_simplebinding_weak_fraction;
    mdl = @(p, x) subfun_simplebinding_weak_fraction_std(p, x);
elseif expmnt == "affinities" && md=="simpleweak" && metric=="fluo"
    mdlo = @subfun_simplebinding_weak_fluo;
    mdl = @(p, x) subfun_simplebinding_weak_fluo_std(p, x);
elseif expmnt == "phases" && md=="simpleweak" && metric=="fraction"
    mdlo = @subfun_simplebinding_weak_fraction_phases;
    mdl = @(p, x) subfun_simplebinding_weak_fraction_std_phases(p, x);
elseif expmnt == "phases" && md=="simpleweak" && metric=="fluo"
    mdlo = @subfun_simplebinding_weak_fluo_phases;
    mdl = @(p, x) subfun_simplebinding_weak_fluo_std_phases(p, x);
end

[b,~,res,~,~,~, J] = lsqcurvefit(mdlo,p0,X,Y, lb, ub, optimoptions);


CI = nlparci(b,res,'jacobian',J);
MSE = mean(res.^2);

CovB = inv(J'*J).*MSE;

covfig = figure;
cv = @(x, y) sqrt(abs(x)) ./ sqrt((y'*y));
imagesc(cv(CovB, b));
colorbar;
ylabel('parameter 1')
xlabel('parameter 2')
title('Covariance matrix of the parameters');

if displayFigures
    xx = (0:1:max(X(:,1)))';
    xxx = repmat(xx, nSets, 1);
    [Ypred,delta] = nlpredci(mdl,xxx,b,res,'Jacobian',full(J));
    yl = Ypred - delta;
    yu = Ypred + delta;
    
    
    figure(1);
    til = tiledlayout(1, nSets);
    
    dsid2 = [];
    for k = 1:nSets
        dsid2 = [dsid2; k*ones(length(xx), 1)];
    end
    
    X2 = [repmat(xx, nSets, 1), dsid2];
    
    if expmnt == "affinities" && md=="simpleweak" && metric=="fraction"
        yfit2 = subfun_simplebinding_weak_fraction(b, X2);
    elseif expmnt == "affinities" && md=="simpleweak" && metric=="fluo"
        yfit2 = subfun_simplebinding_weak_fluo(b, X2);
    elseif expmnt == "phases" && md=="simpleweak" && metric=="fraction"
        yfit2 = subfun_simplebinding_weak_fraction_phases(b, X2);
    elseif expmnt == "phases" && md=="simpleweak" && metric=="fluo"
        yfit2 = subfun_simplebinding_weak_fluo_phases(b, X2);
    end
    
    if expmnt == "affinities"
            vartheta = 'KD = ';
            consttheta = ' \omega'' = ';
    elseif expmnt == "phases"
        consttheta = 'KD = ';
        vartheta = ' \omega'' = ';
    end
    
    for k = 1:nSets
        
        yy = yfit2(X2(:, 2)==k);
        yyl = yl(X2(:, 2)==k);
        yyu = yu(X2(:, 2)==k);
        nexttile;
        plot(xo{k}, yo{k}, 'o-', xx, yy, '-', xx,yyl, '--r', xx, yyu, '--r');
        %     ylim([0, max(yfit2)*1.1]);
        %     xlim([0, max(xx)]);
        xlim([0, 3500]);
        if metric=="fraction"
            title({enhancers{k}, [vartheta, num2str(round2(b(k+1))), ' (', num2str(round2(CI(k+1, 1))), ' ', num2str(round2(CI(k+1, 2))) ' )'],...
                [consttheta, num2str(round2(b(1))), ' (', num2str(round2(CI(1, 1))), ' ', num2str(round2(CI(1, 2))), ' )']})
            ylim([0, 1]);
        elseif metric=="fluo"
            title({enhancers{k}, [vartheta, num2str(round2(b(k+1))), ' (', num2str(round2(CI(k+1, 1))), ' ', num2str(round2(CI(k+1, 2))) ' )'],...
                [consttheta, num2str(round2(b(1))), ' (', num2str(round2(CI(1, 1))), ' ', num2str(round2(CI(1, 2))), ' )'],...
                [' amp = ', num2str(round2(b(end-1))), ' (', num2str(round2(CI(end-1, 1))), ' ', num2str(round2(CI(end-1, 2))), ' )'],...
                [' off = ', num2str(round2(b(end))), ' (', num2str(round2(CI(end, 1))), ' ', num2str(round2(CI(end, 2))), ' )'],...
                })
            ylim([0, y_max]);
        end
        
    end
    title(til, 'global fit, gradient')
end
end



function yfit = subfun_simplebinding(params,X)

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
n = params(nSets + 2);
omegaDP = params(nSets + 3);
offset = params(nSets + 4);

yfit = amplitude.*((n+(x./KD(dsid)).*n.*omegaDP)./(1+x./KD(dsid)+n+(x./KD(dsid)).*n.*omegaDP))+offset;

end


function yfit = subfun_hill(params,X)

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
n = params(nSets + 2);
offset = params(nSets + 3);

yfit = amplitude.*(((x./KD(dsid)).^n)./(1+((x)./KD(dsid)).^n))+offset;

end


function yfit = subfun_simplebinding_weak(params,X)

%simplebinding in the weak promoter limit.

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
omegaDP = params(nSets + 2);
offset = params(nSets + 3);

yfit = amplitude.*(((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP))+offset;

end

%% simple weak fraction

function yfit = subfun_simplebinding_weak_fraction(params,X)
%simplebinding in the weak promoter limit.
x = X(:,1);        % unpack time from X

X3 = [];
n = 1;
X3(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X3(k, 1) = x(k);
    X3(k, 2) = n;
end
dsid3 = X3(:, 2);
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

omegaDP = params(1);
KD = params(2:nSets+1)';

yfit = (((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP));

end


function yfit = subfun_simplebinding_weak_fraction_std(params, x)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

omegaDP = params(1);
KD = params(2:nSets+1)';

yfit = (((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP));

end

%% simple weak fluo

function yfit = subfun_simplebinding_weak_fluo(params,X)
%simplebinding in the weak promoter limit.
x = X(:,1);        % unpack time from X

X3 = [];
n = 1;
X3(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X3(k, 1) = x(k);
    X3(k, 2) = n;
end
dsid3 = X3(:, 2);
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

omegaDP = params(1);
KD = params(2:nSets+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
yfit = amp.*(((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP)) + offset;

end


function yfit = subfun_simplebinding_weak_fluo_std(params, x)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

omegaDP = params(1);
KD = params(2:nSets+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
yfit = amp.*(((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP)) + offset;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% phases

function yfit = subfun_simplebinding_weak_fraction_phases(params,X)
%simplebinding in the weak promoter limit.
x = X(:,1);        % unpack time from X

X3 = [];
n = 1;
X3(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X3(k, 1) = x(k);
    X3(k, 2) = n;
end
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

KD = params(1);
omegaDP = params(2:nSets+1)';

yfit = (((x./KD).*omegaDP(dsid))./(1 + (x./KD) +((x./KD).*omegaDP(dsid))));

end


function yfit = subfun_simplebinding_weak_fraction_std_phases(params, x)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

KD = params(1);
omegaDP = params(2:nSets+1)';

yfit = (((x./KD).*omegaDP(dsid))./(1+ (x./KD)+ ((x./KD).*omegaDP(dsid))));

end

%% simple weak fluo

function yfit = subfun_simplebinding_weak_fluo_phases(params,X)
%simplebinding in the weak promoter limit.
x = X(:,1);        % unpack time from X

X3 = [];
n = 1;
X3(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X3(k, 1) = x(k);
    X3(k, 2) = n;
end
dsid3 = X3(:, 2);
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

KD = params(1);
omegaDP = params(2:nSets+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
yfit = amp.*(((x./KD).*omegaDP(dsid))./(1 + (x./KD) + ((x./KD).*omegaDP(dsid)))) + offset;

end


function yfit = subfun_simplebinding_weak_fluo_std_phases(params, x)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

KD = params(1);
omegaDP = params(2:nSets+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
yfit = amp.*(((x./KD).*omegaDP(dsid))./(1+ (x./KD) +((x./KD).*omegaDP(dsid)))) + offset;

end







% %% example function if i want to use different functional forms for
% different datasets
% function yfit = subfun(params,X)
% T = X(:,1);        % unpack time from X
% dsid = X(:,2);     % unpack dataset id from X
% A0 = params(1);    % same A0 for all datasets
% A1 = params(2:4)'; % different A1 for each dataset
% tau = params(5);   % same tau
% yfit = zeros(size(T));
% % Find separate groups and call the F function for each group
% idx = (dsid==1);
% yfit(idx) = F1(X(idx),[A0 A1(1) tau]);
% idx = (dsid==2);
% yfit(idx) = F2(X(idx),[A0 A1(2) tau]);
% idx = (dsid==3);
% yfit(idx) = F3(X(idx),[A0 A1(3) tau]);
% end
% % Below are the three functions for the three groups
% function y = F1(x,b)
% y = b(1) + b(2)*exp(-x/b(3));
% end
% function y = F2(x,b)
% y = b(1) + b(2)*sin(x/b(3));
% end
% function y = F3(x,b)
% y = b(1) + b(2)*cos(x/(2*b(3)));
% end