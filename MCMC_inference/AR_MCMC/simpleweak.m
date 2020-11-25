function yfit = simpleweak(x, params, options)
%simplebinding in the weak promoter limit.

if ~isfield(options, 'onedsid')
    options.onedsid = false;
end
if ~isfield(options, 'noOff')
    options.noOff = false;
end
if ~isfield(options, 'fraction')
    options.fraction = false;
end
if ~isfield(options, 'dimer')
    options.dimer = false;
end
if ~isfield(options, 'expmnt')
    options.expmnt = "affinities";
end


if isstruct(x)
    if isfield(x, 'X')
        data = x;
        x = data.X(:, 1);
    elseif isfield(x, 'ydata')
        x = x.ydata(:, 1);
    end
end

if ~options.onedsid
    dsid = zeros(length(x), 1);
    n = 1;
    dsid(1) = 1;
    for k = 2:length(x)
        if x(k) < x(k-1)
            n = n + 1;
        end
        dsid(k) = n;
    end
else
    dsid = ones(length(x), 1);
end

params = params(:)'; %need a row vec

if ~options.noOff
    offset = params(max(dsid)+3);
else
    offset = 0;
end

if ~options.fraction
    amp = params(max(dsid)+2);
else
    amp = 1;
    offset = 0;
end

if options.expmnt=="affinities"
    omegaDP = params(1);
    KD = params(2:max(dsid)+1)';
    KD = KD(dsid);
elseif options.expmnt == "phases"
    KD = params(1);
    omegaDP = params(2:max(dsid)+1)';
    omegaDP = omegaDP(dsid);
elseif expmnt == "phaff"
    nph = 3;
    naff = 7;
    omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
    omegaDP = omegaDP(dsid);
    KD = [params(1:naff), params(1)*ones(1, nph-1)]';
    KD = KD(dsid);
    if ~options.fraction
        amp = params(nph+naff+1);
        offset = 0;
    end
end


if options.dimer
    k = params(end); %dimerization kd
    x = (1/2).*(k+2.*x+(-1).*k.^(1/2).*(k+4.*x).^(1/2)); %concentration of dorsal dimers
end

sb = @(x, p) p{1}.*(((x./p{2}).*p{3})./(1+ (x./p{2})+ ((x./p{2}).*p{3}))) + p{4};
%sb inputs are 1. dorsal 2. {amp, KD, omegaDP, offset}. Any of the elements
%of 2 are allowed to be column vectors.
yfit = sb(x, {amp; KD; omegaDP; offset});

end