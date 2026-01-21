% est_irf_lp.m
% -------------------------------------------------------------------------
% Local Projection (Jorda 2005) IRFs with HAC bands

% Original Code: From Replication Code of Bayesian Local Projection(2025),
% Written by Miranda (2014) IRFlocalProj.m 
% 
% 

% Modified by Liu LI(2025) for suitable with the analysis

function [data, vmed, vars, datadraws, y, x] = est_irf_lp( ...
    nirf, bound, ndraws, nconst, prior, begin_date, end_date, ...
    esty1, estq1, esty2, estq2, nlags, choose_vars, samp, DATA_MASTER, labels_vars, varargin)


% --- 1) Use the same data pipeline as est_irf_fev ---
[y, x, vars] = data_choosevars_temp(begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
[T, n] = size(y);

% --- 2) Parse options ---
p = inputParser;
p.FunctionName = 'est_irf_lp';
addParameter(p, 'identification', 'CHOL', @(s)ischar(s) || isstring(s));
addParameter(p, 'shockVar', 1, @(z)isnumeric(z) && isscalar(z) && z>=1);
addParameter(p, 'shockSize', 1, @(z)isnumeric(z) && isscalar(z));
addParameter(p, 'bandsCoverage', 68, @(z)isnumeric(z) && isscalar(z) && z>0 && z<100);
addParameter(p, 'scale100', true, @(z)islogical(z) && isscalar(z));
addParameter(p, 'dates', [], @(z)isnumeric(z) || isdatetime(z) || isempty(z));
addParameter(p, 'instrument', struct(), @(z)isstruct(z));
addParameter(p, 'selectedInstrument', '', @(z)ischar(z) || isstring(z));
parse(p, varargin{:});
opt = p.Results;

iScheme = upper(string(opt.identification));
shockVar = opt.shockVar;
shockSize = opt.shockSize;
nL = nlags;
nH = nirf;

% z critical value for two-sided bands
sLevel = abs(norminv((1 - opt.bandsCoverage/100)/2, 0, 1));

% --- 3) Compute LP IRFs + bands in (nH+1) x n format ---
[irf, irf_l, irf_u, ~] = lp_irf_hac(y, nL, nH, iScheme, shockVar, shockSize, sLevel, opt);

% --- 4) Pack outputs to match est_irf_fev convention ---
data = zeros(nH+1, n, 3);
data(:,:,1) = irf;
data(:,:,2) = irf_l;
data(:,:,3) = irf_u;

datadraws = zeros(nH+1, n, 1);
datadraws(:,:,1) = irf;

vmed = NaN(n,1);  % LP version does not compute FEV shares

if opt.scale100
    data      = data * 100;
    datadraws = datadraws * 100;
end

end


% =========================================================================
% Core LP routine (HAC bands), returns horizon-by-variable matrices
% =========================================================================
function [irf, irf_l, irf_u, B0] = lp_irf_hac(y, nL, nH, iScheme, shockVar, shockSize, sLevel, opt)

[T, n] = size(y);

irf   = NaN(nH+1, n);
irf_l = NaN(nH+1, n);
irf_u = NaN(nH+1, n);

% identify B0 using h=1 regression residuals 
B0 = [];

for h = 1:nH

    % Build lag matrix for this horizon:
    % Dependent: y_{t} with t = nL+h ... T
    % Regressors: [1, y_{t-h}, y_{t-h-1}, ..., y_{t-h-(nL-1)}]
    nT = T - (nL + h) + 1;
    if nT <= 0
        error('Not enough observations: T=%d, nL=%d, h=%d.', T, nL, h);
    end

    YhLag = NaN(nT, n*nL);
    for j = h:(nL + h - 1)
        % regressor block corresponds to y_{t-j}
        YhLag(:, n*(j-h)+1:n*(j-h+1)) = y(nL + h - j : end - j, :);
    end
    Yh = y(nL + h : end, :);

    X = [ones(nT,1) YhLag];

    % OLS local projection
    projCoeffs = X \ Yh;              % (1+n*nL) x n, equations in columns
    irfCoeffs  = projCoeffs(2:n+1,:); % coefficients on the first block (y_{t-h}); size n x n

    % residuals
    u = Yh - X * projCoeffs;          % nT x n

    % Identification (once, at h=1)
    if h == 1
        SigmaU = cov(u);

        switch iScheme
            case "CHOL"
                % lower-triangular impact matrix with unit diagonal normalization
                L  = chol(SigmaU, 'lower');
                B0 = L * diag(1 ./ diag(L));   % normalize diag to ones

            case "PSVAR"
                % Requires opt.dates, opt.instrument, opt.selectedInstrument, and ProxySVARidentification.m in path.
                B0 = psvar_ident(u, opt, nL);

            otherwise
                error('Unknown identification scheme: %s. Use CHOL or PSVAR.', iScheme);
        end

        % Impact response (h=0)
        impactVec = shockSize * (B0(:, shockVar));   % n x 1
        irf(1,:)  = impactVec.';                     % 1 x n

        % Set impact bands equal to point estimate (avoid NaNs)
        irf_l(1,:) = irf(1,:);
        irf_u(1,:) = irf(1,:);
    end

    % Response at horizon h: (impact' * irfCoeffs) gives 1 x n
    impactRow   = (shockSize * B0(:, shockVar)).';
    irf(h+1,:)  = impactRow * irfCoeffs;

    % -------------------------------
    % HAC error bands (Hamilton-style)
    % -------------------------------
    % Center residuals
    u = bsxfun(@minus, u, mean(u,1));

    % HAC covariance estimator for residuals (no 1/nT scaling, following IRFlocalProj.m)
    SigmaHAC = u' * u;   % n x n
    nwLags   = nL + h + 1;

    nwWeights = (nwLags + 1 - (1:nwLags)) ./ (nwLags + 1);
    for j = 1:nwLags
        if j >= nT
            break; % not enough obs for this lag
        end
        Gammaj  = (u(j+1:nT,:)' * u(1:nT-j,:));
        SigmaHAC = SigmaHAC + nwWeights(j) * (Gammaj + Gammaj');
    end

    % Q = (X'X)^{-1}
    Q = inv(X' * X);

    % Approximate block-diagonal variance of coefficients (each equation separately),
    % keep only diagonal residual variances (as in IRFlocalProj.m)
    SigmaB = kron(diag(diag(SigmaHAC)) / max(nT - nL, 1), Q);

    for eq = 1:n
        % Extract variance of the "first block" coefficients (y_{t-h}) for equation eq
        blockSize = n*nL + 1;      % constant + n*nL regressors
        row0 = blockSize*(eq-1);

        idx = (row0 + 2) : (row0 + 1 + n);   % coefficients on first block (size n)
        sigmaB = SigmaB(idx, idx);

        % delta method variance for response = impactRow * beta(:,eq)
        OmegaVar = impactRow * sigmaB * impactRow.'; % scalar
        OmegaSE  = sqrt(max(OmegaVar, 0));
        halfWidth = OmegaSE * sLevel;

        irf_u(h+1, eq) = irf(h+1, eq) + halfWidth;
        irf_l(h+1, eq) = irf(h+1, eq) - halfWidth;
    end
end

end


