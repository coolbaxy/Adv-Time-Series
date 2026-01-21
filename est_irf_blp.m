% est_irf_blp.m
% -------------------------------------------------------------------------
% Bayesian Local Projections (BLP) - self-contained estimator
% -------------------------------------------------------------------------
% Purpose:
%   Compute impulse responses via Bayesian Local Projections with conjugate
%   Matrix-Normal/ Inverse-Wishart posterior sampling, using the SAME
%   interface/output convention as est_irf_fev.m and the SAME data layout as
%   est_irf_lp.m (impact at h=0 included).

% Written by Miranda&Ricco 2015 IRFbayesianLocalProj.m 

% Modified by Liu LI(2025) for suitable with the analysis

% -------------------------------------------------------------------------

function [data, vmed, vars, datadraws, y, x] = est_irf_blp_conjugate( ...
    nirf, bound, ndraws, nconst, prior, begin_date, end_date, ...
    esty1, estq1, esty2, estq2, nlags, choose_vars, samp, DATA_MASTER, labels_vars, varargin)

% keep interface compatibility (some inputs unused here)
%#ok<*NASGU>

% ---- 1) Data: identical pipeline as est_irf_fev / est_irf_lp
[y, x, vars] = data_choosevars_temp(begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
[T, n] = size(y);

% ---- 2) Options
p = inputParser;
p.FunctionName = 'est_irf_blp_conjugate';
addParameter(p, 'identification', 'CHOL', @(s)ischar(s) || isstring(s));
addParameter(p, 'shockVar', 1, @(z)isnumeric(z) && isscalar(z) && z>=1);
addParameter(p, 'shockSize', 1, @(z)isnumeric(z) && isscalar(z));
addParameter(p, 'scale100', true, @(z)islogical(z) && isscalar(z));
addParameter(p, 'tau', 10, @(z)isnumeric(z) && isscalar(z) && z>0);
addParameter(p, 'nu0', [], @(z)isnumeric(z) && (isempty(z) || (isscalar(z) && z>0)));
addParameter(p, 'S0', [], @(z)isnumeric(z) || isempty(z));
addParameter(p, 'M0', [], @(z)isnumeric(z) || isempty(z));
addParameter(p, 'V0', [], @(z)isnumeric(z) || isempty(z));
addParameter(p, 'dates', [], @(z)isnumeric(z) || isdatetime(z) || isempty(z));
addParameter(p, 'instrument', struct(), @(z)isstruct(z));
addParameter(p, 'selectedInstrument', '', @(s)ischar(s) || isstring(s));
addParameter(p, 'seed', [], @(z)isnumeric(z) && (isempty(z) || isscalar(z)));
parse(p, varargin{:});
opt = p.Results;

iScheme  = upper(string(opt.identification));
shockVar = opt.shockVar;
shockSz  = opt.shockSize;

if shockVar > n
    error('shockVar=%d exceeds number of variables n=%d.', shockVar, n);
end

if ~isempty(opt.seed)
    rng(opt.seed);
end

nL = nlags;
nH = nirf;

% regression dimension: intercept + n*nL regressors
k = 1 + n*nL;

% Prior parameters
tau = opt.tau;

if isempty(opt.nu0)
    nu0 = n + 2;          % weakly-informative default
else
    nu0 = opt.nu0;
end

if isempty(opt.S0)
    S0 = eye(n);          % default
else
    S0 = opt.S0;
end

if isempty(opt.M0)
    M0 = zeros(k, n);     % default
else
    M0 = opt.M0;
end

if isempty(opt.V0)
    V0 = (tau^2) * eye(k); % default ridge-like prior on coefficients
else
    V0 = opt.V0;
end

% Basic checks
if ~isequal(size(S0), [n n])
    error('S0 must be n x n.');
end
if ~isequal(size(M0), [k n])
    error('M0 must be k x n where k=1+n*nL.');
end
if ~isequal(size(V0), [k k])
    error('V0 must be k x k where k=1+n*nL.');
end

% Pre-allocate draws: (nH+1) x n x ndraws
datadraws = NaN(nH+1, n, ndraws);

% ---- 3) Precompute X_h and Y_h for each horizon to speed up
Xcell = cell(nH,1);
Ycell = cell(nH,1);

for h = 1:nH
    nT = T - (nL + h) + 1;
    if nT <= 0
        error('Not enough observations: T=%d, nL=%d, h=%d.', T, nL, h);
    end

    YhLag = NaN(nT, n*nL);
    for j = h:(nL + h - 1)
        YhLag(:, n*(j-h)+1:n*(j-h+1)) = y(nL + h - j : end - j, :);
    end

    Xcell{h} = [ones(nT,1) YhLag];     % nT x k
    Ycell{h} = y(nL + h : end, :);     % nT x n
end

% ---- 4) Posterior draws + IRFs
for draw = 1:ndraws

    % identify B0 from horizon h=1 draw
    impactRow = [];

    for h = 1:nH

        Xh = Xcell{h};
        Yh = Ycell{h};
        nT = size(Xh,1);

        % Posterior: Vn, Mn
        iV0 = inv_pd(V0);

        A = iV0 + (Xh' * Xh);             % k x k
        Vn = inv_pd(A);                    % k x k
        Mn = Vn * (iV0 * M0 + Xh' * Yh);    % k x n

        % Posterior IW scale
        E  = Yh - Xh * Mn;
        Sn = S0 + (E' * E) + (Mn - M0)' * iV0 * (Mn - M0);
        nun = nu0 + nT;

        % Draw Sigma_h ~ IW(Sn, nun)
        Sigma = iwishrnd_local(Sn, nun);

        % Draw B_h | Sigma ~ MN(Mn, Vn, Sigma)
        LV = chol_pd(Vn);
        LS = chol_pd(Sigma);
        Z  = randn(k, n);
        Bh = Mn + LV * Z * LS';

        % coefficients on y_{t-h} block
        irfCoeffs = Bh(2:1+n, :);  % n x n

        % Identify B0 only once (at h=1)
         if h == 1
            switch iScheme
                case "CHOL"
                    L  = chol_pd(Sigma);
                    B0 = L * diag(1 ./ diag(L));         % normalize diag to ones
                case "PSVAR"
                    B0 = psvar_ident_from_u(Yh, Xh, Bh, opt, nL);
                otherwise
                    error('Unknown identification scheme: %s. Use CHOL or PSVAR.', iScheme);
            end

            impactRow = (shockSz * B0(:, shockVar)).';   % 1 x n
            datadraws(1, :, draw) = impactRow;           % h=0
        end

        datadraws(h+1, :, draw) = impactRow * irfCoeffs; % h>=1
    end
end

% ---- 5) Summaries: median + quantile bands
data = zeros(nH+1, n, 3);
data(:,:,1) = median(datadraws, 3);

qL = 100 * bound;
qU = 100 * (1 - bound);

data(:,:,2) = prctile(datadraws, qL, 3);
data(:,:,3) = prctile(datadraws, qU, 3);

% ---- 6) FEV placeholder
vmed = NaN(n,1);

% ---- 7) Scaling
if opt.scale100
    data      = data * 100;
    datadraws = datadraws * 100;
end

end




% =========================================================================
% Numerics: PD cholesky (lower) with small jitter if needed
% =========================================================================
function L = chol_pd(A)
A = (A + A')/2;
jitter = 0;
maxTries = 5;
for t = 0:maxTries
    try
        if jitter == 0
            L = chol(A, 'lower');
        else
            L = chol(A + jitter*eye(size(A)), 'lower');
        end
        return
    catch
        jitter = max(1e-12, 10^(-10+t));
    end
end
error('chol_pd failed: matrix not positive definite even after jitter.');
end

function iA = inv_pd(A)
L = chol_pd(A);
iA = L'\(L\eye(size(A)));
end


% =========================================================================
% Inverse-Wishart draw WITHOUT toolboxes (Bartlett decomposition)
% Sigma ~ IW(S, nu)
% =========================================================================
function Sigma = iwishrnd_local(S, nu)

p = size(S,1);
if nu < p
    error('IW degrees of freedom nu must be >= dimension p.');
end

Sinv = inv_pd(S);
L = chol_pd(Sinv); % lower so that Sinv = L*L'

A = zeros(p);
for i = 1:p
    df = nu - i + 1;
    A(i,i) = sqrt(sum(randn(df,1).^2)); % sqrt(chi2)
    for j = 1:i-1
        A(i,j) = randn();
    end
end

W = L * (A*A') * L';
Sigma = inv_pd(W);

end
