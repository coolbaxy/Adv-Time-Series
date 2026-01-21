%niw.m

%procedure designed to generate drawings of parameters from a normal-inverted
%wishart distribution (or normal-inverted gamma if vmat is 1x1)
%
%format:		{betadraw,wish} = niw(betahat,sxx,sinv,df);
%
%inputs:		betahat: coefficient estimates -- xxx*(x'y), where xxx = inv(x'x)
%       		sxx: cholesky decomposition of xxx -- chol(xxx)'
%       		sinv: cholesky decomposition of inv(vmat) -- chol(inv(vmat)), where vmat = (1/df)*(res'res),
%                                                                               and res = y - xbetahat
%               df: degrees of freedom 
%
%output:		drawings of beta and drawings of vcv matrix from wishart 

function [BETADRAW,WISH] = niw(BETAHAT,SXX,SINV,DF,RANDMATRIX,RANDVEC);

NVARS = size(SINV,1);
NPARAMS = size(BETAHAT,1);
RANW = RANDMATRIX/sqrt(DF); %randn(NVARS,DF)/sqrt(DF);
RANTR = RANW'*SINV;
WISH = inv(RANTR'*RANTR);
SWISH = (chol(WISH))';
RANC = RANDVEC; %randn(NPARAMS*NVARS,1);

V=zeros(NVARS*NPARAMS,NVARS*NPARAMS);
for i=1:NVARS;
    for j=1:NVARS;
        V(1+(i-1)*NPARAMS:i*NPARAMS,1+(j-1)*NPARAMS:j*NPARAMS) = SWISH(i,j) * SXX(1+(i-1)*NPARAMS:i*NPARAMS,1+(i-1)*NPARAMS:i*NPARAMS);
    end
end
        
SHOCK = V * RANC;
BETADRAW = vec(BETAHAT) + SHOCK;