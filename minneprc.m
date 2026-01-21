% Adjusted from Andre Kurmann

% minneprc.m
%
% Produces diagonal elements of inv(H) and inv(H)*bprior of Minnesota prior
% to compute posterior mean of reduced form regression coefficients.
%
% Theory:   for a single T x 1 regression y = X*b + u with prior distribution
%           b ~ N(bprior,M), the posterior mean is
%
%           bpost = inv[inv(M) + X'X] * [X'y + inv(M)*bprior]
% 
%           for a multivariate T x n regression Y = X*b + U
%           where X is T x np and b is np x n, with stacked prior
%           vec(b) ~ N(bprior,H) and H diagonal (!!), the posterior mean is
%
%           vec(bpost) = inv[inv(H) + kron(I,X'X)] * [vec[X'Y] + inv(H)*bprior]
%
% Inputs:   Y = LHS data of regressions 
%           nlag = number of lags in VAR
%           quarterly = 1 if quarterly estimation; = 0 if monthly
%           const = 1 if regressions include constant; = 0 otherwise 
%           lev = nvar x 1 vector  with elements
%                       = 1 if estimation in levels [unit root prior]
%                       = 0 otherwise  
%           prior: weight on prior on lags with
%                       1 = standard prior
%                       inf = no prior
%                       < 1 =  tighter prior
%
% Outputs:  hm = diagonal elements of inv(H)
%           bm = inv(H)*bprior            


function [hm,bm] = minneprc(y,x,nlag,quarterly,const,lev,prior);

[t,nvar] = size(y);      %number of variables
k = const+nvar*nlag; %number of regressors per equation   

%values for tightness paramters
lambda0 = prior;
lambda1 = .2;       % = gamma in Hamilton pg 361
lambda2 = .5;       % = w
lambda3 = 1;        %
lambda4 = 1e5;      % prior on intercept



%create monthly lag decay to match 1/q decay in quarterly data where q = quarters
if quarterly == 1; 
    ld=seqa(1,1,nlag).^(-lambda3);    %produces a regularly sequence of lag decay [note ^-mu[4]] */
else
	j=ceil(nlag/3)^(-lambda3);   % last quarter [rounded up] eg. l2=13=>xx2=5
	b=0;
	if nlag>1;
		b= ( log(1)-log(j) ) / (1-nlag);
    end
	a=exp(-b);
	ld=a*exp(b*seqa(1,1,nlag));  % Tao's lag decay to match 13th lag 
end

%add other parameters and take the squared inverse
%this is H^-1 for own lags
ld = (lambda0*lambda1*lambda2*ld).^(-2);  %note: the weight lambda2 is taken
                                          %off from the own lag coeffs
                                          %later on (line 107)

%Scale factors from univariate OLS AR's to compute stdevs of priors of lags
%of other variables
s = zeros(nvar,1);
i=1;

%compute variance of residuals of univariate AR regressions
for i=1:nvar;
   %y = Y(nlag+1:rows(Y),i);
   %t = rows(y);
   %x = ones(t,1);
   if const==1;
     xi=x(:,k);  
   elseif const==2;
   xi=x(:,(k-1):k);
   end
   for j=1:nlag; 
      xi=[xi x(:,i+(j-1)*nvar)];
       %x=[x Y(nlag+1-j:rows(Y)-j,i)];  
   end
   bsh = inv(xi'*xi)*xi'*y(:,i);
   u = y(:,i) - xi*bsh;
   s(i) = (u'*u)/t;
end

%prepare construction of inv(H) for each equation and add prior for intercept
if lambda4>0;
    if const==1;
        H = [kron(ld,s); ((lambda0*lambda4).^(-2))];
    elseif const==2;
        H = [kron(ld,s); ((lambda0*lambda4).^(-2));((lambda0*lambda4).^(-2)) ];
    else
        H = kron(ld,s);
    end

elseif lambda4==0;
    if const
        H = [kron(ld,s); 0];
    else
        H = kron(ld,s);
    end
end

%construct output vectors
hm = zeros(k*nvar,1);
bm = zeros(k*nvar,1);

% stack the H and distinguish own lag vs other 
for i=1:nvar;
    hadd = H;      %NOTE: H/s(i) MEANS WE NORMALIZE BY VARIANCE OF DEPENDENT VARIABLE 
    %for j=0:(nlag-1);
        for j=0:(nlag-1);
       hadd(i+(j*nvar)) = (lambda2^2)*H(i+(j*nvar));
    end
    hm(k*(i-1)+1:k*i,1) = hadd;       %stacking for each equation

    badd = zeros(k,1);
    
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXX%
    if lev(i)==1;   % if lev=1, put unit root prior on own first lag; otherwise prior = 0 
       bm(k*(i-1)+i,1)=hm(k*(i-1)+i)*1; %note: identical results than with Chris' formula
    end
end