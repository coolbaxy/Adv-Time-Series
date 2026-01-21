% Estimate strucutral VAR, get IRFs and FEV
% Adjusted codes of Andre Kurmann, AER 2010 (Kurmann and Otrok)

% Modified by Liu LI (2025)

function [data, vmed, vars,datadraws,y,x] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);

[y, x, vars, lev]=data_choosevars_temp(begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars); 
[T,nvars]=size(y);
%--------------------------------------------------------------------
% Step 3: Estimate VAR, extract shocks, FEV shares explained, and IRFs 
%---------------------------------------------------------------------
if nconst==1;
   const = ones(T,1); 
   x = [x const];
   ncoeffs = nvars*nlags+1;
else
   ncoeffs = nvars*nlags; 
end;

%ncoeffs + 1 x nvars matrix of regression coefficients
b       = x\y;       %OLS coefficients: ncoeffs x nvars
res     = y - x*b;          %T x nvar matrix of residuals
vmat    = (1/T)*(res'*res);
stds    = sqrt(diag(vmat));

% ==============================================================================================================================================================================================
% ==============================================================================================================================================================================================
%Bayesian estimation with Minnesota prior
     

[hm,bm]     = minneprc(y,x,nlags,1,nconst,lev,prior); %Minnesota prior
xxx         = inv(kron(eye(nvars),x'*x) + diagrv(eye(size(hm,1)),hm));
bb          = xxx * (vec(x'*y) + bm);       %stacked posterior means of regression coeffs: vec(bpost) = inv[inv(H) + kron(I,X'X)] * [vec[X'Y] + inv(H)*bprior]
b           = reshape(bb,ncoeffs,nvars);    %posterior coefficient matrix 
res         = y - x*b;                      %T x nvar matrix of residuals
vmat        = (1/T)*(res'*res);

sxx         = chol(xxx)';
sinv        = chol(inv(vmat));
datadraws   =[];
vdraws      =[];
randn('state',0);
randmatrix  = randn(nvars,T,ndraws);
randvec     = randn(nvars*ncoeffs,ndraws);

for j=1:ndraws;
    [bbdraw,vmatdraw] = niw(b,sxx,sinv,T,randmatrix(:,:,j),randvec(:,j)); %drawing from posterior coefficient matrix
	bdraw=reshape(bbdraw,ncoeffs,nvars);
	res = y - x*bdraw;
	%compute results for each draw         
	[output,v]=patent_news(bdraw,res,vmatdraw,nvars,nlags,nirf);
    datadraws(:,:,j) = output;
	vdraws(:,j) = v;
end

%compute median and confidence intervals   
datadraws   = sort(datadraws,3);
low         = round(ndraws*bound); 
high        = round(ndraws*(1-bound));
datalow     = datadraws(:,:,low);
datamed     = median(datadraws,3);
datahigh    = datadraws(:,:,high);
data(:,:,1) = datamed; 
data(:,:,2) = datalow; 
data(:,:,3) = datahigh;
vmed        = median(vdraws,2);   
data        = data*100;
 
