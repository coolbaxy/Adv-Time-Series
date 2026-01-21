
% extracts current Patent-Based News shock using Cholesky decomposition
%
% Adjusted the original function by Andre Kurmann, Federal Reserve Board; last modified December 2012 
%-------------------------------------------------------------------------


%function [output,v]=current_tfp(b,res,vmat,nvars,nlags,nimp,slope,use_yields)

function [output,v]=patent_news(b,res,vmat,nvars,nlags,nimp)
%companion matrix of demeaned VAR (=> don't include constant term)
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=b(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

% Extract impulse vectors 
%------------------------

%unique lower triangular matrix = Cholesky decomp of vmat
Gtilde = chol(vmat)';
gamma = zeros(nvars,1); gamma(1) = 1;  %'selection vector' of same format as in sims extraction


% gamma_proba = zeros(nvars,1); gamma_proba(1)=1;
% gamma_proba1 = zeros(nvars,1); gamma_proba1(3)=1;
%Note: TFP must be ordered first
alpha = Gtilde*gamma;  %extracts first column of cholesky

% Back out time series of structural shock 
v = gamma'*inv(Gtilde)*res';    %1xT vector of structural shocks
v = v';
   

    % Compute impulse responses--------------------------------------------------------------------------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    U1=[alpha; zeros(nvars*nlags-nvars,1)];
    n = size(U1,1); 
   
        for k=1:nimp;
            Zk1(k,:)=(M^(k-1)*U1)';
        end
  
    
    impulse1(:,:)=Zk1(:,1:nvars);

    if Zk1(1,1)<0 %Zk1(1,1) < 0   
       alpha = - alpha; 
       v = -v;
       impulse1 = -impulse1;

    end
% 
% %     
    % Compute fraction of VD explained at different horizons ------------------------------------------------------------------------------------------------------------------------------------------
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    B = zeros(nvars,nvars,nimp+1);
    for l=0:nimp;
        C=M^l;
        B(:,:,l+1) = C(1:nvars,1:nvars);
    end
    sigmak=B(:,:,1)*vmat*B(:,:,1)';
    hh1=B(:,:,1)*alpha*(B(:,:,1)*alpha)';
    vardec1(1,:)=(diag(hh1./sigmak))';

    for k=1:nimp-1;
        sigmak=sigmak+B(:,:,k+1)*vmat*B(:,:,k+1)';
        hh1=hh1+B(:,:,k+1)*alpha*(B(:,:,k+1)*alpha)';
        vardec1(k+1,:)=(diag(hh1./sigmak))';
 end      
    
    output=[vardec1 impulse1];

