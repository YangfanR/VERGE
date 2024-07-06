function [Ln Cs C Cproject] = varGPpartial(CsOld,Cold,diffsq,diffsqs,I,Ismall,m,n,position,lastrhoPosition,candrhoPosition,sigcinv,Js,r)
% Builds new covariance matrices from old by removing old rho_k and
% replacing with new rho_k.  Used on sequential scan sampling of (gamma_k,
% rho_k) to save computation:
% INPUT:
% CsOld, COld:  Previous m x m, CsOld, and m x n, COld.  We recall that m <
%      n random data points are chosen to build CsOld and compute 
%      inverse with at lower m x m dimesion.
% diffsq, diffsqs:  Vectors of unique squared differences in covariate values
%                   used to construct Cs, C and Cproject.
% I, Ismall: indices to allow expansion from unique values to full values
%             for diffsq and diffsqs.
% position: "k" in rho_k
% lastrhoPosition, candrhoPosition: Previous and new values for rho_k
% sigcinv, r:  additional parameters needed to build full GP covariance
%              matrix, Ln
% OUTPUT:
% term : quadratic term of likelihood for GP model with z marginalized out.
% Ln : full n x n covariance matrix with z marginalized out.
% Cs, C, Cproject:  GP covariance matrices
Im = eye(m); Jm = 0.3*Js^2*Im; %Jpartial = Jm;
Cs = Rusym1TPartial(CsOld,diffsqs,Ismall,m,m,position,lastrhoPosition,candrhoPosition,sigcinv);% + Jpartial*Im*((candrhoPosition > lastrhoPosition) + 0.3*(lastrhoPosition > candrhoPosition)); %m x m psample cov matrix
if(max(Cs,[],"all") > 1.0e+11)
    Jm = (max(Cs,[],"all")/1.0e+11) * eye(m);
end
Cs = Cs + Jm;
Gs = chol(Cs,'lower');
C = Rusym1TPartial(Cold,diffsq,I,m,n,position,lastrhoPosition,candrhoPosition,sigcinv); % m x n psample - sample cov matrix
U = Gs\C;
Cproject = U'*U;
Ln = 1/r*eye(n) + Cproject;
