function [R] = Rusym1TPartial(Rold,diffsq,J,n1,n2,position,lastrhoPosition,candrhoPosition,sigcinv)
% Builds new GP covariance matrix from old one by removing term for old
% value of rho (latrhoPosition) and replacing it with new value of 
% rho (candrhoPosition).  This is used in sequential scan sampling of
% (gamma_k, rho_k).
%
% logDeltaRho = log(candrhoPosition/lastrhoPosition);
% deltaVec = exp(diffsq(:,position)*logDeltaRho);
deltaVec = exp(diffsq(:,position)*(log(candrhoPosition/lastrhoPosition)));
deltaVec = deltaVec(J);
% Delta = reshape(deltaVec,n2,n1)';
% Told = Rold - sigcinv;
% R = Told.*Delta + sigcinv;
R = (Rold-sigcinv).*reshape(deltaVec,n2,n1)' + sigcinv;