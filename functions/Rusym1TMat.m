function [R] = Rusym1TMat(diffsq,J,n1,n2,rho,sigcinv,lambdazinv)
% Builds GP covariance from X (through vector of squared differences,
% diffsq), and covariance parameters (rho, sigc, lambdaz).  These are 
% input in their inverted values for faster computation.
selr = rho~=1;
diffsq = diffsq(:,selr); %stacked row matrix formulation
rho = rho(selr);
rVec = sigcinv + lambdazinv*exp(diffsq*log(rho)');
rVec = rVec(J);
% rVec = exp(diffsq*log(rho(selr))');
R = reshape(rVec,n2,n1)';
% R = exp(R);
