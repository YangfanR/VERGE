function [Cs C Cproject Ln] = covComp(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r)
% Computes covariance matrix, Ln, (n x n) reflecting marginalization over
% GP variate, z.  This function doesn't invert Ln to create likelihood term
% because not needed when used. 
Jm = Js^2*eye(m);  
Cs = Rusym1TMat(diffsqs,Ismall,m,m,rho,sigcinv,lambdazinv) + Jm; %m x m psample cov matrix
Gs = chol(Cs,'lower');
C = Rusym1TMat(diffsq,I,m,n,rho,sigcinv,lambdazinv); % m x n psample - sample cov matrix
U = Gs\C;
Cproject = U'*U;
Ln = 1/r*eye(n) + Cproject;

