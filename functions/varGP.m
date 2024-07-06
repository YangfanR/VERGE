function [Ln Cs C Cproject] = varGP(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r)
% Compute GP covariance matrix and resultant quadratic term for posterior
% evaluation
Jm = Js^2*eye(m);
Cs = Rusym1TMat(diffsqs,Ismall,m,m,rho,sigcinv,lambdazinv); % + Jm; %m x m psample cov matrix
if(max(Cs,[],"all") > 1.0e+11)
    Jm = (max(Cs,[],"all")/1.0e+11) * eye(m);
end
% if(cond(Cs) < 1.0e-09 || cond(Cs) > 1.0e+09)
%     Cs = Cs + Jsup^2*eye(m);
% end;
Cs = Cs + Jm;
Gs = chol(Cs,'lower');
C = Rusym1TMat(diffsq,I,m,n,rho,sigcinv,lambdazinv); % m x n psample - sample cov matrix
U = Gs\C;
Cproject = U'*U;
Ln = 1/r*eye(n) + Cproject;
% nut = Cs + r*(C*C') + Jm; % Ln^-1 = (1/J^2)*eye(n) - (1/J^4)*C'[Cs + (1/J^2)*C*C']^-1*C from Woodbury
% % % if(cond(nut) < 1.0e-09 || cond(nut) > 1.0e+09)
% % %     nut = nut + Jsup^2*eye(m);
% % % end;
% V = C'*(nut\C);
% Lninv = r*eye(n) - r^2*(V); % nut is not symmetric positive definite
% term = 0.5*y'*Lninv*y;