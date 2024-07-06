function Ln = quadTerm(Cs,C,Cproject,Js,Jsup,m,n,y,r)
% Compute quadratic 'term' in likelihood and covariance matrix (with z
% marginalized out) for posterior computation on error precision, r.
% This function differs from varGP in that it doesn't re-compute the 
% covariance matrices, C and Cs, but just uses the old values to save
% computation.
% Jm = Js^2*eye(m);
% nut = Cs + r*(C*C') + Jm; % Ln^-1 = (1/J^2)*eye(n) - (1/J^4)*C'[Cs + (1/J^2)*C*C']^-1*C from Woodbury
% V = C'*(nut\C);
% Lninv = r*eye(n) - r^2*(V); % nut is not symmetric positive definite
% term = 0.5*(y')*Lninv*(y);
Ln = 1/r*eye(n) + Cproject +Jsup^2*eye(n);
