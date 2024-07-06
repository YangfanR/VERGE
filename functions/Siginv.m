function Sig_inv = Siginv(taus, sum_term, Js)
    
n = size(sum_term,1);
Jm = Js^2*eye(n);
if(max(sum_term,[],"all") > 1.0e+11)
    Jm = (max(sum_term,[],"all")/1.0e+11) * eye(n);
end
Cs = sum_term + Jm;
Gs = chol(Cs,'lower');
n = size(Cs, 1);
U = Gs\eye(n);
Cinv = U'*U;
V = Cinv + 1/taus*eye(n);%  + Jm; % Ln^-1 = (1/J^2)*eye(n) - (1/J^4)*C'[Cs + (1/J^2)*C*C']^-1*C from Woodbury
Sig_inv = 1/taus*eye(n) - (1/taus)^2*inv(V); 