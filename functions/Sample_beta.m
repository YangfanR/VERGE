function [Beta_j] = Sample_beta(index, X, y, alpha_0, Beta, CJ, gamma, taus)

n = size(X,1);
Sigma_inv = 1/taus * eye(n);
% Beta_j = zeros(n, 1);
Js = 0.065;
i = index;
X_j = diag(X(:,i));
if(ndims(CJ) == 3)
    C_j = CJ(:,:,i);
else
    C_j = CJ;
end
if(gamma(i) == 0)
    Sigma_j = C_j;
    Mu_j = zeros(n,1);
    Beta_j = mvnrnd(Mu_j, Sigma_j);
else  
    % Sigma_j = inv(inv(C_j) + gamma_j' * X_j' * Sigma_inv * X_j * gamma_j);
    % Mu_j = Sigma_j * gamma_j' * X_j' * Sigma_inv * (y - alpha_0 - sum(X.*Beta.*gamma,2) + X(:,i).* (gamma(i) * Beta(:,i)));
    Sigma_j = inv(inv(C_j) + X_j' * Sigma_inv * X_j); 
    if(any(eig(Sigma_j)<0))
        % Sigma_j = Sigma_j + Js^2*eye(n);
        Beta_j = zeros(n,1);
    else
        Mu_j = Sigma_j * Sigma_inv * X_j' * (y - alpha_0 - sum(X.*Beta.*gamma,2) + X_j * Beta(:,i));
        Beta_j = mvnrnd(Mu_j, Sigma_j);
    end
end
end
