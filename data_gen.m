function[X, Z, y, Beta_p, G, Zi] = data_gen(seed, n, K, P, P_t, a, zeta_p)

% Generate data 

% Sample covariates
Z = -1 + 2.*rand(n,K);

% Generate graph
num_tfs = 10;
g_per_tf = (P/num_tfs) - 1;

% Set up Sigma matrix
Sigma = eye(P);
tf_indices = 1:(g_per_tf + 1):P;

for tf = 1:num_tfs
    for gene = 1:g_per_tf
        % Correlation of primary node to subsidiary nodes is 0.7
        Sigma(tf_indices(tf), tf_indices(tf) + gene) = 0.7;
        Sigma(tf_indices(tf) + gene, tf_indices(tf)) = 0.7;
        
        % Correlation among subsidiary nodes are 0.7^2
        for gene2 = (gene + 1):g_per_tf
            Sigma(tf_indices(tf) + gene, tf_indices(tf) + gene2) = 0.7^2;
            Sigma(tf_indices(tf) + gene2, tf_indices(tf) + gene) = 0.7^2;
        end
    end
end

% True graph structure
Omega_true = inv(Sigma);
G = abs(Omega_true) > 0.001;
X = mvnrnd(zeros(1,P),Sigma,n);
X = X - mean(X,1);

% Set up true coefficients
Beta_p = zeros(n,P_t);
Zi = randi([1,K], 1, P_t-1);

Beta_p(:,1) = 0.3;
Beta_p(:,2) = 2 * sin(pi * Z(:,Zi(1)));
Beta_p(:,3) = 2 * Z(:,Zi(2)).^2 - 1;
Beta_p(:,4) = -2 * Z(:,Zi(3));
Beta_p(:,5) = 2 * cos(pi * Z(:,Zi(4)));
Beta_p(:,6) = -2 * exp(-((Z(:,Zi(5)) - 0.3).^2) / (2 * 0.3^2)) - 3 * exp(-((Z(:,Zi(5)) + 0.5).^2) / (2 * 0.3^2));

Beta_t = zeros(n, P);
Beta_t(:,1:P_t) = Beta_p;

% True regression
gamma_t = [ones(1, P_t), zeros(1, P - P_t)];

taus_t = 1;
alpha_0_t = 0;
epsilon = mvnrnd(zeros(n,1), diag(taus_t));
y = alpha_0_t + sum(X.*(Beta_t.*gamma_t),2) + epsilon;

datafile = "data_n" + n + "_P" + P + "_K" + K + "_a" + a + "_pi"+ zeta_p + "_" + seed + ".mat";
save(datafile,"X","Z","y", "Beta_p", "G", "Zi");
