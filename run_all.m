function [] = run_all(seed)

rng(seed + 1000, "twister");

n = 200;
P = 60;
P_t = 6; %number of true predictors
K = 3;

% Priors on predictor selection
a = log(0.1);
b = 0.5;

% Priors on covariate selection
zeta_p = 0.5;

nmc = 10;
burnin = 10;

% hyperparameters
az = 1; %0.5; %shape hyperparameter for exponential weight (lambda_z)
bz = 1; %0.5; %scale hyperparameter for expential weight (lambda_z) 
ac = 1; %0.1; %shape hyperparameter for covariance intercep (lambda_a)
bc = 1; %0.1; %scale hyperparameter for covariance intercept (lambda_a)
ar = 1; %shape hyperparameter for error precision (r)
br = .01; %scale hyperparameter for error precision (r) 
arho = 1.0; % hyperparmeter for beta prior on rho_{k}
brho = 1.0; % hyperparmeter for beta prior on rho_{k}
a_0 = 1; % shape for gamma prior on taus
b_0 = 1; % scale for gamma prior on taus

h = 100^2;
v0 = 0.05^2;

[X, Z, y, Beta_p, G, Zi] = data_gen(seed, n, K, P, P_t, a, zeta_p);
MCMC_all(seed, X, y, Z, G, az, bz, ac, bc, ar, br, arho, brho, a_0, b_0, h, v0, a, b, zeta_p, nmc, burnin);

