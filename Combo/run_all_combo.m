function [] = run_all(seed)

rng(seed + 100, "twister");

y = readmatrix('Y.csv');
y = normalize(y, "center");

X = readmatrix('X.csv');
X = X(2:end, :);
X = normalize(X);

Z = readmatrix('Z.csv');
Z = Z(:,[2, 28, 30]);
Z(Z == 2,1) = -1;

n = size(y,1);
P = size(X,2);
K = size(Z,2);

% Priors on predictor selection
a = log(0.21);
b = 0.5;

% Priors on covariate selection
zeta_p = 0.5;

nmc = 100000;
burnin = 100000;

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

h = 50^2;
v0 = 0.1^2;

MCMC_all(seed, X, y, Z, az, bz, ac, bc, ar, br, arho, brho, a_0, b_0, h, v0, a, b, zeta_p, nmc, burnin);

