function[] = MCMC_all(seed, X, y, Z, G, az, bz, ac, bc, ar, br, arho, brho, a_0, b_0, h, v0, a, b, zeta_p, nmc, burnin)
% MCMC algorithm
% INPUT:
%    X: predictors (n * P)
%    y: response (n * 1)
%    Z: covariates (n * K)
%    ar,br,az,bz,ac,bc: (shape,rate) hyperparameters on gamma priors for precision parameters
%    (arho,brho): hyperparameters for beta prior on rho_k
%    (a_0, b_0): hyperparameters for gamma priors on taus
%    h, v_0 (v_1 = h*v_0): variance on continuous spike-and-slab for graph selection
%    a, b: scalar hyperparameters for MRF prior
%    zeta_p: mean of bernoulli prior on gamma_k
%    nmc, burin: number of sampling iterations and number of burnin
%
% Saved file:
%    gamma_save: posterior samples for gamma
%    rho_save: posterior samples for rho
%    sigc_save: posterior samples for lambda_a
%    lambdaz_save: posterior samples for lambda_z
%    r_save: posterior samples for r
%    gammak_save: posterior samples for gamma_k
%    taus_save: posterior samples for taus
%    Beta_s: mean of sampled Beta
%    adj_avg: mean of estimated graph
%    PPI_pred: PPIs for predictor selection
%    PPI_cov: PPIs for covariate selection
%    elpasedTime: total running time


tic
n = size(X, 1);
P = size(X, 2);
K = size(Z, 2);

v1 = h * v0;
lambda_g = 1;
pii = 2/(P-1);
summary_only = 0;

psample = 1;

Js = 0.065;

% Standardize data for GP
Z = normalize(Z,"range");
% Z = normalize(Z);
[Zs indz] = project(Z,psample);

[diffsq J m n] = difference(Zs,Z);
[diffsqs Jsmall] = difference(Zs,Zs);

lastgammak_all = zeros(P,K);
lastrho_all = ones(P,K);
lastsigc_all =ones(1,P) * (ac/bc);
lastlambdaz_all = ones(1,P) * (az/bz);
lastr_all = ones(1,P) * (ar/br);

% Standardize data so that S(i,i) = n 
S = X' * X;
[n, P] = size(X);
S = corrcoef(S) * n;

% Initial guess for Sigma, precision matrix, and adjacency matrix
Sig = eye(P);
Omega = inv(Sig);
adj = eye(P);

V0 = v0 * ones(P);
V1 = v1 * ones(P);

tau = V1;
ind_noi_all = zeros(P-1, P);

for i = 1:P
    if i==1
        ind_noi = [2:P]';
    elseif i==P
        ind_noi = [1:P-1]';
    else
        ind_noi = [1:i-1,i+1:P]';
    end
    ind_noi_all(:,i) = ind_noi;
    
end

pii_RB = zeros(P);
pii_mat = zeros(P);

taus = 1;
alpha_0 = 0;
gamma = zeros(1,P);

CJ = zeros(n,n,P);
Beta = zeros(n,P);
for j = 1:P
    CJ(:,:,j) = ones(n) + Js^2*eye(n);
    Beta(:,j) = mvnrnd(zeros(1,n), CJ(:,:,j));
end

% Always keep variable selections
alpha_save = zeros(1,nmc);
taus_save = zeros(1,nmc);
gamma_save = zeros(P, nmc);
Beta_save = zeros(n,P,nmc);
rho_save =zeros(P,K,nmc);
gammak_save =zeros(P,K,nmc);
zeta_save =zeros(P,K,nmc);
sigc_save = zeros(P,nmc);
lambdaz_save = zeros(P,nmc);
r_save = zeros(P,nmc);
gammak_temp = zeros(P,K,nmc + burnin);

% Record some diagnostic info
full_gamma_save = zeros(P, burnin + nmc);
node_degrees = zeros(P, burnin + nmc);

% Keep track of info to compute acceptance rates
n_gamma_prop = 0;
n_gamma_accept = 0;
n_add_prop = 0;
n_add_accept = 0;
n_remove_prop = 0;
n_remove_accept = 0;
n_keep_prop = 0;
n_keep_accept = 0;

% Measuring the acceptance rates for rho and nu in our Metropolis steps
% when gamma_k = 1
proposeRhoBetween = zeros(1,P);
proposeRhoWithin = zeros(1,P);
countRhoBetween = zeros(1,P);
countRhoWithin = zeros(1,P);
countgammak = zeros(1,P);

% 'D' is tuning parameter for number of times to re-sample rho_{k}
% for 'within' step when gamma_{k} = 1, meaning x_{k} is in the current
% model.
D = 1;
% 'w' is a tuning parameter to generate multiple candidates for rho_{k}
% from which one is randomly selected to improve mixing.
wr = 20;

% Tuning parameters for updating, zeta (p x 1), the mean for the bernoulli
% proposal for gamma_{k}|gamma_{-k}, ...
lambda = .1;
zetaProposal_all = 0.5*ones(P,K);
updateCount = ones(1,P);
lastUpdateNum = ones(1,P);
updateInterval = 20;
firstUpdate = updateInterval;
stepIncrement = 0.3;
iter_add = ones(1,P);

adj_save = zeros(P, P, nmc);


% MCMC sampling
for iter = 1: burnin + nmc

    % Number of currently included variables
    p_gamma = sum(gamma);

    % UPDATE alpha_0
    % sig_alpha = 1/(n/taus + 1/(h_0*taus));
    % mu_alpha = sum(y - sum(X.*Beta.*gamma ,2),1)/taus*sig_alpha;
    % alpha_0 = normrnd(mu_alpha, sqrt(sig_alpha));
    alpha_0 = 0;

    % UPDATE taus
    sum_term = zeros(n);
    for sum_i = 1:P
        sum_term = sum_term + gamma(sum_i) * diag(X(:,sum_i)) * CJ(:,:,sum_i) * diag(X(:,sum_i));
    end
    tune = 1;
    candtaus = gamrnd(2,taus/2);
    prior_ratio = -(a_0 + 1)*log(candtaus) - b_0/candtaus + ((a_0 + 1)*log(taus) + b_0/taus);
    Sig_y = sum_term + taus * eye(n);
    Sig_y_inv = Siginv(taus, sum_term, Js);
    Sig_y_cand = sum_term + candtaus * eye(n);
    Sig_y_cand_inv = Siginv(candtaus, sum_term, Js);
    likelihood_ratio = - 0.5*log(det(Sig_y_cand)) - 0.5 * y' * Sig_y_cand_inv * y + 0.5*log(det(Sig_y)) + 0.5 * y' * Sig_y_inv * y;
    alpha = prior_ratio + likelihood_ratio - (tune+1)*log(candtaus) - tune*taus/(candtaus) + (tune+1)*log(taus) - tune*candtaus/(taus);
    if (log(rand(1)) < alpha && candtaus > 0.025)
        taus = candtaus;
    end

    % UPDATE gamma
    % Select an entry at random to toggle
    change_index = randsample(P, 1);
    % fprintf('Change index = %d\n', change_index);
    
    n_gamma_prop = n_gamma_prop + 1;
    
    % Toggle value
    gamma_prop = gamma;
    gamma_prop(change_index) = abs(gamma(change_index) - 1);
    if(gamma_prop(change_index) == 0)
        rp = binornd(1,0.5);
        gamma_prop(change_index) = rp;
    end

    % Keep track of number of add vs. remove moves proposed
    if (gamma(change_index) == 0) 
        n_add_prop = n_add_prop + 1;
    elseif(gamma_prop(change_index) == 0)
        n_remove_prop = n_remove_prop + 1;
    else
        n_keep_prop = n_keep_prop + 1;
    end
    
    i = change_index;
    lastrho = lastrho_all(i,:);
    lastsigc = lastsigc_all(i);
    lastlambdaz = lastlambdaz_all(i);
    lastgammak = lastgammak_all(i,:);
    zetaProposal = zetaProposal_all(i,:);
    lastr = lastr_all(i);
    last = [lastrho lastsigc lastr lastlambdaz lastgammak];
    [flast CsOld COld] = fpbeta(last,7,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,ar,br,az,bz,ac,bc,zeta_p,arho,brho);
    flast = flast + sum(gamma) * a + b * (gamma * (adj - eye(P)) * gamma');

    if (gamma_prop(change_index) == 1)
        j = randsample(K, 1);
        candgammak = lastgammak;
        candrho = lastrho;
        component = my_bernoulli(lambda);
        candgammak(j) = (1-component)*my_bernoulli(zetaProposal(j)) + component*my_bernoulli(0.5);
        candsigc = gamrnd(tune,lastsigc/tune);
        candlambdaz = gamrnd(tune,lastlambdaz/tune);
        candr = gamrnd(tune,lastr/tune);
      
        if(candgammak(j)==1)
            proposeRhoBetween(i) = proposeRhoBetween(i) + 1;
            candrho(j) = rand(1);
            cand = [candrho candsigc candr candlambdaz candgammak];
            [fcand, CsNew, CNew] = fpbeta(cand,7,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma_prop,CJ,i,ar,br,az,bz,ac,bc,zeta_p,arho,brho);
            fcand = fcand + sum(gamma_prop) * a + b * (gamma_prop * (adj - eye(P)) * gamma_prop');
            qlastnew = (tune-1)*log(candlambdaz) - tune*candlambdaz/(lastlambdaz) + (tune-1)*log(candr) - tune*candr/(lastr) + (tune-1)*log(candsigc) - tune*candsigc/(lastsigc);
            qnewlast = (tune-1)*log(lastlambdaz) - tune*lastlambdaz/(candlambdaz) + (tune-1)*log(lastr) - tune*lastr/(candr) + (tune-1)*log(lastsigc) - tune*lastsigc/(candsigc);
            alpha = fcand + qnewlast  - flast - qlastnew;
            if (log(rand(1)) < alpha)
                lastgammak(j) = candgammak(j);
                lastrho(j) = candrho(j);
                lastsigc = candsigc;
                lastlambdaz = candlambdaz;
                lastr = candr;
                CsOld = CsNew; COld = CNew;
                countgammak(i) = countgammak(i) + 1;
                countRhoBetween(i) = countRhoBetween(i) + 1;
                Gs = chol(CsOld,'lower');
                U = Gs\COld;
                Cproject = U'*U;
                % CJ(:,:,i) = Cproject + Js^2*eye(n);
                CJ(:,:,i) = Cproject + 1/lastr*eye(n);
                if gamma(change_index) == gamma_prop(change_index)
                    n_keep_accept = n_keep_accept + 1;
                else
                    n_add_accept = n_add_accept + 1;
                    p_gamma = p_gamma + 1;
                end
                gamma = gamma_prop;
                n_gamma_accept = n_gamma_accept + 1;
            end
        % elseif(lastgammak(j) ~= 0) 
        elseif(candgammak(j)==0)
            candrho(j) = 1;
            cand = [candrho candsigc candr candlambdaz candgammak];
            [fcand, CsNew, CNew] = fpbeta(cand,7,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma_prop,CJ,i,ar,br,az,bz,ac,bc,zeta_p,arho,brho);
            fcand = fcand + sum(gamma_prop) * a + b * (gamma_prop * (adj - eye(P)) * gamma_prop');
            qlastnew = (tune-1)*log(candlambdaz) - tune*candlambdaz/(lastlambdaz) + (tune-1)*log(candr) - tune*candr/(lastr) + (tune-1)*log(candsigc) - tune*candsigc/(lastsigc);
            qnewlast = (tune-1)*log(lastlambdaz) - tune*lastlambdaz/(candlambdaz) + (tune-1)*log(lastr) - tune*lastr/(candr) + (tune-1)*log(lastsigc) - tune*lastsigc/(candsigc);
            alpha = fcand + qnewlast  - flast - qlastnew;
            if (log(rand(1)) < alpha)
                lastgammak(j) = candgammak(j);
                lastrho(j) = candrho(j);
                lastsigc = candsigc;
                lastlambdaz = candlambdaz;
                lastr = candr;
                CsOld = CsNew; COld = CNew;
                countgammak(i) = countgammak(i) + 1;
                Gs = chol(CsOld,'lower');
                U = Gs\COld;
                Cproject = U'*U;
                % CJ(:,:,i) = Cproject + Js^2*eye(n);
                CJ(:,:,i) = Cproject + 1/lastr*eye(n);
                if gamma(change_index) == gamma_prop(change_index)
                    n_keep_accept = n_keep_accept + 1;
                else
                    n_add_accept = n_add_accept + 1;
                    p_gamma = p_gamma + 1;
                end
                gamma = gamma_prop;
            end
        end
        if(gamma(change_index)==1)
            last = [lastrho lastsigc lastr lastlambdaz lastgammak];
            [fGamRholast CsOld COld] = fpbeta(last,1,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],[],[],[],[],zeta_p,arho,brho);
            for j = 1:K
                candgammak = lastgammak;
                candrho = lastrho;
                component = my_bernoulli(lambda);
                candgammak(j) = (1-component)*my_bernoulli(zetaProposal(j)) + component*my_bernoulli(0.5);
                if(candgammak(j)==1)
                    proposeRhoBetween(i) = proposeRhoBetween(i) + 1;
                    candrho(j) = rand(1);
                    cand = [candrho lastsigc lastr lastlambdaz candgammak];
                    [fGamRhocand CsNew CNew] = fpbeta(cand,5,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],[],[],[],[],zeta_p,arho,brho,CsOld,COld,[],j,lastrho(j));
                    alpha = fGamRhocand  - fGamRholast;
                    if (log(rand(1)) < alpha)
                        lastgammak(j) = candgammak(j);
                        lastrho(j) = candrho(j);
                        fGamRholast = fGamRhocand; 
                        CsOld = CsNew; COld = CNew;
                        countgammak(i) = countgammak(i) + 1;
                        countRhoBetween(i) = countRhoBetween(i) + 1;
                    end
                elseif(lastgammak(j) ~= 0) % candgammak(j)==0
                    candrho(j) = 1;
                    cand = [candrho lastsigc lastr lastlambdaz candgammak];
                    [fGamRhocand CsNew CNew] = fpbeta(cand,5,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],[],[],[],[],zeta_p,arho,brho,CsOld,COld,[],j,lastrho(j));
                    alpha = fGamRhocand  - fGamRholast;
                    if (log(rand(1)) < alpha)
                        lastgammak(j) = candgammak(j);
                        lastrho(j) = candrho(j);
                        fGamRholast = fGamRhocand; 
                        CsOld = CsNew; COld = CNew;
                        countgammak(i) = countgammak(i) + 1;
                    end 
                end
                % Within-model move resampling of rho_{k}|gamma_    
                if(lastgammak(j)==1)
                    for d = 1:D
                    proposeRhoWithin(i) = proposeRhoWithin(i) + 1;
                    rhoProposal = rand([1,wr]);
                    candrho(j) = randsample(rhoProposal,1);
                      candrho(j) = rand(1);
                    cand = [candrho lastsigc lastr lastlambdaz lastgammak];
                    [fGamRhocand CsNew CNew]= fpbeta(cand,5,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],[],[],[],[],zeta_p,arho,brho,CsOld,COld,[],j,lastrho(j));
                    alpha = fGamRhocand - fGamRholast;
                        if (log(rand(1)) < alpha)
                                lastrho(j) = candrho(j);
                                fGamRholast = fGamRhocand; 
                                CsOld = CsNew; COld = CNew;
                                countRhoWithin(i) = countRhoWithin(i) + 1;
                        end
                    end
                end 
            end 
            lastgammak_all(i,:) = lastgammak;
            lastrho_all(i,:) = lastrho;  

            %updating sigc from full conditional 
            last = [lastrho lastsigc lastr lastlambdaz lastgammak]; 
            fsigclast = fpbeta(last,4,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],[],[],ac,bc,[]);
            tune = 1;
            candsigc = gamrnd(tune,lastsigc/tune);
            cand = [lastrho candsigc lastr lastlambdaz lastgammak]; 
            fsigccand = fpbeta(cand,4,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],[],[],ac,bc,[]);
            qlastnew = (tune-1)*log(candsigc) - tune*candsigc/(lastsigc);
            qnewlast = (tune-1)*log(lastsigc) - tune*lastsigc/(candsigc);
            alpha = (fsigccand+qnewlast-fsigclast-qlastnew); 
            if (log(rand(1)) < alpha)
                lastsigc = candsigc;
            end
            lastsigc_all(i) = lastsigc;

            %updating lambdaZ from full conditional  
            last = [lastrho lastsigc lastr lastlambdaz lastgammak]; 
            [flambdazlast Cs C Cproject] = fpbeta(last,3,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],az,bz,[],[],[]);
            tune = 1;
            candlambdaz = gamrnd(tune,lastlambdaz/tune);
            cand = [lastrho lastsigc lastr candlambdaz lastgammak]; 
            [flambdazcand Csc Cc Cprojectc] = fpbeta(cand,3,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],az,bz,[],[],[]);
            qlastnew = (tune-1)*log(candlambdaz) - tune*candlambdaz/(lastlambdaz);
            qnewlast = (tune-1)*log(lastlambdaz) - tune*lastlambdaz/(candlambdaz);
            alpha = (flambdazcand+qnewlast-flambdazlast-qlastnew); 
            if (log(rand(1)) < alpha) 
                lastlambdaz = candlambdaz;
                Cs = Csc; C = Cc; Cproject = Cprojectc;
            end
            lastlambdaz_all(i) = lastlambdaz;

            %updating r from full conditional 
            last = [lastrho lastsigc lastr lastlambdaz lastgammak]; 
            frlast = fpbeta(last,2,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,ar,br,[],[],[],[],[],[],[],Cs,C,Cproject);
            tune = 1;
            candr = gamrnd(tune,lastr/tune);
            cand = [lastrho lastsigc candr lastlambdaz lastgammak]; 
            frcand = fpbeta(cand,2,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,ar,br,[],[],[],[],[],[],[],Cs,C,Cproject);
            qlastnew = (tune-1)*log(candr) - tune*candr/(lastr);
            qnewlast = (tune-1)*log(lastr) - tune*lastr/(candr);
            alpha = (frcand+qnewlast-frlast-qlastnew); 
            if (log(rand(1)) < alpha)
                lastr = candr;
            end

            last = [lastrho lastsigc lastr lastlambdaz lastgammak];
            [fGamRholast CsOld COld] = fpbeta(last,1,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma,CJ,i,[],[],[],[],[],[],zeta_p,arho,brho);
            Gs = chol(CsOld,'lower');
            U = Gs\COld;
            Cproject = U'*U;
            CJ(:,:,i) = Cproject + 1/lastr*eye(n);
            Beta(:,i) = Sample_beta(i, X, y, alpha_0, Beta, CJ, gamma, taus);
        end
        gammak_temp(i,:,iter_add(i)) = lastgammak;
        testInd = (updateCount(i)==1);
        updateNum = updateInterval*updateCount(i);
        if(iter_add(i) == updateNum*(1-testInd) + firstUpdate*testInd)
            updateCount(i) = updateCount(i) + 1;
            sumGammak = sum(gammak_temp(i,:,lastUpdateNum(i):iter_add(i)),3);
            iterations = iter_add(i) - lastUpdateNum(i) + 1;
            ergAvgZeta = sumGammak/(iterations);
            step = stepIncrement/updateCount(i);
            zetaProposal = zetaProposal + step*(ergAvgZeta - zetaProposal);
            lastUpdateNum(i) = updateNum;
        end
        iter_add(i) = iter_add(i) + 1;
        zetaProposal_all(i,:) = zetaProposal;
        % disp(zetaProposal)
    else
        [fcand, CsNew, CNew] = fpbeta(last,7,diffsq,diffsqs,J,Jsmall,m,n,y,X,alpha_0,taus,gamma_prop,CJ,i,ar,br,az,bz,ac,bc,zeta_p,arho,brho);
        fcand = fcand + sum(gamma_prop) * a + b * (gamma_prop * (adj - eye(P)) * gamma_prop');
        alpha = fcand  - flast;
        
        if (log(rand(1)) < alpha)
            gamma(change_index) = 0;
            p_gamma = p_gamma - 1;
            n_remove_accept = n_remove_accept + 1;
            n_gamma_accept = n_gamma_accept + 1;
            lastgammak = zeros(1,K);
            lastrho = ones(1,K);
            lastsigc = gamrnd(ac,1/bc);
            lastlambdaz = gamrnd(az,1/bz);
            lastr = gamrnd(ar, 1/br);
        end
    end

    lastgammak_all(i,:) = lastgammak;
    lastrho_all(i,:) = lastrho;
    lastsigc_all(i) = lastsigc;
    lastlambdaz_all(i) = lastlambdaz;
    lastr_all(i) =  lastr;
    zetaProposal_all(i,:) = zetaProposal;


    for i = 1:P
        ind_noi = ind_noi_all(:,i);
        tau_temp = tau(ind_noi,i);
        
        Sig11 = Sig(ind_noi,ind_noi);
        Sig12 = Sig(ind_noi,i);
        
        invC11 = Sig11 - Sig12*Sig12'/Sig(i,i);
        
        Ci = (S(i,i)+lambda_g)*invC11+diag(1./tau_temp);
        Ci = (Ci+Ci')./2;
        Ci_chol = chol(Ci);
        mu_i = -Ci_chol\(Ci_chol'\S(ind_noi,i));
        eta = mu_i+ Ci_chol\randn(P-1,1);

        Omega(ind_noi,i) = eta;
        Omega(i,ind_noi) = eta;
        
        a_gam = 0.5*n+1;
        b_gam = (S(i,i)+lambda_g)*0.5;
        gam = gamrnd(a_gam,1/b_gam);
        
        c = eta'*invC11*eta;
        Omega(i,i) = gam+c;
        
        %% Below updating Covariance matrix according to one-column change of precision matrix
        invC11eta = invC11*eta;
        
        Sig(ind_noi,ind_noi) = invC11+invC11eta*invC11eta'/gam;
        Sig12 = -invC11eta/gam;
        Sig(ind_noi,i) = Sig12;
        Sig(i,ind_noi) = Sig12';
        Sig(i,i) = 1/gam;  
        
        v0 = V0(ind_noi,i);
        v1 = V1(ind_noi,i);
        
        w1 = -0.5*log(v0) -0.5*eta.^2./v0+log(1-pii);
        w2 = -0.5*log(v1) -0.5*eta.^2./v1+log(pii);
        
        w_max = max([w1,w2],[],2);
        
        w = exp(w2-w_max)./sum(exp([w1,w2]-repmat(w_max,1,2)),2);
        
        z = (rand(P-1,1)<w);
        
        v = v0;
        v(z) = v1(z);
        
        pii_mat(ind_noi,i) = w;
        
        tau(ind_noi,i) = v;
        tau(i,ind_noi) = v;
        
        adj(ind_noi,i) = z;
        adj(i,ind_noi) = z; 
    end
    % adj = G + eye(P);
    
    if iter > burnin
        alpha_save(iter-burnin) = alpha_0;
        taus_save(iter-burnin) = taus;
        Beta_save(:,:,iter-burnin) = Beta;
        gamma_save(:, iter-burnin) = gamma;
        rho_save(:,:,iter-burnin) = lastrho_all;
        gammak_save(:,:,iter-burnin) = lastgammak_all;
        zeta_save(:,:,iter-burnin) = zetaProposal_all;
        sigc_save(:,iter-burnin) = lastsigc_all;
        lambdaz_save(:,iter-burnin) = lastlambdaz_all;
        r_save(:,iter-burnin) = lastr_all;
        
        pii_RB = pii_RB + pii_mat/nmc;

        adj_save(:, :, iter-burnin) = adj;
    end

    % Print out info every 100 iterations
    if mod(iter, 1) == 0
        fprintf('Iteration = %d\n', iter);
        fprintf('Number of included variables = %d\n', sum(gamma));
        fprintf('Number of add variable moves proposed %d and accepted %d\n', n_add_prop, n_add_accept);
        fprintf('Number of remove variable moves proposed %d and accepted %d\n', n_remove_prop, n_remove_accept);
        fprintf('Number of keep variable moves proposed %d and accepted %d\n', n_keep_prop, n_keep_accept);
        fprintf('Number of included edges %d \n\n', (sum(sum(adj)) - P) / 2);
    end

    % iter = iter + 1;
    
    full_gamma_save(:, iter) = gamma;
    node_degrees(:, iter) = sum(adj, 2) - 1;

end

PPI_pred = sum(gamma_save,2)/(nmc);
p_sel = sum(PPI_pred > 0.5);
index_sel = find(PPI_pred > 0.5);
PPI_cov = zeros(p_sel,K);

for j = 1:p_sel
    i = index_sel(j);
    PPI_cov(j,:) = sum(gammak_save(i,:,gamma_save(i,:) == 1),3)/sum(gamma_save(i,:) == 1);
end

adj_avg = sum(adj_save,3)/nmc;

Beta_s = zeros(n, P);
for j = 1:P
    Beta_s(:,j) = mean(Beta_save(:,j,gamma_save(j,:) == 1), 3);
end

elapsedTime = toc;

filename = "results_n" + n + "_P" + P + "_K" + K + "_a" + a + "_pi" + zeta_p + "_" + seed + ".mat";
save(filename, "gamma_save", "rho_save", "sigc_save", "lambdaz_save", "r_save", "gammak_save", "taus_save", "Beta_s", "adj_avg", "PPI_pred", "PPI_cov", "elapsedTime");
