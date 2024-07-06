n = 200;
P = 60;
P_t = 6;
K = 3;
a = log(0.1);
zeta_p = 0.5;
n_t = 50;
nmc = 30000;
Jit = .01;

results = "results_n" + n + "_P" + P + "_K" + K + "_a" + a + "_pi"+ zeta_p + "_" + iteration + ".mat";
data = "data_n" + n + "_P" + P + "_K" + K + "_a" + log(0.1) + "_pi"+ 0.5 + "_" + iteration + ".mat";
testfile = "data_test_n" + n + "_P" + P + "_K" + K + "_a" + log(0.1) + "_pi"+ 0.5 + "_" + iteration + ".mat";
load(testfile);
load(data);
load(results);

Z = normalize(Z, "range");
Z_t = normalize(Z_t,"range");

keep_pred = find(PPI_pred > 0.5);
n_pred = length(keep_pred);
Z_list = cell(n_pred,1);

PPI_cov_temp = zeros(P, K);
for i = 1:size(PPI_cov,1)
    PPI_cov_sorted = sort(PPI_cov(i,:), 'descend');
    fdr_all = cumsum(1-PPI_cov_sorted)./(1:K);
    t = PPI_cov_sorted(find(fdr_all <= 0.1,1, 'last'));
    if(isempty(t))
        t = 1.1;
    end
    PPI_cov(i,PPI_cov(i,:) >= t) = 1;
    PPI_cov(i, PPI_cov(i,:) ~= 1) = 0; 
end
PPI_cov_temp(keep_pred,:) = PPI_cov;
for i = 1:n_pred
    keep = find(PPI_cov(i,:) == 1);
    Ztrain = Z(:,keep);
    Ztest = Z_t(:,keep);
    [difftesq Ite numte numtr] = difference(Ztest,Ztrain);
    [difftrsq Itr] = difference(Ztrain,Ztrain);
    Z_list{i} = {difftesq, Ite, numte, numtr, difftrsq, Itr};
end

Beta = Beta_s(:, keep_pred);

inds1 = find(gamma_save(keep_pred,1:nmc) == 1);
[ind_cov1, ind_cov2] = find(PPI_cov_temp == 1);
combined_condition = true(nmc, 1);
for k = 1:length(ind_cov1)
    current_condition = squeeze(gammak_save(ind_cov1(k), ind_cov2(k), 1:nmc)) == 1;
    combined_condition = combined_condition & current_condition;
end

inds2 = find(combined_condition);
inds = intersect(inds1, inds2);
ypred = 0;
for ind = 1:length(inds)
    i = inds(ind);
    
    gamma = gamma_save(:,i);
    gammak = gammak_save(:,:,i);
    r = r_save(:,i);
    lambdaz = lambdaz_save(:,i);
    rho = rho_save(:,:,i);
    sigc = sigc_save(:,i);
    taus = taus_save(i);
    Beta_t = zeros(n_t, n_pred);
    for p_i = 1:n_pred
        j = keep_pred(p_i);
        r_j = r(j);
        lambdazinv_j = 1/lambdaz(j);
        sigcinv_j = 1/sigc(j);
        rho_j = rho(j,PPI_cov_temp(j,:) == 1);
        difftesq = Z_list{p_i}{1};
        Ite = Z_list{p_i}{2};
        numte = Z_list{p_i}{3};
        numtr = Z_list{p_i}{4};
        difftrsq = Z_list{p_i}{5};
        Itr = Z_list{p_i}{6};
        Ltetr = Rusym1TMat(difftesq,Ite,numte,numtr,rho_j,sigcinv_j,lambdazinv_j);
        Ltrtr = (1/r_j+Jit)*eye(numtr)+Rusym1TMat(difftrsq,Itr,numtr,numtr,rho_j,sigcinv_j,lambdazinv_j);
        Gtrtr = chol(Ltrtr,'lower');
        Utr = Gtrtr\Ltetr';
        vtr = Gtrtr\Beta(:,p_i);    
        Beta_t(:,p_i) = Utr'*vtr;
    end
    ypred = ypred + sum(Beta_t .* X_t(:,keep_pred),2);
end

ypred_m = ypred/length(inds);
pmse = mean((y_t - ypred_m).^2);
fprintf('pmse = %d\n', pmse);


