function [p Cs C Cproject]= fpbeta(parms,type,diffsq,diffsqs,I,Ismall,m,n,y,X,alpha_0,taus,gamma,CJ,i,ar,br,az,bz,ac,bc,zeta,arho,brho,Cs,C,Cproject,position,lastRhoPosition)
% Returns parameter log-posterior densities for Metropolis-Hastings moves
% in a Gibbs step
% Input:
%   parms : vector of all model parameters constructed in metGibbs.
%           [rho, sigc, r, lambdaz, gamma]
%   type:  indicates which parameter among above group to be sampled
%   diffsq, diffsqs : Covariate data in squared difference form,
%                (X1(i,:)-X2(j,:)).^2, to build GP covariance matrix
%                 to construct likelihood.  Posteriors are proportional to
%                 likelihood x prior.  'diffsq' contains only n(n-1)/2
%                 unique values.  'diffsqs' is of length m(m-1)/2 from
%                 smaller set of m randomly selected points used for
%                 faster inverstion
%   I, Ismall: Indices for diffsq and diffsqs to reinflate from unique
%              to full set of values 
%   m, n    : see above.
%   y       : Response (nsample,1)
%   ar,br,az,bz,ac,bc: (shape,rate) hyperparameters on gamma priors
%                       for precision parameters
%  zeta:   mean of bernoulli prior on gamma_k
%  (arho,brho): hyperparameters for beta prior on rho_k
%  Cs,C, Cproject : Previously built values of covariance matrices, m x m, 
%        Cs (used  for inversion) and m x n, C, from which the
%        n x n GP covariance Cproject may be constructed.
%        These old values are able to be used
%        for sampling error precion, 'r'.  They are modified by varGPpartial.m
%        for the gibbs scan on (gamma_k, rho_k) to avoid full re-build.
%  position:  which 'k' when sampling (gamma_k, rho_k)
%  lastRhoPosition: previously sampled value for rho_k used by 
%              varGPpartial to build new Cs, C by adjusting old ones.
% Output:
%   p       : value of conditional posterior for sampled parameter
%   Cs, C   : new covariance matrices
%   Savitsky, Vannucci (2008 - 2011)


% Calculating the covariance matrix, R

% Definining size indices and (rho,nu,lambdaz,r) to support calculation of R
p = size(X,2);
k = size(diffsq,2);
rho = parms(1:k);
sigc = parms(k+1);
sigcinv = 1/sigc;
r = parms(k+2);
lambdaz = parms(k+3);
lambdazinv = 1/lambdaz;
gammak = parms(k+4:2*k+3);

Js = 0.065; %Jitter on m x m projection covariance matrix
Jsup = 0;

switch type
    case 1 %joint posterior for rho & gamma
        [Ln Cs C] = varGP(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r);
        CJ_temp = CJ;
        CJ_temp(:,:,i) = Ln;
        sum_term = zeros(n);
        for sum_i = 1:p
            sum_term = sum_term + gamma(sum_i) * diag(X(:,sum_i)) * CJ_temp(:,:,sum_i) * diag(X(:,sum_i));
        end
        Sig = sum_term + taus * eye(n);
        Sig_inv = Siginv(taus, sum_term, Js);
        like_all = -0.5*log(det(Sig)) - 0.5 * y' * Sig_inv * y;
        kgam = length(gammak(gammak==1));
        piGamRho = kgam*log(zeta) + (k-kgam)*log((1-zeta)) + sum((arho-1)*log(rho(rho<1))+(brho-1)*log(1-rho(rho<1)));
        p = like_all + piGamRho;
    case 2 %posterior for r
        if r < 1 || r > 130
            p = -99e99;
        else
            % [term Ln] = quadTerm(Cs,C,Cproject,Js,Jsup,m,n,bi,r);
            CJ_temp = CJ;
            CJ_temp(:,:,i) = 1/r*eye(n) + Cproject;
            sum_term = zeros(n);
            for sum_i = 1:p
                sum_term = sum_term + gamma(sum_i) * diag(X(:,sum_i)) * CJ_temp(:,:,sum_i) * diag(X(:,sum_i));
            end
            Sig = sum_term + taus * eye(n);
            Sig_inv = Siginv(taus, sum_term, Js);
            like_all = -0.5*log(det(Sig)) - 0.5 * y' * Sig_inv * y;
            p= like_all + (ar-1)*log(r) -(br*r);
        end
    case 3 %posterior for lambdaz
        if lambdaz < .001 || lambdaz > 125
            Cs = eye(m);
            C = eye(n);
            Cproject = eye(n);
            % [Cs C Cproject] = covComp(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r);
            p = -99e99;
        else
            [Ln Cs C Cproject] = varGP(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r);
            CJ_temp = CJ;
            CJ_temp(:,:,i) = Ln;
            sum_term = zeros(n);
            for sum_i = 1:p
                sum_term = sum_term + gamma(sum_i) * diag(X(:,sum_i)) * CJ_temp(:,:,sum_i) * diag(X(:,sum_i));
            end
            Sig = sum_term + taus * eye(n);
            Sig_inv = Siginv(taus, sum_term, Js);
            like_all = -0.5*log(det(Sig)) - 0.5 * y' * Sig_inv * y;
            p = like_all + (az-1)*log(lambdaz) -(bz*lambdaz);
        end
    case 4 %posterior for sigc
        if sigc < .001 || sigc > 125
            Cs = eye(m);
            C = eye(n);
            Cproject = eye(n);
            % [Cs C Cproject] = covComp(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r);
            p = -99e99;
        else
            Ln = varGP(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r);
            CJ_temp = CJ;
            CJ_temp(:,:,i) = Ln;
            sum_term = zeros(n);
            for sum_i = 1:p
                sum_term = sum_term + gamma(sum_i) * diag(X(:,sum_i)) * CJ_temp(:,:,sum_i) * diag(X(:,sum_i));
            end
            Sig = sum_term + taus * eye(n);
            Sig_inv = Siginv(taus, sum_term, Js);
            like_all = -0.5*log(det(Sig)) - 0.5 * y' * Sig_inv * y;
            p= like_all + (ac-1)*log(sigc) -(bc*sigc); 
        end
    case 5 %joint posterior for rho & gamma - changing one position, j, at a time
        [Ln,~, ~] = varGPpartial(Cs,C,diffsq,diffsqs,I,Ismall,m,n,position,lastRhoPosition,rho(position),sigcinv,Js,r);
        CJ_temp = CJ;
        CJ_temp(:,:,i) = Ln;
        sum_term = zeros(n);
        for sum_i = 1:p
            sum_term = sum_term + gamma(sum_i) * diag(X(:,sum_i)) * CJ_temp(:,:,sum_i) * diag(X(:,sum_i));
        end
        Sig = sum_term + taus * eye(n);
        Sig_inv = Siginv(taus, sum_term, Js);
        like_all = -0.5*log(det(Sig)) - 0.5 * y' * Sig_inv * y;
        [Ln Cs C] = varGPpartial(Cs,C,diffsq,diffsqs,I,Ismall,m,n,position,lastRhoPosition,rho(position),sigcinv,Js,r);
        kgam = length(gammak(gammak==1));
        piGamRho = kgam*log(zeta) + (k-kgam)*log((1-zeta)) + sum((arho-1)*log(rho(rho<1))+(brho-1)*log(1-rho(rho<1)));
%         piGam = gamma*log(zeta)' + (1-gamma)*log((1-zeta))';
        p = like_all + piGamRho;
    case 7
        if lambdaz < .001 || lambdaz > 125 || sigc < .001 || sigc > 125 || r < 1 || r > 130
            % [Cs C Cproject] = covComp(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r);
            Cs = eye(m);
            C = eye(n);
            p = -99e99;
        else
            [Ln Cs C Cproject] = varGP(diffsq,diffsqs,I,Ismall,m,n,rho,sigcinv,lambdazinv,Js,r);
            CJ_temp = CJ;
            CJ_temp(:,:,i) = Ln;
            sum_term = zeros(n);
            for sum_i = 1:p
                sum_term = sum_term + gamma(sum_i) * diag(X(:,sum_i)) * CJ_temp(:,:,sum_i) * diag(X(:,sum_i));
            end
            Sig = sum_term + taus * eye(n);
            Sig_inv = Siginv(taus, sum_term, Js);
            like = -0.5*log(det(Sig)) - 0.5 * y' * Sig_inv * y;
            pr= (ar-1)*log(r) -(br*r);
            kgam = length(gammak(gammak==1));
            piGamRho = kgam*log(zeta) + (k-kgam)*log((1-zeta)) + sum((arho-1)*log(rho(rho<1))+(brho-1)*log(1-rho(rho<1)));
            plambdaz = (az-1)*log(lambdaz) -(bz*lambdaz);
            psigc = (ac-1)*log(sigc) -(bc*sigc);
            p = like + piGamRho + plambdaz + psigc + pr; 
        end
end