function [te Xtr Xte ytr yitr yite myi syi] = transform(type,Xi,yi,numte)
% Given a data matrix, Xi, response vector, yi, and number of observations
% to set aside for validation, this function standardizes X
% according to 1 of 4 options chosen by the user, conducts a z-score 
% transformation of the response, y, and then selects a random
% test set, returning training and test set data entities after
% standardization.
% Input:
%  type:  1 = centers the columns of X; 2 = transforms the columns of X so
%            that all values live in [0,1]; 3 = transforms columns of X to
%            gaussian cdf values (by first z-score transforming); 4 =
%            untransformed
%  Xi: Un-standardized n x p design matrix
%  yi: Un-standardized n x 1 response
%  numte: number of observations to set aside for test/validation
% Output:
%  te: row index of X values selected for test/validation
%  Xtr: standardized training design matrix
%  Xte: standardized test design matrix
%  ytr: normalized response (z-score)
%  yitr,yite: unnormalized response split by training and test
%  myi, syi:  mean and std deviance of unnormalized training response
%  yite: unnormalized response and associated mean and std deviation
% Savitsky, Vanucci (2008 - 2011)

[N L] = size(Xi);

% TRANSFORM {Xi}
switch type
    case 1 % Center Xi
        xbar = mean(Xi);
        Xbar = repmat(xbar,N,1);
        X = Xi - Xbar;

    case 2 % Transform X to [0,1]
        X = zeros(N,L);
        for j = 1:L
            bottom = min(Xi(:,j));
            top = max(Xi(:,j));
            X(:,j)=(Xi(:,j)- bottom)/(top-bottom);
        end;
        
    case 3 % Standardize X
        X = zscore(Xi);
        X = normcdf(X);
        
    case 4 % Untransformed
        X=Xi;
end;

%SPLIT {Xi,X,yi,y} INTO TRAIN/TEST SETS

% Randomly divide sample between train and test
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
indice = randperm(N);
te = indice(1:numte);
tr = indice(numte+1:N);
Xte = X(te,:);
yite = yi(te);
Xtr = X(tr,:);
yitr = yi(tr);

% Normalize training data response
ytr = zscore(yitr);
myi = mean(yitr);
syi = std(yitr);





