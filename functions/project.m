function [Xtrs indx] = project(Xtr,psample)
% Given a the training design matrix, Xtr, and psample, a percent subset of the data
% to use for inverse computation of the GP covariance matrix, C, this
% function randomly selects a psample subset of the records, using latin
% hypercube sampling to ensure values are drawn from across the quantiles
% of the data.  The sampled design matrix is Xtrs and the index of training
% cases sampled is returned in indx
% Savitsky/Vannucci (2008 - 2011)

% Extract size of training set
numtr = size(Xtr,1);

% Randomly select sub-set of points to use for "m < numtr" dimension covariance
% matrix.

m = ceil(psample*numtr); % m < numtr determines size of cov matrix
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
% indx = randperm(numtr);
% s = indx(1:m);
% Xtrs = Xtr(s,:);

indx = round(lhsu(1,numtr,m));
[inds,ind] = sort(indx);
unsorted = 1:length(indx);
newInd(ind) = unsorted;
for i = 2:length(inds)
    if inds(i)==inds(i-1) 
        inds(i) = inds(i) + 1;
    end;
end;
indx = inds(newInd);
        
Xtrs = Xtr(indx,:);