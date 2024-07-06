function [diffsq J n1 n2] = difference(X1,X2)
% Computes n1*n2 x p matrix of euclidean distances between cell values
% of matrices X1 (n1 x p), X2 (n2 x p), where each observation x_j is p x 1
% Savitsky and Vannucci (2008 - 2011)
%
[n1 p] = size(X1);
n2 = size(X2,1);

% diff = zeros(n1,n2,p);
% diffsq = zeros(n1,n2,p);
% for i = 1:n1
%     for j = 1:n2
%         diff(i,j,:) = X1(i,:)-X2(j,:);
%         diffsq(i,j,:) = diff(i,j,:).^2;
%     end;
% end;

diffsqAll = zeros((n1*n2),p);
count = 1;
for i = 1:n1
    for j = 1:n2
    diffsqAll(count,:) = (X1(i,:)-X2(j,:)).^2;
    count = count + 1;
    end;
end;

[diffsq I J] = unique(diffsqAll,'rows');