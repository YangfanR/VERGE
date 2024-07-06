function y = logdet(A)
% log(det(A)) where A is positive-definite.
% This is faster and more stable than using log(det(A)).

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.
if(any(eig(A) < 0))
    A = A + 10e3 * eye(size(A,1));
end
U = chol(A);
y = 2*sum(log(diag(U)));
