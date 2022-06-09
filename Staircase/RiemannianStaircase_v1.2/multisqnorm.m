function sqnorms = multisqnorm(A)
% Given a 3D array A of size nxmxk, returns a column vector of length k
% such that the i-th entry is the squared Frobenius norm of the i-th slice
% of A, sqnorms(i) = norm(A(:, :, i), 'fro')^2.
%
% Nicolas Boumal, UCLouvain, March 4, 2014.

    sqnorms = sum(reshape(A, [size(A, 1)*size(A, 2), size(A, 3)]).^2)';

end
