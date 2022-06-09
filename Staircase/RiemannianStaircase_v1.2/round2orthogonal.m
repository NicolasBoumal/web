function Q = round2orthogonal(Y, d)
% function Q = round2orthogonal(Y, d)
% 
% Given an nxp matrix Y and a dimension d such that n = md with m integer,
% returns Q, a 3D array of size dxdxm such that each slice Q(:, :, i) of Q
% is an orthogonal matrix. Q is a projection of the slices of the matrix Y
% to the space of tuples of m orthogonal matrices. See the code for
% details.
%
% See also: linearcost huberlosscost
%
% Nicolas Boumal, UCLouvain, March 5, 2014.

    % Determine the dimensions of the problem.
    [n, ~] = size(Y);
    m = round(n/d);
    assert(round(m*d) == round(n), 'Y must have m*d rows, with m integer.');
    manifold = stiefelstackedfactory(m, d, d);
    
    % Extract an approximation of Y to a matrix with d <= p columns.
    [U, S, ~] = svds(Y, d);
    Yd2 = U*S;
    
    % Reshape the approximated Y to a 3D array of square slices.
    Yd = manifold.to3D(Yd2);
    
    % Project each slice to the orthogonal group O(d).
    Q = zeros(d, d, m);
    for i = 1 : m
        [u, ~, v] = svd(Yd(:, :, i));
        Q(:, :, i) = u*v';
    end
    
end



% % % This part is optional: it further moves from Q such that the new
% % % model Ynew*Ynew' is as close as possible to the given model Y*Y', in
% % % Frobenius norm. If the given model has rank d, this step does
% % % nothing.
% % problem = linearcost(@(U) -Y*((Y'*U)/(m*n)), m, d, d);
% % Y0 = manifold.to2D(Q);
% % Ynew = trustregions(problem, Y0);
% % Q = manifold.to3D(Ynew);

