function [X, info] = linearcost_sedumi(C, m, d)
% function [X, info] = linearcost_sedumi(C, m, d)
%
% Uses the SeDuMi algorithm to solve the problem
%
%   min Trace(CX)
% 
%     such that X is symmetric positive semidefinite of size nxn (n = m*d)
%     and the m diagonal blocks of X of size dxd each equal eye(d).
%
%     C must be a real matrix of size nxn.
%     Only the symmetric part of C matters.
% 
% The output 'info' is a structure of information returned by SeDuMi.
% 
% Nicolas Boumal, UCLouvain, May 12, 2014.

    n = m*d;
    
    % Setting up the problem in SeDuMi format.
    % The equality constraints are Ax = b, with x = X(:).
    % The cost is c'*x, with c = C(:).

    c = C(:);

    % Build a mask for the constraints, to identify indices
    Z = sparse([], [], [], n, n, n*d);
    for i = 1 : m
        ii = (i-1)*d + (1:d);
        Z(ii, ii) = ones(d); %#ok<SPRIX>
    end
    indices = find(Z);
    n_constraints = length(indices);
    A = sparse(1:n_constraints, indices, ones(n_constraints, 1), ...
               n_constraints, n^2, n_constraints);

    Id = eye(d);
    b = repmat(Id(:), [m, 1]);

    % We further ask that X be symmetric, positive semidefinite.
    K.s = n;
    
    % Define options (see help sdplr)
    pars = struct();
    pars.eps = 1e-8;
    pars.errors = 1;
    
    % Solve with SeDuMi
    [x, ~, info] = sedumi(A, b, c, K, pars);
    X = reshape(x, [n, n]);
    
end
