function [Y, info] = linearcost_sdplr(C, m, d, p)
% function [Y, info] = linearcost_sdplr(C, m, d)
% function [Y, info] = linearcost_sdplr(C, m, d, p)
%
% Uses the SDPLR algorithm to solve the problem
%
%   min Trace(CX)
% 
%     such that X is symmetric positive semidefinite of size nxn (n = m*d)
%     and the m diagonal blocks of X of size dxd each equal eye(d).
%
%     C must be a real matrix of size nxn.
%     Only the symmetric part of C matters.
%
% If p is specified, than the extra constraint rank(X) <= p is added.
% If so, there is no optimality guarantee, but it holds that if the
% returned Y is rank deficient, than it corresponds to a global optimum.
% You may check for rank deficiency by computing rank(Y'*Y) < p or by
% comparing cond(Y) to a large number (say, cond(Y) > 1e5 for example).
% 
% The optimal X is returned in factored form: Y is a matrix of size nxp for
% some p, such that X = Y*Y'.
% 
% The output 'info' is a structure of information returned by SDPLR.
% 
% Nicolas Boumal, UCLouvain, May 12, 2014.

    if ~exist('p', 'var') || isempty(p)
        p = [];
    end

    n = m*d;
    assert(all(size(C) == n), 'C must be square of size md x md.');
    
    % Setting up the problem in SDPLR format (same as SeDuMi format).
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
    pars.feastol = 1e-6;
    pars.centol = 1e-1;
    pars.dir = 1; % 1 or 2
    pars.penfac = 2.0;
    pars.reduce = 0;
    pars.limit = 3600;
    pars.printlevel = 1;
    pars.forcerank = p;
    pars.soln_factored = 1;
    
    % Solve with SDPLR. The actual computed solution is X = Y*Y'.
    [r, ~, info] = sdplr(A, b, c, K, pars);
    Y = r{1};
    
end
