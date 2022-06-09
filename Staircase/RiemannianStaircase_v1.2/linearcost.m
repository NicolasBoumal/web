function problem = linearcost(C, m, d, p)
% function problem = linearcost(C, m, d, p)
%
% Returns a Manopt problem structure corresponding to the problem
%
%   min Trace(Y'*C*Y)
% 
%     such that Y has size n x p (n = m*d) and is divided in m
%     blocks of size d x p each. Call them Yi. Each Y_i is
%     constrained to be orthonormal: Y_i*Y_i' = eye(d).
%
%     C must be a real SYMMETRIC matrix of size nxn, or a function handle
%     such that C(V) = C*V, for some real SYMMETRIC matrix C of size nxn
%     and any real matrix V with n rows. If C is not a symmetric matrix,
%     the algorithm will stagnate without converging.
%
% This corresponds to constraining the matrix X = Y*Y' to be symmetric,
% positive semidefinite with rank(X) <= p and the m diagonal blocks of X of
% size dxd each equal eye(d).
%
% If, when solving this problem, the returned Y is a rank-deficient local
% optimizer of this problem, than X = Y*Y' is a global optimizer of the
% semidefinite program which consists in minimizing Trace(C*X) such that X
% is as described above. Numerically, one can check for rank deficiency by
% checking rank(Y'*Y) < p (not rank(Y*Y')). You may also compute cond(Y)
% and compare it against a chosen large number (the larger cond(Y), the
% closer Y is to rank deficiency).
%
% To solve the problem, call one of Manopt's solvers. For example:
%
% Y = trustregions(problem);
%
% will use Riemannian Trust-Regions with a random initial guess and default
% options. See www.manopt.org for more help.
%
% See also linearcost_staircase huberlosscost prepare_data data2matrix round2orthogonal
%
% Nicolas Boumal, UCLouvain, March 4, 2014.

    % The geometry is a product of m Stiefel manifolds where each
    % orthonormal matrix has size dxp (d unit-norm, orthogonal vectors of
    % length p).
    manifold = stiefelstackedfactory(m, d, p);
    
    problem.M = manifold;
    
    if ~isa(C, 'function_handle')
        assert(all(size(C) == [m*d, m*d]), ...
               'C must be an m*d x m*d real matrix.');
        % Let's make sure C is symmetric.
        C = (C+C')/2;
        Cfun = @(Y) C*Y;
    else
        % Here, we cannot enforce symmetri of C. Let's hope the user read
        % the documentation :).
        Cfun = C;
    end

    % The computation of C*Y can be reused when computing the cost or the
    % gradient at a given point Y. This function makes sure the store
    % structure contains the product CY.
    function store = prepare(Y, store)
        if ~isfield(store, 'CY')
            store.CY = Cfun(Y);
        end
    end

    % The cost is trace(Y'*C*Y).
    problem.cost = @cost;
    function [f, store] = cost(Y, store)
        store = prepare(Y, store);
        CY = store.CY;
        f = Y(:)'*CY(:);
    end

    % The Euclidean (classical) gradient w.r.t. Y is 2*C*Y.
    % For the cost function, it does not matter whether C is symmetric or
    % not: only the symmetric part will play a role. For the gradient and
    % the Hessian though, it is critical that C be symmetric.
    problem.egrad = @egrad;
    function [G, store] = egrad(Y, store)
        store = prepare(Y, store);
        G = 2*store.CY;
    end

    % The Euclidean (classical) Hessian w.r.t. Y along the tangent vector
    % Ydot is 2*C*Ydot (that is, the derivative of the gradient 2*C*Y).
    problem.ehess = @ehess;
    function [H, store] = ehess(Y, Ydot, store) %#ok<INUSL>
        H = 2*Cfun(Ydot);
    end

end
