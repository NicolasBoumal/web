function problem = pseudohuberlosscost(data, epsilon, p)
% function problem = pseudohuberlosscost(data, epsilon, p)
%
% Returns a Manopt problem structure corresponding to the problem
%
%   min Sum_{i~j} pseudohuber(norm(C_ij Y_j - Y_i, 'fro'))
% 
%     such that Y has size n x p (n = m*d) and is divided in m slices of
%     size d x p each, called Y_i, i = 1 : m. Each Y_i is constrained to be
%     orthonormal: Y_i*Y_i' = eye(d). The sum is over the edges of a graph
%     described in the data structure. The pseudo Huber-loss function is a
%     smoothed version of the absolute value function, with smoothing
%     parameter epsilon. The closer epsilon is to zero, the closer
%     pseudohuber looks like the absolute value function, but the harder
%     the problem becomes from a numerical point of view. It is typically
%     better to start with a reasonable value of epsilon (say, .1) and to
%     decrease that value progressively, using the previous estimator as
%     initial guess.
%
%     The data structure describes the blocks C_ij and the graph on which
%     they are defined. Such a structure is obtained by calling
%     prepare_data.
%
% This corresponds to constraining the matrix X = Y*Y' to be symmetric,
% positive semidefinite with rank(X) <= p and the m diagonal blocks of X of
% size dxd each equal eye(d).
%
% Note that the cost function is strictly concave in X. Ignoring the rank
% constraint, the search space for X is convex; thus, the local optimizers
% are extreme points of the search space, and those are known to have small
% rank. Unfortunately, the problem as a whole is nonconvex, hence we cannot
% guarantee that the returned solution is optimal.
%
% If, when solving this problem, the returned Y is a rank-deficient
% second-order critical point of this problem (that is: the gradient is
% (almost) zero and the Hessian is (almost) positive semidefinite), than
% X = Y*Y' is a KKT point of the problem which consists in minimizing the
% same cost function such that X is as described above (without rank
% constraint). Numerically, one can check for rank deficiency by checking
% rank(Y'*Y) < p (not rank(Y*Y')). You may also compute cond(Y) and compare
% it against a chosen large number (the larger cond(Y), the closer Y is to
% rank deficiency).
%
% To solve the problem, call one of Manopt's solvers. For example:
%
%    Y = trustregions(problem);
%
% will use Riemannian Trust-Regions with a random initial guess and default
% options. See www.manopt.org for help.
%
% See also pseudohuberlosscost_staircase prepare_data linearcost round2orthogonal
%
% Nicolas Boumal, UCLouvain, March 4, 2014.
% Modified on August 14, 2014.
    
    % Extract data from the synchronization problem structure.
    % Such a data structure is obtained from prepare_data.
    %
    m = data.m; % number of orthonormal matrices to find
    d = data.d; % dimension of the ambient space
    k = data.k; % the number of measurements available
    C = data.C; % data on the edges (relative measurements C_ij)
    I = data.I; % edges of the graph
    J = data.J;
    maskI = data.maskI;
    maskJ = data.maskJ;
    
    manifold = stiefelstackedfactory(m, d, p);
    problem.M = manifold;
    
    % Pseudo Huber loss: a smoothed version of the absolute value function.
    % Note that the expected input is the square of x.
    function y = phi(x_squared)
        y = sqrt(epsilon^2 + x_squared) - epsilon;
    end

    % Derivative of phi at x, divided by x.
    function y = dphi_x(x_squared)
        y = 1./sqrt(epsilon^2 + x_squared);
    end

    % Cost function. Remember that Y is given as a 2D array which is
    % implicitly "sliced" in as many orthonormal matrices as there are
    % diagonal blocks in the constraints. The function manifold.to3D()
    % transforms the 2D array in a 3D array to ease access to these slices.
    % The function manifold.to2D() reverts that change.
    problem.cost = @cost;
    function f = cost(Y)
        
        Y3 = manifold.to3D(Y);
        
        Yi = Y3(:, :, I);
        Yj = Y3(:, :, J);
        residues = multiprod(C, Yj) - Yi;
        residues_sqnorms = multisqnorm(residues);
        f = sum(phi(residues_sqnorms)) / k;
        
    end

    % Classical (Euclidean) gradient of the cost function. It will be
    % transformed into the Riemannian gradient automatically by Manopt.
    problem.egrad = @egrad;
    function G = egrad(Y)
        
        Y3 = manifold.to3D(Y);
        G3 = zeros(size(Y3));
        
        Yi = Y3(:, :, I);
        Yj = Y3(:, :, J);
        residues = multiprod(C, Yj) - Yi;
        residues_sqnorms = multisqnorm(residues);
        
        weights = dphi_x(residues_sqnorms);
        weighted_residues = multiscale(weights, residues);
        alternate_weighted_residues = ...
                              multiprod(multitransp(C), weighted_residues);
        
        % We write the code for the gradient in this way to avoid looping
        % over k elements, seen as k may be of order m^2 and Matlab doesn't
        % like big loops. An equivalent but much slower code is given
        % below, for readability.
        for k1 = 1 : d
            for k2 = 1 : p
                G3(k1, k2, :) = - maskI * (squeeze(weighted_residues(k1, k2, :))) ...
                                + maskJ * (squeeze(alternate_weighted_residues(k1, k2, :)));
            end
        end
        
        G = manifold.to2D(G3);
        G = G / k;
        
    end

end


% Slow code kept here to help with readability: this is equivalent
% to what the code in egrad with maskI and maskJ does.
% for kk = 1 : k
%     i = I(kk);
%     j = J(kk);
%     G3(:, :, i) = G3(:, :, i) -           weighted_residues(:, :, kk);
%     G3(:, :, j) = G3(:, :, j) + alternate_weighted_residues(:, :, kk);
% end
