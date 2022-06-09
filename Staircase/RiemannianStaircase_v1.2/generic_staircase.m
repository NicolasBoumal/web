function [Y, problem] = generic_staircase(genproblemfun, pp, varargin)
% function [Y, problem] = generic_staircase(genproblemfun, pp)
% function [Y, problem] = generic_staircase(genproblemfun, pp, Y0)
% function [Y, problem] = generic_staircase(genproblemfun, pp, Y0, options)
%
% The first input genproblemfun is a function handle that, given a single
% input p, returns a Manopt problem structure corresponding to a problem of
% this sort:
%
%   min f(Y*Y')
% 
%     such that Y has size n x p (n = m*d) and is divided in m
%     blocks of size d x p each. Call them Yi. Each Yi is
%     constrained to be orthonormal: Yi*Yi' = eye(d).
%
% For example:
% genproblemfun = @(p) linearcost(C, m, d, p);
% genproblemfun = @(p) huberlosscost(synchroproblem, epsilon, p);
% etc.
%
% The vector pp is a sorted (increasing) list of integers representing
% the values of p to try (p is a bound on the rank of X = Y*Y'). The values
% in p must be >= d and <= n. Typically, the smallest value should be at
% least d+1. The staircase algorithm implemented here will use Manopt to
% find a local optimizer of the problem with the first value p = pp(1). If
% the returned matrix Y is rank deficient, then X = Y*Y' is a solution of
% the associated problem without rank constraint and Y is returned. If not,
% p is set to pp(2) and we iterate.
%
% It is normal that it takes some tuning to determine a good rank sequence.
% If you are feeling lucky, try first with d+1. Another reasonable choice
% (as first or as second p) is the square roort of the number of
% constraints (besides semidefiniteness).
%
% If Y0 is specified, it must be a feasible point for the optimization
% problem with p <= pp(1). It will be used as first initial guess.
% Alternatively, Y0 may be a 3D array (3D matrix) of size d x d x m with
% each slice Y0(:, :, i) an orthogonal matrix. Set to [] to use a default
% value.
%
% If options is specified, it must be a structure of options that will be
% passed to the Manopt solver in use (trustregions). Set to [] to ignore.
%
% The output 'problem' is the Manopt problem structure from which Y was
% obtained by optimization. Thus, Y is a point on the manifold problem.M.
%
% See www.manopt.org for help regarding the Manopt toolbox for optimization
% on manifolds.
%
% See also: linearcost huberlosscost
%
% Nicolas Boumal, UCLouvain, May 13, 2014.

    assert(issorted(pp), 'Values in pp should be sorted increasingly.');

    if numel(varargin) >= 1 && ~isempty(varargin{1})
        Y = varargin{1};
        if ndims(Y) == 2
            assert(size(Y, 2) <= pp(1), ...
                   'The initial guess Y must have at most pp(1) columns.');
        elseif ndims(Y) == 3    
            % Convert the 3D array to a 2D matrix
            d = size(Y, 1);
            assert(size(Y, 2) == d, ...
                   ['If Y is a 3D array, its first and second ' ...
                    'dimension must be equal.']);
            problem = genproblemfun(d);
            Y = problem.M.to2D(Y);
        end
    else
        Y = [];
    end
    
    if numel(varargin) >= 2 && ~isempty(varargin{2})
        options = varargin{2};
    else
        options = struct();
    end
    
    % For each rank to be checked
    for p = pp
        
        % Generate a Manopt problem structure corresponding to the
        % optimization problem min f(Y*Y') such that Y has orthonormal
        % slices and rank(Y*Y') <= p.
        problem = genproblemfun(p);
        
        % If we already had a Y from a previous iteration, warm start from
        % there by expanding it to a feasible point on the new, larger
        % dimensional search space: simply add zero columns. The issue is
        % that this new point will be a saddle point. Thus, we escape from
        % this saddle using a small random vector, which we retract.
        if ~isempty(Y)
            manifold = problem.M;
            Y = [Y zeros(size(Y, 1), p-size(Y, 2))]; %#ok<AGROW>
            escape_vector = 1e-6*sqrt(manifold.dim())*manifold.randvec(Y);
            Y = manifold.retr(Y, escape_vector);
        end
        
        % Optimize to attain a local optimizer of this problem.
        Y = trustregions(problem, Y, options);
        
        % If the attained local optimizer is (numerically) rank deficient,
        % we may stop, because Y*Y' is a global optimizer of the rank
        % unconstrained problem.
        if cond(Y) >= 1e+5
            break;
        end
        
    end

end
