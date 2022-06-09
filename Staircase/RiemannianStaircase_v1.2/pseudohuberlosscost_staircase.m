function [Y, problem] = pseudohuberlosscost_staircase(data, epsilon, pp, varargin)
% function [Y, problem] = pseudohuberlosscost_staircase(data, epsilon)
% function [Y, problem] = pseudohuberlosscost_staircase(data, epsilon, pp)
% function [Y, problem] = pseudohuberlosscost_staircase(data, epsilon, pp, Y0)
% function [Y, problem] = pseudohuberlosscost_staircase(data, epsilon, pp, Y0, options)
%
% Helper function to call generic_staircase with a pseudohuberlosscost
% driver. The optional parameters pp, Y0 and options may be skipped by
% passing the empty matrix [] in their place.
%
% See also: pseudohuberlosscost generic_staircase
%
% Nicolas Boumal, UCLouvain, May 15, 2014.

    % Default choice of pp.
    m = data.m;
    d = data.d;
    if ~exist('pp', 'var') || isempty(pp)
        pataki = floor(max(roots([1, 1, -m*d*(d+1)-2])));
        pp = unique([ d+1, pataki:pataki:m*d, m*d ]);
    end
    
    % We know that Manopt will use an approximate Hessian for this problem.
    state = warning('off', 'manopt:getHessian:approx');

    [Y, problem] = generic_staircase(...
                            @(p) pseudohuberlosscost(data, epsilon, p), ...
                            pp, varargin{:} );
    
    % Restore the warning state
    warning(state);

end
