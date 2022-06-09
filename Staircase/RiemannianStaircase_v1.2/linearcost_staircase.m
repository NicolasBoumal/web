function Y = linearcost_staircase(C, m, d, pp, varargin)
% function Y = linearcost_staircase(C, m, d)
% function Y = linearcost_staircase(C, m, d, pp)
% function Y = linearcost_staircase(C, m, d, pp, Y0)
% function Y = linearcost_staircase(C, m, d, pp, Y0, options)
%
% Helper function to call generic_staircase with a linearcost driver.
% The optional parameters pp, Y0 and options may be skipped by passing the
% empty matrix [] in their place.
%
% Nicolas Boumal, UCLouvain, May 14, 2014.

    % Default choice of pp.
    if ~exist('pp', 'var') || isempty(pp)
        pataki = floor(max(roots([1, 1, -m*d*(d+1)-2])));
        pp = unique([ d+1, pataki:pataki:m*d, m*d ]);
    end

    Y = generic_staircase(@(p) linearcost(C, m, d, p), pp, varargin{:});

end
