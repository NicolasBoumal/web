function data = prepare_data(m, d, C, I, J)
% function data = prepare_data(m, d, C)
% function data = prepare_data(m, d, C, I, J)
%
% Helper function for use with huberlosscost to prepare a type of
% synchronization problem.
%
% If your data is in the format (C, I, J) and you would like to have it in
% the format (C), call this function, then pass the result to data2matrix.
%
% INPUTS:
%
%  m, d : Dimensions of the problem. m is the number of nodes in the graph.
%         d is the size of the dxd data matrices associated to the edges.
%
%  C :    If only C is given as input, it must be a matrix of size nxn,
%         where n = m*d. C is assumed symmetric and only its upper
%         triangular part is read. C is assumed to be partitioned in blocks
%         of size dxd. The diagonal blocks are ignored. The (i, j)th block
%         of size dxd is then considered as the data associated to the edge
%         (i, j) in the synchronization problem. In the data structure, C
%         will be transformed to a 3D array as described below: zero blocks
%         are ignored, nonzero blocks are stored in slices and their
%         position is recored in the vectors I and J.
% 
%  C, I, J :
%
%         C is d-by-d-by-k real 3D array (matrix with three dimensions),
%         such that the slice C(:, :, t) is a dxd matrix of data
%         corresponding to the edge (I(t), J(t)). It is implicitly assumed
%         that the data associated to edge (J(t), I(t)) is C(:, :, t)'.
%         I and J are length-k vectors of indices between 1 and m, defining
%         the edges of a graph: (I(t), J(t)) is an edge of the graph for
%         each t. The graph is assumed undirected: if the edge (4, 9) is
%         specified, do not specify the edge (9, 4).
%
% OUTPUTS:
%
%  data : a structure containing all the given information plus some
%         precomputed data.
%
% See also: huberlosscost data2matrix
%
% Nicolas Boumal, UCLouvain, May 15, 2014
    
    data = struct();
    
    if nargin == 5
    
        assert(isvector(I) && isempty(setdiff(I, 1:m)), ...
               'I must be a vector of indices in 1:m.');
        assert(isvector(J) && isempty(setdiff(J, 1:m)), ...
               'J must be a vector of indices in 1:m.');

        k = length(I);
        assert(length(J) == k, 'I and J must be vectors of same length');

        assert(all(size(C) == [d, d, k]), ...
         'C must be a 3D array of size d by d by k, where k = length(I).');
        
        data.C = C;
    
    elseif nargin == 3
        
        n = m*d;
        assert(all(size(C) == [n, n]), ...
               'If C is passed without I and J, it must have size md x md');
        
        % Do a first pass on the data to record where the nonzero blocks of
        % C are. When we know that, we can preallocate for I, J and CC.
        counter = 0;
        for j = 2 : m
            for i = 1 : (j-1)
                ii = (i-1)*d + (1:d);
                jj = (j-1)*d + (1:d);
                Cij = C(ii, jj);
                if norm(Cij) ~= 0
                    counter = counter + 1;
                end
            end
        end
           
        I = zeros(counter, 1);
        J = zeros(counter, 1);
        CC = zeros(d, d, counter);
        counter = 1;
        for j = 2 : m
            for i = 1 : (j-1)
                ii = (i-1)*d + (1:d);
                jj = (j-1)*d + (1:d);
                Cij = C(ii, jj);
                if norm(Cij) ~= 0
                    I(counter) = i;
                    J(counter) = j;
                    CC(:, :, counter) = Cij;
                    counter = counter + 1;
                end
            end
        end
        
        data.C = CC;
        
        k = length(I);
        
    else
        error('prepare_data accepts either 3 or 5 inputs.');
    end

    data.m = m;
    data.d = d;
    data.k = k;
    data.I = I(:);
    data.J = J(:);
    
    % These matrices are used, for example in huberlosscost, to speed up
    % the computation of certain gradients.
    [maskI, maskJ] = computemasks(m, I, J);
    data.maskI = maskI;
    data.maskJ = maskJ;
    
end

function [maskI, maskJ] = computemasks(m, I, J)

    M = length(I);
    maskI = sparse(I, 1:M, ones(M, 1), m, M, M);
    maskJ = sparse(J, 1:M, ones(M, 1), m, M, M);
    
end
