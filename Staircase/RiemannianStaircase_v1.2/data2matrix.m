function C = data2matrix(data)
% function C = data2matrix(data)
% 
% Given a data structure as obtained from the function prepare_data,
% returns a matrix C of size md x md, full or sparse depending on the
% amount of data to hold, such that C is structured by dxd blocks and
% calling prepare_data(m, d, C) would return the same data structure as the
% one received in input (up to ordering of the edges I, J and the slices
% C). The dimensions m and d are as specified in the data structure.
%
% Important: it is assumed that, in the data structure, if an edge (i, j)
% is specified, then the edge (j, i) is /not/ specified. If this condition
% is not fulfilled, one of the data matrices in data.C will be used for the
% blocks (i, j) and (j, i), but it is not specified which one will be used.
%
% See also: prepare_data linearcost
%
% Nicolas Boumal, UCLouvain, May 15, 2014.

    CC = data.C;
    I = data.I;
    J = data.J;
    m = data.m;
    d = data.d;
    k = data.k;
    
    n = m*d;

    % Depending on whether the graph is almost complete or not, allocate C
    % as a full or as a sparse matrix.
    if k >= nchoosek(m, 2) * 0.5
        C = zeros(n, n);
    else
        C = spalloc(n, n, 2*d*d*k);
    end
    
    for t = 1 : k
        i = I(t);
        j = J(t);
        ii = (i-1)*d + (1:d);
        jj = (j-1)*d + (1:d);
        C(ii, jj) = CC(:, :, t);
    end
    C = C + C';

end
