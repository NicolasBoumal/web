function Q = randperms(n, N)
% Returns a 3D matrix of size nxnxN such that each slice of size nxn is
% a permutation.

    I = eye(n);
    Q = zeros(n, n, N);
    for i = 1 : N
        p = randperm(n);
        Q(:, :, i) = I(1:n, p);
    end

end
