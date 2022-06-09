function Q = randortho(n, N)
% Returns a 3D matrix of size nxnxN such that each slice of size nxn is
% orthogonal.

    Q = zeros(n, n, N);
    for i = 1 : N
        [q, r] = qr(randn(n));
        Q(:, :, i) = q*diag(sign(diag(r)));
    end

end