function MSE = mse(R1, R2)
% Given two 3D arrays R1 and R2 of same size dxdxm, computes the Mean
% Squared Error separating R1 and R2, computed as
%
%  min_Q sum_i norm( R1(:, :, i)*Q - R2(:, :, i), 'fro')^2, Q orthogonal.
%

    assert(all(size(R1) == size(R2)));
    assert(size(R1, 1) == size(R1, 2));
    d = size(R1, 1);
    P = mean(multiprod(multitransp(R1), R2), 3);
    MSE = max(0, 2*(d-sum(svd(P))));

end
