% Example using the Riemannian staircase method on a random synchronization
% problem.
%
% See http://perso.uclouvain.be/nicolas.boumal/ and http://www.manopt.org.
% Contact: nicolasboumal@gmail.com
%
% Nicolas Boumal, UCLouvain, May 17, 2014.

clear all;
close all;
clc;

% First things first: you must have Manopt on your Matlab path. Follow the
% installation instructions on www.manopt.org. As of writing this, you
% just need to go to the manopt folder and execute "importmanopt.m".

% We synchronize m orthogonal transformations in O(d), the orthogonal group
% containing matrices of size dxd. In particular, O(1) = {-1, +1}, so you
% can use this code to solve the Max-Cut SDP relaxation too.
m = 250;
d = 3;
n = m*d;

% Here are the target matrices, picked uniformly at random, arranged in a
% 3D array of size d x d x m. Each slice is one of the orthogonal matrices.
Rtrue = randortho(d, m);

% Re-arrange them in a matrix of size m*d x d = n x d, that is: simply
% stack the slices on top of each other. The functions multitransp and
% others (multiprod, multitrace, ...) are part of the Manopt package and
% make it easier (and faster) to work with 3D arrays.
Rtrue_stacked = reshape(multitransp(Rtrue), [d, n])';

% Generate a symmetric noise matrix W. This is Wigner noise, with entries
% above the diagonal being iid Gaussians with zero mean and standard
% deviation sigma. The entries on the dxd diagonal blocks are irrelevant.
sigma = 0.1;
W = sigma*randn(n);
W = (W+W')/2;

% This is our data matrix: a signal of rank d corresponding to the correct
% relative orthogonal transformations, superimposed with noise. C is
% symmetric, nxn.
C = Rtrue_stacked * Rtrue_stacked' + W;

% The problem is the following: the dxd blocks Cij correspond to noisy
% measurements of the relative rotations Ri*Rj'. Based on all these
% measurements, we wish to estimate the Ri's.

% We will call Manopt solvers with this optional options structure.
% Set the verbosity to 0 to mute the solvers. Other possible values are 1
% and 2 for increasingly much information.
options.verbosity = 0;

% Give this data to the Riemannian staircase solver. It will minimize the
% trace of -C*X under the constraint that X is positive semidefinite
% and the diagonal blocks of X are eye(d). This means we are maximizing the
% trace of C*X. X is returned in factored form: X = Y*Y'. We hope that the
% rank of Y will be d, in which case X is rank d and we have actually
% computed the maximum likelihood estimator of our problem (provided Y is a
% second-order critical point of the Riemannian problem, which is typically
% the case). If the rank of the optimal X is not d, the Riemannian
% staircase will return a solution anyway, but it will have higher rank.
% Y is an nxp matrix, for some value p >= d (typically, d+1). Each dxp
% slice of Y is orthonormal.
%
% Two important points:
%
%   (1) C must be symmetric. This is particularly important if C is not
%       known as a matrix but is known as a function handle.
%
%   (2) The staircase minimizes. Since we want to maximize, we pass -C. If
%       the sign is omitted, results typically make no sense and take more
%       time to compute. Also: it is good practice to normalize C so that
%       the expected optimal cost value will be on the order of 1 in
%       magnitude. In this way, the tolerance on the gradient norm will
%       always make sense.
%
% Nota bene: the empty matrices represent omitted parameters: default
% values will be used. If we do not need to pass the options structure, we
% may call this function with just the three first inputs.
Y = linearcost_staircase(-C/(n*m), m, d, [], [], options);

% We now project the orthonormal slices of Y to the orthogonal group O(d).
% Rhat is our estimator of the orthogonal transformations Rtrue. Rhat is a
% 3D array of the same size as Rtrue. If the rank of Y is d, which we
% hope, then there is no loss going from Y to Rhat: we just rearrange Y so
% that all the information is in the d first columns and we drop the last
% ones.
Rhat = round2orthogonal(Y, d);

% Compute the distance between Rtrue and Rhat (the error). Notice that we
% can only estimate Rtrue up to a global orthogonal transformation since we
% only have relative information. The function mse (mean squared error)
% takes that into account and computes the minimal distance between Rtrue
% and Rhat, after these two are best aligned with each other.
MSE = mse(Rtrue, Rhat);

fprintf('For a noise level sigma = %g, linearcost attained an MSE of %g.\n', sigma, MSE);
if cond(Y'*Y) >= 1e8
    fprintf('Y is rank-deficient.\n');
else
    fprintf('Y is not rank deficient.\n');
end

% The following spectral method often performs very well for this specific
% problem and is very simple to implement. It is fast and the SDP method
% above does not perform significantly better on this vanila version of
% the problem. The advantages of the approach above are:
%  1) We now have much better control over the cost function (spectral
%     methods are by nature least-squares methods): see
%     pseudohuberlosscost_staircase for example, for robust versions.
%  2) When the solution of the SDP is rank d, we actually know that we have
%     the maximum likelihood estimator. The eigenvector method does not
%     attain the MLE.
%  3) For small m (say, 3 or 5), the SDP does outperform the spectral
%     method markedly.
[V, D] = eigs(C, d);
Reig = round2orthogonal(V*sqrt(D), d);
fprintf('The eigenvector method achieves an MSE of %g.\n', mse(Rtrue, Reig));
% Nota bene: for graphs with nodes of different degrees, this is not the
% correct formulation of the eigenvector method. See this paper:
% Amit Singer, Applied and Computational Harmonic Analysis, 2011,
% Angular synchronization by eigenvectors and semidefinite programming.

% Let's replace some of the measurements by outliers: they will be random
% orthogonal matrices (instead of Ri*Rj') with noise. Thus, these outliers
% look exactly like legit measurements, except they try to confuse us.
outlier_fraction = 0.50;
[II, JJ] = erdosrenyi(m, outlier_fraction);
n_outliers = length(II);
fprintf('\nNow replacing about %.2g%% of the data with outliers... ', ...
                                                     100*outlier_fraction);
for k = 1 : n_outliers
    i = II(k);
    j = JJ(k);
    ii = (i-1)*d + (1:d);
    jj = (j-1)*d + (1:d);
    Cij = randortho(d, 1);
    C(ii, jj) = Cij;
    C(jj, ii) = Cij';
end
actual_fraction = n_outliers / nchoosek(m, 2);
fprintf('done, with actually %.2g%% outliers.\n\n', 100*actual_fraction);


% Compute the SDP solution again on the corrupted data
Y = linearcost_staircase(-C/(n*m), m, d, [], [], options);
Rhat = round2orthogonal(Y, d);
fprintf('MSE reached with the linear cost (least squares): %g\n', mse(Rtrue, Rhat));

% And same for the eigenvector method
[V, D] = eigs(C, d);
Reig = round2orthogonal(V*sqrt(D), d);
fprintf('MSE reached with the eigenvector method: %g\n', mse(Rtrue, Reig));

% We now minimize a nonlinear, concave cost function with the corrupted
% data. Because of concavity, we cannot guarantee global optimality
% anymore. Nevertheless, concavity promotes low-rank solutions and we
% observe excellent performance in practice. We first obtain a data
% structure containing all the necessary information.
data = prepare_data(m, d, C);
% We obtain an initial guess with the linear cost.
Y = linearcost_staircase(-C/(n*m), m, d, [], [], options);
Rlud = round2orthogonal(Y, d);
% Then, for decreasing values of the smoothing parameter epsilon, we solve
% the problem with the pseudo Huber loss cost, using the previously
% obtained solution as initial guess. This is particularly interesting when
% sigma is zero: this cost function seems then to attain exact recovery
% even for a large fraction of outliers.
% An alternative, very similar but convex function is described in
% Wang & Singer 2012, Exact and Stable Recovery of Rotations for Robust
% Synchronization. They show that with their cost, they can guarantee
% exact recovery for up to 50% outliers. Unfortunately, we found that the
% algorithm here described is rather slow on that function.
for epsilon = [1, .1, .01, .001]
    Y = pseudohuberlosscost_staircase(data, epsilon, [], Rlud, options);
    Rlud = round2orthogonal(Y, d);
    fprintf('MSE reached with the pseudo Huber loss cost (epsilon = %g): %g\n', epsilon, mse(Rtrue, Rlud));
    fprintf('Is the found solution Y rank-deficient? %d\n', cond(Y'*Y) > 1e8);
end
