%% MATLAB main for the exercise RGD on Stiefel
% Output : This code shows the results of RGD with constant step size or LS 
% applied on the map X \in \St(n,p) \mapsto Tr(X'AX) for a random 
% symmetric matrix A \in \R^{n \times n}.
% The algorithm also display the final point (X, the final function 
% value Tr(X'AX) and the number of iterations of RGD. 

% Initialization
close all;
clear;
clc; 

% Creation of random symetric matrix and starting point for RGD
n = 10; % Size of Stiefel (St(n,p))
p = 4;
A = rand(n,n); % get a random matrix
A = (A+A') / 2; % make it symmetric
X0 = randn(n, p); % get a random point
[X0,~]=qr(X0,0); % put it on Stiefel
LS = true; % Line search


% Step size for RGD. The choice is based on the behaviour Gradient Descent
% in an euclidean space.
alpha = 1/(2*norm(A, 'Fro'));

% Applying Riemannian Gradient Descent
[X,iterates,conv] = RGDLS(A, X0, LS, alpha, 1000);