%% MATLAB main for the exercise Riemannian Hessian on Stiefel
% Output : The error in Frobenius norm between the exact hessian and 
% the finit difference approximation with time step 
% t=10^(-1),10^(-2),10^(-4),10^(-8)

% Initialization
close all;
clear;
clc;


% Creation of random symetric matrix
A = normrnd(3,10,[5,5]);
A = Sym(A);

% Creation of random point on Stiefel
X = normrnd(3,10,[5,3]);
[X,~]=qr(X,0);

% Creation of random tangent vector in T_X St(n,p)
U = normrnd(3,10,[5,3]);
U = Proj(X,U);

% Check Hessian
fprintf('Error of finite difference approximation with time step t=10^(-1): \n')
t = 10^(-1);
err = CheckH(A,X,U,t);
disp(err)

fprintf('Error of finite difference approximation with time step t=10^(-2): \n')
t = 10^(-2);
err = CheckH(A,X,U,t);
disp(err)

fprintf('Error of finite difference approximation with time step t=10^(-4): \n')
t = 10^(-4);
err = CheckH(A,X,U,t);
disp(err)

fprintf('Error of finite difference approximation with time step t=10^(-8): \n')
t = 10^(-8);
err = CheckH(A,X,U,t);
disp(err)



