%% MATLAB main for the exercise Fréchet mean on the sphere
% Output : This code shows the used data for the computation, the results of 
% checkgradient applied to the implementation of the cost function and
% the gradient of the cost function, the results of RGD, RCG and RTR 
% applied to this data with a random initial point and the north pole 
% as initial point and the spectrum of the Hessian at the final point of
% RTR.

% ATTENTION : Parts of this code use ManOpt. You need to have ManOpt to
% run the code.

% Run different parts of the exercise. With 'generatingdata' one can
% generate the data, with 'checkgrad' one can ran checkgradient, with
% InitialRand one can run RGD, RCG and RTR with a random initial point 
% on the current data and with InitialPole one can run RGD, RCG and RTR
% with the initialization at the north pole on the current data.
% InitialRand and InitialPole show also the spectrum of the Hessian at the 
% final point of RTR
% ATTENTION : Do not run InitialRand and InitialPole without data.
generatingdata = true;
checkgrad = false;
InitialRand = true;
InitialPole = false;

% Initialization
clc;
close all;

%% Generating data
% Generating points x_1,...,x_n and parameteres for f.
% The points are computed by first taking random normal points, then
% multiplying the first entry by a coefficent m and finily finily dividing
% every vector by its norm.

if generatingdata

f.n = 50; % number of points used
m = 10; % coefficent m
f.d = 3; % dimension
f.X = randn(f.d,f.n);
f.X(1,:) = m * abs(f.X(1,:));
for i=1:f.n
    f.X(:,i) = f.X(:,i)/norm(f.X(:,i));
end
fprintf('Data points:\n')
disp(f.X)

end

%% Check gradient

% Check gradient with ManOpt
problem.M = spherefactory(f.d);
problem.cost = @(x) cost(x,f);
problem.grad = @(x) gradcost(x,f);

if checkgrad
checkgradient(problem)
end

%% Runninf different optimization algorithm
%% Random Initial Point

if InitialRand

% Initialization by ManOpt
fprintf('Random Initalization on the sphere:\n')
x =steepestdescent(problem);
fprintf('Final point RGD\n')
disp(x)
x = conjugategradient(problem);
fprintf('Final point RCG\n')
disp(x)
x = trustregions(problem);
fprintf('Final point RTR\n')
disp(x)

% Spectrum of Hessian at optimum
hessianspectrum(problem,x)

end

%% Initial Point North Pole

if InitialPole

% Initialization at north pole
fprintf('Initialization at north pole:\n')
x0 = zeros(f.d,1);
x0(1,1) = 1;
x =steepestdescent(problem,x0);
fprintf('Final point RGD\n')
disp(x)
x = conjugategradient(problem,x0);
fprintf('Final point RCG\n')
disp(x)
x = trustregions(problem,x0);
fprintf('Final point RTR\n')
disp(x)

% Spectrum of Hessian at optimum
hessianspectrum(problem,x)

end

