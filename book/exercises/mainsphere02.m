%% MATLAB main for the exercise Rayleigh quotient on the sphere
% Output : This code shows the values of the Rayleigh quotient of a
% randomly generated symetric matrix A and the steps for Riemannian
% Gradient Descent as red points on the sphere. The algorithm also displays
% the final point x, the final function value x'Ax and the number of
% iterations.

% Initialization
close all;
clear;
clf; 
clc; 

% Creation of random symetric matrix and starting point for RGD
A = randn(3); % get a random matrix
A = 0.5*(A + A'); % making it symmetric 
x0 = randn(3, 1); % get a random point
x0 = x0/norm(x0); % normalize to get it on the sphere

% Constant step size for RGD.
% One can show that this works, but it's likely to be slow.
alpha = 1/(2*norm(A, 'fro'));

% Applying Riemannian Gradient Descent
[x, iterates] = RGDsphere02(A, x0, alpha, 100);

% Plotting the sphere and the RGD steps
f = @(x) x'*A*x;
resolution = 51;
C = zeros(resolution);
[X,Y,Z] = sphere(resolution);

v = zeros(3, 1);
for i = 1:resolution
    for j = 1:resolution
        v(1) = X(i, j); 
        v(2) = Y(i, j);
        v(3) = Z(i, j); 
        C(i,j) = f(v); 
    end 
end

surf(X,Y,Z,C);
colorbar('eastoutside');
axis equal;
hold on;
% Adding the iterates of the algorithm as red points on the sphere
plot3(iterates(1, :), iterates(2, :), iterates(3, :), ...
      '.', 'MarkerSize', 20, 'Color', 'r');
