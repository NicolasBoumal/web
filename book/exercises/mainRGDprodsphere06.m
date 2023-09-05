%% MATLAB main for the exercise RGD on product spheres
% Output : This code shows the values of the map (x,y) \ampsto x'My of a
% randomly generated matrix M on the square [0,2 \pi]^2.
% The algorithm also display the final point (x,y), the final function 
% value x'My and the number of iterations of RGD with a constant step size.

% Initialization
close all;
clear;
clc; 

% Creation of random symetric matrix and starting point for RGD
M = randn(2); % get a random matrix
x0 = randn(2, 1); % get a random point
x0 = x0/norm(x0); % normalize to get it on the sphere
y0 = randn(2, 1); % get a random point
y0 = y0/norm(y0); % normalize to get it on the sphere

% Step size for RGD. The choice is based on the behaviour Gradient Descent
% in an euclidean space.
alpha = 1/(sqrt(2)*norm(M, 'Fro'));

% Applying Riemannian Gradient Descent
[x, iteratesX,iteratesY,conv] = RGD(M, x0, y0, alpha, 100);

% Plotting the function on the square [0,2 \pi]^2
f = @(x,y) [cos(x) sin(x)]*M*[cos(y); sin(y)];
[X,Y] = meshgrid(0:(pi/20):2*pi,0:(pi/20):2*pi);
Z = zeros(size(X));

for i = 1:41
   for j = 1:41
       Z(i,j) = f(X(i,j),Y(i,j)); 
   end 
end

surf(X,Y,Z);
xlabel('x')
ylabel('y')
zlabel('f')
colorbar
