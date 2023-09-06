%% Implementation of Riemannian Gradient Descent for 
%% the exercise Rayleigh quotient on the sphere.
% Input : - Symmetric matrix A
%         - Starting position x0
%         - Fixed step length alpha
%         - Maximum number of itrations
% Output : - Final position x
%          - List of iterates
%          - conv, which is true if it did converge false if not
% The algorithm stops if difference between x'Ax and \lambda_min(A) is less
% than \epsilon = 10^(-6).
 
function [x, iterates,conv] = RGDsphere02(A, x0, alpha, maxiter)

% Initalization
conv = false;
x = x0;
iterates = zeros(length(x), maxiter);
epsilon = 1e-6;
mineigA = min(eig(A));
iterates(:, 1) = x;

% RGD loop
for k = 1:(maxiter-1)
    Ax = A*x; 
    v = -alpha*2*(Ax - x'*Ax*x); % v = -alpha * Proj_xk(gradf(xk))
    x = (x+v)/norm(x+v); % x_(k+1) = R_xk(v)
    iterates(:, k+1) = x;

    % Control of convergence
    if abs(mineigA - x'*A*x) < epsilon
        conv = true;
        iterates = iterates(:,1:(k+1));
        fprintf('Riemannian Gradient Descent did converge\n')
        fprintf('Final position:\n')
        disp(x)
        fprintf('Final function value:\n')
        disp(x'*A*x)
        fprintf('Number of iterations:\n')
        disp(k)
        break;
    end
end

if conv == false
    fprintf('Riemannian Gradient Descent did not converge\n')
    fprintf('Final position:\n')
    disp(x)
    fprintf('Final function value:\n')
    disp(x'*A*x)
end

end
