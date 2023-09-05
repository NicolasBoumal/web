%% Implementation of Riemannian Gradient Descent for 
%% the exercise RGD on products of spheres.
% Input : - Matrix M \in \R^{m \times n}
%         - Starting position (x0,y0)
%         - Fixed step length alpha
%         - Maximum number of itrations
% Output : - Final position (x,y)
%          - List of iterates for x and y (iteratesX,iteratesY)
%          - conv, which is true if it did converge false if not
% The algorithm stops if difference between x'My and \sigma_max(M) is less
% than \epsilon = 10^(-6).
 
function [x, iteratesX,iteratesY,conv] = RGD(M, x0, y0, alpha, maxiter)

% Initalization
conv = false;
x = x0;
y = y0;
iteratesX = zeros(length(x), maxiter);
iteratesY = zeros(length(y), maxiter);
epsilon = 10^(-6);
maxsingM = norm(M,2);
iteratesX(:, 1) = x;
iteratesY(:, 1) = y;

% RGD loop
for k = 1:(maxiter-1)
    My = M*y;
    Mtx = M'*x;
    f = x'*M*y;
    % v = -alpha * Proj_xk(gradf(xk))
    vx = alpha * (My - f*x);
    vy = alpha * (Mtx - f*y);
    % x_(k+1) = R_xk(v)
    x = (x+vx)/norm(x+vx);
    y = (y+vy)/norm(y+vy);
    iteratesX(:, k+1) = x;
    iteratesY(:, k+1) = y;

    % Control of convergence
    if abs(maxsingM - x'*M*y) < epsilon
        conv = true;
        iteratesX = iteratesX(:,1:(k+1));
        iteratesY = iteratesY(:,1:(k+1));
        fprintf('Riemannian Gradient Descent did converge\n')
        fprintf('Final position x:\n')
        disp(x)
        fprintf('Angle x:\n')
        gamma = acos(x(1,1));
        if x(2,1) < 0
            gamma = 2*pi - gamma;
        end
        disp(gamma)
        fprintf('Final position y:\n')
        disp(y)
        fprintf('Angle y:\n')
        gamma = acos(y(1,1));
        if y(2,1) < 0
            gamma = 2*pi - gamma;
        end
        disp(gamma)
        fprintf('Final function value:\n')
        disp(x'*M*y)
        fprintf('Number of iterations:\n')
        disp(k)
        break;
    end
end

if conv == false
    fprintf('Riemannian Gradient Descent did not converge\n')
    fprintf('Final position x:\n')
        disp(x)
        fprintf('Angle x:\n')
        gamma = acos(x(1,1));
        if x(2,1) < 0
            gamma = 2*pi - gamma;
        end
        disp(gamma)
        fprintf('Final position y:\n')
        disp(y)
        fprintf('Angle y:\n')
        gamma = acos(y(1,1));
        if y(2,1) < 0
            gamma = 2*pi - gamma;
        end
        disp(gamma)
        fprintf('Final function value:\n')
        disp(x'*M*y)
end

end