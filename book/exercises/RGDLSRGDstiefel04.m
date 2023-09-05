%% Implementation of Riemannian Gradient Descent for 
%% the exercise RGD on Stiefel.
% Input : - Symmetric matrix A \in \R^{n \times n}
%         - Starting position X0 \in St(n,p)
%         - LS true/false
%         - Fixed step length alpha
%         - Maximum number of itrations
% Output : - Final position X
%          - List of iterates for X
%          - conv, which is true if it did converge false if not
% The algorithm stops if ||\grad f(X)|| < 10^(-4).
 
function [X,iterates,conv] = RGDLS(A, X0, LS, alpha, maxiter)

% Initalization
conv = false;
X = X0;
s = size(X0);
iterates = zeros(s(1,1),s(1,2),maxiter);
epsilon = 10^(-4);
iterates(:,:,1) = X;

% RGD loop
for k = 1:(maxiter-1)
    AX = A*X;
    gradfX = 2 * (AX - X*X'*AX);
    if LS
        r = 10^(-4);
        tau = 1/2;
        alpha = 1 / (2 * norm(A,'fro'));
        V = -alpha * gradfX;
        [U,~,W] = svd(X+V,0);
        Xtry = U*W';

        while (trace(X'*A*X) - trace(Xtry'*A*Xtry)) < (r*alpha*norm(gradfX,'fro')^2)
            alpha = tau*alpha;
        end
        X = Xtry;
    else
        V = -alpha * gradfX;
        [U,~,W] = svd(X+V,0);
        X = U*W';
    end

    % Control of convergence
    if norm(gradfX,'fro') < epsilon
        conv = true;
        fprintf('Riemannian Gradient Descent did converge\n')
        if LS
            fprintf('Line search was used\n')
        else
            fprintf('Line search was not used\n')
        end
        fprintf('Final position X:\n')
        disp(X)
        fprintf('Final function value:\n')
        disp(trace(X'*A*X))
        fprintf('Number of iterations:\n')
        disp(k)
        break;
    end
end

if conv == false
    fprintf('Riemannian Gradient Descent did not converge\n')
    if LS
        fprintf('Line search was used\n')
    else
        fprintf('Line search was not used\n')
    end
    fprintf('Final position X:\n')
    disp(X)
    fprintf('Final function value:\n')
    disp(trace(X'*A*X))
end

end