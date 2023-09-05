%% Truncated Conjugat Gradient
% This is the code for the truncated Conjugat Gradient, which is a solver for the subproblem
% in RTR.          
% Input : - fhandle, which takes as input (x,v) in TM and computes function 
%           value and Riemannian gradient at x and Riemannian hessian at x.
%           in direction v.
%         - A point x.
%         - A redius delta.
%         - params.
%
% Output : - Approximative solution "s" to the subproblem.

function v = tCG(fhandle,x,delta,params)

if params.dispintstep
fprintf('Subproblem solver tCG\n')
end

% Initialization
v = zeros(size(x));
numiter = 0;
c.f = false; c.g = true; c.H = false;
[~,b,~] = fhandle(x,zeros(size(x)),c);
b = -b;
r = b;
p = r;

for i=1:params.max_iter_subprob

    numiter = numiter + 1;

    % Computation of hessian at x in direction p
    c.f = false; c.g = false; c.H = true;
    [~,~,H] = fhandle(x,p,c);

    % Computation of (g,H_x(g))_x
    scal = params.scal(x,p,H);

    % Computation of alpha
    alpha = params.norm(x,r)^2 / scal;

    % Computation attempt next v
    v_plus = v + alpha * p;

    % Checking for extrem values 
    if (scal <= 0) || (params.norm(x,v_plus) >= delta)
        D = params.scal(x,v,p)^2 + params.norm(x,p)^2 * (delta^2 - params.norm(x,v)^2);
        t = (- params.scal(x,v,p) + sqrt(D)) / params.norm(x,p)^2;

        % output v + tp
        v = v + t * p;
        break;
    end

    % Computation next v
    v = v_plus;

    r_min = r;
    r = r - alpha * H;

    % Checking for convergence

    if params.norm(x,r)<= params.norm(x,b) * min(params.norm(x,b),0.1)
        break;
    end

    beta = params.norm(x,r)^2 / params.norm(x,r_min)^2;
    p = r + beta * p;
end

if params.dispintstep
fprintf('Number of iteration in tCG: ')
disp(numiter)
end

end