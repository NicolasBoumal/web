%% Riemannian trust regions
% Input : - fhandle, which takes as input (x,v) in TM and computes function 
%           value and Riemannian gradient at x and Riemannian hessian at x
%           in direction v.
%         - subproblem, which is a solver for the subproblem like Cauchy
%           step or tCG. It takes as input fhandle, x, delta and params and gives
%           the approximat solution "s" to the subproblem
%         - retraction, which takes as input (x,v) in TM and gives R_x(v)
%           for a given retraction R.
%         - model, which takes as input (x,v) in TM and params and computes the model
%           value at x in direction v. This model is 
%           m_x(v) = f(x) + (gradf(x),v)_x + 0.5 * (H(x)[v],v)_x.
%         - params, which containes different parameters for RTR.
%
% Output : - Approximative solution x
%          - Final cost value f
%          - Final gradient value g
%          - Final number of iterations iter

function [x,f,g,iter] = TrustRegions(fhandle,subprob,retraction,model,params)

% Initialization of parameters
x = params.xinit;
delta = params.dinit;
delta_bar = params.dbar;
rho_min = params.rhom;
epsilon = params.epsilon;
max_iter = params.max_iter;
iter = 1;

% Loop for riemannian trust region
fprintf('Begin of loop for trust regions\n')
for k = 1:max_iter

    % Computation of function value and Riemannian gradient at x 
    p.f = true; p.g = true; p.H = false;
    [f,g,~] = fhandle(x,0,p);

    % Printing actual values
    if params.dispintstep
        fprintf("Iteration number ")
        disp(iter)
        fprintf("Value of the cost\n")
        disp(f);
        fprintf("Value of the gradient norm\n")
        disp(norm(g));
    end
    
    % Checking for convergence of the algorithm
    if norm(g) < epsilon
        fprintf("The algorithm did converge the results are the following\n")
        fprintf("Final cost value\n")
        disp(f);
        fprintf("Final gradient norm\n")
        disp(norm(g));
        fprintf("Final position\n")
        disp(x);
        fprintf("Number of total iteration\n")
        disp(iter)
        break;
    end
    
    % Solving the Subproblem
    s = subprob(x,delta);
    x_plus = retraction(x,s);
    p.f = true; p.g = false; p.H = false;
    [fint,~,~] = fhandle(x_plus,zeros(size(x)),p);
    f_plus = f - fint;
    m_plus =  model(x,zeros(size(x))) - model(x,s);
    
    % !! Could be problematic if both f_plus and m_plus are small !!
    rho = f_plus / m_plus;
    
    % Compution next step
    if rho > rho_min
        x = x_plus;
    end
    
    % Computation next delta
    if rho < 1/4
        delta = delta/4;
    elseif rho > 3/4 && abs(norm(s)-delta) < epsilon
        delta = min(2*delta,delta_bar);
    end
    
    % Next iter
    iter = iter + 1;

end

% Output of the results if the algorithm didn't converge
if norm(g) > epsilon
    fprintf("The algorithm didn't converge the results are the following\n")
    fprintf("Final cost value\n")
    disp(f);
    fprintf("Final gradient norm\n")
    disp(norm(g));
    fprintf("Final position\n")
    disp(x);
    fprintf("Number of total iteration\n")
    disp(max_iter)
end

end