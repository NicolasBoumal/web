%% Cauchy Step
% This is the code for the Cauchy step, which is a solver for the subproblem
% in RTR.          
% Input : - fhandle, which takes as input (x,v) in TM and computes function 
%           value and Riemannian gradient at x and Riemannian hessian at x.
%           in direction v.
%         - A point x.
%         - A redius delta.
%         - params.
%
% Output : - Approximative solution "s" to the subproblem.

function s = CauchyStep(fhandle,x,delta,params)

if params.dispintstep
fprintf('Subproblem solver Cauchy step\n')
end

% Computation of gradient at x
p.f = false; p.g = true; p.H = false;
[~,g,~] = fhandle(x,zeros(size(x)),p);

% Computation of hessian at x in direction g
p.f = false; p.g = false; p.H = true;
[~,~,H] = fhandle(x,g,p);

% Computation of (g,H_x(g))_x and norm(g)
scal = params.scal(x,g,H);
norm = params.norm(x,g);

% Computation Cauchy step
t = delta / norm;
if scal > 0
    t = min(norm^2 / scal,delta/norm);
end

s = -t*g;

end