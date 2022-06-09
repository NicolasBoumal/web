function [X, optval, status] = linearcost_cvx(C, m, d)
% function [X, optval, status] = linearcost_cvx(C, m, d)
%
% Uses CVX (with the solver set up for CVX) to solve the problem
%
%   min Trace(CX)
% 
%     such that X is symmetric positive semidefinite of size nxn (n = m*d)
%     and the m diagonal blocks of X of size dxd each equal eye(d).
%
%     C must be a real matrix of size nxn.
%     Only the symmetric part of C matters.
%
% To choose between SeDuMi, SDPT3 or other CVX solvers, call for example:
% cvx_solver sedumi;
% cvx_solver sdpt3;
% before calling this function. More generally, this function uses the
% current CVX options. So you can also choose precision by calling:
% cvx_precision ...
% where ... is either low, medium, default, high or best.
% See http://cvxr.com/cvx/doc/solver.html.
% 
% Nicolas Boumal, UCLouvain, May 13, 2014.

    n = m*d; %#ok<NASGU>
    
    cvx_begin
    
        variable X(n, n) symmetric;
        
        minimize( C(:)'*X(:) ); %#ok<NODEF>
        
        for i = 1 : m
            
            ii = (i-1)*d + (1:d);
            X(ii, ii) == eye(d); %#ok<EQEFF>
            
        end
        
        X == semidefinite(n);
        
    cvx_end
    
    optval = cvx_optval;
    status = cvx_status;
    
end
