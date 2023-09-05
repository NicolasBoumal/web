%% This is the main for the exercise Implementation of RTR
% Output : Results of Riemannian Trust Regions with CauchyStep or tCG for
% the Rayleigh quotient on the sphere.

clear;
clc;
close all;
format long;

% Setup of the parameteres for the problem
% min f on the d-1 dimensional sphere with f(x) = x'Ax.
% We see the sphere as a Riemannian subamnigold of R^d.
% We have that gradf(x) = Proj_x(2Ax) = 2 (A- x'Ax)x 
% with Proj_x(v) = (I - xx')v.
% Riemannian Hessianf(x)[v] = 2((I-xx')Av - x'Axv)
params.d = 10;
params.xinit = randn(params.d,1);
params.xinit = params.xinit / norm(params.xinit);
params.dbar = sqrt(2);
params.dinit = params.dbar / 8;
params.epsilon = 10^(-6);
params.max_iter = 5000;
params.max_iter_subprob = 1000;
params.rhom = 0.1;
params.scal = @(x,u,v) u'*v;
params.norm = @(x,u) sqrt(params.scal(x,u,u));
params.dispintstep = false;

% Function handle
A = randn(params.d,params.d);
A = A'*A;
fhandle = @(x,v,p) RayQuot(A,x,v,p);

% Subproblem
subprob = @(x,delta) CauchyStep(fhandle,x,delta,params);

% Retraction
retraction = @(x,v) (x + v) / norm(x + v);

% Model
model = @(x,v) modelQ(x,v,fhandle,params);

% Run trust regions
TrustRegions(fhandle,subprob,retraction,model,params);

% Disp minimum of f(x) = x'Ax on the d-1 dimensional sphere
fprintf('Disp minimum of f on the d-1 dimensional sphere')
min(eig(A))

