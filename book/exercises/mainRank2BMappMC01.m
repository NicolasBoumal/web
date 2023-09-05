%% This is the main for the exercise Rank-2 Burer-Monteiro approach to Max-Cut
% Output : Results of Riemannian Trust Regions with CauchyStep or tCG for
% the Rank-2 Burer-Monteiro approach to Max-Cut.

clear;
clc;
close all;
format long;

% Setup of the parameteres for the problem
params.n = 10;
params.p = 2;
params.xinit = randn(params.n,params.p);
for i = 1:params.n
    params.xinit(i,:) = params.xinit(i,:) / norm(params.xinit(i,:));
end
params.dbar = sqrt(2);
params.dinit = params.dbar / 8;
params.epsilon = 10^(-5);
params.max_iter = 10000;
params.max_iter_subprob = 1000;
params.rhom = 0.1;
params.scal = @(x,u,v) u(:)'*v(:);
params.norm = @(x,u) sqrt(params.scal(x,u,u));
params.dispintstep = false;

% Function handle
C = randn(params.n,params.n);
C = (C + C') / 2;
fhandle = @(Y,V,p) fmaxcut(C,Y,V,p,params);

% Subproblem
subprob = @(x,delta) tCG(fhandle,x,delta,params);

% Retraction
retraction = @(Y,V) RetobMan(Y,V,params);

% Model
model = @(x,v) modelQ(x,v,fhandle,params);

% Run trust regions
TrustRegions(fhandle,subprob,retraction,model,params);
