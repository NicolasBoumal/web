%% MATLAB main for the exercise The Riemannian Newton method
% Output : The different parts of the exercise "The Riemannian Newton method"
% ATTENTION : Make sure to carefully read the code to understand the
% different computations.
% ATTENTION : Parts of this code use ManOpt. You need to have ManOpt to
% run the code. Also look up how ManOpt handels tangent vectors for SO(d).

% Initialization
close all;
clear;
clc;
format long

% Run different parts of the exercise
% ATTENTION : Do not run everything at ones. There are different variable
% names reused in the different parts. Running everything at ones may not
% work.
question4 = false;
checkGradHess = false;
question6 = false;
question7 = false;
question8 = false;

% Frob functions
frob = @(U,V) U(:)'*V(:);
frobsq = @(U) U(:)'*U(:);

%% Question 4 generation of an basis by projection
% In order to analyse if we get a basis we comput two tensors T containing
% the random created matrices and PT containing the projected onto the 
% tangent space of Q random matrices, i.e., the projection of the matrices
% in T. After that we vectorize the matrices and put the vectors of T in MT
% and the vectorize the matrices of PT in MPT. The rank of this two
% matrices is then the dimension of the space spaned by the matrices of T,
% respectively PT. If both have the value d(d-1)/2 we have a basis for the
% tangent space.

if question4 == true

fprintf('Run question 4\n')

% Parameters
d = 10;
n = d*(d-1) / 2;
T = randn(d,d,n);
PT = zeros(d,d,n);
[Q,~] = qr(randn(d,d));
Q = Q^2;

% Tensor of projected matrices
for i=1:n
    PT(:,:,i) = Proj(Q,T(:,:,i));
end

% Check that the matrices are in the tangent space
QoPT = mumult(PT,Q',1);
s = zeros(n,1);
for i = 1:n
    s(i,1) = norm(QoPT(d,d,i) + QoPT(d,d,i)');
end
fprintf('Check that the matrices are in the tangent space\n')
fprintf('Should be close to zero')
norm(s)

% Check linear independence
MT = zeros(d*d,n);
MPT = zeros(d*d,n);
for i = 1:n
MT(:,i) = reshape(T(:,:,i),[],1);
MPT(:,i) = reshape(PT(:,:,i),[],1);
end

fprintf('numerical ranks of MT:')
rank(MT)
fprintf('numerical ranks of MPT:')
rank(MPT)

end

%% Check GradHess with Manopt
% This uses ManOpt to check that we have the correct gradient and hessian
% for f.

if checkGradHess == true

fprintf('Run checkGradHess\n')

% Parameters
d = 5;
dim = d*(d-1)/2;
A = randn(d,d);
A = A' + A;
B = randn(d,d);
B = B' + B;

% Verification
problem.M = rotationsfactory(d, 1);
problem.cost = @(Q) frobsq(A*Q-Q*B);
problem.grad = @(Q) Q' * (- 2*Q*(Q'*A*Q*B-B*Q'*A*Q)); 
problem.hess = @(Q,U) Q' * Proj(Q,2*(Q*B*Q'*A*Q*U+Q*B*(Q*U)'*A*Q+Q*U*B*Q'*A*Q - A*Q*U*B));

fprintf('Check Gradient:\n')
figure
checkgradient(problem)
fprintf('Check Hessian:\n')
figure
checkhessian(problem)

end

%% Question 6
% Solves the Newton system using the backslash symbol from Matlab. The basis, which 
% is used to rewrite the Newton system as matrix vector equation is the one
% from question 5.
% The matrix Q \in SO(d) in the Newton system, which is solved is compued by
% the trustregions method from ManOpt. The reason why Q isn't an arbitary
% (or random) matrix in SO(d) is that Hess f(Q)[X] = grad f(Q) must be
% solvable. This is the case if Hess f(Q) is positive definite, which
% should be the case close to a minimizer.

if question6 == true

fprintf('Run question 6\n')

% Parameters
d = 50;
dim = d*(d-1)/2;
A = randn(d,d);
A = A' + A;
B = randn(d,d);
B = B' + B;

problem.M = rotationsfactory(d, 1);
problem.cost = @(Q) frobsq(A*Q-Q*B);
problem.grad = @(Q) Q' * (- 2*Q*(Q'*A*Q*B-B*Q'*A*Q)); 
problem.hess = @(Q,U) Q' * Proj(Q,2*(Q*B*Q'*A*Q*U+Q*B*(Q*U)'*A*Q+Q*U*B*Q'*A*Q - A*Q*U*B));
Q = trustregions(problem);

% Creatin of orthonormal basis, using the results of question 5
tic
U = zeros(d,d,dim);
iter = 1;
for j=2:d
    for i=1:(j-1)
        U(i,j,iter) = 1;
        U(j,i,iter) = -1;
        iter = iter + 1;
    end
end
U = mumult(U,Q,1);
U = U / sqrt(2);
t = toc;

% Check basis
QoU = mumult(U,Q',1);
s = zeros(dim,1);
for i = 1:dim
    s(i,1) = norm(QoU(d,d,i) + QoU(d,d,i)');
end
fprintf('Check that the matrices are in the tangent space\n')
fprintf('Should be close to zero\n')
norm(s)

for i = 1:dim
MU(:,i) = reshape(U(:,:,i),[],1);
end
fprintf('Numerical ranks of MT should equal dim(SO(d))')
rank(MU)

% Computing matrix H and vector g from theory
tic;
G = Grad(A,B,Q);
g = zeros(dim,1);
H = zeros(dim,dim);
for i=1:dim
    g(i,1) = frob(U(:,:,i),G);
end

for j=1:dim
    He = Hess(A,B,Q,U(:,:,j));
    for i=1:dim
        H(i,j) = frob(U(:,:,i),He);
    end
end

% Solving Newton system
fprintf('Solution of Newton system for good choosen Q in SO(d)\n')
fprintf('By matrix vector equation using backslash\n')

% Computing of solution to the Newton system
x = H \ (-g);
X = zeros(d,d);
for i=1:dim
    X = x(i)*U(:,:,i) + X;
end
t = t + toc;

% Control of outcome 
G = Grad(A,B,Q);
H = Hess(A,B,Q,X);
fprintf('Difference between Hess(Q)[X]+Grad(Q)\n')
fprintf('Should be close to zero\n')
sqrt(frobsq(H+G))

% Computation time
fprintf('Computation time\n')
disp(t)

end

%% Question 7
% Solves the Newton system using Gradient Descent.
% The matrix Q \in SO(d) in the Newton system, which is solved is compued by
% the trustregions method from ManOpt. The reason why Q isn't an arbitary
% (or random) matrix in SO(d) is that Hess f(Q)[X] = grad f(Q) must be
% solvable and the GD algorithm solving this equation requires Hess f(Q) to
% be symmetric positive definite in order to work. We should have that
% Hess f(Q) is positive definite close to a minimizer.

if question7 == true

fprintf('Run question 7\n')

% parameters
d = 10;
dim = d*(d-1)/2;
A = randn(d,d);
A = A' + A;
B = randn(d,d);
B = B' + B;
tol = 10^(-4);
iter = 1000;

problem.M = rotationsfactory(d, 1);
problem.cost = @(Q) frobsq(A*Q-Q*B);
problem.grad = @(Q) Q' * (- 2*Q*(Q'*A*Q*B-B*Q'*A*Q)); 
problem.hess = @(Q,U) Q' * Proj(Q,2*(Q*B*Q'*A*Q*U+Q*B*(Q*U)'*A*Q+Q*U*B*Q'*A*Q - A*Q*U*B));
Q = trustregions(problem);

% Solving Newton system
fprintf('Solution of Newton system for good choosen Q in SO(d)\n')
fprintf('By Gradient Descent\n')

% !!! ATTENTION algorithm may not work if g.H isn't s.p.d. !!!
g.b = - Grad(A,B,Q);
g.H = @(U) Hess(A,B,Q,U);
tic
X = GD(g,tol,iter);
t = toc;

% Control of outcome 
G = Grad(A,B,Q);
H = Hess(A,B,Q,X);
fprintf('Difference between Hess(Q)[X]+Grad(Q)\n')
fprintf('Should be close to zero\n')
sqrt(frobsq(H+G))

% Computation time
fprintf('Computation time\n')
disp(t)

end

%% Question 8
% Solves the Newton system using Conjugate Gradient.
% The matrix Q \in SO(d) in the Newton system, which is solved is compued by
% the trustregions method from ManOpt. The reason why Q isn't an arbitary
% (or random) matrix in SO(d) is that Hess f(Q)[X] = grad f(Q) must be
% solvable and the CG algorithm solving this equation requires Hess f(Q) to
% be symmetric positive definite in order to work. We should have that
% Hess f(Q) is positive definite close to a minimizer.

if question8 == true

fprintf('Run question 8\n')

% Parameters
d = 50;
dim = d*(d-1)/2;
A = randn(d,d);
A = A' + A;
B = randn(d,d);
B = B' + B;
tol = 10^(-4);
iter = 1000;

problem.M = rotationsfactory(d, 1);
problem.cost = @(Q) frobsq(A*Q-Q*B);
problem.grad = @(Q) Q' * (- 2*Q*(Q'*A*Q*B-B*Q'*A*Q)); 
problem.hess = @(Q,U) Q' * Proj(Q,2*(Q*B*Q'*A*Q*U+Q*B*(Q*U)'*A*Q+Q*U*B*Q'*A*Q - A*Q*U*B));
Q = trustregions(problem);

% Solving Newton system
fprintf('Solution of Newton system for good choosen Q in SO(d)\n')
fprintf('By Conjugate Gradient\n')

% !!! ATTENTION algorithm may not work if g.H isn't s.p.d. !!!
g.b = - Grad(A,B,Q);
g.H = @(U) Hess(A,B,Q,U);
tic
X = CG(g,tol,iter);
t = toc;

% Control of outcome 
G = Grad(A,B,Q);
H = Hess(A,B,Q,X);
fprintf('Difference between Hess(Q)[X]+Grad(Q)\n')
fprintf('Should be close to zero\n')
sqrt(frobsq(H+G))

% Computation time
fprintf('Computation time\n')
disp(t)

end
