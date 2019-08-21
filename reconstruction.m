function [ans,error] = reconstruction(y,A,A_adj,info,varargin)
% Image reconstruction algorithm based on FISTA/FISTA+ADMM
%   solve argmin 1/2*||y-Ax||^2+tau*psi(x)

%   Input: y: the observed sensor data
%          A & A_adj: the function handler of measurement A and the 
%                     adjoint operator of A
%          info: parameters for FISTA algorithm
%               info.a: stepsize
%               info.tau: penalty parameter
%               info.radius: the Schatten norm: 1, 2 or Inf
%               info.Nx: the amount of rows of original image
%               info.Ny: the amount of columns of original image
%               info.maxIter: the maximal iteration number
%               info.tol: the tolerance
%          varargin: sigma: the Gaussian filter's standard deviation
%
%   Output: ans: the reconstructed image
%           error: The normstep at every iteration


% initialization
a = info.a;
radius = info.radius;
tau = info.tau;
m = info.Nx;
n = info.Ny;
type = [];
error = [];

% obtain the number of sigma
if nargin > 4
    sigma = cell2mat(varargin);
    if size(sigma,2) == 1
        type = 'FISTA';
    elseif size(sigma,2) == 2
        type = 'FADMM';
    end
end

x_n = 0.01*ones(m,n);
u_n = x_n;
t = 1;
tol = info.tol;
maxIter = info.maxIter;

StopCond = false;
iter = 0;
while StopCond ~= true
    z = u_n + a^(-1) * A_adj(y - A(u_n));
    x_n_l = x_n;
    switch type
        case 'FISTA'
            x_n = denoise(z,tau/a,radius,sigma);
            
        case 'FADMM'
            x_n = ADMM(z,sigma(1),sigma(2),info.tau1,info.tau2,radius,info.rho,a);
    end
    t_l = t;
    t = (1 + sqrt(1+4*t_l*t_l))/2;
    u_n = x_n + (t_l - 1) / t * (x_n - x_n_l);
    
    iter = iter + 1;
    
    normStep = norm(x_n - x_n_l)/norm(x_n);
    error = [error,normStep];
    
    if iter >= maxIter
        break
    end
    if normStep <= tol
        break
    end   
end
x_n(x_n < 0) = 0;
ans = x_n;
end
