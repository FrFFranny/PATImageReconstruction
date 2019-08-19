function [x_Min] = ADMM(z,sig_1,sig_2,tau1,tau2,radius,rho,a)
% solve min_x f(x)+f(z) subject to x = z = [z1; z2]

% Input:
%   sig_1 & sig_2: two different sigma for gaussian filter
%   tau1 & tau2: two penalty parameters for ||HGx||_(s.p)
%   radius: the radius for schattern norm
%   rho: the penalty parameters for ADMM
%   a: stepsize
% Output: the optimal x

%obtain the size of image
[m,n] = size(z); 

% The function handler to use Gaussian filter
z1 = @(x) imgaussfilt(x, sig_1);
z2 = @(x) imgaussfilt(x, sig_2);

% initialization
x_k = z(:);
z1_k = zeros(m*n,1);
z2_k = zeros(m*n,1);
u1_k = zeros(m*n,1);
u2_k = zeros(m*n,1);
w1 = rho*(x_k-z1_k);
w2 = rho*(x_k-z2_k);

tol = 1e-4;
maxIter = 50;
stopCon = false;
iter = 0;

% create the matrix to update the x_k
i = 1:m*n;
j = i;
v = ones(m*n,1);
temp1 = sparse(i,j,v);
A = [temp1;sqrt(rho)*temp1;sqrt(rho)*temp1];
b = @(z1,z2,u1,u2) [z(:); sqrt(rho)*(z1-u1);sqrt(rho)*(z2-u2)];

while (~stopCon && iter<=maxIter)
    % assign the last value
    x_k_l = x_k;
    z1_k_l = z1_k;
    z2_k_l = z2_k;
    u1_k_l = u1_k;
    u2_k_l = u2_k;
    
    % function handler for calculating the vector b based on last z1,z2,u1,u2
    b_vec = b(z1_k_l,z2_k_l,u1_k_l,u2_k_l);
    
    % find the minimal value of F
    x_k = (A'*A)\A'*b_vec;
    
    % reshape the image vector to matrix
    temp1 = reshape(x_k+u1_k_l,m,n); 
    temp2 = reshape(x_k+u2_k_l,m,n);
    
    % find the minimal value of z1,z2
    z1_k = denoise(temp1,tau1/(a*rho),radius,sig_1);
    z2_k = denoise(temp2,tau2/(a*rho),radius,sig_2);
    z1_k = z1_k(:); %vectorize
    z2_k = z2_k(:); %vectorize
    
    % update the u1 u2
    residual1 = x_k - z1_k;
    residual2 = x_k - z2_k;
    u1_k = u1_k_l + residual1;
    u2_k = u2_k_l + residual2;
    
    % update the dual variables
    w1 = w1 + rho*residual1;
    w2 = w2 + rho*residual2;
    
    % obtain the residuals
    r_k = [residual1;residual2]; 
    % dual residual
    s_k = -rho*[z1_k-z1_k_l;z2_k-z2_k_l];
    norm_r_k = norm(r_k);
    norm_s_k = norm(s_k);
    
    iter = iter+1;
    
    if iter > maxIter
        break
    end
    
    % Tolerence
    tol_primal = sqrt(2*length(x_k)) * tol + tol * max(norm(x_k),norm(-[z1_k;z2_k]));
    tol_dual = sqrt(length(r_k)) * tol + tol * norm([w1;w2]);
    % stop when residuals converge
    stopCon = (norm_r_k<tol_primal && norm_s_k<tol_dual);
    
end
    

x_Min = reshape(x_k,m,n);

end
