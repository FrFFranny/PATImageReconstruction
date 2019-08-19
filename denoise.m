function [ans] = denoise(z,tau,radius,sigma)
% To solve the problem
%     argmin 1/2*||x-z||^2+tau*||HGx||_1,p for any p>=1  
%
%  Input:
%       z: the noise image
%       tau: the penalty parameter
%       radius: the Schatten norm 1, 2 or Inf
%       sigma: the standard deviation of Gaussian filter
% Output: ans: the optimal x


% initialization 
stopCond = false;
t = 1;
[m,n]=size(z);
maxIter = 20;

% obtain the number of pixels
N = m*n;

omega = zeros(2,2,N);
psi = zeros(2,2,N);
iter = 0;
while stopCond ~= true
    v = z - tau * adjmulti(imgaussfilt(psi,sigma),m,n);
    omega_l = omega;
    temp = imgaussfilt(projection(v),sigma);
    omega = normalization(psi + 1/(64*tau) * Hessian(temp),radius);
    t_l = t;
    t = (1 + sqrt(1 + 4*t_l*t_l))/2;
    psi = omega + (t_l - 1) / t * (omega - omega_l);
    iter = iter + 1;
    if iter >= maxIter
        break;
    end
end
v = z - tau * adjmulti(imgaussfilt(omega_l,sigma),m,n);
ans = projection(v);
end

function [ans] = normalization(omega,radius)
% Calculate the Schatten norm projections for q = 1, 2, Inf
[~,~,N]=size(omega);
ans=[];
switch radius
    case '1'
        for i=1:N
            [U,S,V]= svd(omega(:,:,i));
            sigma = diag(S);
            if sigma(1)<=1-sigma(2)
                gamma = 0;
            elseif sigma(1)>1-sigma(2)
                gamma = sigma(1)-1;
            else
                gamma = (sigma(1)+sigma(2)-1)/2;
            end
            S_new = max(sigma-gamma,zeros(size(sigma)));
            ans(:,:,i) = U*diag(S_new)*V';
        end
        
    case '2'
        for i=1:N
            FroNorm = norm(omega(:,:,i),'fro');
            if FroNorm > 1
                ans(:,:,i) = omega(:,:,i)/FroNorm;
            else
                ans(:,:,i) = omega(:,:,i);
            end
        end
        
        
    case 'Inf'
        for i=1:N
            [U,S,V]= svd(omega(:,:,i));
            sigma = min(diag(S),ones(size(diag(S))));
            ans(:,:,i) = U*diag(sigma)*V';
        end
end
end

function [Hx] = projection(Hx)
% Calculation the projection onto the convex set C = {x = R^N|x_n belongs to [0,1] for any n = 1, ..., N}
Hx(Hx<0) = 0;
Hx(Hx>1) = 1;
end

function [ans] = adjmulti(Y,m,n)
% Calculate the adjoiont operators that correspond to backward operators
% with Neumann boundary conditions

[~,~,N]=size(Y);
% obtain the vectors for every entries of the 2X2 matrix Y_n where n = 1,...,N
Y11 = Y(1,1,:);
Y12 = Y(1,2,:);
Y21 = Y(2,1,:);
Y22 = Y(2,2,:);

% reshape the vector to a (Nx X Ny) matrix
for i=1:N
    col = mod(i,n);
    if col == 0
        col = n;
    end
    row = fix(i/n)+1;
    if row == m+1
        row = m;
    end
    Y11_reshape(row,col)=Y11(i);
    Y12_reshape(row,col)=Y12(i);
    Y21_reshape(row,col)=Y21(i);
    Y22_reshape(row,col)=Y22(i);
end
% Get the second order backward difference operators
[Ny,Nx] = size(Y11_reshape);
dydy = [Y11_reshape(2,:)-Y11_reshape(1,:); Y11_reshape(2,:)-Y11_reshape(1,:); -diff(-diff(Y11_reshape))];
dxdx = [Y22_reshape(:,2)-Y22_reshape(:,1),Y22_reshape(:,2)-Y22_reshape(:,1), -diff(-diff(Y22_reshape,1,2),1,2)];
dxdy1 = [zeros(1,Nx-1);-diff(-diff(Y12_reshape),1,2)];
dxdy1 = [zeros(Ny,1),dxdy1];
dxdy2 = [zeros(1,Nx-1);-diff(-diff(Y21_reshape),1,2)];
dxdy2 = [zeros(Ny,1),dxdy2];
ans = dydy + dxdx + dxdy1 + dxdy2;
end
