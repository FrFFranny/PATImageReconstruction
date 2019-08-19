function Hx = Hessian(image)
% Calculate the Hessian matrix for a 2D image
[Nx,Ny] = size(image);
dxdx = [diff(diff(image));image(Nx-1,:)-image(Nx,:);image(Nx-1,:)-image(Nx,:)];
dydy = [diff(diff(image,1,2),1,2),image(:,Ny-1)-image(:,Ny),image(:,Ny-1)-image(:,Ny)];
dxdy = [diff(diff(image),1,2),zeros(Nx-1,1)];
dxdy = [dxdy;zeros(1,Ny)];
for i = 1:Nx
    for j = 1:Ny
        n = (i-1) * Ny + j;
        Hx(:,:,n) = [dxdx(i,j), dxdy(i,j); dxdy(i,j), dydy(i,j)];
    end
end
end

