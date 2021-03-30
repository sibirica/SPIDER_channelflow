function [U, V, W] = pullVelocity(x, y, z, t)

authkey = 'edu.jhu.meneveau-hiSmxkae';
dataset = 'Channel';
Lag4 = 'Lag4';          % 4th order Lagrangian interpolation in space
PCHIPInt = 'PCHIP';     % Piecewise cubic Hermit interpolation in time

nx = length(x);
ny = length(y);
nz = length(z);
nt = length(t);
npoints = nx*ny*nz;

[X, Y, Z] = ndgrid(x, y, z);
points = zeros(3, npoints);
points(1, :) = X(:)';
points(2, :) = Y(:)';
points(3, :) = Z(:)';

U = zeros(nx, ny, nz, nt);
V = zeros(nx, ny, nz, nt);
W = zeros(nx, ny, nz, nt);

for n = 1:nt

    time = t(n);
    result3 = getVelocity(authkey, dataset, time, Lag4, PCHIPInt, npoints, points);
    
    U(:, :, :, n) = reshape(result3(1, :), nx, ny, nz);
    V(:, :, :, n) = reshape(result3(2, :), nx, ny, nz);
    W(:, :, :, n) = reshape(result3(3, :), nx, ny, nz);
    
end

end