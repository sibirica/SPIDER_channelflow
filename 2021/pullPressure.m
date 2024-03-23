function [p] = pullPressure(x, y, z, t)

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

p = zeros(nx, ny, nz, nt);

for n = 1:nt

    time = t(n);
    result = getPressure(authkey, dataset, time, Lag4, PCHIPInt, npoints, points);
    
    p(:, :, :, n) = reshape(result, nx, ny, nz);
    
end

end