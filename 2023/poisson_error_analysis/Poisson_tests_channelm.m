%
% Turbmat - a Matlab library for the JHU Turbulence Database Cluster
%   
% Sample code, part of Turbmat
%

%
% Written by:
%  
% Perry Johnson
% The Johns Hopkins University
% Department of Mechanical Engineering
% pjohns86@jhu.edu, johnson.perry.l@gmail.com
%

% 2017: Modified by Zhao Wu and Rohit Ravoori, adding getInvariant 
% 2019: Modified by Zhao Wu for getThreshold function

%
% This file is part of Turbmat.
% 
% Turbmat is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
% 
% Turbmat is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with Turbmat.  If not, see <http://www.gnu.org/licenses/>.
%

clear all;
close all;

authkey = 'edu.jhu.pha.turbulence.testing-201406';
dataset = 'channel';

% ---- Temporal Interpolation Options ----
NoTInt   = 'None' ; % No temporal interpolation
PCHIPInt = 'PCHIP'; % Piecewise cubic Hermit interpolation in time

% ---- Spatial Interpolation Flags for getVelocity & getVelocityAndPressure ----
NoSInt = 'None'; % No spatial interpolation
Lag4   = 'Lag4'; % 4th order Lagrangian interpolation in space
Lag6   = 'Lag6'; % 6th order Lagrangian interpolation in space
Lag8   = 'Lag8'; % 8th order Lagrangian interpolation in space

% ---- Spatial Differentiation & Interpolation Flags for getVelocityGradient & getPressureGradient ----
FD4NoInt = 'None_Fd4' ; % 4th order finite differential scheme for grid values, no spatial interpolation
FD6NoInt = 'None_Fd6' ; % 6th order finite differential scheme for grid values, no spatial interpolation
FD8NoInt = 'None_Fd8' ; % 8th order finite differential scheme for grid values, no spatial interpolation
FD4Lag4  = 'Fd4Lag4'  ; % 4th order finite differential scheme for grid values, 4th order Lagrangian interpolation in space

% ---- Spline interpolation and differentiation Flags for getVelocity,
% getPressure, getVelocityGradient, getPressureGradient,
% getVelocityHessian, getPressureHessian
M1Q4   = 'M1Q4'; % Splines with smoothness 1 (3rd order) over 4 data points. Not applicable for Hessian.
M2Q8   = 'M2Q8'; % Splines with smoothness 2 (5th order) over 8 data points.
M2Q14  = 'M2Q14'; % Splines with smoothness 2 (5th order) over 14 data points.

 % fprintf(1,'%3i: duxdx=%13.6e, duxdy=%13.6e, duxdz=%13.6e, ', p, result9(1,p), result9(2,p), result9(3,p));
 % fprintf(1,'duydx=%13.6e, duydy=%13.6e, duydz=%13.6e, ', result9(4,p), result9(5,p), result9(6,p));
 % fprintf(1,'duzdx=%13.6e, duzdy=%13.6e, duzdz=%13.6e\n', result9(7,p), result9(8,p), result9(9,p)); 
 % fprintf(1,'%3i: d2pdxdx=%13.6e, d2pdxdy=%13.6e, d2pdxdz=%13.6e, ', p, result6(1,p), result6(2,p), result6(3,p));
 % fprintf(1,'d2pdydy=%13.6e, d2pdydz=%13.6e, d2pdzdz=%13.6e\n', result6(4,p), result6(5,p), result6(6,p));

 s1=0;
 s2=0;
 s3=0;
 s4=0;

 % generate 20 planes comparisons at random times and (x,y,z) locations
for k = 1:20
time = 25.*rand;  % random time in interval
spacing = 8 * 0.001;  % spacing of 8 x viscous scales 
xoff =  8*pi*rand; % random x position of plane
yoff =  2.*rand-1.0; % random y-position of plane
zoff =  3*pi*rand; % random x position of plane
nx = 25;
nz = nx;
npoints = nx*nz;
%
points = zeros(3,npoints);
result1  = zeros(1,npoints);
result3  = zeros(3,npoints);
result6  = zeros(6,npoints);
result9  = zeros(9,npoints);
%
A = zeros(3,3,npoints);
Q = zeros(1,npoints);
Laplp = zeros(1,npoints);
%
% Create grid:
x = linspace(0, (nx-1)*spacing, nx) + xoff;
z = linspace(0, (nz-1)*spacing, nz) + zoff;
[X Z] = meshgrid(x, z);
points(1,:) = X(:)';
points(3,:) = Z(:)';
points(2,:) = yoff;
%    
fprintf('\nRequesting velocity, gradients, pressure and its Hessian at (%ix%i) points \n',nx,nz);
%
result9 = getVelocityGradient(authkey, dataset, time, M2Q14, PCHIPInt , npoints, points);
result3 = getVelocity(authkey, dataset, time, M2Q14, PCHIPInt , npoints, points);
result6 = getPressureHessian(authkey, dataset, time, M2Q14, PCHIPInt , npoints, points);
result1 = getPressure(authkey, dataset, time, M2Q14, PCHIPInt , npoints, points);
%
A(1,1,:) = result9(1,:);
A(1,2,:) = result9(2,:);
A(1,3,:) = result9(3,:);
A(2,1,:) = result9(4,:);
A(2,2,:) = result9(5,:);
A(2,3,:) = result9(6,:);
A(3,1,:) = result9(7,:);
A(3,2,:) = result9(8,:);
A(3,3,:) = result9(9,:);
%
% Compute -Trace(A^2) as LHS of Poisson equation (note, this is 2Q): 
for p = 1:npoints
 Q(1,p) = 0;
  for i=1:3
    for j=1:3
        Q(1,p) = Q(1,p) - A(i,j,p)*A(j,i,p);
    end
  end
end
%
% 
% Evaluate Laplacian of pressure
Laplp(1,:) = result6(1,:)+result6(4,:)+result6(6,:); 
%
% accumulate sum of square errors:
s1=s1+sum((Laplp(1,:)-Q(1,:)).^2);
s2=s2+sum(Laplp(1,:).^2);
s3=s3+sum((A(1,1,:)+A(2,2,:)+A(3,3,:)).^2);
s4=s4+sum(A(1,1,:).^2);
%
% reshape for plotting:
Y1 = reshape(Q, nz, nx);
Y2 = reshape(Laplp, nz, nx);
Y3 = reshape(result1, nz, nx);
Y4 = reshape(squeeze(A(1,1,:)+A(2,2,:)), nz, nx);
Y5 = - reshape(squeeze(A(3,3,:)), nz, nx);
Y6 = reshape(squeeze(result3(1,:)), nz, nx);
%
lim1 = max(max(Q),max(Laplp));
lim2 = max(max(squeeze(A(3,3,:))));
%
t1 = tiledlayout(2,2);
t1.TileSpacing = 'compact';
t1.Padding = 'compact';
txt = ['y=',num2str(yoff)];
sgtitle(txt);
% Plot -2Q contours
nexttile
contourf(X, Z, Y1, 50, 'LineStyle', 'none','LineColor','none');
set(gca, 'FontSize', 11)
title('-2Q','FontSize', 11);
xlabel('x', 'FontSize', 11);
ylabel('z', 'FontSize', 11);
caxis([-lim1 lim1]);
colorbar('FontSize', 11);
axis([min(x) max(x) min(z) max(z)]);
pbaspect([1 1 1])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');
%
% Plot Laplacian of pressure:
nexttile
contourf(X, Z, Y2, 50, 'LineStyle', 'none','LineColor','none');
set(gca, 'FontSize', 11)
title('Lapl(p)','FontSize', 11);
xlabel('x', 'FontSize', 11);
ylabel('z', 'FontSize', 11);
caxis([-lim1 lim1]);
colorbar('FontSize', 11);
axis([min(x) max(x) min(z) max(z)]);
pbaspect([1 1 1])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');
%
% Plot error:
nexttile
contourf(X, Z, (Y1-Y2), 50, 'LineStyle', 'none');
set(gca, 'FontSize', 11)
title('(Lapl(p)+2Q)','FontSize', 11);
xlabel('x', 'FontSize', 11);
ylabel('z', 'FontSize', 11);
caxis([-lim1 lim1]);
colorbar('FontSize', 11);
axis([min(x) max(x) min(z) max(z)]);
pbaspect([1 1 1])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on'); 
%
% Plot pressure:
nexttile
contourf(X, Z, Y3, 50, 'LineStyle', 'none');
set(gca, 'FontSize', 11)
title('pressure','FontSize', 11);
xlabel('x', 'FontSize', 11);
ylabel('z', 'FontSize', 11);
colorbar('FontSize', 11);
axis([min(x) max(x) min(z) max(z)]);
pbaspect([1 1 1])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');
%
exportgraphics(t1,'plot.pdf','Append',true)
%
t2 = tiledlayout(2,2);
t2.TileSpacing = 'compact';
t2.Padding = 'compact';
txt = ['y=',num2str(yoff)];
sgtitle(txt);
% Plot dudx+dvdy:
nexttile
contourf(X, Z, Y4, 50, 'LineStyle', 'none');
set(gca, 'FontSize', 11)
title('dudx+dvdy','FontSize', 11);
xlabel('x', 'FontSize', 11);
ylabel('z', 'FontSize', 11);
caxis([-lim2 lim2]);
colorbar('FontSize', 11);
axis([min(x) max(x) min(z) max(z)]);
pbaspect([1 1 1])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on'); 

% Plot -dwdz:
nexttile
contourf(X, Z, Y5, 50, 'LineStyle', 'none');
set(gca, 'FontSize', 11)
title('-dwdz','FontSize', 11);
xlabel('x', 'FontSize', 11);
ylabel('z', 'FontSize', 11);
caxis([-lim2 lim2]);
colorbar('FontSize', 11);
axis([min(x) max(x) min(z) max(z)]);
pbaspect([1 1 1])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on'); 
 
% Plot dudx+dvdy+dwdz:
nexttile
contourf(X, Z, (Y4-Y5), 50, 'LineStyle', 'none');
set(gca, 'FontSize', 11)
title('dudx+dvdy+dwdz','FontSize', 11);
xlabel('x', 'FontSize', 11);
ylabel('z', 'FontSize', 11);
caxis([-lim2 lim2]);
colorbar('FontSize', 11);
axis([min(x) max(x) min(z) max(z)]);
pbaspect([1 1 1])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on'); 

% Plot u:
nexttile
contourf(X, Z, Y6, 50, 'LineStyle', 'none');
set(gca, 'FontSize', 11)
title('u','FontSize', 11);
xlabel('x', 'FontSize', 11);
ylabel('z', 'FontSize', 11);
colorbar('FontSize', 11);
axis([min(x) max(x) min(z) max(z)]);
pbaspect([1 1 1])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on'); 

exportgraphics(t2,'plot.pdf','Append',true)

end
% ouput averaged rms of errors compared to rms of terms if uncorrelated: 
disp(' ');
fprintf('relative error Poisson Eq: = %f',(s1/(2*s2))^0.5 );
disp(' ');
disp(' ');
fprintf('relative error continuity Eq: = %f', (s3/(3*s4))^0.5);
disp(' ');
 
