function [P, P0, U, V, W, x, y, z, t] = read_h5_channelflow_data()
  %{
  PURPOSE:
  I used the web interface at http://turbulence.idies.jhu.edu/cutout/jobs
  to pull pressure and velocity data. The resultant h5 files are stored in
  the h5_cutouts folder.

  This function reads the data from those h5 files and puts them into a
  format usable by our scripts.

  INPUT:
  None at the moment

  OUTPUT:
  P - pressure
  U,V,W - components of velocity
  x,y,z - grids in each dimension 
  %}

  pressure_file = "h5_cutouts/pressure_middle.h5";
  velocity_file = "h5_cutouts/velocity_middle.h5";
  
  % I pulled data from points 128 to 192 in every dimension
  %
  Lx = 65;
  Ly = 65;
  Lz = 65;
  Lt = 64;
  
  P   = zeros(   Lx, Ly, Lz, Lt);
  vel = zeros(3, Lx, Ly, Lz, Lt);
  
  for t = 1:Lt
    %tt = 128 + t;
    tt = 2000 + t;
    P(:,:,:,t) = h5read(pressure_file, sprintf("/Pressure_%04d", tt));
    
    vel(:,:,:,:,t) = h5read(velocity_file, sprintf("/Velocity_%04d", tt));
  end
  
  %Read grid as well
  x = h5read( velocity_file, '/xcoor');
  y = h5read( velocity_file, '/ycoor');
  z = h5read( velocity_file, '/zcoor');
  
  U = squeeze( vel(1,:,:,:,:) );
  V = squeeze( vel(2,:,:,:,:) );
  W = squeeze( vel(3,:,:,:,:) );
  
  U = U - 0.45; %since simulation was done in comoving frame

  P = double(P);
  U = double(U);
  V = double(V); 
  W = double(W);
  
  dt = 0.0065; %JHU timestep pulled from the README
  t = dt * (1:size(U,4)); %make a time coordinate array
  
  
  [X,~,~,~] =ndgrid(x,y,z,t);
  P0 = P; %save before subtracting X 
  P = P - 0.0025*X;
  
  %{
  %interpolate data onto uniform y 
  y_uniform = linspace( min(y), max(y), Ly );
  tic
  for i = 1:Lx
      for k = 1:Lz
          for t = 1:Lt
              method = 'spline';
              P(i,:,k,t) = interp1( y, P(i,:,k,t), y_uniform, method );
              U(i,:,k,t) = interp1( y, U(i,:,k,t), y_uniform, method );
              V(i,:,k,t) = interp1( y, V(i,:,k,t), y_uniform, method );
              W(i,:,k,t) = interp1( y, W(i,:,k,t), y_uniform, method );
          end
      end
  end
  toc
  %}
end