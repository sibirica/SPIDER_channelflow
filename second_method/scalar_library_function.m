function [G, labels, scales] = scalar_library_function( noise, envelope_power, size_vec, location, scales_override, seed )


addpath('functions/');
addpath('../')

%JHU data load-in

location
if location == "edge"
  [p, p0, U, V, W, x, y, z, t] = read_h5_channelflow_data_edge();
end

if location == "center"
  [p, p0, U, V, W, x, y, z, t] = read_h5_channelflow_data_center();
end


%seed random number generation
rng(0);

grid = cell(4,1);
grid{1} = x;
grid{2} = y;
grid{3} = z;
grid{4} = t;

p = p - mean(p,'all'); %subtract mean pressure

U_mean = mean( sqrt(U.^2 + V.^2 + W.^2),    'all' );
std_vec= [std(  U, 0, 'all' );
          std(  V, 0, 'all' );
          std(  W, 0, 'all' );];
U_std  = max(std_vec);

P_mean = mean( abs(p0), 'all' );
P_std  = std(  p0, 0,   'all' );

p_mean = P_mean;
p_std  = P_std;

%Compute derivatives
U_x = diff_dim( U, x, 1);
U_y = diff_dim( U, y, 2);
U_z = diff_dim( U, z, 3);
V_x = diff_dim( V, x, 1);
V_y = diff_dim( V, y, 2);
V_z = diff_dim( V, z, 3);
W_x = diff_dim( W, x, 1);
W_y = diff_dim( W, y, 2);
W_z = diff_dim( W, z, 3);
p_x = diff_dim( p, x, 1);
p_y = diff_dim( p, y, 2);
p_z = diff_dim( p, z, 3);

U_t = diff_dim( U, t, 4);
V_t = diff_dim( V, t, 4);
W_t = diff_dim( W, t, 4);
p_t = diff_dim( p, t, 4);

%compute length and time scales for fluid flow and pressure separately
Lu = U_std/sqrt( mean(U_x.^2 + U_y.^2 + U_z.^2 + V_x.^2 + V_y.^2 + V_z.^2 + W_x.^2 + W_y.^2 + W_z.^2, 'all') );
Lp = P_std/sqrt( mean(p_x.^2 + p_y.^2 + p_z.^2, 'all') );
Tu = U_std/sqrt( mean(U_t.^2 + V_t.^2 + W_t.^2, 'all') );
Tp = P_std/sqrt( mean(p_t.^2, 'all') );

num_lib = 6; %number of library terms
num_windows = 256;
nw = num_windows;
dof = 1; %3 degrees of freedom for vector equation, 1 for scalar

%Use a high order polynomial envelope
dimension = 4;
pol = envelope_pol( envelope_power, dimension );

G = zeros( dof*num_windows, num_lib );
labels = cell(num_lib, 1);

corners  = zeros(4, num_windows);
size_vec = round(size_vec); %make sure it's integers

uniform_noise = @() 2*rand( size(U) ) - 1;

U = U + noise * std(U,0,'all') * uniform_noise();
V = V + noise * std(V,0,'all') * uniform_noise();
W = W + noise * std(W,0,'all') * uniform_noise();
p = p + noise * std(p,0,'all') * uniform_noise();

rng(seed);
L = 1;
Lx = size(U,1); Ly = size(U,2); Lz = size(U,3); Lt = size(U,4);
for i=1:num_windows
  for j=1:4
    corners(j,i) = randi( size(U,j) - size_vec(j) - 2*L ) + L;
  end
end


%recompute derivatives since noise was added
U_x = diff_dim( U, x, 1);
U_y = diff_dim( U, y, 2);
U_z = diff_dim( U, z, 3);
V_x = diff_dim( V, x, 1);
V_y = diff_dim( V, y, 2);
V_z = diff_dim( V, z, 3);
W_x = diff_dim( W, x, 1);
W_y = diff_dim( W, y, 2);
W_z = diff_dim( W, z, 3);
p_x = diff_dim( p, x, 1);
p_y = diff_dim( p, y, 2);
p_z = diff_dim( p, z, 3);


U_t = diff_dim( U, t, 4);
V_t = diff_dim( V, t, 4);
W_t = diff_dim( W, t, 4);
p_t = diff_dim( p, t, 4);



a=1; %library index we will increment as we add terms



labels{a} = "p";
G(:,a)    = SPIDER_integrate( p, [], grid, corners, size_vec, pol );   
scales(a) = P_mean;
ideal(a)  = 0;
a = a+1;


labels{a} = "\nabla_i u_i";
G(:,a)    = SPIDER_integrate( U, [1], grid, corners, size_vec, pol ) + ...
            SPIDER_integrate( V, [2], grid, corners, size_vec, pol ) + ...
            SPIDER_integrate( W, [3], grid, corners, size_vec, pol );
scales(a) = U_std/Lu; %U_std/dx;
ideal(a)  = 0;
a = a+1;

labels{a} = "\partial_t p";
G(:,a)    = SPIDER_integrate( p, [4], grid, corners, size_vec, pol );     
scales(a) = P_std/Tp; %P_std/dt;
ideal(a)  = 0;
a = a+1;


labels{a} = "p^2";
G(:,a)    = SPIDER_integrate( p.^2, [], grid, corners, size_vec, pol );   
scales(a) = P_mean^2;
ideal(a)  = 0;
a = a+1;

labels{a} = "E";
E = (U.^2 + V.^2 + W.^2)/2;
G(:,a)    = SPIDER_integrate( E, [], grid, corners, size_vec, pol );   
scales(a) = U_mean^2;
ideal(a)  = 0;
a = a+1;


labels{a} = "u_i \nabla_i p";
G(:,a)    = SPIDER_integrate( U.*p_x + V.*p_y + W.*p_z, [], grid, corners, size_vec, pol );   
scales(a) = U_mean*P_std/Lp;
ideal(a)  = 0;
a = a+1;

labels{a} = "p^3";
G(:,a)    = SPIDER_integrate( p.^3, [], grid, corners, size_vec, pol );   
scales(a) = P_mean^2;
ideal(a)  = 0;
a = a+1;

labels{a} = "\partial_t p^2";
G(:,a)    = SPIDER_integrate( p.^2, [4], grid, corners, size_vec, pol );     
scales(a) = P_mean*P_std/Tp;
ideal(a)  = 0;
a = a+1;


labels{a} = "\partial_t^2 p";
G(:,a)    = SPIDER_integrate( p, [4,4], grid, corners, size_vec, pol );     
scales(a) = P_std/Tp;
ideal(a)  = 0;
a = a+1;

labels{a} = "p \nabla_i u_i";
G(:,a)    = SPIDER_integrate( (U_x + V_y + W_z).*p, [], grid, corners, size_vec, pol );   
scales(a) = U_std*P_mean/Lu;
ideal(a)  = 0;
a = a+1;

labels{a} = "E p";
G(:,a)    = SPIDER_integrate( E.*p, [], grid, corners, size_vec, pol );   
scales(a) = U_mean.*U_mean*P_mean;
ideal(a)  = 0;
a = a+1;


labels{a} = "\partial_t E";
G(:,a)    = SPIDER_integrate( E, [4], grid, corners, size_vec, pol );
scales(a) = U_mean.*U_std/Tu;
ideal(a)  = 0;
a = a+1;


%Source term of PP
labels{a} = "\nabla_i \nabla_j (u_i u_j)";
G(:,a) =   SPIDER_integrate( U.*U, [1,1], grid, corners, size_vec, pol ) ...
       +   SPIDER_integrate( V.*V, [2,2], grid, corners, size_vec, pol ) ...
       +   SPIDER_integrate( W.*W, [3,3], grid, corners, size_vec, pol ) ...
       + 2*SPIDER_integrate( U.*V, [1,2], grid, corners, size_vec, pol ) ...
       + 2*SPIDER_integrate( U.*W, [1,3], grid, corners, size_vec, pol ) ...
       + 2*SPIDER_integrate( V.*W, [2,3], grid, corners, size_vec, pol );
scales(a) = U_std^2/Lu^2;
ideal(a) = 1;
a = a+1;

labels{a} = "\nabla_i(u_i E)";
G(:,a)    =  SPIDER_integrate( U.*E, [1], grid, corners, size_vec, pol ) ...
          +  SPIDER_integrate( V.*E, [2], grid, corners, size_vec, pol ) ...
          +  SPIDER_integrate( W.*E, [3], grid, corners, size_vec, pol );
scales(a) = U_mean.^2*U_std/Lu;
ideal(a)  = 0;
a = a+1;

labels{a} = "\nabla^2 E";
G(:,a) =   SPIDER_integrate( E, [1,1], grid, corners, size_vec, pol ) ...
       +   SPIDER_integrate( E, [2,2], grid, corners, size_vec, pol ) ...
       +   SPIDER_integrate( E, [3,3], grid, corners, size_vec, pol );
scales(a) = U_std*U_mean/Lu^2;
ideal(a) = 1;
a = a+1;

labels{a} = "(\nabla_i u_j)(\nabla_i u_j)";
G(:,a) =   SPIDER_integrate( U_x.^2 + U_y.^2 + U_z.^2 + V_x.^2 + V_y.^2 + V_z.^2 + W_x.^2 + W_y.^2 + W_z.^2, [], grid, corners, size_vec, pol );
scales(a) = U_std^2/Lu^2;
ideal(a) = 1;
a = a+1;

labels{a} = "\nabla^2 p";
G(:,a)    = SPIDER_integrate( p0, [1,1], grid, corners, size_vec, pol ) + ...
            SPIDER_integrate( p0, [2,2], grid, corners, size_vec, pol ) + ...
            SPIDER_integrate( p0, [3,3], grid, corners, size_vec, pol );      
scales(a) = P_std/Lp^2;
ideal(a)  = 1;
a = a+1;

%{
labels{a} = "1";
G(:,a)    = SPIDER_integrate( 0*p+1, [], grid, corners, size_vec, pol );
scales(a) = 1;
ideal(a)  = 0;
a = a+1;
%}

labels{a} = "u^4";
G(:,a)    = SPIDER_integrate( 4 * E.^2, [], grid, corners, size_vec, pol );
scales(a) = U_mean^4;
ideal(a)  = 0;
a = a+1;

%Lastly, normalize by the integral of unity to standardize things
norm_vec = SPIDER_integrate( ones(size(U)), [], grid, corners, size_vec, pol );      

if numel( scales_override ) ~= 0
   scales = scales_override; 
end

G = G./norm_vec; % domain-wise normalization
G = G./scales;   %   term-wise normalization
end






function [ L_scale, T_scale ] = generate_length_and_time_scales(U, V, W, grid)
  dimension = 4; %4D data: x,y,z,t
  
  corners = ones(dimension,1);
  size_vec= size(U);
  pol     = envelope_pol(4, dimension);
  
  int_dvaldt = SPIDER_integrate( W, [4], grid, corners, size_vec, pol );
  int_val    = SPIDER_integrate( W,  [], grid, corners, size_vec, pol );

  T_scale = abs(int_val/int_dvaldt);
  L_scale = mean( sqrt(U.^2 + V.^2 + W.^2), 'all' )*T_scale;
end