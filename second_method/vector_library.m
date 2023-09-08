
function [G, labels, scales] = scalar_library_function( noise, envelope_power )

%JHU data
tic
[p, p0, U, V, W, x, y, z, t] = read_h5_channelflow_data();
toc

%{
[p, U, V, W, x,y,z,t] = load_in_Channelflow_data();
dt = t(2) - t(1);
%}

%{
[p, U, V, W, x,y,z,t] = load_interpolated_JHU();
dt = t(2) - t(1);
%}

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

%alternative names
p_mean = P_mean;
p_std  = P_std;

dx = x(2) - x(1);
dt = t(2) - t(1);


num_lib = 5; %25; %purposely lower than the "true" 28
num_windows = 100; %25*num_lib;
nw = num_windows;
dof = 3; %3 degrees of freedom for vector equation, 1 for scalar

%Use a high order polynomial envelope
dimension = 4;
%envelope_power = 5; %found 5 to be best for JHU data
pol = envelope_pol( envelope_power, dimension );

G = zeros( dof*num_windows, num_lib );
labels = cell(num_lib, 1);

corners  = zeros(4, num_windows);
%size_vec = [27 32 33 38]; %directly from SPIDER paper
size_vec = [35, 35, 35, 35];
%size_vec = [25, 25, 25, 25];
size_vec = round(size_vec); %make sure it's integers

U = U + noise * std(U,0,'all');
V = V + noise * std(V,0,'all');
W = W + noise * std(W,0,'all');
p = p + noise * std(p,0,'all');


seed = 1;
rng(seed);
L = 0;
Lx = size(U,1); Ly = size(U,2); Lz = size(U,3); Lt = size(U,4);
for i=1:num_windows
  for j=1:4
    corners(j,i) = randi( size(U,j) - size_vec(j) - 2*L ) + L;
  end
end

U_x = diff_dim( U, x, 1);
V_y = diff_dim( V, y, 2);
W_z = diff_dim( W, z, 3);

U_t = diff_dim( U, t, 4);
V_t = diff_dim( V, t, 4);
W_t = diff_dim( W, t, 4);
p_t = diff_dim( p, t, 4);

a=1; %library index we will increment as we add terms

labels{a} = "u_i";
G(:,a)    = [SPIDER_integrate( U, [], grid, corners, size_vec, pol );
             SPIDER_integrate( V, [], grid, corners, size_vec, pol );
             SPIDER_integrate( W, [], grid, corners, size_vec, pol );];
scales(a) = U_mean;
ideal(a)  = 0;
a = a+1;


labels{a} = "\partial_t u_i";
G(:,a)    = [SPIDER_integrate( U, [4], grid, corners, size_vec, pol );
             SPIDER_integrate( V, [4], grid, corners, size_vec, pol );
             SPIDER_integrate( W, [4], grid, corners, size_vec, pol );];
scales(a) = U_std/dt;
ideal(a)  = 1;
a = a+1;


labels{a} = "\nabla_i p";
G(:,a)    = [SPIDER_integrate( p, [1], grid, corners, size_vec, pol );
             SPIDER_integrate( p, [2], grid, corners, size_vec, pol );
             SPIDER_integrate( p, [3], grid, corners, size_vec, pol );];
scales(a) = P_std/dx;
ideal(a)  = 1;
a = a+1;


labels{a} = "p u_i";
G(:,a)    = [SPIDER_integrate( p.*U, [], grid, corners, size_vec, pol );
             SPIDER_integrate( p.*V, [], grid, corners, size_vec, pol );
             SPIDER_integrate( p.*W, [], grid, corners, size_vec, pol );];
scales(a) = p_mean*U_mean;
ideal(a)  = 0;
a = a+1;

%Replacing (u dot \nabla)u
labels{a} = "\nabla_j( u_i u_j)";
G(:,a)    = [SPIDER_integrate( U.*U, [1], grid, corners, size_vec, pol ) + SPIDER_integrate( U.*V, [2], grid, corners, size_vec, pol ) + SPIDER_integrate( U.*W, [3], grid, corners, size_vec, pol );
             SPIDER_integrate( V.*U, [1], grid, corners, size_vec, pol ) + SPIDER_integrate( V.*V, [2], grid, corners, size_vec, pol ) + SPIDER_integrate( V.*W, [3], grid, corners, size_vec, pol );
             SPIDER_integrate( W.*U, [1], grid, corners, size_vec, pol ) + SPIDER_integrate( W.*V, [2], grid, corners, size_vec, pol ) + SPIDER_integrate( W.*W, [3], grid, corners, size_vec, pol );];
scales(a) = U_std*U_mean/dx;
ideal(a)  = 1;
a = a+1;


labels{a} = "\nabla^2 u_i";
G(:,a)    = [SPIDER_integrate( U, [1,1], grid, corners, size_vec, pol ) + SPIDER_integrate( U, [2,2], grid, corners, size_vec, pol ) + SPIDER_integrate( U, [3,3], grid, corners, size_vec, pol );
             SPIDER_integrate( V, [1,1], grid, corners, size_vec, pol ) + SPIDER_integrate( V, [2,2], grid, corners, size_vec, pol ) + SPIDER_integrate( V, [3,3], grid, corners, size_vec, pol );
             SPIDER_integrate( W, [1,1], grid, corners, size_vec, pol ) + SPIDER_integrate( W, [2,2], grid, corners, size_vec, pol ) + SPIDER_integrate( W, [3,3], grid, corners, size_vec, pol );];
scales(a) = U_std/dx/dx;
ideal(a)  = 1;
a = a+1;


labels{a} = "\partial_t^2 u_i";
G(:,a)    = [SPIDER_integrate( U, [4,4], grid, corners, size_vec, pol );
             SPIDER_integrate( V, [4,4], grid, corners, size_vec, pol );
             SPIDER_integrate( W, [4,4], grid, corners, size_vec, pol );];
scales(a) = U_std/dx/dx;
ideal(a)  = 1;
a = a+1;


labels{a} = "u^2 u_i";
U_sq = U.^2 + V.^2 + W.^2;
G(:,a)    = [SPIDER_integrate( U_sq.*U, [], grid, corners, size_vec, pol );
             SPIDER_integrate( U_sq.*V, [], grid, corners, size_vec, pol );
             SPIDER_integrate( U_sq.*W, [], grid, corners, size_vec, pol );];
scales(a) = U_mean^3;
ideal(a)  = 1;
a = a+1;


labels{a} = "p^2 u_i";
G(:,a)    = [SPIDER_integrate( p.*p.*U, [], grid, corners, size_vec, pol );
             SPIDER_integrate( p.*p.*V, [], grid, corners, size_vec, pol );
             SPIDER_integrate( p.*p.*W, [], grid, corners, size_vec, pol );];
scales(a) = U_mean*p_mean^2;
ideal(a)  = 1;
a = a+1;


labels{a} = "\nabla_i p^2";
G(:,a)    = [SPIDER_integrate( p.*p, [1], grid, corners, size_vec, pol );
             SPIDER_integrate( p.*p, [2], grid, corners, size_vec, pol );
             SPIDER_integrate( p.*p, [3], grid, corners, size_vec, pol );];
scales(a) = p_mean.*p_std/dx;
ideal(a)  = 1;
a = a+1;

labels{a} = "(\nabla_j u_j) u_i";
div = U_x + V_y + W_z;
G(:,a)    = [SPIDER_integrate( div.*U, [], grid, corners, size_vec, pol );
             SPIDER_integrate( div.*V, [], grid, corners, size_vec, pol );
             SPIDER_integrate( div.*W, [], grid, corners, size_vec, pol );];
scales(a) = U_mean*U_std/dx;
ideal(a)  = 1;
a = a+1;

labels{a} = "\nabla_i u^2";
G(:,a)    = [SPIDER_integrate( U_sq, [1], grid, corners, size_vec, pol );
             SPIDER_integrate( U_sq, [2], grid, corners, size_vec, pol );
             SPIDER_integrate( U_sq, [3], grid, corners, size_vec, pol );];
scales(a) = U_mean.*U_std/dx;
ideal(a)  = 1;
a = a+1;


labels{a} = "\nabla_i \nabla_j u_j";
G(:,a)    = [SPIDER_integrate( U, [1,1], grid, corners, size_vec, pol ) + SPIDER_integrate( V, [2,1], grid, corners, size_vec, pol ) + SPIDER_integrate( W, [3,1], grid, corners, size_vec, pol );
             SPIDER_integrate( U, [1,2], grid, corners, size_vec, pol ) + SPIDER_integrate( V, [2,2], grid, corners, size_vec, pol ) + SPIDER_integrate( W, [3,2], grid, corners, size_vec, pol );
             SPIDER_integrate( U, [1,3], grid, corners, size_vec, pol ) + SPIDER_integrate( V, [2,3], grid, corners, size_vec, pol ) + SPIDER_integrate( W, [3,3], grid, corners, size_vec, pol );];
scales(a) = U_std/dx/dx;
ideal(a)  = 1;
a = a+1;

labels{a} = "p \partial_t u_i";
G(:,a)    = [SPIDER_integrate( p.*U_t, [], grid, corners, size_vec, pol );
             SPIDER_integrate( p.*V_t, [], grid, corners, size_vec, pol );
             SPIDER_integrate( p.*W_t, [], grid, corners, size_vec, pol );];
scales(a) = p_mean.*U_std/dt;
ideal(a)  = 1;
a = a+1;

labels{a} = "u_i \partial_t p";
G(:,a)    = [SPIDER_integrate( p_t.*U, [], grid, corners, size_vec, pol );
             SPIDER_integrate( p_t.*V, [], grid, corners, size_vec, pol );
             SPIDER_integrate( p_t.*W, [], grid, corners, size_vec, pol );];
scales(a) = p_std.*U_mean/dt;
ideal(a)  = 1;
a = a+1;




labels{a} = "e_x";
G(:,a)    = [SPIDER_integrate( 0*p+1, [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*p,   [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*p,   [], grid, corners, size_vec, pol );];
scales(a) = 1;
ideal(a)  = 0;
a = a+1;

labels{a} = "e_y";
G(:,a)    = [SPIDER_integrate( 0*p,   [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*p+1, [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*p,   [], grid, corners, size_vec, pol );];
scales(a) = 1;
ideal(a)  = 0;
a = a+1;

labels{a} = "e_z";
G(:,a)    = [SPIDER_integrate( 0*p,   [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*p,   [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*p+1, [], grid, corners, size_vec, pol );];
scales(a) = 1;
ideal(a)  = 0;
a = a+1;


%Lastly, normalize by the integral of unity to standardize things
norm_vec = SPIDER_integrate( ones(size(U)), [], grid, corners, size_vec, pol );      
norm_vec = repmat( norm_vec, [3,1] );

G = G./norm_vec;
G = G./scales;
return

end