%%%% Weak-Form Parameter Estimation %%%%
%% ----------------------------------------------------------------------
%                             PREAMBLE
% -----------------------------------------------------------------------
%{
    INPUTS 
    ------
    term_names
        * string vector containing keys for terms to include in library
        * CHECK glossaryNS.m
    N_d             
        * number of integration domains
    N_h = [Q R S] 
        * max number of weight function harmonics. S >= 1, odd; Q, R >= 0
    opts (struct)
        * optional inputs: includes...
        filename        
            * string for path to -v7.3 .mat file containing saved data        
        mode
            * which equation to identify: NS, div, or BC
        domain_loc
            * location of domains to use: center or edge
        D = [Hx, Hy, Ht]
            * size of integration domains with units {nd, nd, s}
        env_exp = [alpha, beta, gamma]
            * envelope function exponents. alpha, beta >= 3; gamma >= 1. 
        noise           
            * amount of Gaussian noise to apply to data after loading
        seed           
            * random seed for selecting domain locations and noise
        track           
            * boolean to track library filling in cmd window

    OUTPUTS
    -------
    Q   
        * matrix containing integral evaluations of gov eq terms    
    P   
        * 4 x N_d matrix containing integration domain locations
    H
        * size of domains
    char_sizes
        * characteristic sizes of each term
    valid_single
        * vector representing which terms can appear alone in the model
    lib_names
        * names of library terms
%}
%%
function [Q, P, H, char_sizes, valid_single, lib_names] = ...
         SPIDR_NS(term_names, N_d, N_h, opts)
%% ----------------------------------------------------------------------
%                          DEFINE LIBRARY
% -----------------------------------------------------------------------

[lib_names, krnl_info, tex_names] = defineLibrary(term_names, opts.mode);

N = length(lib_names);

%% ----------------------------------------------------------------------
%                       INITIALIZE VARIABLES
% -----------------------------------------------------------------------
% Read options
if isfield(opts, 'mode')
     mode = opts.mode;
else
     mode = "NS";
end
if isfield(opts, 'file')
     file = opts.file;
else
     file = "check_this.mat";
end
if isfield(opts, 'env_exp')
     env_exp = opts.env_exp;
else
     env_exp = [6 6 6 6];
end
if isfield(opts, 'domain_loc')
     domain_loc = opts.domain_loc;
else
     domain_loc = "edge";
end
if isfield(opts, 'track')
     track = opts.track;
else
     track = 1;
end
if isfield(opts, 'noise')
     noise = opts.noise;
else
     noise = 0;
end
if isfield(opts, 'seed') % random seed for choosing domains/noise
     seed = opts.seed;
else
     seed = 1;
end

% Full channel grid parameters
Lx = 8*pi;
Ly = 2;
Lz = 3*pi;
Lt = 26;

Nx = 2048;
Ny = 512;
Nz = 1536;
Nt = 4000;

dx = 8*pi/Nx;
dy = 2/Ny;
dz = 3*pi/Nz;
dt = 0.0065;

% Full single-dimension grids
x_full = linspace(0, 8*pi - dx, Nx);
y_full = load('y.txt');
if mode == "BC"
    y_full = y_full(1:3); % only keep the locations of the bottom 3 points
end
z_full = linspace(0, 3*pi - dz, Nz);
t_full = linspace(0, 26 - dt, Nt);

% Rough size of integration domains
Hx = 1/3;

if domain_loc == "edge"
    Hy = 0.02; % for NS
else
    Hy = 1/5; % for Euler
end

Hz = 1/5;
Ht = 1/4;

H = [Hx, Hy, Hz, Ht];

% Initialize random seed for reproducibility
rng(seed);

% Sampling scheme
P = zeros(4, N_d);
P(1, :) = (8*pi - 1.1*Hx)*rand(1, N_d);
if domain_loc == "edge"
    P(2, :) = -1;  % attach domain to the wall
elseif domain_loc == "center"
    P(2, :) = -0.1; % center domain about middle
else
    P(2, :) = (2 - 1.1*Hy)*rand(1, N_d) - 1;
end
P(3, :) = (3*pi - 1.1*Hz)*rand(1, N_d);
P(4, :) = (26 - 1.1*Ht)*rand(1, N_d);

% Store envelope exponents
var.alpha = env_exp(1);
var.beta = env_exp(2);
var.gamma = env_exp(3);
var.mu = env_exp(4);

% Initialize Library
K = N_d*(N_h(1) + 1)*(N_h(2) + 1)*(N_h(3) + 1)*(N_h(4) + 1); % number of library rows
if mode == "BC"
    Q = {};
    Q.t = zeros(2*K, N);
    Q.n = zeros(K, N);
elseif mode == "NS"
    Q = zeros(3*K, N);
elseif mode == "div"
    Q = zeros(K, N);
end

char_sizes = zeros(N, 1);
if mode == "BC"
    valid_single = {};
    valid_single.n = zeros(N, 1);
    valid_single.t = zeros(N, 1);
else
    valid_single = zeros(N, 1);
end

%% ----------------------------------------------------------------------
%                           CONSTRUCT LIBRARY
% -----------------------------------------------------------------------
k = 0;

% Iterate through integration domains
for np = 1:N_d  

    % Starting corner of integration domain
    x_min = P(1, np);
    y_min = P(2, np);
    z_min = P(3, np);
    t_min = P(4, np);

    % Create local grids
    x = x_full(x_full >= x_min & x_full <= x_min + Hx);
    if mode == "BC"
        y = y_full(1);
    else
        y = y_full(y_full >= y_min & y_full <= y_min + Hy);
    end
    z = z_full(z_full >= z_min & z_full <= z_min + Hz);
    t = t_full(t_full >= t_min & t_full <= t_min + Ht);

    % Create derivative conversion factors
    S_x = 2/(x(end) - x(1));
    S_y = 2/(y(end) - y(1));
    S_z = 2/(z(end) - z(1));
    S_t = 2/(t(end) - t(1));

    % Create weight function grids
    var.x = linspace(-1, 1, length(x));
    if mode == "BC"
        var.y = y;
    else
        var.y = 2*(y - mean([y(1), y(end)]))/(y(end) - y(1));
    end
    var.z = linspace(-1, 1, length(z));
    var.t = linspace(-1, 1, length(t));
    
    % Load velocity fields on integration domain    
    tic
    
    download_data = 1;
    file_loc = sprintf("%s_%s.mat",[file np]);
    if exist(file_loc, 'file') % don't need to load data
        download_data = 0;
    end

    if download_data
        if track
            disp(['Loading the velocities on domain # ', num2str(np), ' from JHU ...'])
        end
        if mode == "BC"
            [U, V, W] = pullVelocity(x, y_full, z, t);
        else
            [U, V, W] = pullVelocity(x, y, z, t);
        end
        save(file_loc, 'U', 'V', 'W');
    else
        if track
            disp(['Loading the velocities on domain # ', num2str(np), ' from disk!'])
        end
        load(file_loc, 'U', 'V', 'W');
    end
    
    U = U+noise*std(U(:))*randn(size(U));
    V = V+noise*std(V(:))*randn(size(V));
    W = W+noise*std(W(:))*randn(size(W));
    
    file_loc_p = sprintf('%s_p_%s.mat',[file np]);
    if mode ~= "BC"
        if download_data
            p = pullPressure(x, y, z, t);
            save(file_loc_p, 'p');
        else
            load(file_loc_p, 'p');
        end
        p = p+noise*std(p(:))*randn(size(p));
    end
    toc

    if np == 1 % compute pieces of characteristic term size
        sz = size(U)
        norm_u = sqrt(U.*2+V.*2+W.*2);
        mean_u = mean(abs(norm_u(:)))
        std_u = sqrt(std(U(:))^2+std(V(:))^2+std(W(:))^2)
        if mode ~= "BC"
            mean_p = mean(abs(p(:)))
            std_p = std(p(:))
            dx = norm([1/S_x 1/S_y 1/S_z])/sqrt(3)
        else
            dx = norm([1/S_x 1/S_z])/sqrt(2)
        end
        dt = 1/S_t
    end
    
    tic
    if mode == "BC"
        % find first normal derivatives
        dU = diff_dim(U, y_full, 2);
        dV = diff_dim(V, y_full, 2);
        dW = diff_dim(W, y_full, 2);
        dyyU = diff_dim(dU, y_full, 2);
        dyyV = diff_dim(dV, y_full, 2);
        dyyW = diff_dim(dW, y_full, 2);
    
        % keep only the boundary terms
        U = U(:, 1, :, :); V = V(:, 1, :, :); W = W(:, 1, :, :);
        dU = dU(:, 1, :, :); dV = dV(:, 1, :, :); dW = dW(:, 1, :, :);
        dyyU = dyyU(:, 1, :, :); dyyV = dyyV(:, 1, :, :); dyyW = dyyW(:, 1, :, :);
    else
        % find first derivatives
        dxp = diff_dim(p, x, 1);
        dyp = diff_dim(p, y, 2);
        dzp = diff_dim(p, z, 3);
        dtU = diff_dim(U, t, 4);
        dtV = diff_dim(V, t, 4);
        dtW = diff_dim(W, t, 4);
        dtp = diff_dim(p, t, 4);
        dxU = diff_dim(U, x, 1);
        dyV = diff_dim(V, y, 2);
        dzW = diff_dim(W, z, 3);
        %dxxp = diff_dim(diff_dim(p, x, 1), x, 1);
        %dyyp = diff_dim(diff_dim(p, y, 2), y, 2);
        %dzzp = diff_dim(diff_dim(p, z, 3), z, 3);
        %dxxU = diff_dim(diff_dim(U, x, 1), x, 1);
        %dyyU = diff_dim(diff_dim(U, y, 2), y, 2);
        %dzzU = diff_dim(diff_dim(U, z, 3), z, 3);
        %dxxV = diff_dim(diff_dim(V, x, 1), x, 1);
        %dyyV = diff_dim(diff_dim(V, y, 2), y, 2);
        %dzzV = diff_dim(diff_dim(V, z, 3), z, 3);
        %dxxW = diff_dim(diff_dim(W, x, 1), x, 1);
        %dyyW = diff_dim(diff_dim(W, y, 2), y, 2);
        %dzzW = diff_dim(diff_dim(W, z, 3), z, 3);
    end

    if track 
        disp('Computing integrals for various weight fns ...')
    end

    % Iterate through weight function modes
    for q = 0:N_h(1)
    for r = 0:N_h(2)
    for s = 0:N_h(3)
    for o = 0:N_h(4)
                
            k = k + 1;
            
            % Put wave numbers into var struct for use in weight_full()
            var.q = q;
            var.r = r;
            var.s = s;
            var.o = o;

            % Iterate through library terms
            for n = 1:N
                krnl_struct = krnl_info(lib_names(n));
                if np+q+s+r+o == 1
                    char_sizes(n) = eval(krnl_struct.size_string);
                end
                if mode == "BC"
                    % Iterate over tangential terms
                    if (np+q+s+r+o == 1)
                         valid_single.t(n) = 1;
                         valid_single.n(n) = 1;
                    end
                    for d=1:2
                        % Compute kernel in pieces
                        krnl = 0;
                        for n_piece = 1:length(krnl_struct.t_pieces)
                            if isfield(krnl_struct, 't_skip')
                                valid_single.t(n) = 0;
                                break
                            end
                            w{1} = 0; w{2} = 0;
                            w{d} = weight_3d(krnl_struct.t_modes(n_piece, :), var);
                            wx = w{1}; wz = w{2};
                            krnl = krnl + eval(krnl_struct.t_pieces(n_piece));
                        end                                

                        % Compute library element
                        if all(krnl==0)
                           integral_val = 0; 
                        else
                           integral_val = trapz(var.x, trapz(var.z, ...
                           trapz(var.t, krnl, 4), 3), 1);
                        end
                        Q.t(2*k+d-2,n) = integral_val;
                    end
                    % Evaluate normal terms
                    krnl = 0;
                    for n_piece = 1:length(krnl_struct.n_pieces)
                        if isfield(krnl_struct, 'n_skip')
                            valid_single.n(n) = 0;
                            break
                        end
                        wy = weight_3d(krnl_struct.n_modes(n_piece, :), var);
                        krnl = krnl + eval(krnl_struct.n_pieces(n_piece));
                    end                                

                    % Compute library element
                    if all(krnl==0)
                       integral_val = 0; 
                    else
                       integral_val = trapz(var.x, trapz(var.z, ...
                       trapz(var.t, krnl, 4), 3), 1);
                    end
                    Q.n(k, n) = integral_val;
                elseif mode == "div"
                    if (np+q+s+r+o == 1) && (any(contains(krnl_struct.pieces, "d"))+any(krnl_struct.modes(:))) 
                        valid_single(n) = 1;
                    end
                    % Compute kernel in pieces
                    krnl = 0;
                    
                    for n_piece = 1:length(krnl_struct.pieces)                                             
                        w = weight_full(krnl_struct.modes(n_piece, :), var);                                        
                        krnl = krnl + eval(krnl_struct.pieces(n_piece));
                    end                                

                    % Compute library element
                    integral_val = trapz(var.x, trapz(var.y, trapz(var.z, ...
                                   trapz(var.t, krnl, 4), 3), 2), 1);
                    %disp(integral_val)
                    Q(k,n) = integral_val;
                elseif mode == "NS"
                % Iterate over dimensions of vector weight
                    if (np+q+s+r+o == 1) && (any(contains(krnl_struct.pieces, "d"))+any(krnl_struct.modes(:)))
                        valid_single(n) = 1;
                    end
                    for d=1:3
                        krnl = 0;
                        for n_piece = 1:length(krnl_struct.pieces)                                             
                            w = weight_full(krnl_struct.modes(n_piece, :), var);                                        
                            krnl = krnl + eval(krnl_struct.pieces(n_piece));
                        end                                

                        % Compute library element
                        integral_val = trapz(var.x, trapz(var.y, trapz(var.z, ...
                                       trapz(var.t, krnl, 4), 3), 2), 1);
                        Q(3*k+d-3,n) = integral_val;
                        % Move to terms corresponding to next component of A
                        krnl_struct = permute_krnl_struct(krnl_struct);
                    end
                end

            end % n
   
    end % o
    end % s
    end % r
    end % q
    
    toc
end % np

end

function [new_string] = permute_dims(old_string)

new_string = regexprep(old_string, {'U', 'V', 'W', 'x', 'y', 'z'}, {'@', '#', '`', '%', '~', '&'});
new_string = regexprep(new_string, {'@', '#', '`', '%', '~', '&'}, {'V', 'W', 'U', 'y', 'z', 'x'});

end

function [new_array] = permute_array(old_array)

new_array = [old_array(3), old_array(1:2), old_array(4)];

end

function [new_struct] = permute_krnl_struct(krnl_struct)

for n_piece = 1:length(krnl_struct.pieces)                                             
    new_struct.modes(n_piece, :) = permute_array(krnl_struct.modes(n_piece, :));                                      
    new_struct.pieces(n_piece) = permute_dims(krnl_struct.pieces(n_piece));
end

end
