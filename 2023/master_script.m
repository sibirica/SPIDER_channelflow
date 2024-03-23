%{
Written by Matthew Golden. Questions can be directed to
mgolden30@gatech.edu

PURPOSE:
This script should contain sections for producing all results from the paper. 
This carries out model discovery for 3+1D fluid turbulence

Velocity and pressure data will have to be downloaded manually and placed in 
the h5_cutouts folder.
%}

clear

addpath('functions/');

%global parameters: these are used by all routines
gamma = 1.3;             %sparsification parameter
envelope_power = 8;       %parameter of weight function w = (1-x^2)^envelope_power.
size_vec = 32*ones(1,4 ); %number of points to use in each spacetime subdomain

%seeds for randomly selecting spacetime domains
seed_train = 1;
seed_test  = 5432112345;





%% SCALAR LIBRARY CENTER OF CHANNEL
noise = 0.0;

location = "center";
[G,  labels, scales]  = scalar_library_function( noise, envelope_power, size_vec, location, [], seed_train );

%Find div(u) = 0
[center_div_res_ave, res_std, center_div_cs_ave, cs_std] = compact_regression( G, labels );
n = select_relation( center_div_res_ave, gamma );
print_table_entry( center_div_cs_ave, cs_std, G, n, scales, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
[G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");

%Energy equation
[center_energy_res_ave, res_std, center_energy_cs_ave, cs_std] = compact_regression( G, labels );
n = select_relation( center_energy_res_ave, gamma );
print_table_entry( center_energy_cs_ave, cs_std, G, n, scales, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t E");

%Pressure Poisson
[center_pp_res_ave, res_std, center_pp_cs_ave, cs_std] = compact_regression( G, labels );
n = select_relation( center_pp_res_ave, gamma );
print_table_entry( center_pp_cs_ave, cs_std, G, n, scales, labels );

%Pressure Poisson 3-term
[center_pp_res_ave, res_std, center_pp_cs_ave, cs_std] = compact_regression( G, labels );
n=3;
print_table_entry( center_pp_cs_ave, cs_std, G, n, scales, labels );

%print entry for 2-term PP as well
n = 2;
print_table_entry( center_pp_cs_ave, cs_std, G, n, scales, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");


%% Edge of Channel

location = "edge";
[G, labels, scales] = scalar_library_function( noise, envelope_power, size_vec, location, [], seed_train );

%Find divu
[edge_div_res_ave, res_std, edge_div_cs_ave, cs_std] = compact_regression( G, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
[G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");

%Pressure Poisson
[edge_pp_res_ave, res_std, edge_pp_cs_ave, cs_std] = compact_regression( G, labels );
n = select_relation( edge_pp_res_ave, gamma );
print_table_entry( edge_pp_cs_ave, cs_std, G, n, scales, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");


%Energy equation
[edge_energy_res_ave, res_std, edge_energy_cs_ave, cs_std] = compact_regression( G, labels );
n = select_relation( edge_energy_res_ave, gamma );
print_table_entry( edge_energy_cs_ave, cs_std, G, n, scales, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t u^2");


%% Vector stuff

location = "center";
[center_G, labels, center_scales] = vector_library_function( noise, envelope_power, location, size_vec, seed_train );
[center_NS_res_ave, res_std, center_NS_cs_ave, center_NS_cs_std] = compact_regression( center_G, labels );
n = select_relation( center_NS_res_ave, gamma );
print_table_entry( center_NS_cs_ave, center_NS_cs_std, center_G, n, center_scales, labels );

location = "edge";
[edge_G, labels, edge_scales] = vector_library_function( noise, envelope_power, location, size_vec, seed_train );
[edge_NS_res_ave, res_std, edge_NS_cs_ave, edge_NS_cs_std] = compact_regression( edge_G, labels );
n = select_relation( edge_NS_res_ave, gamma );
print_table_entry( edge_NS_cs_ave, edge_NS_cs_std, edge_G, 4, edge_scales, labels );


%% Make plots

plot_sparsity_curve( center_div_res_ave, edge_div_res_ave, 0 )
%ylim([1e-6 1e-3]); yticks([1e-6, 1e-5, 1e-4, 1e-3]);
drawnow;
saveas( gcf, 'figures/div.pdf' );

plot_sparsity_curve( center_energy_res_ave, edge_energy_res_ave, 1)
%yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0]); 
drawnow;
saveas( gcf, 'figures/energy.pdf' );


plot_sparsity_curve( center_pp_res_ave, edge_pp_res_ave, 0 )
%yticks([ 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0]); drawnow;
drawnow;
saveas( gcf, 'figures/pressure.pdf' );


plot_sparsity_curve( center_NS_res_ave, edge_NS_res_ave, 0 )
drawnow;
saveas( gcf, 'figures/NS.pdf' );

return;





%% Noisy Energy Equation (center of channel)

%specify noise levels we want to add
%noise_levels = [0, 0.01, 0.1];
%location = "center";


%noise_levels = [0, 0.1, 0.5];
noise_levels = [0.1];
location = "center";


% specify the model you are aiming for, in the order you want terms to
% appear in the table. Energy equation is these five terms
target = {"\partial_t E", "u_i \nabla_i p", "\nabla_i(u_i E)", "(\nabla_i u_j)(\nabla_i u_j)", "\nabla^2 E"};


nn = numel(noise_levels); %num noise
nt = numel(target);       %num target

%use three-letter names for the rows of the table
csa = zeros(nn, nt);
unc = zeros(nn, nt);
chi = zeros(nn, nt);


for i = 1:numel(noise_levels)
  fprintf( "iteration %d: noise = %f\n", i, noise_levels(i) );
  noise = noise_levels(i);
  [G,  labels, scales]  = scalar_library_function( noise, envelope_power, size_vec, location, [], seed_train );

  [G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");
  %[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t u^2");
  [G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");

  [res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );
  %gamma = 3;
  n = select_relation( res_ave, gamma );

  print_table_entry( cs_ave, cs_std, G, n, scales, labels );

  [ csa(i,:), unc(i,:), chi(i,:) ] = compute_table_entry( cs_ave, cs_std, G, n, scales, labels, target  )
end

generate_table( csa, unc, chi, target, noise_levels, "tables/energy_coeffs_" + location + ".txt" );

%save('tables/energy_' + location + "_residuals.mat", "res_ave" );


%% Noisy Pressure-Poisson table

%Change between center and edge
location = "edge";

if location == "center"
  %noise levels for center of domain
  noise_levels = [0, 0.1, 0.2];
  % specify the model you are aiming for, in the order you want terms to
  target = {"\nabla^2 p", "\nabla_i \nabla_j (u_i u_j)", "1"};
end

if location == "edge"
  %noise levels for edge of domain
  noise_levels = [0, 1, 5];
  % specify the model you are aiming for, in the order you want terms to
  target = {"\nabla^2 p", "\nabla_i \nabla_j (u_i u_j)"};
end


nn = numel(noise_levels); %num noise
nt = numel(target);       %num target

%use three-letter names for the rows of the table
csa = zeros(nn, nt);
unc = zeros(nn, nt);
chi = zeros(nn, nt);

for i = 1:numel(noise_levels)
  fprintf( "iteration %d: noise = %f\n", i, noise_levels(i) );
  noise = noise_levels(i);

  [G,  labels, scales]  = scalar_library_function( noise, envelope_power, size_vec, location, [], seed_train );

  [G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t E");
  
  % remove these algebraic terms for the noisy case
  % otherwise not fitting physics, fitting noise
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "u^4");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "p^3");

  %[G, labels, scales] = remove_term_by_name( G, labels, scales, "p");
  %[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");

  [res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );
  n = select_relation( res_ave, gamma );

  %force n = 2 to get PP
  %n = 2;

  print_table_entry( cs_ave, cs_std, G, n, scales, labels );

  [ csa(i,:), unc(i,:), chi(i,:) ] = compute_table_entry( cs_ave, cs_std, G, n, scales, labels, target  )
end

generate_table( csa, unc, chi, target, noise_levels, "tables/PP_coeffs_" +location + ".txt" );
save('tables/PP_' + location + "_residuals.mat", "res_ave" );


%% Make plot for energy of noisy data

rc = load("tables/energy_center_residuals.mat", "res_ave");
re = load("tables/energy_edge_residuals.mat",   "res_ave");

res_c = rc.res_ave;
res_e = re.res_ave;

plot_sparsity_curve( res_c, res_e, 1 )
ylim([1e-3 1]);
yticks([ 1e-3 1e-2 1e-1 1e0])
drawnow;
saveas( gcf, 'figures/noisy_energy.pdf' );


%% Make plot for PP of noisy data

rc = load("tables/PP_center_residuals.mat", "res_ave");
re = load("tables/PP_edge_residuals.mat",   "res_ave");

res_c = rc.res_ave;
res_e = re.res_ave;

plot_sparsity_curve( res_c, res_e, 1 )
ylim([1e-2 1]);
yticks([1e-2 1e-1 1e0])
drawnow;
saveas( gcf, 'figures/noisy_PP.pdf' );




%% Make table for noisy data
%This part isnt automated, manually change noise and locations
noise = 0.15; %0.5 for 50%, 1.0 for 100%
location = "center";
size_vec = 32*ones(4,1);
envelope_power = 8;
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec );

[res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

n = select_relation( res_ave, gamma );

% c is mean coefficient value
% s is standard deviation
% x is the relative size of each term
[c,s,x] = get_NS_statistics( G, cs_ave, cs_std, scales, labels, n ); %last arg is number of terms in model

fprintf('%.6e\t', c); fprintf("\n");
fprintf('%.6e\t', s); fprintf("\n");
fprintf('%.6e\t', x); fprintf("\n");





%% 10% noise

noise = 0.1;

location = "center";
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec, seed_train );
[center_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

location = "edge";
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec, seed_train );
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t^2 u_i");
[edge_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

plot_sparsity_curve( center_NS_res_ave, edge_NS_res_ave, 1 )
ylim([1e-4, 1]);
yticks([1e-4 1e-2 1e-0]);
drawnow
saveas( gcf, 'figures/NS_noisy_10.pdf' );


%% 100% noise
noise = 1; 

location = "center";
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec, seed_train );
[center_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

location = "edge";
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec, seed_train );
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t^2 u_i");
[edge_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );


plot_sparsity_curve( center_NS_res_ave, edge_NS_res_ave, 1 )
ylim([5e-3, 1]);
yticks([1e-2 1e-1 1e-0])
drawnow;
saveas( gcf, 'figures/NS_noisy_100.pdf' );


gamma = 1.25;
n = select_relation( center_NS_res_ave, gamma );
n = select_relation(   edge_NS_res_ave, gamma );


%% A loop over beta

envelope_powers = 2:1:16;
size_vec = 32*[1,1,1,1];

res_center = 0*envelope_powers;
res_edge   = 0*envelope_powers;


for i = 1:numel(envelope_powers)
  i
  noise = 0;
  envelope_power = envelope_powers(i);

  location = "center";
  [G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec );
  [center_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels ); 

  location = "edge";
  [G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec );
  [edge_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

  res_center(i) = center_NS_res_ave(1,4);
  res_edge(i)   = edge_NS_res_ave(1,4);
end
save('beta_scaling_NS.mat');

%%
load('beta_scaling_NS.mat');

ms = 100;
scatter(envelope_powers, res_edge, ms, 's', 'filled', 'MarkerFaceColor', 'black' )
hold on
  scatter(envelope_powers, res_center, ms, 's', 'MarkerEdgeColor', 'black', 'LineWidth', 1 );
hold off
set( gca, 'yscale', 'log' );

drawnow
legend({'edge','center'},'location', 'north');
xlabel('\beta');
ylabel('|Gc|');

fs = 16; %fontsize
xlabel( "$$\beta$$", "Interpreter", "latex", "FontSize", fs);
xlim([min(envelope_powers) max(envelope_powers)]);
%changing ylabel to r instead of ||Qc||
ylabel('$$r$$', 'Interpreter', 'latex', 'FontSize', fs, 'Rotation', 0);
%ylabel("$$\|Q{\bf c}\|$$", "Interpreter", "latex", "FontSize", fs);
set(gcf, 'color', 'w');
xlim([2-0.5 14+0.5]);
xticks(2:14);
xticklabels( {"2","","","","","","8","","","","","","14"} )
set(gca,'XTickLabelRotation',0);

ax = gca;
ax.FontSize = fs;
drawnow;

saveas( gcf, 'figures/beta_scaling_NS.pdf');
%export_fig('./figures2/beta_scaling_NS.pdf')

%% Another loop over beta, but for pressure Poisson now


envelope_powers = 2:1:16;
size_vec = 32*[1,1,1,1];

res_center = 0*envelope_powers;
res_edge   = 0*envelope_powers;


for i = 1:numel(envelope_powers)
  i
  noise = 0;
  envelope_power = envelope_powers(i);

  location = "center";
  [G, labels, scales] = scalar_library_function( noise, envelope_power, size_vec, location, [] );
  [G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t u^2");
  [center_pp_res_ave, res_std, center_pp_cs_ave, cs_std] = compact_regression( G, labels );

  location = "edge";
  [G, labels, scales] = scalar_library_function( noise, envelope_power, size_vec, location, [] );
  [G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t u^2");
  [edge_pp_res_ave, res_std, edge_pp_cs_ave, cs_std] = compact_regression( G, labels );

  labels(   edge_pp_cs_ave(:,2) ~= 0 )
  labels( center_pp_cs_ave(:,2) ~= 0 )

  res_center(i) = center_pp_res_ave(1,4);
  res_edge(i)   = edge_pp_res_ave(1,4);
end

save('beta_scaling_pp.mat');


%%

load('beta_scaling_pp.mat');

ms = 100;
scatter(envelope_powers, res_edge, ms, 's', 'filled', 'MarkerFaceColor', 'black' )
hold on
  scatter(envelope_powers, res_center, ms, 's', 'MarkerEdgeColor', 'black', 'LineWidth', 1 );
hold off
set( gca, 'yscale', 'log' );

%semilogy( envelope_powers, res_center, 'Color', 'black' );
%hold on
%semilogy( envelope_powers, res_edge, '--', 'Color', 'black' );
%hold off
drawnow
legend({'edge', 'center'}, 'location', 'north');
xlabel('\beta');
ylabel('|Gc|');

fs = 16; %fontsize
xlabel( "$$\beta$$", "Interpreter", "latex", "FontSize", fs);
%xlim([min(envelope_powers) max(envelope_powers)]);
drawnow
xlim([2-0.5 14+0.5]);
xticks(2:14);
xticklabels( {"2","","","","","","8","","","","","",14} )
set(gca,'XTickLabelRotation',0);
yl = ylabel('$$r$$', 'Interpreter', 'latex', 'FontSize', fs, 'Rotation', 0);
yl.Position(1) = -0.5; %adjust position of y axis a bit
%ylabel("$$\|Q{\bf c}\|$$", "Interpreter", "latex", "FontSize", fs);
yticks([1e-5 1e-3 1e-1])
ylim([5e-6 1e-1]);
set(gcf, 'color', 'w');
ax = gca;
ax.FontSize = fs;
drawnow;

saveas( gcf, 'figures/beta_scaling_pp.pdf');
%export_fig('./figures2/beta_scaling_pp.pdf')





function plot_sparsity_curve( center_res_ave, edge_res_ave, label_y_axis )
  residual_labels = { "$\|Q{\bf c}\|$", "$|Gc|/|G_{c\neq 0}|$", "$|Gc|/\textrm{max}\{|G c_n|\}$" };
  index = 1; %Which residual we use
  
  ms = 16*8*2; %marker size 
  fs = 20; %Font Size


  %lim = size(center_res_ave,2);
  lim = 10;
  %h = plot( 1:lim, center_res_ave(index,1:lim), 's', 'Color', 'black', 'MarkerSize', ms );
  h = scatter( 1:lim, center_res_ave(index,1:lim), ms, 's', 'MarkerEdgeColor', 'black' );

  hold on
    scatter( 1:lim, edge_res_ave(index,1:lim), ms, 's', 'filled', 'MarkerFaceColor', 'black' );
  hold off
  set(get(h,'Parent'), 'YScale', 'log');

  xlabel('$N$', 'interpreter', 'latex', 'FontSize', fs);
  if label_y_axis == 1
    %ylabel( residual_labels{index}, 'interpreter', 'latex', 'FontSize', fs );
    ylabel( '$r\,\,\,\,\,$', 'interpreter', 'latex', 'FontSize', fs, "rotation", 0 );
  end
  xticks(1:10);
  xticklabels({'','2','','','','6', '', '', '', '10'});
  set(gca,'XTickLabelRotation',0);
  
  xlim( [0.5, lim+0.5] );
  ylim( [1e-6, 5]);
  set(gcf,'color','w');
  pbaspect([2 1 1]) 
  drawnow

  %Try fixing the ticks font size
  ax = gca;
  ax.FontSize = fs;
  drawnow

  yticks([1e-5, 1e-3, 1e-1]);
  if label_y_axis == 0
    %yticks([]);
  end
end

function [c,s,x] = get_NS_statistics( G, cs, std, scales, labels, num )
  I = find( cs(:,num) );
  labels(I)

  c = cs(:,num)./scales';
  s = std(:,num)./scales';
  x = vecnorm( G.*cs(:,num)');

  c = c(I); %restrict to nonzero terms
  s = s(I);
  x = x(I);

  normal = max(abs(c));
  c = c/normal;
  s = s/normal;

  normal2 = max(abs(x));
  x = x/normal2;

  if(num == 4)
    %reorder advection and pressure gradient
    c = c([1 3 2 4]);
    s = s([1 3 2 4]);
    x = x([1 3 2 4]);
  end

   if(num == 3)
    %reorder advection and pressure gradient
    c = c([1 3 2 ]);
    s = s([1 3 2 ]);
    x = x([1 3 2 ]);
  end
end

function [G, labels, scales] = remove_term( G, labels, scales, n)
  G(:,n)    = [];
  scales(n) = [];
  labels(n) = [];
end

function [G, labels, scales] = remove_term_by_name( G, labels, scales, name)
  for i = 1:numel(labels)
    if name == labels{i}
      n = i;
    end
  end

  if exist("n")
    [G, labels, scales] = remove_term( G, labels, scales, n);
  end
end


function [res_ave, res_std, cs_ave, cs_std, Gs] = compact_regression( G, labels )
  number_of_subsamples = 128;             %Do parameter estimation 100 times to estimate uncertainty
  fraction_of_windows_in_subsample = 1/2; %use half of windows randomly each time we estimate parameters
  res1 = @(G,c) norm(G*c)/norm(c);
  res2 = @(G,c) norm(G*c)/norm(c)/norm(G);
  res3 = @(G,c) norm(G*c)/max( vecnorm(G*diag(c)) );
  residual_functions = {res1, res2, res3}; %residuals to compute during sparsification
  index = 1; %this picks the residual
  res_func_for_sparsification = residual_functions{index};
  
  
  starting_model = ones( numel(labels), 1 );
  [res_ave, res_std, cs_ave, cs_std, Gs] = reverse_regression3( G, starting_model, ...
                                                             number_of_subsamples, ...
                                                             fraction_of_windows_in_subsample, ...
                                                             residual_functions, ...
                                                             res_func_for_sparsification );
end

function n = select_relation( res_ave, gamma )
  %{
  PURPOSE:
  Decide how many terms an optimal relation ought to have.
  %}

  res_ave = res_ave(1,:); %restrict to |Gc|
  assert( res_ave(end-1)/res_ave(end) < gamma )

  for n = numel(res_ave):-1:1
    if n == 1
      return; %single term model
    end

    ratio = res_ave(n-1)/res_ave(n);
    %[n,ratio] %for debugging
    if ratio > gamma
      return;
    end
  end
end


function print_table_entry( cs_ave, cs_std, G, n, scales, labels  )
  %{
  PURPOSE: 
  The paper has a few tables specifying mean coefficients, uncertainties,
  and characteristic sizes.

  INPUT:
  c - cs_ave and cs_std are outputs of regression.
  n - the number of terms you have selected
  scales - scales vector from library integration

  OUTPUT:
  none. This prints to Command Window
  %}

  fprintf("\n%d term model:\n", n);

  c = cs_ave(:,n)./scales';
  normalization = max(abs(c));
  c = c/normalization;

  s = cs_std(:,n)./scales'/normalization;

  chi = vecnorm(G*diag(cs_ave(:,n)));
  chi = chi./max(chi);

  idx = find(c~=0);
  for j = 1:numel(idx)
    i = idx(j);
    fprintf("%s : c_n = %.8f : s_n = %.8f : chi_n = %.8f\n", labels{i}, c(i), s(i), chi(i));
  end

  eta = norm(G*cs_ave(:,n))/max(vecnorm(G*diag(cs_ave(:,n))));
  fprintf("eta = %e\n", eta );

  fprintf("\n", n);
end


function [cs, un, ch] = compute_table_entry( cs_ave, cs_std, G, n, scales, labels, target  )
  %{
  PURPOSE: 
  The paper has a few tables specifying mean coefficients, uncertainties,
  and characteristic sizes.

  INPUT:
  c - cs_ave and cs_std are outputs of regression.
  n - the number of terms you have selected
  scales - scales vector from library integration

  OUTPUT:
  cs  - average coefficient values
  un  - uncertainty
  ch  - relative size
  %}

  fprintf("\n%d term model:\n", n);

  c = cs_ave(:,n)./scales';
  normalization = max(abs(c));
  c = c/normalization;

  s = cs_std(:,n)./scales'/normalization;

  chi = vecnorm(G*diag(cs_ave(:,n)));
  chi = chi./max(chi);

  %idx is a boolean vector indicating the nonzero terms
  idx = find(c~=0);
  
  %subsample to nonzero terms
  l   = labels(idx)
  c   = c(idx);
  u   = s(idx);
  chi = chi(idx);

  %These are the output arrays to fill
  nt = numel(target);
  cs = zeros(nt,1);
  un = zeros(nt,1);
  ch = zeros(nt,1);

  %Now shape it in the desired target format. Reorder terms
  for i = 1:numel(l)
    %find the corresponding index of the target
    j = 1;
    while(true)
      if( j > nt )
        % No match was found in the target. Error out
        error("The term " + l(i) + " is not in the target model.");
      end 

      [i,j]
      if( l{i} == target{j} )
        break;
      end
      j = j+1;
    end

    %you found the target index
    cs(j) = c(i);
    un(j) = u(i);
    ch(j) = chi(i);
  end
end

function str = print_uncertainty( u )
  % Make a latex string for unceratinty.
  % I'm going to round to the first digit and put in scientific notation
  if( u == 0 )
    str = " 0 ";
    return;
  end

  %first nonzero digit in base-10
  digit = -floor(log10(u));
  
  num = round(10^digit * u);
  if( num == 10)
    digit = digit - 1; 
    num = 1;
  end

  str = "$ " + round( 10^digit*u ) + "\!\times\! 10^{-" + digit + "} $ ";
end

function str = print_chi( chi )
  % Make a latex string for characteristic size chi
  %print first two siginificant digits
  if( chi == 1 )
    str = " 1 ";
    return;
  end

  if( chi == 0)
    str = " 0 ";
    return;
  end

  %first nonzero digit in base-10
  digit = -floor(log10(chi));

  str = sprintf("%0."+ (digit+1) + "f ", chi );
end

function str = print_coeff( c, u)
  %u is the corresponding uncertainty so that we can print an intelligent
  %number of digits
  
  %check for integer coefficient
  if( c == round(c) )
    str = "" + round(c);
    return;
  end

  %first nonzero digit in base-10
  digit = -floor(log10(u));

  str = sprintf("%0."+ (digit) + "f ", c );
end


function generate_table( csa, unc, chi, target, noise_levels, filename )
%This magical function creates LaTeX table

nn = size(csa,1);
nt = size(csa,2);

%make sure sign of c is okay
[~, I] = max( abs(csa') );
signs = zeros(nn,1);
for i = 1:nn
  signs(i) = sign(csa(i,I(i)));
end
csa = csa .*signs;

%write the table to a text file systematically

my_file = fopen( filename, "w");

header = " & $\sigma$ ";
for i =1:nt
  header = header + " & $ " + target{i} + " $ ";
end
header = header + "\\ \hline";

%hack to write a string with 'escape characters'
fprintf( my_file, '%s', header );
fprintf( my_file, "\n");

coeff_string = "$\bar{c}_n$";
for i = 1:nn
  coeff_string = coeff_string + "& " + round(noise_levels(i)*100) + "\% ";
  for j = 1:nt
    coeff_string = coeff_string + "& " + print_coeff( csa(i,j), unc(i,j) ); %sprintf("%.6e", csa(i,j)); 
  end
  coeff_string = coeff_string + " \\ ";
end
coeff_string = coeff_string + "\hline ";
fprintf( my_file, '%s', coeff_string );
fprintf( my_file, "\n");


coeff_string = "$s_n$";
for i = 1:nn
  coeff_string = coeff_string + "& " + round(noise_levels(i)*100) + "\% ";
  for j = 1:nt
    coeff_string = coeff_string + "& " + print_uncertainty( unc(i,j) ); 
  end
  coeff_string = coeff_string + " \\ ";
end
coeff_string = coeff_string + "\hline ";
fprintf( my_file, '%s', coeff_string );
fprintf( my_file, "\n");


coeff_string = "$\chi_n$";
for i = 1:nn
  coeff_string = coeff_string + "& " + round(noise_levels(i)*100) + "\% ";
  for j = 1:nt
    coeff_string = coeff_string + "& " + print_chi(chi(i,j)); %sprintf("%.6e", chi(i,j)); 
  end
  coeff_string = coeff_string + " \\ ";
end
coeff_string = coeff_string + "\hline ";
fprintf( my_file, '%s', coeff_string );
fprintf( my_file, "\n");

fclose(my_file);

fprintf("done\n");
end