%{
Written by Matthew Golden. Questions can be directed to
mgolden30@gatech.edu (although I might have a new gatech email soon)

PURPOSE:
This script should contain sections for producing all results from the
channelflow paper. This carries out model discovery for 3+1D fluid turbulence
%}

clear

addpath('functions/');
addpath('export/');

%global parameters: these are used by all routines
gamma = 1.25;            %sparsification parameter
envelope_power = 8;      %parameter of weight function w = (1-x^2)^envelope_power.
size_vec = 32*ones(1,4); %number of points to use in each spacetime subdomain


%% SCALAR LIBRARY CENTER OF CHANNEL
noise = 0.0;

location = "center";
[G,  labels, scales]  = scalar_library_function( noise, envelope_power, size_vec, location, [] );

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
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t u^2");

%Pressure Poisson
[center_pp_res_ave, res_std, center_pp_cs_ave, cs_std] = compact_regression( G, labels );
n = select_relation( center_pp_res_ave, gamma );
print_table_entry( center_pp_cs_ave, cs_std, G, n, scales, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");

return;

%% Edge of Channel

location = "edge";
[G, labels, scales] = scalar_library_function( noise, envelope_power, size_vec, location, [] );

%Find divu
[edge_div_res_ave, res_std, edge_div_cs_ave, cs_std] = compact_regression( G, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
[G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");

%Pressure Poisson
[edge_pp_res_ave, res_std, edge_pp_cs_ave, cs_std] = compact_regression( G, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");

%Energy equation
[edge_energy_res_ave, res_std, edge_energy_cs_ave, cs_std] = compact_regression( G, labels );
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t u^2");


%% Vector stuff
noise = 0;
location = "center";
[center_G, labels, center_scales] = vector_library_function( noise, envelope_power, location, size_vec );
[center_NS_res_ave, res_std, center_NS_cs_ave, center_NS_cs_std] = compact_regression( center_G, labels );

location = "edge";
[edge_G, labels, edge_scales] = vector_library_function( noise, envelope_power, location, size_vec );

[edge_NS_res_ave, res_std, edge_NS_cs_ave, edge_NS_cs_std] = compact_regression( edge_G, labels );


%% Make plots
addpath("export/")

plot_sparsity_curve( center_div_res_ave, edge_div_res_ave, 1 )
%ylim([1e-6 1e-3]); yticks([1e-6, 1e-5, 1e-4, 1e-3]);
drawnow;
export_fig('figures2/div.pdf')

plot_sparsity_curve( center_energy_res_ave, edge_energy_res_ave, 1)
%yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0]); 
drawnow;
export_fig('figures2/energy.pdf')

plot_sparsity_curve( center_pp_res_ave, edge_pp_res_ave, 1 )
%yticks([ 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0]); drawnow;
drawnow;
export_fig('figures2/pressure.pdf')

plot_sparsity_curve( center_NS_res_ave, edge_NS_res_ave, 1 )
%yticks([ 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0]);
drawnow;
export_fig('figures2/NS.pdf')

return;


%% Pareto curve 3.2 from review paper

for noise = [0:10]/10

noise

location = "center";

[center_G, labels, center_scales] = vector_library_function( noise, envelope_power, location, size_vec );
[center_NS_res_ave, res_std, center_NS_cs_ave, center_NS_cs_std] = compact_regression( center_G, labels );

[G2, labels, center_scales] = vector_library_function( noise, envelope_power, location, size_vec );

clf

%tiledlayout(1,2)
%nexttile
semilogy( 1:size(center_G, 2), center_NS_res_ave(1,:), 'LineWidth', 3 );
hold on
  semilogy( 1:size(center_G, 2), vecnorm(G2*center_NS_cs_ave), 'LineWidth', 3 );
hold off
legend( {'train', 'test'} )
ylabel('r');
xlabel('complexity N');
%{
nexttile
semilogy( 1:size(center_G, 2), vecnorm( (center_G-G2)*center_NS_cs_ave), 'LineWidth', 3 ) ;
%}
title("noise = " + noise);
drawnow;
saveas(gcf, "pareto_NS_noise_" + noise + ".png" );

end


%% Pareto curve for Pressure Poisson


for noise = [0:10]/10

noise

location = "center";

location = "center";
[G,  labels, scales]  = scalar_library_function( noise, envelope_power, size_vec, location, [] );
[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
[G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");
%[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t E");

[res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

[G2,  labels, scales]  = scalar_library_function( noise, envelope_power, size_vec, location, [] );
[G2, labels, scales] = remove_term_by_name( G2, labels, scales,   "\nabla_i u_i");
[G2, labels, scales] = remove_term_by_name( G2, labels, scales, "p \nabla_i u_i");
%[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");
[G2, labels, scales] = remove_term_by_name( G2, labels, scales, "\partial_t E");

clf

semilogy( 1:size(G, 2), res_ave(1,:), 'LineWidth', 3 );
hold on
  semilogy( 1:size(G, 2), vecnorm(G2*cs_ave), 'LineWidth', 3 );
hold off
legend( {'train', 'test'} )
ylabel('r');
xlabel('complexity N');
title("noise = " + noise);
drawnow;
saveas(gcf, "pareto_PP_noise_" + noise + ".png" );

end






%% Noisy Pressure-Poisson (center of channel)

noise_levels = [0, 0.01, 0.1]; %for energy equation

for i = 1:numel(noise_levels)
  fprintf( "iteration %d: noise = %f\n", i, noise_levels(i) );
  noise = noise_levels(i);
  location = "center";
  [G,  labels, scales]  = scalar_library_function( noise, envelope_power, size_vec, location, [] );

  [G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla_i u_i");
  [G, labels, scales] = remove_term_by_name( G, labels, scales, "p \nabla_i u_i");
  %[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t u^2");
  [G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");


  %Pressure Poisson
  [center_pp_res_ave, res_std, center_pp_cs_ave, cs_std] = compact_regression( G, labels );
  %gamma = 3;
  n = select_relation( center_pp_res_ave, gamma );
  print_table_entry( center_pp_cs_ave, cs_std, G, n, scales, labels );
  %[G, labels, scales] = remove_term_by_name( G, labels, scales,   "\nabla^2 p");
end






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
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec );
[center_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );


location = "edge";
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec );
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t^2 u_i");
[edge_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

plot_sparsity_curve( center_NS_res_ave, edge_NS_res_ave, 1 )
ylim([1e-4, 1]);
yticks([1e-4 1e-2 1e-0]);
drawnow
saveas( gcf, 'figures2/NS_noisy_10.pdf' );


%% 100% noise
noise = 1; 

location = "center";
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec );
[center_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

location = "edge";
[G, labels, scales] = vector_library_function( noise, envelope_power, location, size_vec );
[G, labels, scales] = remove_term_by_name( G, labels, scales, "\partial_t^2 u_i");
[edge_NS_res_ave, res_std, cs_ave, cs_std] = compact_regression( G, labels );

%%
plot_sparsity_curve( center_NS_res_ave, edge_NS_res_ave, 1 )
ylim([5e-3, 1]);
yticks([1e-2 1e-1 1e-0])
drawnow;
saveas( gcf, 'figures2/NS_noisy_100.pdf' );



gamma = 1.25;
n = select_relation( center_NS_res_ave, gamma )
n = select_relation(   edge_NS_res_ave, gamma )


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

semilogy( envelope_powers, res_center, 'Color', 'black' );
hold on
semilogy( envelope_powers, res_edge, '--', 'Color', 'black' );
hold off
drawnow
legend({'center', 'edge'},'location', 'north');
xlabel('\beta');
ylabel('|Gc|');

fs = 16; %fontsize
xlabel( "$$\beta$$", "Interpreter", "latex", "FontSize", fs);
xlim([min(envelope_powers) max(envelope_powers)]);
ylabel("$$\|Q{\bf c}\|$$", "Interpreter", "latex", "FontSize", fs);
set(gcf, 'color', 'w');
xlim([2 14]);
xticks(2:14);
xticklabels( {"2","","","","","","8","","","","","","14"} )
set(gca,'XTickLabelRotation',0);

ax = gca;
ax.FontSize = fs;
drawnow;
export_fig('./figures2/beta_scaling_NS.pdf')

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

semilogy( envelope_powers, res_center, 'Color', 'black' );
hold on
semilogy( envelope_powers, res_edge, '--', 'Color', 'black' );
hold off
drawnow
legend({'center', 'edge'}, 'location', 'north');
xlabel('\beta');
ylabel('|Gc|');

fs = 16; %fontsize
xlabel( "$$\beta$$", "Interpreter", "latex", "FontSize", fs);
%xlim([min(envelope_powers) max(envelope_powers)]);
drawnow
xlim([2 14]);
xticks(2:14);
xticklabels( {"2","","","","","","8","","","","","",14} )
set(gca,'XTickLabelRotation',0);

ylabel("$$\|Q{\bf c}\|$$", "Interpreter", "latex", "FontSize", fs);
yticks([1e-5 1e-3 1e-1])
set(gcf, 'color', 'w');
ax = gca;
ax.FontSize = fs;
drawnow;

export_fig('./figures2/beta_scaling_pp.pdf')





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
    yticks([]);
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

  fprintf("\n", n);
end

