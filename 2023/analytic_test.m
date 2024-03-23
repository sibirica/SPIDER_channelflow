%{
Let's do an analytic test 

\int_{-1}^1 dx (1-x^2)^n cos(x) dx = \sqrt{\pi}
{}_0\tilde{F}_1(;n+3/2;-1/4) \gamma(n+1)

where 
%}

%Time to finally test my scripts for 1D data

clear;

addpath("functions/");

trials = 8;

N = 32;

min_pow = 2;
max_pow = 14;

noise_level = 1e-8;

M = max_pow - min_pow + 1;
ps = min_pow:max_pow;

estimate = zeros(M, trials);
analytic = zeros(M, 1);
estimate2 = zeros(M, trials);


fprintf('loading data from center...\n');
addpath('functions/');
[~, ~, ~, ~, ~, ~, y, ~, ~] = read_h5_channelflow_data_center();

for trial = 1:trials
  trial
  seed = trial*N;
  rng(seed);

  x = linspace(-1,1,N)';

  noise = noise_level * (2*rand(N,1) - 1);

  f = cos(2*pi*x) + noise;

  grid = {x};
  corners = [1];
  size_vec = [N];

  dimension = 1;



  



  for i = 1:M
    p = ps(i);
    pol = envelope_pol( p, dimension );

    estimate(i,trial) = SPIDER_integrate( f, [], grid, corners, size_vec, pol );
    analytic(i) = pi^(-p)*gamma(p+1)*besselj(p+1/2, 2*pi);
  end

  y = y(1:N); %select same size grid
  y = y - ( max(y) + min(y))/2; %center it
  y = y/max(y);
  grid = {y};

  f = cos(2*pi*y) + noise;


  for i = 1:M
    p = ps(i);
    pol = envelope_pol( p, dimension );

    estimate2(i,trial) = SPIDER_integrate( f, [], grid, corners, size_vec, pol );
  end
end



%% Make figure

clf

err1 = mean( abs(estimate - analytic), 2 );
err2 = mean( abs(estimate2 - analytic), 2 );

ms = 100; %marker size

scatter(ps, err1', ms, 's', 'filled', 'MarkerFaceColor', 'black' )
hold on
  scatter(ps, err2', ms, 's', 'MarkerEdgeColor', 'black', 'LineWidth', 1 );
  %xline(8);
hold off
set(gca,'YScale', 'log');
drawnow

fs = 16; %fontsize
legend({'uniform grid', 'nonuniform grid'}, 'location', 'NorthWest');
xlabel( "$$\beta$$", "Interpreter", "latex", "FontSize", fs);
xlim([min_pow-0.5 max_pow+0.5]);
xticks(min_pow:max_pow)
xticklabels(["2", "", "", "", "", "", "8", "", "", "", "", "", "14", "", ""]);
set(gca,'XTickLabelRotation',0);
yticks([1e-9, 1e-7, 1e-5])
ylim([1e-10 1e-4])
%ylabel("$$|I_\beta - I_{\beta, \textrm{trap}}|$$", "Interpreter", "latex", "FontSize", fs);
ylabel("$$\Delta I$$", "Interpreter", "latex", "FontSize", fs);
set(gcf, 'color', 'w');

ax = gca;
ax.FontSize = fs;
drawnow;

saveas( gcf, './figures/integration_error.pdf' );