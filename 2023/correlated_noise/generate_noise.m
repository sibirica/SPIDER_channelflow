%{
N= 32;
k_max = 4;

noise = correlated_noise( [N,N,N], k_max );
for i = 1:N
  imagesc(squeeze(noise(:,:,i)));
  colorbar();
  clim([-1 1]);
  drawnow;
end
%}

function noise = correlated_noise( Ns, k_max )
  %Produce periodic correlated noise

  k = cell(3,1);
  for i = 1:3
    k{i} = 0:Ns(i)-1; 
    k{i}(k{i}>Ns(i)/2) = k{i}(k{i}>Ns(i)/2) - Ns(i);
  end

  noise =      (2*rand(Ns) - 1) ...
        + 1i * (2*rand(Ns) - 1);

  k_sq = reshape( k{1} / k_max(1), [Ns(1),1,1]).^2 + ...
       + reshape( k{2} / k_max(2), [1,Ns(2),1]).^2 + ...
       + reshape( k{3} / k_max(3), [1,1,Ns(3)]).^2;

  clf; 
  tiledlayout(2,2);
  nexttile

  noise( k_sq < 1 ) = 0;
  noise(1,1,1,1) = 0; %set the mean to zero!
  imagesc(abs(noise(:,:,1)));
  title("Fourier coefficients k_z = 0, unscaled, cutoff = " + k_max(1) );
  axis square
  colorbar();

  noise = real(fftn(noise));
  noise = noise / std(noise, 0, "all");
  nexttile
  imagesc( noise(:,:,1) );
  title("real space noise (z=0)");
  axis square
  colorbar();

  nexttile
  histogram(noise);
  drawnow;
end