function [res_ave, res_std, cs_ave, cs_std, Gs] = reverse_regression3( G, starting_model, ...
                                                                   number_of_subsamples, ...
                                                                   fraction_of_windows_in_subsample, ...
                                                                   residual_functions, ...
                                                                   res_func_for_sparsification )
  %{
  PURPOSE:
  This function starts with a model as sparse as you's like, and both adds
  and removes terms to it to see how the residual changes. This can be used
  to perform sequential thresholding (starting with full library and
  dropping terms one at a time) or reverse regression (start with a coarse
  model and add terms one at a time).

                                                                   
                                                                   
  INPUT: apologies for the long-winded names. I'm trying to make this
         understandable
  
  G - integrated library matrix. Called Q in some papers, but I had a naming
      conflict in the active nematics project.
      size(G) = [number_of_windows, number_of_library_terms]
  
  starting_model - a column vector with nonzero terms you'd like to start
                   with. The coeefificent values do not matter at all, as
                   this function will determine the best values anyways.
                   For example, ones(number_of_library_terms,1) will start
                   with the entire library and perform sequential
                   thresholding.
  
  number_of_subsamples - This tries to estimate uncertainty in coefficients
                         and residuals. number_of_subsamples is the number
                         of times to perform regression to estimate these
                         uncertainties. I usually use ~100.
                                                                   
  fraction_of_windows_in_subsample - I usually use 1/2. To get an idea of
                                     uncertainty, regression will happen
                                     number_of_subsamples times, using this
                                     fraction of windows (selected
                                     randomly)

  residual_functions - a cell of function handles {res1, res2, ...} where 
                       each function should take arguments (G,c) and return
                       a single number.
                                                                   
  res_func_for_sparsification - a function that will be minimized greedily 
                                during sparsification
  
  OUTPUT:
                                                                   
  res_ave - average residuals. 
            size(res_ave) = [ numel(residual_functions), num_library_terms] 
  
  res_std - standard deviation of residuals. Same size as above
  
  cs_ave  - average coefficient values. Columns correspond to models.
  
  cs_std  - standard deviation of coefficient values  
             
                                                                   
                                                                   
  EXAMPLE USAGE:
                                                                   
  starting_model = cs(1,:); %start with a decent coarse model from
                            %combinatoric search
  number_of_subsamples = 100; %Do parameter estimation 100 times to estimate uncertainty
  fraction_of_windows_in_subsample = 1/2; %use half of windows randomly each time we estimate parameters

  res1 = @(G,c) norm(G*c);
  res2 = @(G,c) norm(G*c)/norm(c)/norm(G);
  res3 = @(G,c) norm(G*c)/max( vecnorm(G*diag(c)) );
  residual_functions = {res1, res2, res3}; %residuals to compute during sparsification

  res_func_for_sparsification = res1;



  [res_ave, res_std, cs_ave, cs_std] = reverse_regression2( G, starting_model, ...
                                                                   number_of_subsamples, ...
                                                                   fraction_of_windows_in_subsample, ...
                                                                   residual_functions, ...
                                                                   res_func_for_sparsification );

  %}

  [nw, nl] = size(G); %nw == num_windows, nl == num_library_terms
  nr = numel(residual_functions); %Number of different residuals we want to compute
  
  Gs = cell(nl,1);

  res_ave = zeros( nr, nl);
  res_std = zeros( nr, nl);
  cs_ave  = zeros( nl, nl);
  cs_std  = zeros( nl, nl);


  %For reproducibility, I should seed randsample
  s = RandStream('dsfmt19937');
  %According to documentation, this is the SIMD-Oriented Fast Mersenne Twister.


  
  %1. Start with current number of terms
  idx = (starting_model ~= 0); %These terms are nonzero
  Gr = G( :, idx );            %Restrict G to this subspace (hence Gr)
  
  Gs{ size(Gr,2) } = Gr;
  
  c_temp   = zeros( sum(idx), number_of_subsamples ); %temporary coefficient matrix for doing statistics
  res_temp = zeros( nr,       number_of_subsamples ); %temporary residual matrix for doing statistics
  for i = 1:number_of_subsamples
    G_temp = Gr( randsample(s, nw, round(nw*fraction_of_windows_in_subsample) ), : ); %Further restrict GR by randomly subsampling
    
    [~,~,V] = svd(G_temp, 'econ'); %compute coefficients with SVD (minimizing |G_temp*c|/|c|)
    c = V(:,end);
    
    %Note that the sign of c is arbitrary. To prevent this from ruining
    %statistics, make sure the dot product with the first result is
    %positive.
    if i ~= 1
      c = sign( c'*c_temp(:,1) )*c;
    end
    
    c_temp( :, i ) = c;
    
    %Compute all the residuals
    for j = 1:nr
      res_func      = residual_functions{j};
      res_temp(j,i) = res_func(G_temp, c);
    end
  end
  
  res_ave( :, sum(idx) ) = mean(res_temp, 2);
  res_std( :, sum(idx) ) = std( res_temp, 0, 2);
  cs_ave(idx, sum(idx) ) = mean( c_temp, 2);
  cs_std(idx, sum(idx) ) = std(  c_temp, 0, 2);

  
  %2. Add more terms in order that least increases eta
  initial_size_of_model = sum(idx);
  c_best = starting_model;
  for j=1:(nl - initial_size_of_model)
    idx  = (c_best ~= 0);
    idx2 = (c_best == 0);
    
    missing = find(idx2);
    
    eta_competition = zeros( sum(idx2), 1 );
    for k=1:sum(idx2)
      idx_temp = idx;
      idx_temp(missing(k)) = 1; %Add a missing term to this list of library terms
      
      G_temp = G( :, idx_temp );
      [~,~,V] = svd(G_temp, 'econ');
      c = V(:,end);
    
      eta_competition(k) = res_func_for_sparsification(G_temp,c);
    end
    
    [~, k] = min( eta_competition );
    idx(missing(k)) = 1;
    
    Gr = G(:, idx);
    Gs{ size(Gr,2) } = Gr;


    [~,~,V] = svd(Gr, 'econ');
    c_best = V(:,end);
    
    c_temp   = zeros(sum(idx), number_of_subsamples);
    for i = 1:number_of_subsamples
      G_temp = Gr( randsample(s, nw, round(nw*fraction_of_windows_in_subsample) ), : );
    
      [~,~,V] = svd(G_temp, 'econ');
      c = V(:,end);
      if i ~= 1
        %Make sure we use consistent signs
        c = sign( c'*c_temp(:,1) )*c;
      end
      c_temp( :, i ) = c;

      for k = 1:nr
        res_func      = residual_functions{k};
        res_temp(k,i) = res_func(G_temp, c);
      end
    end
  
    res_ave( :, sum(idx) ) = mean( res_temp, 2);
    res_std( :, sum(idx) ) = std(  res_temp, 0, 2);
    cs_ave(idx, sum(idx) ) = mean( c_temp, 2);
    cs_std(idx, sum(idx) ) = std(  c_temp, 0, 2);
    
    %update c_best
    c_best = cs_ave( :, sum(idx) );
  end
  
    
  %3. Trim the model down
  idx = (starting_model ~= 0);
  initial_size_of_model = sum(idx);
  c_best = starting_model;
  for j=1:(initial_size_of_model-1)
    idx  = (c_best ~= 0);
    idx2 = (c_best == 0);
    
    missing = find(idx); %I am lying for this part. missing is actually the ones not missing
    
    eta_competition = zeros( sum(idx), 1 );
    for k=1:sum(idx)
      idx_temp = idx;
      idx_temp(missing(k)) = 0; %remove term to this list of library terms
      
      G_temp = G( :, idx_temp );
      if( numel(G_temp) == 0 )
        break; 
      end
      
      [~,~,V] = svd(G_temp, 'econ');
      c = V(:,end);
      
      eta_competition(k) = res_func_for_sparsification(G_temp,c);
    end
    [~, k] = min( eta_competition );
    idx(missing(k)) = 0;
    
    Gr = G(:, idx);
    
    Gs{ size(Gr,2) } = Gr;

    c_temp   = zeros(sum(idx), number_of_subsamples);
    for i = 1:number_of_subsamples
      G_temp = Gr( randsample(s, nw, round(nw*fraction_of_windows_in_subsample) ), : );
          
      [~,~,V] = svd(G_temp, 'econ');
      c = V(:,end);
      if i ~= 1
        %Make sure we use consistent signs
        c = sign( c'*c_temp(:,1) )*c;
      end
      c_temp( :, i ) = c;
    
      for k = 1:nr
        res_func      = residual_functions{k};
        res_temp(k,i) = res_func(G_temp, c);
      end
    end
  
    res_ave( :, sum(idx) ) = mean( res_temp, 2);
    res_std( :, sum(idx) ) = std(  res_temp, 0, 2);
    cs_ave(idx, sum(idx) ) = mean( c_temp, 2);
    cs_std(idx, sum(idx) ) = std(  c_temp, 0, 2);
    
    %update c_best
    c_best = cs_ave( :, sum(idx) );
  end
end