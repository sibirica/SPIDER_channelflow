function [Xi, lambda, best_term, lambda1] = SparseReg(Theta, char_sizes, valid_single, opts)

% compute sparse regression on 0 = Theta * Xi

    % read options
    if nargin>3
        if isfield(opts,'threshold')
            threshold = opts.threshold;
        else
            threshold = "multiplicative";
        end
        if isfield(opts,'brute_force')
            brute_force = opts.brute_force;
        else
            brute_force = 1;
        end
        if isfield(opts,'delta')
            delta = opts.delta;
        else
            delta = 1e-15;
            %delta = 0;
        end
        if isfield(opts,'gamma')
            gamma = opts.gamma;
        else
            gamma = 3;
        end
        if isfield(opts,'epsilon') % parameter for error-based method
            epsilon = opts.epsilon;
        else
            epsilon = 1e-2;
        end
        if isfield(opts,'verbose')
            verbose = opts.verbose;
        else
            verbose = 0;
        end
        if isfield(opts,'n_terms')
            n_terms = opts.n_terms;
        else
            n_terms = -1;
        end
    else
        threshold = "threshold";
        brute_force = 1;
        gamma = 3;
        delta = 1e-15;
        %delta = 0;
        verbose = 0;
    end

    [h,w] = size(Theta);
    
    if nargin>2
        for term = 1:w
            Theta(:, term) = Theta(:, term) / char_sizes(term); % renormalize by char. size
        end
    end
    if nargin<4
        valid_single = ones(w, 1);
    end

    [~,Sigma,V] = svd(Theta);
    Xi = V(:,end);
    if verbose
        Sigmas = Sigma(Sigma(:)>0);
        V, log(Sigmas/min(Sigmas))
    end
    lambda = norm(Theta*Xi);
    if verbose
       lambda 
    end

    % find best one-term model as well
    for term = 1:w
        nrm(term) = norm(Theta(:, term))/valid_single(term);
        [lambda1, ind_single] = min(nrm);
        if verbose
            nrm(term)
        end
    end

    smallinds = zeros(w, 1);
    margins = zeros(w, 1); % increases in residual per time step
    lambdas = zeros(w, 1);
    lambdas(1) = lambda;
    if threshold ~= "multiplicative"
        Xis = {}; % record coefficients
    end
    for i=1:min(100, w)
      if threshold~="multiplicative"
         Xis{i} = Xi; 
      end
      if brute_force
         % product of the coefficient and characteristic size of library function
         res_inc = ones(w,1)*Inf;
      end
      for p_ind = 1:w
         if brute_force
             if smallinds(p_ind)==0
                % Try dropping each term
                smallinds_copy = smallinds;
                smallinds_copy(p_ind) = 1;
                Xi_copy = Xi;
                Xi_copy(p_ind) = 0;
                [~,~,V] = svd(Theta(:, smallinds_copy==0));
                Xi_copy(smallinds_copy==0) = V(:,end);
                res_inc(p_ind) = norm(Theta*Xi_copy)/lambda;
             end 
         else
            col = Theta(:, p_ind);
            % project out other columns
            for q_ind = 1:w
                if (p_ind ~= q_ind) && smallinds(q_ind)==0
                   other_col = Theta(:, q_ind);
                   col = col - dot(col, other_col)/norm(other_col)^2*other_col;
                end
            end
            product(p_ind) = norm(Xi(p_ind)*col./sqrt(sum(Theta(:,:).^2,2)));
            %product(p_ind) = norm(Xi(p_ind)*col);
         end
      end
      
      if brute_force
         [Y, I] = min(res_inc);
         margins(i) = Y;
         if verbose
            res_inc
            if threshold ~= "multiplicative"
                i, lambda % for pareto plot
            end
         end
         if (Y<=gamma) || (threshold ~= "multiplicative")
            smallinds(I) = 1;
            if sum(smallinds==0)==1
               break;
            end
            Xi(I) = 0;
            [~,~,V] = svd(Theta(:, smallinds==0));
            Xi(smallinds==0) = V(:,end);
            lambda = norm(Theta*Xi);
            lambdas(i+1) = lambda;
         else
            Y, I
            break
         end
      else
         product(smallinds==1)=Inf;
         [Y, I] = min(product);
         smallinds(I) = 1;
         if sum(smallinds==0)==0
            break;
         end
         if verbose
            product'
         end
         Xi_old = Xi;
         Xi(smallinds==1) = 0;    % set negligible terms to 0
         [~,~,V] = svd(Theta(:,smallinds==0));
         Xi(smallinds==0) = V(:,end);
         lambda_old = lambda;
         lambda = norm(Theta*Xi);
         lambdas(i+1) = lambda;
         margin = lambda/lambda_old;
         if verbose
            lambda, margin 
         end
         margins(i) = margin;
         if (margin > gamma && lambda>delta && threshold=="multiplicative")
            I
            Xi = Xi_old;
            break
         end
      end   
    end
    Xis{i+1} = Xi;
    if threshold=="pareto"
       [Y_mar, I_mar] = max(margins);
       if n_terms>1
          I_mar = length(margins)-n_terms+1;
       end
       if verbose
           margins
           Y_mar, I_mar
       end
       stopping_point = I_mar-1;
       Xi = Xis{I_mar}; %stopping_point
       lambda = norm(Theta*Xi);
    elseif threshold=="error"
       I_sm = max(find(lambdas>epsilon*lambda1, 1, 'first')-1, 1);
       if isempty(I_sm)
          I_sm = 1; 
       end
       Xi = Xis{I_sm}; %stopping_point
       lambda = norm(Theta*Xi);
    end
    
    % now compare single term and sparsified model
    if verbose
        lambda, lambda1
        lambdas
    end
    best_term = ind_single;
    if nargin>2
        Xi = Xi ./ char_sizes; % renormalize by char. size
    end
end