function vals = SPIDER_integrate( data, derivs, grid, corners, size_vec, pol )
  %{
  PURPOSE:
  This function is meant to integrate space-time fields with nice
  polynomial weight functions to be used for model discovery. Here are some
  key points.

  1. Integration by parts is handled ANALYTICALLY
  2. Non-uniform grid spacings are allowed. Integration is approximated
     with the trapezoid rule.
  3. The "data" argument should only have dimensions corresponding to space
     and time. If you have a "velocity" matrix of size [3,128,128,...]
     where the first index is the components, you are going to have to feed
     components in one by one. This is a pain, but you can write a wrapper
     that uses this function as a building block.
  
  
  INPUT:
  data    - the data you need weighted integrals of

  derivs  - [1,1] for instance uses integration by parts to integrate two
            derivatives along the first dimension (could be x or y
            depending on your convention)

  grid    - a cell containing [---points in x---;
             ---points in y---;
             ...]

  corners - a nxnum_windows matrix containing lower corners of your integration
            domains.

  size_vec- [24 24 48] as an example. The size of the cube (in points) you want to integrate over.
  
  pol - a [l, 3] matrix describing the polynomial weight function

 
  OUTPUT:
  vals - a vector of integrated values.
  %}

  n           = size(corners,1);        %dimensionality of the data
  num_windows = size(corners,2);        %number of spacetime windows we want to integrate over
  vals        = zeros(num_windows,1);   %vals are the integrated values
           
  %Do a sanity check and make sure the number of coordinates given for the
  %spacetime cubes matches the dimensionality of the data
  data_size = size(data);
  if numel(data_size) == 2
    %If this is the case, then your data is stored as a simple matrix.
    %If data is 1D, it still is stored as a matrix since MATLAB is built on
    %matrices. Check for one of these dimensions being trivial.
    data_size( data_size == 1 ) = []; %kill trivial dimension
  end
  assert( numel(data_size) == n ); 
  assert( numel(grid)      == n );
  
  
  % subsample_indices is a cell that will be used to subsample the data for 
  % spacetime domains.
  subsample_indices = cell(n, 1);
  for s = 1:n
    subsample_indices{s} = 0:size_vec(s)-1; 
  end
  
  % Do integration by parts to figure out the new polynomial weight.
  poly_sign = 1;
  l = size(pol,1);  %upper bound of degree of the polynomial weight
  derivs0 = derivs; %save a copy since we will destroy the derivs object
  while( numel(derivs) ~= 0 )
    d = derivs(1); %direction to differentiate
      
    poly_sign = -poly_sign; %each integration by parts accumulates a -1
    for ll = 2:l
      pol( ll-1, d ) = (ll-1) * pol( ll, d );
    end
    pol(ll, d) = 0; %Fixing a bug! Don't forget to kill the last term manually!
    
    derivs(1) = []; %Delete this derivative as we have taken it analytically
  end
  
 
  %Now we are ready to start the integration loop
  for i = 1:num_windows
      
    %STEP 1: subsample the data
    specific_subsample = cell(n,1);
    for j = 1:n
      specific_subsample{j} = corners(j,i) + subsample_indices{j}; 
    end
    small_data = data( specific_subsample{:} );
    %^small data now contains the chunk of data to integrate
    
    %STEP 2: Loop over dimension and construct the polynomial weight
    conversion_factor = poly_sign; 
    %^this will account for coordinate changes and the sign change from
    %integration by parts.
    for j = 1:n
      x_full = grid{j};
      x_spec = x_full( specific_subsample{j} );
      % Let x be the grid passed in (assumed to have units), and let y be 
      % Some subset of x translated and rescaled to the interval [-1,1]
      % Our weight functions are naturally functions of y instead of x
      % y = mx+b
      m = 2/(x_spec(end) - x_spec(1));
      b = -(x_spec(end) + x_spec(1))/(x_spec(end) - x_spec(1));
      
      y = m*x_spec + b; %[-1,1] coordinates
      conversion_factor = conversion_factor * m^sum(derivs0 == j);
      %^this accounts for changing derivatives from d/dx = m d/dy
      
      poly_weight_1d = polyval( fliplr( pol(:,j)' ), y ); %evaluate the polynomial weight for this dimension
      poly_weight_1d = reshape( poly_weight_1d, [numel(y), 1]); %make sure its a column vector
      
      %compute grid spacings to left and right
      hr      = circshift(x_spec, -1) - x_spec;
      hr(end) = 0;
      hl      = x_spec - circshift(x_spec, 1);
      hl(1)   = 0;
      hl      = reshape( hl, [numel(y), 1] );
      hr      = reshape( hr, [numel(y), 1] );

      small_data = poly_weight_1d.* small_data; %add polynomial weight
      small_data = (hr+hl)/2     .* small_data; %add trapezoid rule weights
      small_data = sum(small_data);             %sum along first dimension
      small_data = kill_first_dimension( small_data ); %kill first dimension
    end
    vals(i) = small_data * conversion_factor;
  end
end


function data = kill_first_dimension( data )
  %Squeeze almost always does what you want, until you have something of
  %size [1,10]. Squeeze does nothing in this case, so we can just take the
  %transpose.
  data = squeeze(data);
  if(size(data,1) == 1)
    data = data';
  end
end