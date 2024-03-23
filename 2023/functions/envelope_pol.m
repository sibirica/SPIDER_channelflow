function pol = envelope_pol(m,n)
  syms x;
  pol_1d = flip( coeffs( (1-x)^m * (x+1)^m, 'all' ) );
  pol_1d = pol_1d'; %make column vector
  
  pol = repmat( pol_1d, [1,n] );
  
  pol = double(pol);
end