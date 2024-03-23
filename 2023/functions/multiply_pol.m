function pol_out = multiply_pol( pol1, pol2 )
  pols = cell(3,1);
  for i = 1:3
    %apparently conv does polynomial multiplication
    p1 = flip(pol1(:,i));
    p2 = flip(pol2(:,i));
    pols{i} = flip( conv( p1, p2  ));
  end
  
  ls = zeros(3,1);
  for i= 1:3
    ls(i) = numel( pols{i} );
  end
  
  pol_out = zeros( max(ls), 3 );
  for i=1:3
    pol_out(1:ls(i), i) = pols{i};
  end
end