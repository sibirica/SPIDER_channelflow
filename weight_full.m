%% -----------------------------------------------------------------------
%                       FULL 3D WINDOWING FUNCTION
% ------------------------------------------------------------------------
function W = weight_full(k, var)
%{
    Combine base windows and harmonics in 1D, then assemble into 3D
    k = [kx, ky, kz, kt]: order of derivative(s)

    - Take derivative of composite functions using Liebniz Rule
    
    NOTE: can also use weight_harm() instead of weight_Leg()
%}

% X dependence
wx = 0;
for n = 0:k(1)
    p = weight_poly(var.x, var.alpha, k(1) - n); 
    g = weight_Leg(var.x, n, var.q); 
    wx = wx + nchoosek(k(1), n)*p.*g;  
end

% Y dependence
wy = 0;
for n = 0:k(2)
    p = weight_poly(var.y, var.beta, k(2) - n);
    g = weight_Leg(var.y, n, var.r);
    wy = wy + nchoosek(k(2), n)*p.*g;
end

% Z dependence
wz = 0;
for n = 0:k(3)
    p = weight_poly(var.z, var.gamma, k(3) - n);
    g = weight_Leg(var.z, n, var.s);
    wz = wz + nchoosek(k(3), n)*p.*g;
end

% T dependence
wt = 0;
for n = 0:k(4)
    p = weight_poly(var.t, var.mu, k(4) - n);
    g = weight_Leg(var.t, n, var.o);
    wt = wt + nchoosek(k(4), n)*p.*g;
end

% Make 1D components 4D
[wX, wY, wZ, wT] = ndgrid(wx, wy, wz, wt);        

% combine all components
W = wX.*wY.*wZ.*wT;                         

end