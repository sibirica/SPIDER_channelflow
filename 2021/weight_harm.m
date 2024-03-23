% %% -----------------------------------------------------------------------
% %           COMPLEX EXPONENTIAL WEIGHT FUNCTION HARMONICS
% % ------------------------------------------------------------------------
function g = weight_harm(x, k, q)
    
%{
g : k^th derivative of the q^th mode of the complex exponentional
%}
g = (-1i*pi*q)^k*exp(-1i*pi*q*x);
    
end