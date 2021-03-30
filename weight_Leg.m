%% -----------------------------------------------------------------------
%              LEGENDRE WEIGHT FUNCTION HARMONICS
% ------------------------------------------------------------------------
function g = weight_Leg(x,k,q)
    %{
    g : k^th derivative of the q^th Legendre polynomial on x = [-1,1]
    %}
    c = flipud(LegendrePoly(q));
    if k > 0 
        for l = 1:k
            c = (1:length(c)-1)'.*c(2:end); % take derivative
        end
    end
    g = 0;
    for n = 0:length(c)-1
        g = g + c(n+1)*x.^n;
    end

end