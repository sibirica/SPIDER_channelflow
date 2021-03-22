function [dF] = diff_dim(F, y, dim)
% 3-point forward/center/backwards differencing to approximate first derivative

dF = zeros(size(F));

% boundary points
y0 = y(1); y1 = y(2); y2 = y(3);
c0 = (2*y0-y1-y2)/((y0-y1)*(y0-y2));
c1 = (2*y0-y0-y2)/((y1-y0)*(y1-y2));
c2 = (2*y0-y0-y1)/((y2-y0)*(y2-y1));
if dim == 1
    dF(1, :, :, :) = c0*F(1, :, :, :)+c1*F(2, :, :, :)+c2*F(3, :, :, :);
elseif dim == 2
    dF(:, 1, :, :) = c0*F(:, 1, :, :)+c1*F(:, 2, :, :)+c2*F(:, 3, :, :);
elseif dim == 3
    dF(:, :, 1, :) = c0*F(:, :, 1, :)+c1*F(:, :, 2, :)+c2*F(:, :, 3, :);
else
    dF(:, :, :, 1) = c0*F(:, :, :, 1)+c1*F(:, :, :, 2)+c2*F(:, :, :, 3);
end

for i=2:length(y)-1
    y0 = y(i-1); y1 = y(i); y2 = y(i+1);
    c0 = (2*y1-y1-y2)/((y0-y1)*(y0-y2));
    c1 = (2*y1-y0-y2)/((y1-y0)*(y1-y2));
    c2 = (2*y1-y0-y1)/((y2-y0)*(y2-y1));
    if dim == 1
        dF(i, :, :, :) = c0*F(i-1, :, :, :)+c1*F(i, :, :, :)+c2*F(i+1, :, :, :);
    elseif dim == 2
        dF(:, i, :, :) = c0*F(:, i-1, :, :)+c1*F(:, i, :, :)+c2*F(:, i+1, :, :);
    elseif dim == 3
        dF(:, :, i, :) = c0*F(:, :, i-1, :)+c1*F(:, :, i, :)+c2*F(:, :, i+1, :);
    else
        dF(:, :, :, i) = c0*F(:, :, :, i-1)+c1*F(:, :, :, i)+c2*F(:, :, :, i+1);
    end
end

y0 = y(end-2); y1 = y(end-1); y2 = y(end);
c0 = (2*y2-y1-y2)/((y0-y1)*(y0-y2));
c1 = (2*y2-y0-y2)/((y1-y0)*(y1-y2));
c2 = (2*y2-y0-y1)/((y2-y0)*(y2-y1));
if dim == 1
    dF(end, :, :, :) = c0*F(end-2, :, :, :)+c1*F(end-1, :, :, :)+c2*F(end, :, :, :);
elseif dim == 2
    dF(:, end, :, :) = c0*F(:, end-2, :, :)+c1*F(:, end-1, :, :)+c2*F(:, end, :, :);
elseif dim == 3
    dF(:, :, end, :) = c0*F(:, :, end-2, :)+c1*F(:, :, end-1, :)+c2*F(:, :, end, :);
else
    dF(:, :, :, end) = c0*F(:, :, :, end-2)+c1*F(:, :, :, end-1)+c2*F(:, :, :, end);
end

end