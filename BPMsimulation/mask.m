function mask2D = mask(x, y, truncate)

if nargin < 3
    truncate = true;        %default
end
    
mask2D = x(1:end).^2 + y(1:end).^2 <= 1;

if truncate
    mask2D = repmat(mask2D, [1, 1, 2]);    % replicate mask
end

end