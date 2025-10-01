function idx = checkerboard_index(F)
% Correlation value for the normalized grid to that of a checkerboard
% 1  : perfect correlation ( checkerboard )
% 0  : no corelation
% -1 : inverse perfect checkerboard

% F: 2D array (may contain NaN)
[ny, nx] = size(F);

% Build alternating mask: +1 -1 +1 -1 ...
[X, Y] = meshgrid(1:nx, 1:ny);
mask = (-1).^(X + Y);

% Remove global mean or local mean? use anomaly w.r.t global mean:
A = F - mean(F(~isnan(F)),'omitnan');

% Flatten and remove NaNs
valid = ~isnan(A);
v = A(valid);
m = mask(valid);

% normalized dot product
idx = sum(v(:).*m(:)) / sqrt(sum(v(:).^2) * sum(m(:).^2));
% idx in [-1,1]: 1 = perfect matching checkerboard, -1 = inverted checkerboard

end