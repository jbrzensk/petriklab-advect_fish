function score = checkerboard_contrast(F)
% returns RMS of checkerboard-filtered field (higher -> more alternation)
% Handle NaNs by filling using inpaint_nans
Ffill = inpaint_nans(F);

% Laplacian-thing with alternating sign sensitivity
% 3x3 kernel that responds to + - / - + pattern (centered)
kernel = [ 1 -2  1;
          -2  4 -2;
           1 -2  1 ];

% Convolution
convF = conv2(Ffill, kernel, 'same');

% Normalize by local variance, or anything really
score = sqrt(nanmean(convF(:).^2)) / (nanstd(Ffill(:)) + eps);

end