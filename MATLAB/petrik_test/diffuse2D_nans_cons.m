function Zsmooth = diffuse2D_nans_cons(Z, D)
% diffuse2D_nans_cons  Apply Laplacian diffusion to 2D data with NaN mask
%
%   Zsmooth = diffuse2D_nans_cons(Z, D)
%
%   Performs one explicit Laplacian diffusion step:
%       Z_new = Z + D * ∇²Z
%
%   NaNs are treated as solid boundaries (no flux through them).
%   D controls how strong the diffusion is (typical stable range: D ≤ 0.25).
%
%   Inputs:
%     Z : 2D array (with NaNs for invalid cells)
%     D : diffusion coefficient (scalar)
%
%   Output:
%     Zsmooth : diffused field

    if nargin < 2
        D = 0.1;
    end

    tol = 1e-12; % tolerance for comparison after diffusing.

    mask = ~isnan(Z);    % valid data mask, basically where nans are.

    Zfilled = Z;
    Zfilled(~mask) = 0;

    % Laplacian kernel
    % laplacianKernel = [0  1  0;
    %                    1 -4  1;
    %                    0  1  0];
    laplacianKernel = [1   4  1;
                       4 -20  4;
                       1   4  1] / 6;
    %kernel = ones(window);

    % Compute sum of neighbors and count of valid neighbors
    lapZ = conv2(Zfilled, laplacianKernel, 'same');
    lapZ(~mask) = 0;

    % Find Neighbors
    validNeighbors = conv2(double(mask), abs(laplacianKernel), 'same');
    % Avoid dividing by zero
    validNeighbors(validNeighbors == 0) = 1;

    % Normalize Laplacian so that NaN edges are handled properly
    lapZ = lapZ ./ validNeighbors;

    % Compute diffused field
    Zsmooth = Zfilled + D * lapZ;

    % Compute total before and after, renormalize to conserve total mass
    total_before = sum(Z(:), 'omitnan');
    total_after  = sum(Zsmooth(:), 'omitnan');
    
    % Error here is consistently 0.2%. 
    if (abs(total_after-total_before)) > tol
        fprintf('Adjusting smooth vals, sum= %4.8f, diff= %4.8f\n',total_before, total_after-total_before);
        Zsmooth = Zsmooth * (total_before / total_after);
    end

    % Apply mask to keep NaNs
    Zsmooth(~mask) = NaN;
end
