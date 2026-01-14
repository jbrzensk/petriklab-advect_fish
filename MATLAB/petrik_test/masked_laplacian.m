function diff_phi = masked_laplacian(phi, kappa, dx, dy)
% masked_laplacian 2D Laplacian with NaN mask and variable dx, dy
%
% Uses the flux-form of the Laplacian which can handle irregular grids,
% land points (given as NaNs in the data). This sums the diffusive fluxes
% through the available faces ( valid neighbors ). This has the advantage
% of doing all points at the same time ( no looping ).
%
% INPUT
%       phi, dx, dy : 2D arrays (same size), NaN in phi = land
%       kappa       : diffusion coefficient - non-dimensionalized!!
% 
% OUTPUT
%       lap         : Laplacian, NaN on land
%
% This is 
% BY: Jared Brzenski
% Jan 08, 2026
%

    % Error tolerance to renormalize
    tol = 1e-8;

    % Initialize lap with NaNs
    lap = NaN(size(phi));

    % Cell areas
    A = dx .* dy;

    % Shift fields (periodic in x) to get values on faces
    phiE = circshift(phi, [-1  0]);
    phiW = circshift(phi, [ 1  0]);
    phiN = circshift(phi, [ 0 -1]);
    phiS = circshift(phi, [ 0  1]);

    % Grid spacing in each of the four directions
    dxE = circshift(dx, [-1  0]);
    dxW = circshift(dx, [ 1  0]);
    dyN = circshift(dy, [ 0 -1]);
    dyS = circshift(dy, [ 0  1]);

    % Valid neighbors
    vE = ~isnan(phiE);
    vW = ~isnan(phiW);
    vN = ~isnan(phiN);
    vS = ~isnan(phiS);
    % is current location valid?
    vC = ~isnan(phi);

    % Harmonic mean of grid spacing at faces
    dx_e = 2 ./ (1./dx + 1./dxE);
    dx_w = 2 ./ (1./dx + 1./dxW);
    dy_n = 2 ./ (1./dy + 1./dyN);
    dy_s = 2 ./ (1./dy + 1./dyS);

    % Zero out land faces
    dx_e(~vE | ~vC) = NaN;
    dx_w(~vW | ~vC) = NaN;
    dy_n(~vN | ~vC) = NaN;
    dy_s(~vS | ~vC) = NaN;

    % Gradients at faces
    dphi_e = (phiE - phi) ./ dx_e;
    dphi_w = (phi - phiW) ./ dx_w;
    dphi_n = (phiN - phi) ./ dy_n;
    dphi_s = (phi - phiS) ./ dy_s;

    % Remove land contributions
    dphi_e(~vE | ~vC) = 0;
    dphi_w(~vW | ~vC) = 0;
    dphi_n(~vN | ~vC) = 0;
    dphi_s(~vS | ~vC) = 0;

    % Flux divergence ( Lap * phi )
    lap = ( dy .* (dphi_e - dphi_w) + dx .* (dphi_n - dphi_s) ) ./ A;

    % Diffusion ( phi + Lap(phi))
    diff_phi = phi + ( kappa * lap );

    % Check for conservation
    total_before = sum(phi(:), 'omitnan');
    total_after  = sum(diff_phi(:), 'omitnan');

    % Error here is consistently 
    if (abs(total_after-total_before)) > tol
        fprintf('Adjusting mask_laplacian vals, sum= %4.8f, diff= %4.8f\n',total_before, total_after-total_before);
        diff_phi = diff_phi * (total_before / total_after);
    end

end
