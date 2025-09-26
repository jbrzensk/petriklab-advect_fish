%% Test output from 30 days
%
close all; clear all; clc;

load('test_output.mat');

bioLd_smooth = sub_1Dto2D(GRD1,S_Lrg_d(:,31), param);

figure(2); pcolor(bioLd_smooth'); shading interp; title('Ld');colorbar; clim([0 4])
hold on; plot(60, 75,'ro');

load('test_output_allsmooth.mat');
bioLd = sub_1Dto2D(GRD1,S_Lrg_d(:,31), param);

figure(21); pcolor(bioLd'); shading interp; title('Ld no smooth');colorbar;clim([0 4]);
hold on; plot(60, 75,'ro');

%% TryDFT
mask = ~isnan(bioLd);
F = bioLd;

F_filled = bioLd;
F_filled(~mask) = 0;

F_hat = fft2(F_filled);
M_hat = fft2(mask);

F_hat_corrected = F_hat ./ (M_hat + eps );

F_hat_shifted = fftshift(F_hat_corrected);

figure
imagesc(log(abs(F_hat_shifted) + 1))   % +1 avoids log(0)
colormap jet
colorbar
title('2D Fourier Amplitude Spectrum')

[Ny, Nx] = size(F);

kx = (-floor(Nx/2):ceil(Nx/2)-1) / Nx;  % normalized frequency
ky = (-floor(Ny/2):ceil(Ny/2)-1) / Ny;

figure
imagesc(kx, ky, log(abs(F_hat_shifted) + 1))
axis xy
xlabel('k_x (cycles per grid cell)')
ylabel('k_y (cycles per grid cell)')
colorbar
title('2D Spectrum with Frequency Axes')

% Compute power spectrum
P = abs(F_hat_shifted).^2;

% Radial averaging
[KY, KX] = meshgrid(ky, kx);
K = sqrt(KX.^2 + KY.^2);   % radial wavenumber

% Bin into radial shells
nbins = 50;
edges = linspace(0, max(K(:)), nbins);
Pradial = zeros(1, nbins-1);

for i = 1:nbins-1
    mask = (K >= edges(i)) & (K < edges(i+1));
    Pradial(i) = mean(P(mask), 'omitnan');
end

% Plot radial spectrum
k_centers = (edges(1:end-1) + edges(2:end))/2;
figure
loglog(k_centers, Pradial)
xlabel('Wavenumber |k|')
ylabel('Power')
title('Radially Averaged Spectrum')
grid on

checkerboard_ns = checkerboard_index(F);
checkerboard_s = checkerboard_index(bioLd_smooth);

con_ns = checkerboard_contrast(F);
con_s = checkerboard_contrast(bioLd_smooth);

%% Convolution Check
Ffill = inpaint_nans(F);

% 3x3 kernel that responds to + - / - + pattern (centered)
kernel = [ 1 -2  1;
          -2  4 -2;
           1 -2  1 ];   % Laplacian-like with alternating sign sensitivity

convF = conv2(Ffill, kernel, 'same');
figure;
subplot(1,2,1); imagesc(F); colorbar; title('Original field F');
subplot(1,2,2); imagesc(convF); colorbar; title('Checkerboard kernel response (convF)');
colormap('parula');

figure; histogram(convF(:) ./ nanstd(Ffill(:))); title('Histogram of convF normalized by std(F)');
xlabel('convF / std(F)');

[ny,nx] = size(F);
[X,Y] = meshgrid(1:nx,1:ny);
ideal = mod(X+Y,2);         % ±1 checkerboard
%ideal = ideal * 5;%std(Ffill(:));  % scale to have same std as F
conv_ideal = conv2(ideal, kernel, 'same');
score_ideal = sqrt(mean(conv_ideal(:).^2)) / nanstd(ideal(:));
disp(score_ideal);

figure(34); surf(F');
