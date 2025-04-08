%{
This script simulates the detection of thermal residuals from a collimated beam 
by adding a faint temperature excess to a synthetic noisy background. The goal 
is to visually assess the detectability of narrow thermal traces (Bonaque traces) 
under varying noise conditions using spatial binning and normalized profile analysis.

The simulation generates synthetic thermal maps with Gaussian noise, injects a 
beam-induced signal, and processes the result via binning and Z-score profile 
extraction. The output reproduces the core structure of Figure 3 in the associated 
paper, illustrating the effect of increasing noise levels on trace detectability.

NOTE: Due to the use of randomized noise (randn), the exact visual output of Figure 3 
cannot be identically reproduced. However, the qualitative behavior and structure 
of the result remain consistent across runs.

Created by Sergio Bonaque-González, PhD
Postdoctoral Fellow (Viera y Clavijo Program)
Department of Applied Physics
University of La Laguna
sbonaque@ull.edu.es

April, 2025
%}
clearvars; close all; clc; % Clear environment and prepare workspace

%% Parameters
img_size = 3000;                  % Image resolution (pixels)
T_background = 15;               % Background temperature [K]
sigma_noise_all = [0.001, 0.002, 0.003, 0.004];  % Thermal noise std devs
dT_beam = 0.001;                 % Beam-induced temperature increase [K]
beam_width_pix = 1;             % Beam width in pixels
block_size = 20;                % Spatial binning size

%% --- Initialization ---
n_cases = length(sigma_noise_all);
T_binned_all = cell(1, n_cases);
profile_z_all = cell(1, n_cases);
max_z_all = zeros(1, n_cases);
center_binned = round(img_size / 2 / block_size); % Beam center (after binning)

%% --- Create figure canvas ---
figure('Color','w','Position',[100, 100, 1600, 1200]);

for idx = 1:n_cases
    sigma_noise = sigma_noise_all(idx);

    %1. Generate background map + beam
    rng(1);  % For reproducibility
    T_map = T_background + sigma_noise * randn(img_size);

    % Superimpose beam trace (vertical line)
    center = round(img_size/2);
    half_w = floor(beam_width_pix/2);
    T_map(:, center-half_w:center+half_w) = T_map(:, center-half_w:center+half_w) + dT_beam;

    % Optional smoothing for better visualization (instrumental PSF)
    T_map_smooth = imgaussfilt(T_map, 2);

    % High-resolution map (row 1)
    subplot(3, n_cases, idx);
    imagesc(T_map_smooth); axis image off;
    colormap(gca, 'hot'); c = colorbar;
    title(sprintf('\\sigma = %.3f K', sigma_noise));

    %2. Binning to simulate resolution
    M = size(T_map,1);
    N = size(T_map,2);
    M_cropped = floor(M/block_size)*block_size;
    N_cropped = floor(N/block_size)*block_size;
    T_crop = T_map(1:M_cropped, 1:N_cropped);
    T_binned = zeros(M_cropped/block_size, N_cropped/block_size);

    for i = 1:block_size:M_cropped
        for j = 1:block_size:N_cropped
            block = T_crop(i:i+block_size-1, j:j+block_size-1);
            T_binned((i-1)/block_size+1, (j-1)/block_size+1) = mean(block(:));
        end
    end

    T_binned_smooth = imgaussfilt(T_binned, 2);

    subplot(3, n_cases, idx + n_cases);
    imagesc(T_binned_smooth); axis image off;
    colormap(gca, 'hot'); colorbar;
    title('Binned map (20×20)');

    %3. Z-profile (column-averaged)
    profile = mean(T_binned, 1);
    profile_z = (profile - mean(profile)) / std(profile);
    profile_z_all{idx} = profile_z;
    max_z_all(idx) = max(profile_z);

    subplot(3, n_cases, idx + 2*n_cases);
    plot(profile_z, 'k-', 'LineWidth', 1.5); hold on;
    plot(xlim, [3 3], '--r');  % detection threshold
    plot(xlim, [0 0], ':k');
    plot(center_binned, 0, 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(center_binned+1, 2.5, 'True beam', 'Color', 'b', 'FontSize', 8);
    ylim([-4 4]);
    title(sprintf('Normalized profile (max z = %.2f)', max_z_all(idx)));
    xlabel('Column index'); ylabel('Z-score');
end

%% --- Save final figure ---
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 16 12]);
print(gcf, 'Figure_5', '-dpng', '-r600');
print(gcf, 'Figure_5', '-dpdf');
