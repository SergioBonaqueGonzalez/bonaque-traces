%{
This script performs a parametric simulation of thermal signatures induced
by highly collimated electromagnetic beams propagating through different
interstellar structures. It calculates the expected temperature increase
(DeltaT) in various media, based on beam geometry, wavelength, absorption
coefficients, power, and emission duration. It also computes the
signal-to-noise ratio (SNR) by accounting for intrinsic thermal fluctuations
in the ISM and instrumental sensitivity limits.

The output includes DeltaT matrices and SNR maps across multiple configurations,
as well as graphical visualizations for each distance–power pair, enabling
a robust assessment of detectability conditions for beam-induced thermal
residuals (Bonaque traces).

MATLAB R2016 or later recommended (due to graphical functions and array operations)

Created by Sergio Bonaque-González, PhD
Postdoctoral Fellow (Viera y Clavijo Program)
Department of Applied Physics
University of La Laguna
sbonaque@ull.edu.es

April, 2025
%}
clearvars; close all; clc; % Clear environment and prepare workspace
% Physical constants
G = 6.6743e-11;               % [m^3/kg/s^2]
M = 0.35 * 1.989e30;          % [kg]
c = 2.99792458e8;             % [m/s]
AU = 1.496e11;                % [m]
R_sun = 6.96e8;               % [m]
pc_to_m = 3.086e16;           % [m]
k_B = 1.380649e-23;           % Boltzmann constant [J/K].

%% Simulation parameters: distances, divergences, and emitter geometry
L_pc_sel = [16, 75, 230, 920, 15767];  % Selected propagation distances [pc]
b_vals = linspace(1, 5, 5) * R_sun;    % Impact parameters [m], 1–5 solar radii
theta_vals = [1e-7, 1e-8, 1e-9, 1e-10]; % Fixed beam divergences [rad]
D_emit = 27e3;                         % Emitting section width for post-focal beam [m]
wavelengths_m = [200e-9, 500e-9, 1e-6, 300e-6, 1e-2, 3e-2, 0.18, 0.21]; %Selected wavelengths [m]
wavelength_labels = {'200nm','500nm','1um','300um','1cm','3cm','18cm','21cm'};
powers_W = [40, 400, 4000, 4e4, 4e5, 4e6, 4e7, 4e8, 4e9, 4e10, 4e11]; % Selected beam powers [W]
% powers_W = [4e8, 4e9, 4e10, 4e11]; % High powers [W]
% powers_W = [40, 400, 4000, 4e4, 4e5, 4e6, 4e7]; % Non-high powers [W]
threshold_dT = 1e-3; % Minimum DeltaT [K] to display in console (1 mK)
ring_width_ratio = 0.0000388; % Ring thickness as a fraction of solar radius
all_beam_labels = ['Postfocus avg', strcat('Collimated \theta=', {'1e-7','1e-8','1e-9','1e-10'}), '27x27 km beam postfocus', 'Collimated 27x27 km beam' ];

% === USER-DEFINED EMISSION DURATION FOR SIMULATION ===
duration_s = 1e7;  % [s] Adjust to 1e9 for Figure 4, 1e7 for Figures 2 & 5

%% Definition of interstellar structures and absorption properties
structures(1).name = 'Cold molecular cloud core';
structures(1).size_ly = [0.1, 1];
structures(1).n_H = 1e4;
structures(1).T = [10, 20];
structures(1).filling_factor = 1e-4;
structures(1).cooling_time_yr = 1e7;
structures(1).alpha = [2, 1.5, 0.1, 0.005, 1e-4, 5e-5, 1e-5, 1e-5];

structures(2).name = 'Giant molecular cloud';
structures(2).size_ly = [10, 100];
structures(2).n_H = 500;
structures(2).T = [10, 20];
structures(2).filling_factor = 0.01;
structures(2).cooling_time_yr = 1e6;
structures(2).alpha = [0.5, 0.4, 0.05, 0.002, 1e-4, 5e-5, 1e-5, 1e-5];

structures(3).name = 'Cold neutral medium';
structures(3).size_ly = [5, 50];
structures(3).n_H = 40;
structures(3).T = [50, 100];
structures(3).filling_factor = 0.05;
structures(3).cooling_time_yr = 5e5;
structures(3).alpha = [0.2, 0.1, 0.01, 0.0005, 0.0005, 0.0002, 1e-6,1e-6];

structures(4).name = 'Galactic cirrus (diffuse dust)';
structures(4).size_ly = [1, 10];
structures(4).n_H = 1;
structures(4).T = [15, 30];
structures(4).filling_factor = 0.01;
structures(4).cooling_time_yr = 5e4;
structures(4).alpha = [1, 0.8, 0.3, 0.001, 1e-4, 5e-5, 1e-5, 1e-5];

structures(5).name = 'Circumstellar dust';
structures(5).size_ly = [0.001, 0.01];
structures(5).n_H = 1e6;
structures(5).T = [50, 500];
structures(5).filling_factor = 0;
structures(5).cooling_time_yr = 1e3;
structures(5).alpha = [5, 4, 3, 1, 0.01, 0.005, 0.001, 0.001];

%%Save the database if desired
% save('ISM_structures.mat', 'structures', 'wavelengths_m', 'wavelength_labels');

%% Simulation
% Initialization
L_m_sel = L_pc_sel * pc_to_m;
n_powers = length(powers_W);
W_postfocus = NaN(length(b_vals), length(L_pc_sel));  % [pc]
arc_fraction = D_emit / (2 * pi * R_sun);  % Angular fraction of the Einstein ring covered by a 'D_emit' km emitter

% CASE 1: Calculates the post-focal beam envelope for rays grazing a lensing star
% at various impact parameters (b = 1–5 R_sun). Each ray experiences a
% deflection Theta = 4GM/(c²b), resulting in a focal distance f = b/Theta.
% Beyond the focal region, the beam expands linearly with Theta. The width
% at each selected distance is computed as W = 2·tan(Theta)·(L - f), expressed in parsecs.
for i = 1:length(b_vals)
    b = b_vals(i);
    theta = 4 * G * M / (c^2 * b);       % rad
    f = b / theta;                       % [m]
    
    for j = 1:length(L_m_sel)
        d = L_m_sel(j) - f;
        if d > 0
            W_postfocus(i, j) = 2 * tan(theta) * d / pc_to_m;  % [pc]
        else
            W_postfocus(i, j) = NaN;
        end
    end
end

W_postfocus_mean = nanmean(W_postfocus, 1);  % Averaging in b (ignoring NaNs)


% CASE 2: Ideal collimated beams with fixed divergence (no gravitational lensing at destination)
W_collimated = NaN(length(theta_vals), length(L_pc_sel));  % [pc]

for i = 1:length(theta_vals)
    theta = theta_vals(i);
    for j = 1:length(L_m_sel)
        W_collimated(i,j) = 2 * theta * L_m_sel(j) / pc_to_m;  % [pc]
    end
end

% Show results (optional)
disp('Mean post-focus divergence (b = 1–5 R_sun), beam width in pc:')
disp(array2table(W_postfocus_mean, ...
    'VariableNames', strcat('d_', strtrim(cellstr(num2str(L_pc_sel'))))))

disp('Collimated beams (theta = 1e-7 to 1e-10), beam width in pc:')
disp(array2table(W_collimated, ...
    'VariableNames', strcat('d_', strtrim(cellstr(num2str(L_pc_sel')))), ...
    'RowNames', {'1e-7','1e-8','1e-9','1e-10'}))

% CASE 3: Post-focal beam divergence from a finite 27 km emitting section
% This case simulates emission from a 27x27 km region on the Einstein ring.
theta = 4 * G * M / (c^2 * R_sun);     % angular deflection [rad].
f = R_sun / theta;                      % focal length [m]

% Initialize results
W_27km_m = NaN(size(L_m_sel));         % width in meters
W_27km_pc = NaN(size(L_m_sel));        % width in parsecs

% Calculation of the beamwidth at these distances
for k = 1:length(L_m_sel)
    if L_m_sel(k) >= f
        d = L_m_sel(k) - f;
        W = D_emit + 2 * (D_emit / f) * d; % beam width at distance d, expanding linearly from focus
        W_27km_m(k) = W;
        W_27km_pc(k) = W / pc_to_m;
    end
end

% Show results (optional)
fprintf('\nBeam width from section 27x27 km:\n');
for k = 1:length(L_pc_sel)
    fprintf('Distance = %6.0f pc -> Beam width = %.3e m = %.3e pc\n', ...
        L_pc_sel(k), W_27km_m(k), W_27km_pc(k));
end

% CASE 4: Ideal collimated beams with fixed divergence (no gravitational lensing at destination), emitted from a 27×27 km aperture
% This case reuses the divergence geometry from CASE 2 (collimated beams with theta=1e-10), but rescales the power spatially
% according to the ratio of effective emitting area (27×27 km) to that of a full annular ring.
% It represents a collimated segment of the ring, pointing directly at the target, without any post-focal divergence.


fprintf('\n----------------------------------------\n');

% --------------------------------------------
% Energy deposition and temperature rise across all beam–structure–wavelength configurations
% --------------------------------------------
% All beam geometries to be analyzed, expressed as full beamwidths [pc].
% Rows: 1 = post-focus average (grazing rays), 2–5 = fixed divergence beams,
% row 6 = post-focal 27×27 km source, row 7 = collimated 27×27 km ring emitter.
all_beam_widths_pc = [W_postfocus_mean; W_collimated; W_27km_pc; W_collimated(4,:)];  %

% Initialize result matrix log10(DeltaT)
num_beams = length(all_beam_labels);
num_structures = length(structures);
num_distances = length(L_pc_sel);
num_wavelengths = length(wavelengths_m);
dT_log10_matrix = -Inf(num_beams, num_structures * num_wavelengths, n_powers, num_distances);
count = 0;

for s = 1:length(structures)
    struc = structures(s);
    n_H = struc.n_H;                                 % [cm^-3]
    l_ly = max(struc.size_ly);                       % [ly]
    l_pc = l_ly* 0.3066;                             % [pc]
    l_m = l_pc * 9.461e15;                           % [m]
    cooling_time_s = struc.cooling_time_yr * 365.25 * 86400;
    effective_duration_s = min(duration_s, cooling_time_s);
    n_H_m3 = n_H * 1e6;                              % [m^-3]
    
    for w = 1:length(wavelengths_m)
        col_idx = (s-1)*num_wavelengths + w;         % Linear index for matrix column corresponding to this structure–wavelength pair
        alpha_pc = struc.alpha(w);                   % [pc^-1]
        f_abs = 1 - exp(-alpha_pc * l_pc);           % absorbed fraction
        
        for b = 1:size(all_beam_widths_pc, 1)
            for d = 1:length(L_pc_sel)
                R_pc = all_beam_widths_pc(b, d)/2;
                R_m = R_pc * 3.086e16;
                % Select beam area depending on the geometry:
                % - Cases 1–5: only a narrow radial slice carries power (ring_width_ratio)
                % - Case 6: full beam post-focal (limited angular segment)
                % - Case 7: collimated beam from 27 km ring (limited angular segment)
                if b < 6
                    A_m2 = pi * R_m^2 * ring_width_ratio;  % Cases 1-5: only a radial fraction of the ring carries energy
                elseif b == 6
                    A_m2 = pi * R_m^2 * arc_fraction;      % Case 6: post-focus, angular fraction of the ring.
                elseif b == 7
                    A_m2 = pi * R_m^2 * arc_fraction;      % Caso 7: collimated, angular fraction of the ring
                end
                
                % Calculate number of hydrogen atoms within beam volume, and resulting thermal capacity
                N = A_m2 * l_m * n_H_m3;
                C = N * (3/2) * k_B;
                
                for i = 1:n_powers
                    P = powers_W(i);
                    E_total = P * effective_duration_s;
                    E_abs = f_abs * E_total;
                    dT = E_abs / C;
                    dT_log10_matrix(b, col_idx, i, d) = floor(log10(dT)); % Store log10(?T) in discrete form (floor) for clearer visualization; change if needed
                    if dT > threshold_dT
                        % If temperature rise exceeds display threshold, print case summary
                        fprintf('%.3e K | %s | Struct: %-30s | Lambda = %-6s | P = %.1e W | d = %.0f pc\n', ...
                            dT, all_beam_labels{b}, struc.name, wavelength_labels{w}, P, L_pc_sel(d));
                        count = count+1;
                        
                    end
                end
            end
        end
    end
end
% Prepare column labels: one per structure–wavelength pair
col_labels = {};
for s = 1:num_structures
    for w = 1:num_wavelengths
        col_labels{end+1} = [structures(s).name ' ' wavelength_labels{w}]; %#ok<AGROW>
    end
end


% Plot heatmaps of log10(DeltaT) for each combination of distance and power
for d = 1:num_distances
    for i = 1:n_powers
        figure('Color','w','Position', [100, 100, 1600, 400]);
        imagesc(dT_log10_matrix(:,:,i,d), [-8, 2]);
        colormap(jet); colorbar;
        caxis([-8 2]);
        
        % Axes and labels
        set(gca, 'XTick', 1:length(col_labels), 'XTickLabel', col_labels, 'XTickLabelRotation', 90);
        set(gca, 'YTick', 1:length(all_beam_labels), 'YTickLabel', all_beam_labels);
        title(sprintf('log_{10}(\\DeltaT [K]) at %g pc, %g s, %g W', L_pc_sel(d), duration_s, powers_W(i)), 'FontSize', 12);
        
        % Overlay integer log10(DeltaT) values inside each cell
        for r = 1:size(dT_log10_matrix,1)
            for c = 1:size(dT_log10_matrix,2)
                val = dT_log10_matrix(r,c,i,d);
                if ~isinf(val)
                    text(c, r, sprintf('%d', val), ...
                        'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'k');
                end
            end
        end
    end
end

fprintf('\nTotal number of detectable cases (DeltaT > %.1e K): %d\n', threshold_dT, count);


%% Generate de figures of the paper
% Optionally close all figures except those relevant to the paper
disp('Note: All figures are being closed to display only the paper figures.');
disp('If this behavior is not desired, comment out the "close all" command below (Line 291)');
close all;
%% Figure 2
%% --- Final: Show and save selected figure for 230 pc, 1e7 s, 4e6 W ---
% Identify indices
if duration_s == 1e7
    
    target_d = find(L_pc_sel == 230);
    target_p = find(powers_W == 4e6);
    
    if isempty(target_d) || isempty(target_p)
        error('Target distance or power not found in predefined lists.');
    end
    
    % Re-plot only the selected heatmap
    fig2 = figure('Color','w','Position', [100, 100, 1600, 400]);
    imagesc(dT_log10_matrix(:,:,target_p,target_d), [-8, 2]);
    colormap(jet); colorbar;
    caxis([-8 2]);
    
    % Axes and labels
    set(gca, 'XTick', 1:length(col_labels), 'XTickLabel', col_labels, ...
        'XTickLabelRotation', 90, ...
        'FontSize', 10);
    set(gca, 'YTick', 1:length(all_beam_labels), 'YTickLabel', all_beam_labels);
    title(sprintf('log_{10}(\\DeltaT [K]) at %g pc, %g s, %g W', ...
        L_pc_sel(target_d), duration_s, powers_W(target_p)), 'FontSize', 14);
    
    % Overlay integer values
    for r = 1:size(dT_log10_matrix,1)
        for c = 1:size(dT_log10_matrix,2)
            val = dT_log10_matrix(r,c,target_p,target_d);
            if ~isinf(val)
                text(c, r, sprintf('%d', val), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'k');
            end
        end
    end
    
    
    % Save figure as high-resolution PNG and vector PDF
    print(fig2, 'Figure_2', '-dpng', '-r600');
    print(fig2, 'Figure_2', '-dpdf');
else
    warning(['The current emission duration (duration_s = %.1e s) differs from the canonical value of 1e7 s.\n' ...
        'This may lead to a mismatch with the published version of Figure 2.'], duration_s);
end

%% Figure 4
if duration_s ~= 1e9
    warning(['duration_s = %.1e. This does not match the expected 1e9 s for Figure 4. ' ...
        'Figure_4_Top and Figure_4_Bottom will NOT be saved.'], duration_s);
else
    %% Save top panel (high power)
    fig_top = figure('Color','w','Position', [100, 100, 1600, 400]);
    imagesc(dT_log10_matrix(:,:,find(powers_W==4e10),find(L_pc_sel==75)), [-8, 2]);
    colormap(jet); colorbar; caxis([-8 2]);
    set(gca, 'XTick', 1:length(col_labels), 'XTickLabel', col_labels, ...
        'XTickLabelRotation', 90, 'YTick', 1:length(all_beam_labels), ...
        'YTickLabel', all_beam_labels, 'FontSize', 10);
    title('log_{10}(\DeltaT [K]) at 75 pc, 1e9 s, 4e10 W', 'FontSize', 14);
    
    for r = 1:size(dT_log10_matrix,1)
        for c = 1:size(dT_log10_matrix,2)
            val = dT_log10_matrix(r,c,find(powers_W==4e10),find(L_pc_sel==75));
            if ~isinf(val)
                text(c, r, sprintf('%d', val), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'k');
            end
        end
    end
    print(fig_top, 'Figure_4_Top', '-dpng', '-r600');
    print(fig_top, 'Figure_4_Top', '-dpdf');
    
    %% Save bottom panel (low power)
    fig_bot = figure('Color','w','Position', [100, 100, 1600, 400]);
    imagesc(dT_log10_matrix(:,:,find(powers_W==40000),find(L_pc_sel==75)), [-8, 2]);
    colormap(jet); colorbar; caxis([-8 2]);
    set(gca, 'XTick', 1:length(col_labels), 'XTickLabel', col_labels, ...
        'XTickLabelRotation', 90, 'YTick', 1:length(all_beam_labels), ...
        'YTickLabel', all_beam_labels, 'FontSize', 10);
    title('log_{10}(\DeltaT [K]) at 75 pc, 1e9 s, 40000 W', 'FontSize', 14);
    
    for r = 1:size(dT_log10_matrix,1)
        for c = 1:size(dT_log10_matrix,2)
            val = dT_log10_matrix(r,c,find(powers_W==40000),find(L_pc_sel==75));
            if ~isinf(val)
                text(c, r, sprintf('%d', val), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'k');
            end
        end
    end
    print(fig_bot, 'Figure_4_Bottom', '-dpng', '-r600');
    print(fig_bot, 'Figure_4_Bottom', '-dpdf');
end

%% Figure 5
if duration_s ~= 1e7
    warning(['duration_s = %.1e. This does not match the expected 1e7 s for Figure 5. ' ...
        'Figure_5 will NOT be saved.'], duration_s);
else
    %% Save Figure 4
    fig5 = figure('Color','w','Position', [100, 100, 1600, 400]);
    idx_power = find(powers_W == 4e10);
    idx_dist = find(L_pc_sel == 15767);
    
    imagesc(dT_log10_matrix(:,:,idx_power,idx_dist), [-8, 2]);
    colormap(jet); colorbar; caxis([-8 2]);
    set(gca, 'XTick', 1:length(col_labels), 'XTickLabel', col_labels, ...
        'XTickLabelRotation', 90, 'YTick', 1:length(all_beam_labels), ...
        'YTickLabel', all_beam_labels, 'FontSize', 10);
    title('log_{10}(\DeltaT [K]) at 15767 pc, 1e7 s, 4e10 W', 'FontSize', 14);
    
    % Overlay integer log10(DeltaT) values
    for r = 1:size(dT_log10_matrix,1)
        for c = 1:size(dT_log10_matrix,2)
            val = dT_log10_matrix(r,c,idx_power,idx_dist);
            if ~isinf(val)
                text(c, r, sprintf('%d', val), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'k');
            end
        end
    end
    
    print(fig5, 'Figure_5', '-dpng', '-r600');
    print(fig5, 'Figure_5', '-dpdf');
end
