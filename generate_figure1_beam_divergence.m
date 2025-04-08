%{
Beam Divergence after Gravitational Focusing and in Collimated Transmission

This script computes the post-focal divergence of beams focused by a 
stellar-mass object (e.g., red dwarf star), using the weak-field 
approximation for general relativistic deflection. It compares the 
resulting beam width at various distances with that of purely collimated 
emission with fixed angular divergences.

The script produces Figure 1 of the paper:
"Microthermal Residuals from Directed Electromagnetic Emissions as a New Class of Technosignature" 
by Sergio Bonaque-González (2025), submitted to publication on April 2025.

The figure consists of 8 subpanels, showing the evolution of beam width 
across four spatial regimes (local, regional, long-range, and galaxy-scale),
for both gravitationally-focused and collimated configurations.

Created by Sergio Bonaque-González, PhD
Postdoctoral Fellow (Viera y Clavijo Program)
Department of Applied Physics
University of La Laguna
sbonaque@ull.edu.es

April, 2025
%}

% Beam divergence after passing near a star, with reference divergence lines
clear; clc; close all % Clear environment and prepare workspace

% Physical constants
G = 6.6743e-11;               % [m^3/kg/s^2]
M = 0.35*1.989e30;            % solar mass [kg]
c = 2.99792458e8;             % [m/s]
AU = 1.496e11;                % [m]
R_sun = 6.96e8;               % [m]

%% Simulation parameters:
b_vals = linspace(1, 5, 5) * R_sun; % Impact parameters: multiples of the solar radius
L_pc = linspace(1, 30000, 2000);     % X axis: distance in [pc]
L_m = L_pc * 3.086e16;               % [m]
theta_vals = [1e-7, 1e-8, 1e-9, 1e-10];  % Fixed beam divergences [rad]
baseline_theta = 1e-7; % for comparison purposes [rad]
regions = [153, 307, 1534, 30000];

%% Simulation
% Initialization
W = NaN(length(b_vals), length(L_m));  % NaN to hide pre-focus zone

% CASE 1: Calculates the post-focal beam envelope for rays grazing a lensing star 
% at various impact parameters (b = 1–5 R_sun). Each ray experiences a 
% deflection Theta = 4GM/(c²b), resulting in a focal distance f = b/Theta.
% Beyond the focal region, the beam expands linearly with Theta. The width 
% at each selected distance is computed as W = 2·tan(Theta)·(L - f), expressed in parsecs.
for i = 1:length(b_vals)
    b = b_vals(i);
    theta = 4 * G * M / (c^2 * b);       % angular deflection [rad]
    f = b / theta;                       % focal length [m]

    for j = 1:length(L_m)
        if L_m(j) >= f
            d = L_m(j) - f;
            W(i,j) = 2 * tan(theta) * d; % post-focus divergence
        end
    end
end

% Plot the individual figures
figure('Color','w'); hold on;
colores = lines(length(b_vals)+1);

% Curves for each b
for i = 1:length(b_vals)
    f_au = b_vals(i)^2 * c^2 / (4 * G * M) / AU;
    label = sprintf('b = %.1f R_\\odot, f = %.1f AU', b_vals(i)/R_sun, f_au);
    plot(L_pc, W(i,:) / 3.086e16, 'LineWidth', 2, ...
         'Color', colores(i,:), 'DisplayName', label);
end

ss_width_pc = ones(size(W(i,:)));
ss_width_pc = ss_width_pc*(100 * 1.496e11 / 3.086e16);
plot(L_pc, ss_width_pc, 'LineWidth', 3, 'Color', colores(i+1,:), 'LineStyle', '--','DisplayName', 'Solar system width');

% Baseline (comparative)
W_div = 2 * baseline_theta * L_m; 
W_div_pc = W_div / 3.086e16;
W_div_pc(~isfinite(W_div_pc)) = NaN;
plot(L_pc, W_div_pc, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, ...
     'DisplayName', '\theta = 10^{-7} rad');
set(gca, 'YScale', 'log', ...
         'FontSize', 14, ...
         'FontName', 'Helvetica', ...
         'LineWidth', 1.2, ...
         'Box', 'off', ...
         'TickDir', 'out');
xlabel('Distance along focal line (pc)', 'FontSize', 16);
ylabel('Beam width (pc)', 'FontSize', 16);
title('Beam divergence - Local communication', 'FontSize', 16);

legend('show', 'Location','northwest');
xlim([0 regions(1)]);
grid on;

%2nd figure
figure('Color','w'); hold on;
% Curves for each b
for i = 1:length(b_vals)
    f_au = b_vals(i)^2 * c^2 / (4 * G * M) / AU;
    label = sprintf('b = %.1f R_\\odot, f = %.1f AU', b_vals(i)/R_sun, f_au);
    plot(L_pc, W(i,:) / 3.086e16, 'LineWidth', 2, ...
         'Color', colores(i,:), 'DisplayName', label);
end
plot(L_pc, ss_width_pc, 'LineWidth', 3, 'Color', colores(i+1,:), 'LineStyle', '--','DisplayName', 'Solar system width');

%Baseline
plot(L_pc, W_div_pc, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, ...
     'DisplayName', '\theta = 10^{-7} rad');
set(gca, 'YScale', 'log', ...
         'FontSize', 14, ...
         'FontName', 'Helvetica', ...
         'LineWidth', 1.2, ...
         'Box', 'off', ...
         'TickDir', 'out');

xlabel('Distance along focal line (pc)', 'FontSize', 16);
ylabel('Beam width (pc)', 'FontSize', 16);
title('Beam divergence - Regional communication', 'FontSize', 16);

legend('show', 'Location','northwest');
xlim([regions(1) regions(2)]);
grid on;

%3rd figure
figure('Color','w'); hold on;
% Curves for each b
for i = 1:length(b_vals)
    f_au = b_vals(i)^2 * c^2 / (4 * G * M) / AU;
    label = sprintf('b = %.1f R_\\odot, f = %.1f AU', b_vals(i)/R_sun, f_au);
    plot(L_pc, W(i,:) / 3.086e16, 'LineWidth', 2, ...
         'Color', colores(i,:), 'DisplayName', label);
end
plot(L_pc, ss_width_pc, 'LineWidth', 3, 'Color', colores(i+1,:), 'LineStyle', '--','DisplayName', 'Solar system width');

% Baseline
plot(L_pc, W_div_pc, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, ...
     'DisplayName', '\theta = 10^{-7} rad');

set(gca, 'YScale', 'log', ...
         'FontSize', 14, ...
         'FontName', 'Helvetica', ...
         'LineWidth', 1.2, ...
         'Box', 'off', ...
         'TickDir', 'out');

xlabel('Distance along focal line (pc)', 'FontSize', 16);
ylabel('Beam width (pc)', 'FontSize', 16);
title('Beam divergence - Long-range communication', 'FontSize', 16);

legend('show', 'Location','northwest');
xlim([regions(2) regions(3)]);
grid on;

%4th figure
figure('Color','w'); hold on;
% Curves for each b
for i = 1:length(b_vals)
    f_au = b_vals(i)^2 * c^2 / (4 * G * M) / AU;
    label = sprintf('b = %.1f R_\\odot, f = %.1f AU', b_vals(i)/R_sun, f_au);
    plot(L_pc, W(i,:) / 3.086e16, 'LineWidth', 2, ...
         'Color', colores(i,:), 'DisplayName', label);
end
plot(L_pc, ss_width_pc, 'LineWidth', 3, 'Color', colores(i+1,:), 'LineStyle', '--','DisplayName', 'Solar system width');

% Baseline
plot(L_pc, W_div_pc, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, ...
     'DisplayName', '\theta = 10^{-7} rad');
set(gca, 'YScale', 'log', ...
         'FontSize', 14, ...
         'FontName', 'Helvetica', ...
         'LineWidth', 1.2, ...
         'Box', 'off', ...
         'TickDir', 'out');

xlabel('Distance along focal line (pc)', 'FontSize', 16);
ylabel('Beam width (pc)', 'FontSize', 16);
title('Beam divergence - Galaxy-scale communication', 'FontSize', 16);

legend('show', 'Location','northwest');
xlim([regions(3)  regions(4)]);
grid on;

%% CASE 2: Ideal collimated beams with fixed divergence (no gravitational lensing at destination)
colores_col = lines(length(theta_vals)+1);
%1st figure

figure('Color','w'); hold on;
for i = 1:length(theta_vals)
    W_col = 2 * theta_vals(i) * L_m;
    W_col_pc = W_col / 3.086e16;
    W_col_pc(~isfinite(W_col_pc)) = NaN;
    W_col_pc = max(W_col_pc, eps);  % para evitar ceros en escala log
    label = sprintf('\\theta = 10^{%d} rad', round(log10(theta_vals(i))));
    plot(L_pc, W_col_pc, 'LineWidth', 2, 'Color', colores_col(i,:), 'DisplayName', label);
end

plot(L_pc, ss_width_pc, 'LineWidth', 3, 'Color', colores_col(i+1,:), 'LineStyle', '--','DisplayName', 'Solar system width');
set(gca, 'YScale', 'log', ...
         'FontSize', 14, ...
         'FontName', 'Helvetica', ...
         'LineWidth', 1.2, ...
         'Box', 'off', ...
         'TickDir', 'out');

xlabel('Distance (pc)', 'FontSize', 16);
ylabel('Beam width (pc)', 'FontSize', 16);
title('Beam divergence - Local communication', 'FontSize', 16);

legend('show', 'Location','northwest');
xlim([0 regions(1)]);
grid on;

% 2nd figure
figure('Color','w'); hold on;
for i = 1:length(theta_vals)
    W_col = 2 * theta_vals(i) * L_m;
    W_col_pc = W_col / 3.086e16;
    W_col_pc(~isfinite(W_col_pc)) = NaN;
    W_col_pc = max(W_col_pc, eps);  % para evitar ceros en escala log
    label = sprintf('\\theta = 10^{%d} rad', round(log10(theta_vals(i))));
    plot(L_pc, W_col_pc, 'LineWidth', 2, 'Color', colores_col(i,:), 'DisplayName', label);
end

plot(L_pc, ss_width_pc, 'LineWidth', 3, 'Color', colores_col(i+1,:), 'LineStyle', '--','DisplayName', 'Solar system width');
set(gca, 'YScale', 'log', ...
         'FontSize', 14, ...
         'FontName', 'Helvetica', ...
         'LineWidth', 1.2, ...
         'Box', 'off', ...
         'TickDir', 'out');

xlabel('Distance (pc)', 'FontSize', 16);
ylabel('Beam width (pc)', 'FontSize', 16);
title('Beam divergence - Regional communication', 'FontSize', 16);

legend('show', 'Location','northwest');
xlim([regions(1) regions(2)]);
grid on;

% 3rd figure
figure('Color','w'); hold on;
for i = 1:length(theta_vals)
    W_col = 2 * theta_vals(i) * L_m;
    W_col_pc = W_col / 3.086e16;
    W_col_pc(~isfinite(W_col_pc)) = NaN;
    W_col_pc = max(W_col_pc, eps);  % para evitar ceros en escala log
    label = sprintf('\\theta = 10^{%d} rad', round(log10(theta_vals(i))));
    plot(L_pc, W_col_pc, 'LineWidth', 2, 'Color', colores_col(i,:), 'DisplayName', label);
end

plot(L_pc, ss_width_pc, 'LineWidth', 3, 'Color', colores_col(i+1,:), 'LineStyle', '--','DisplayName', 'Solar system width');
set(gca, 'YScale', 'log', ...
         'FontSize', 14, ...
         'FontName', 'Helvetica', ...
         'LineWidth', 1.2, ...
         'Box', 'off', ...
         'TickDir', 'out');

xlabel('Distance (pc)', 'FontSize', 16);
ylabel('Beam width (pc)', 'FontSize', 16);
title('Beam divergence - Long-range communication', 'FontSize', 16);

legend('show', 'Location','northwest');
xlim([regions(2) regions(3)]);
grid on;

%4th figure
figure('Color','w'); hold on;
for i = 1:length(theta_vals)
    W_col = 2 * theta_vals(i) * L_m;
    W_col_pc = W_col / 3.086e16;
    W_col_pc(~isfinite(W_col_pc)) = NaN;
    W_col_pc = max(W_col_pc, eps);  % para evitar ceros en escala log
    label = sprintf('\\theta = 10^{%d} rad', round(log10(theta_vals(i))));
    plot(L_pc, W_col_pc, 'LineWidth', 2, 'Color', colores_col(i,:), 'DisplayName', label);
end

plot(L_pc, ss_width_pc, 'LineWidth', 3, 'Color', colores_col(i+1,:), 'LineStyle', '--','DisplayName', 'Solar system width');
set(gca, 'YScale', 'log', ...
         'FontSize', 14, ...
         'FontName', 'Helvetica', ...
         'LineWidth', 1.2, ...
         'Box', 'off', ...
         'TickDir', 'out');

xlabel('Distance (pc)', 'FontSize', 16);
ylabel('Beam width (pc)', 'FontSize', 16);
title('Beam divergence - Galaxy-scale communication', 'FontSize', 16);

legend('show', 'Location','northwest');
xlim([regions(3) regions(4)]);
grid on;


%% Figure as it appears in the paper
fig = figure('Color','w','Position',[100 100 1800 1000]);
titles = {
    'Post-focus divergence - Local'
    'Post-focus divergence - Regional'
    'Post-focus divergence - Long-range'
    'Post-focus divergence - Galaxy-scale'
    'Collimated beam - Local'
    'Collimated beam - Regional'
    'Collimated beam - Long-range'
    'Collimated beam - Galaxy-scale'
};

xlims = [0 regions(1); regions(1) regions(2); regions(2) regions(3); regions(3) regions(4)];

% Plot the 4 post-focus subgraphs
for idx = 1:4
    subplot(2,4,idx); hold on;
    for i = 1:length(b_vals)
        plot(L_pc, W(i,:) / 3.086e16, 'LineWidth', 2, 'Color', colores(i,:));
    end
    plot(L_pc, ss_width_pc, 'LineWidth', 2.5, 'Color', colores(i+1,:), 'LineStyle', '--');
    plot(L_pc, W_div_pc, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

    xlim(xlims(idx,:));
    set(gca, 'YScale', 'log', 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
    xlabel('Distance (pc)');
    ylabel('Beam width (pc)');
    title(titles{idx}, 'FontSize', 12);
    grid on;
end

% Plot the 4 collimated subgraphs
for idx = 1:4
    subplot(2,4,idx+4); hold on;
    for i = 1:length(theta_vals)
        W_col = 2 * theta_vals(i) * L_m;
        W_col_pc = W_col / 3.086e16;
        W_col_pc = max(W_col_pc, eps);
        plot(L_pc, W_col_pc, 'LineWidth', 2, 'Color', colores_col(i,:));
    end
    plot(L_pc, ss_width_pc, 'LineWidth', 2.5, 'Color', colores_col(end,:), 'LineStyle', '--');

    xlim(xlims(idx,:));
    set(gca, 'YScale', 'log', 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
    xlabel('Distance (pc)');
    ylabel('Beam width (pc)');
    title(titles{idx+4}, 'FontSize', 12);
    grid on;
end

% -------------------------
% Add common legend for the first 4 subplots (top)
subplot(2,4,2); 
h1 = findall(gca, 'Type', 'Line');
legend([h1(end:-1:1)], {'b = 1 R_\odot','b = 2 R_\odot','b = 3 R_\odot','b = 4 R_\odot','b = 5 R_\odot','Solar system width','\theta = 10^{-7} rad'}, ...
       'Position',[0.15 0.94 0.7 0.05], 'Orientation','horizontal', 'Box','off');

% -------------------------
% Common legend for the 4 lower subplots (collimated)
subplot(2,4,6); 
h2 = findall(gca, 'Type', 'Line');
legend([h2(end:-1:1)], {'\theta = 10^{-7} rad','\theta = 10^{-8} rad','\theta = 10^{-9} rad','\theta = 10^{-10} rad','Solar system width'}, ...
       'Position',[0.15 0.00 0.7 0.05], 'Orientation','horizontal', 'Box','off');

%% Save as high-resolution PNG
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 16 8]);  % Size in inches
print(fig, 'Figure_1', '-dpng', '-r600');  % PNG  600 dpi

%% Save as vectorial
print(fig, 'Figure_1', '-dpdf', '-painters');  % Alta calidad para impresión
