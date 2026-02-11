%% meniscus_bottom_comparison.m
% Thesis-ready comparison figure:
%   (A) 3D concave meniscus triangulated surface colored by interface power flux
%   (B) 2D bottom-plane triangulated footprints colored by delivered irradiance
%
% Required CSV (from tile_overlap.cpp):
%   - triangle_mesh_meniscus.csv
%   - triangle_mesh_bottom.csv
%
% Optional legacy CSV (not used here): field_map_cpp.csv, triangle_power_*.csv

clear; clc; close all;

%% ===================== CONFIG =====================
cfg.fileMeniscus = 'triangle_mesh_meniscus.csv';
cfg.fileBottom   = 'triangle_mesh_bottom.csv';

cfg.powerModeMeniscus = 'Sdotn'; % 'Sdotn' [W/m^2] OR 'PinDensity' [W/m^2]
cfg.unitLabel = 'mW/cm^2';       % display unit for colorbars
cfg.toDisplay = 0.1;             % W/m^2 -> mW/cm^2

cfg.fontSize = 12;
cfg.edgeAlpha3D = 0.12;
cfg.edgeAlpha2D = 0.08;
cfg.edgeColor = [0.2 0.2 0.2];
cfg.cmap = turbo(256);

% To emphasize small differences in weak-loss media, use robust percentile limits
cfg.cLimPercentile = [2 98];   % [low high] percentiles
cfg.sameCLim = true;           % true => same color limits for both panels

cfg.exportFigure = false;
cfg.exportName = 'thesis_meniscus_vs_bottom_power.png';
cfg.exportDPI = 350;
%% ==================================================

% ---------- Load ----------
Tm = readtable(cfg.fileMeniscus);
Tb = readtable(cfg.fileBottom);

reqM = {'tri_id','v0x_m','v0y_m','v0z_m','v1x_m','v1y_m','v1z_m','v2x_m','v2y_m','v2z_m','area_face_m2','Sdotn_Wm2','Pin_W','R','T','A'};
reqB = {'tri_id','p0x_m','p0y_m','p1x_m','p1y_m','p2x_m','p2y_m','foot_area_m2','att','Ponplane_W','I_Wm2'};
check_columns(Tm, reqM, cfg.fileMeniscus);
check_columns(Tb, reqB, cfg.fileBottom);

% Keep only triangles with physical transmission to bottom in Tb
Tb = Tb(Tb.Ponplane_W > 0 & isfinite(Tb.Ponplane_W), :);
if isempty(Tb)
    error('No bottom triangles with Ponplane_W > 0 found.');
end

% Join to map R/T/A and other meniscus values if needed
T = innerjoin(Tm, Tb, 'Keys', 'tri_id');
if isempty(T)
    error('No common tri_id between meniscus and bottom mesh CSV files.');
end

%% ---------- Build triangulations ----------
% Meniscus mesh (3D): 3 unique vertices per triangle -> explicit triangulation
[Vm, Fm] = build_tri_mesh_3d(T.v0x_m, T.v0y_m, T.v0z_m, T.v1x_m, T.v1y_m, T.v1z_m, T.v2x_m, T.v2y_m, T.v2z_m);

% Bottom mesh (2D in xy): triangulated footprints
[Vb, Fb] = build_tri_mesh_2d(T.p0x_m, T.p0y_m, T.p1x_m, T.p1y_m, T.p2x_m, T.p2y_m);

%% ---------- Triangle scalar fields ----------
% Meniscus scalar to show on 3D surface
switch lower(cfg.powerModeMeniscus)
    case 'sdotn'
        Ctri_m = T.Sdotn_Wm2;                  % interface power flux density [W/m^2]
        mLabel = 'Meniscus interface flux S\cdot n';
    case 'pindensity'
        Ctri_m = T.Pin_W ./ max(T.area_face_m2, eps); % [W/m^2]
        mLabel = 'Meniscus Pin / area';
    otherwise
        error('Unknown cfg.powerModeMeniscus: %s', cfg.powerModeMeniscus);
end

% Bottom scalar (delivered irradiance)
Ctri_b = T.I_Wm2;                              % [W/m^2]
bLabel = 'Bottom delivered irradiance';

% Convert display units (mW/cm^2)
Ctri_m_disp = Ctri_m * cfg.toDisplay;
Ctri_b_disp = Ctri_b * cfg.toDisplay;

% Color limits (robust)
if cfg.sameCLim
    cAll = [Ctri_m_disp; Ctri_b_disp];
    clim = prctile(cAll, cfg.cLimPercentile);
else
    climM = prctile(Ctri_m_disp, cfg.cLimPercentile);
    climB = prctile(Ctri_b_disp, cfg.cLimPercentile);
end

%% ---------- Figure: only the requested comparison ----------
fig = figure('Color','w','Position',[80 80 1300 580]);
tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact'); %#ok<NASGU>

% --- Panel A: 3D meniscus, triangulated, colored by interface power ---
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');

% Face colors per triangle -> repeat for each face
trisurf(Fm, Vm(:,1)*1e3, Vm(:,2)*1e3, Vm(:,3)*1e3, ...
    'FaceColor','flat', ...
    'FaceVertexCData', repelem(Ctri_m_disp, 3), ...
    'EdgeColor', cfg.edgeColor, 'EdgeAlpha', cfg.edgeAlpha3D, 'LineWidth', 0.15);

view(ax1, 36, 22);
axis(ax1,'equal'); grid(ax1,'on');
xlabel(ax1, 'x (mm)'); ylabel(ax1, 'y (mm)'); zlabel(ax1, 'z (mm)');
title(ax1, sprintf('Meniscus 3D triangulated (%s)', mLabel));
colormap(ax1, cfg.cmap);
cb1 = colorbar(ax1); cb1.Label.String = sprintf('%s (%s)', mLabel, cfg.unitLabel);
if cfg.sameCLim, caxis(ax1, clim); else, caxis(ax1, climM); end
set(ax1,'FontSize',cfg.fontSize,'LineWidth',1);

% --- Panel B: Bottom plane triangulated footprints, colored by delivered I ---
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');

% 2D triangulation (z=0 just for patch rendering)
Vbz = [Vb(:,1)*1e3, Vb(:,2)*1e3, zeros(size(Vb,1),1)];
trisurf(Fb, Vbz(:,1), Vbz(:,2), Vbz(:,3), ...
    'FaceColor','flat', ...
    'FaceVertexCData', repelem(Ctri_b_disp, 3), ...
    'EdgeColor', cfg.edgeColor, 'EdgeAlpha', cfg.edgeAlpha2D, 'LineWidth', 0.1);

view(ax2, 2);
axis(ax2,'equal'); grid(ax2,'on');
xlabel(ax2, 'x (mm)'); ylabel(ax2, 'y (mm)');
title(ax2, sprintf('Bottom plane triangulated (%s)', bLabel));
colormap(ax2, cfg.cmap);
cb2 = colorbar(ax2); cb2.Label.String = sprintf('%s (%s)', bLabel, cfg.unitLabel);
if cfg.sameCLim, caxis(ax2, clim); else, caxis(ax2, climB); end
set(ax2,'FontSize',cfg.fontSize,'LineWidth',1);

% Global summary annotation
PinTot = nansum(T.Pin_W);
PoutTot = nansum(T.Ponplane_W);
etaG = PoutTot / max(PinTot, eps);
lossG = 1 - etaG;

sgtitle(sprintf('Power comparison meniscus \rightarrow bottom   |   \eta_{global}=%.4f   (loss=%.2f%%)', ...
    etaG, 100*lossG), 'FontWeight','bold');

if cfg.exportFigure
    exportgraphics(fig, cfg.exportName, 'Resolution', cfg.exportDPI);
end

fprintf('\n=== Thesis Comparison Summary ===\n');
fprintf('Triangles in comparison: %d\n', height(T));
fprintf('Pin,total       = %.6g W\n', PinTot);
fprintf('Ponplane,total  = %.6g W\n', PoutTot);
fprintf('eta_global      = %.6f\n', etaG);
fprintf('global loss     = %.3f %%\n', 100*lossG);

%% ---------------- local functions ----------------
function check_columns(T, requiredCols, fileName)
missing = setdiff(requiredCols, T.Properties.VariableNames);
if ~isempty(missing)
    error('File %s is missing required columns: %s', fileName, strjoin(missing, ', '));
end
end

function [V,F] = build_tri_mesh_3d(v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z)
% Build explicit triangulation with no vertex sharing (robust/simple)
n = numel(v0x);
V = zeros(3*n,3);
F = zeros(n,3);
for i = 1:n
    b = 3*(i-1);
    V(b+1,:) = [v0x(i), v0y(i), v0z(i)];
    V(b+2,:) = [v1x(i), v1y(i), v1z(i)];
    V(b+3,:) = [v2x(i), v2y(i), v2z(i)];
    F(i,:) = [b+1, b+2, b+3];
end
end

function [V,F] = build_tri_mesh_2d(p0x,p0y,p1x,p1y,p2x,p2y)
n = numel(p0x);
V = zeros(3*n,2);
F = zeros(n,3);
for i = 1:n
    b = 3*(i-1);
    V(b+1,:) = [p0x(i), p0y(i)];
    V(b+2,:) = [p1x(i), p1y(i)];
    V(b+3,:) = [p2x(i), p2y(i)];
    F(i,:) = [b+1, b+2, b+3];
end
end
