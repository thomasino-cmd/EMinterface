%% meniscus_bottom_comparison.m
% Thesis-ready comparison:
%   (A) 3D concave meniscus triangulated surface with power heatmap
%   (B) 2D bottom-plane triangulated footprints with delivered irradiance heatmap
%
% Required CSV (from tile_overlap.cpp):
%   - triangle_mesh_meniscus.csv
%   - triangle_mesh_bottom.csv

clear; clc; close all;

%% ===================== CONFIG =====================
cfg.dataDir = 'D:\cose random\università\3anno\CAMPI ELETTROMAGNETICI\TESI\outputV2\testA\150';
%cfg.dataDir = 'C:\Users\thoma\source\repos\EMinterface\EMinterface';
cfg.fileMeniscus = 'triangle_mesh_meniscus.csv';
cfg.fileBottom   = 'triangle_mesh_bottom.csv';

cfg.powerModeMeniscus = 'Sdotn'; % 'Sdotn' [W/m^2] OR 'PinDensity' [W/m^2]
cfg.unitLabel = 'mW/cm^2';
cfg.toDisplay = 0.1;             % W/m^2 -> mW/cm^2

cfg.fontSize = 12;
cfg.edgeAlpha3D = 0.12;
cfg.edgeAlpha2D = 0.08;
cfg.edgeColor = [0.2 0.2 0.2];
cfg.cmap = turbo(256);

% To emphasize small differences in weak-loss media
cfg.cLimPercentile = [2 98];
cfg.sameCLim = true;

% Optional: keep only positive delivered power triangles on bottom panel
cfg.filterBottomPositive = false;

% --- AUTO EXPORT SETTINGS ---
cfg.exportFigures = true;      % <--- metti false se non vuoi salvare
cfg.exportDir = cfg.dataDir;   % <--- automatico: stessa cartella dei CSV
cfg.exportDPI = 350;

% nomi file (senza path): verranno salvati in cfg.exportDir
cfg.exportMeniscusName = 'meniscus_3D_power.png';
cfg.exportBottomName   = 'bottom_2D_irradiance.png';

% Global-loss reference metric:
%  - 'meniscus_to_bottom': loss relative to power that reached meniscus
%  - 'source_to_bottom'  : loss relative to total source power emitted
cfg.globalLossMode = 'source_to_bottom';

% Source power model (used when globalLossMode='source_to_bottom').
% If sourcePowerOverride_W is provided, it has priority.
cfg.sourcePowerOverride_W = [];

% Default source model coherent with tile_overlap.cpp demo parameters
% (uniform plane-wave beam over disk of radius Rf).
cfg.sourceBeamRadius_m = 5.64e-3;      % Rf in tile_overlap.cpp
cfg.sourceEamp_Vpm = 2744.9;           % |A_inc| in tile_overlap.cpp
%% ==================================================

% ---------- Load ----------
pathMeniscus = fullfile(cfg.dataDir, cfg.fileMeniscus);
pathBottom   = fullfile(cfg.dataDir, cfg.fileBottom);

if ~isfile(pathMeniscus)
    error('Unable to find meniscus CSV: %s', pathMeniscus);
end
if ~isfile(pathBottom)
    error('Unable to find bottom CSV: %s', pathBottom);
end

Tm = readtable(pathMeniscus);
Tb = readtable(pathBottom);

reqM = {'tri_id','v0x_m','v0y_m','v0z_m','v1x_m','v1y_m','v1z_m','v2x_m','v2y_m','v2z_m','area_face_m2','Sdotn_Wm2','Pin_W','R','T','A'};
reqB = {'tri_id','valid','p0x_m','p0y_m','p1x_m','p1y_m','p2x_m','p2y_m','foot_area_m2','att','Ponplane_W','I_Wm2','status_code','status_label'};
check_columns(Tm, reqM, pathMeniscus);
check_columns(Tb, reqB, pathBottom);

if isempty(Tm)
    error('Meniscus table is empty.');
end
if isempty(Tb)
    error('Bottom table is empty.');
end

% Diagnostics: why triangles are invalid/excluded in CPP
if any(~Tb.valid)
    fprintf('\nBottom exclusion diagnostics (from CSV):\n');
    [g,labels] = findgroups(string(Tb.status_label));
    counts = splitapply(@numel, Tb.status_label, g);
    for ii = 1:numel(labels)
        fprintf('  %-28s : %d\n', labels(ii), counts(ii));
    end
end

TbAll   = Tb;                               % output completo simulazione
TbValid = TbAll(TbAll.valid==1, :);         % triangoli accettati dal C++
TbPlot  = TbValid;                          % ciò che visualizzi

if cfg.filterBottomPositive
    TbPlot = TbPlot(TbPlot.Ponplane_W > 0 & isfinite(TbPlot.Ponplane_W), :);
end

if isempty(TbPlot)
    error('No valid bottom triangles to plot (valid==1).');
end

%% ---------- Meniscus triangulation (ALL triangles, no exclusion) ----------
[Vm, Fm] = build_tri_mesh_3d(Tm.v0x_m, Tm.v0y_m, Tm.v0z_m, ...
                             Tm.v1x_m, Tm.v1y_m, Tm.v1z_m, ...
                             Tm.v2x_m, Tm.v2y_m, Tm.v2z_m);

switch lower(cfg.powerModeMeniscus)
    case 'sdotn'
        Ctri_m = Tm.Sdotn_Wm2;
        mLabel = 'Meniscus interface flux S\cdot n';
    case 'pindensity'
        Ctri_m = Tm.Pin_W ./ max(Tm.area_face_m2, eps);
        mLabel = 'Meniscus Pin / area';
    otherwise
        error('Unknown cfg.powerModeMeniscus: %s', cfg.powerModeMeniscus);
end
Ctri_m_disp = Ctri_m * cfg.toDisplay;

%% ---------- Bottom triangulation (all available bottom triangles) ----------
[Vb, Fb] = build_tri_mesh_2d(TbPlot.p0x_m, TbPlot.p0y_m, TbPlot.p1x_m, TbPlot.p1y_m, TbPlot.p2x_m, TbPlot.p2y_m);
Ctri_b = TbPlot.I_Wm2;
Ctri_b_disp = Ctri_b * cfg.toDisplay;
bLabel = 'Bottom delivered irradiance';

%% ---------- Color limits (robust) ----------
if cfg.sameCLim
    cAll = [Ctri_m_disp; Ctri_b_disp];
    cAll = cAll(isfinite(cAll));
    clim = prctile(cAll, cfg.cLimPercentile);
else
    climM = prctile(Ctri_m_disp(isfinite(Ctri_m_disp)), cfg.cLimPercentile);
    climB = prctile(Ctri_b_disp(isfinite(Ctri_b_disp)), cfg.cLimPercentile);
end

%% ---------- (Optional) on-screen combined figure ----------
fig = figure('Color','w','Position',[80 80 1300 580]);
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% --- Panel A: 3D meniscus ---
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');
trisurf(Fm, Vm(:,1)*1e3, Vm(:,2)*1e3, Vm(:,3)*1e3, ...
    'FaceColor','flat', ...
    'FaceVertexCData', repelem(Ctri_m_disp, 3), ...
    'EdgeColor', cfg.edgeColor, 'EdgeAlpha', cfg.edgeAlpha3D, 'LineWidth', 0.15);

view(ax1, 36, 22);
axis(ax1,'equal'); grid(ax1,'on');
% xlabel(ax1, 'x (mm)'); ylabel(ax1, 'y (mm)'); zlabel(ax1, 'z (mm)');
colormap(ax1, cfg.cmap);
cb1 = colorbar(ax1); cb1.Label.String = sprintf('%s (%s)', mLabel, cfg.unitLabel);
if cfg.sameCLim, caxis(ax1, clim); else, caxis(ax1, climM); end
set(ax1,'FontSize',cfg.fontSize,'LineWidth',1);

% --- Panel B: 2D bottom plane ---
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');
Vbz = [Vb(:,1)*1e3, Vb(:,2)*1e3, zeros(size(Vb,1),1)];
trisurf(Fb, Vbz(:,1), Vbz(:,2), Vbz(:,3), ...
    'FaceColor','flat', ...
    'FaceVertexCData', repelem(Ctri_b_disp, 3), ...
    'EdgeColor', cfg.edgeColor, 'EdgeAlpha', cfg.edgeAlpha2D, 'LineWidth', 0.1);

view(ax2, 2);
axis(ax2,'equal'); grid(ax2,'on');
% xlabel(ax2, 'x (mm)'); ylabel(ax2, 'y (mm)');
colormap(ax2, cfg.cmap);
cb2 = colorbar(ax2); cb2.Label.String = sprintf('%s (%s)', bLabel, cfg.unitLabel);
if cfg.sameCLim, caxis(ax2, clim); else, caxis(ax2, climB); end
set(ax2,'FontSize',cfg.fontSize,'LineWidth',1);

%% ---------- Export TWO separate figures automatically ----------
if cfg.exportFigures
    if ~exist(cfg.exportDir, 'dir')
        error('Export directory does not exist: %s', cfg.exportDir);
    end

    % Figure 1: Meniscus only
    figM = figure('Color','w','Position',[120 120 760 640]);
    axM = axes(figM); hold(axM,'on'); box(axM,'on');
    trisurf(Fm, Vm(:,1)*1e3, Vm(:,2)*1e3, Vm(:,3)*1e3, ...
        'FaceColor','flat', ...
        'FaceVertexCData', repelem(Ctri_m_disp, 3), ...
        'EdgeColor', cfg.edgeColor, 'EdgeAlpha', cfg.edgeAlpha3D, 'LineWidth', 0.15);
    view(axM, 36, 22);
    axis(axM,'equal'); grid(axM,'on');
    % xlabel(axM, 'x (mm)'); ylabel(axM, 'y (mm)'); zlabel(axM, 'z (mm)');
    colormap(axM, cfg.cmap);
    cbM = colorbar(axM); cbM.Label.String = sprintf('%s (%s)', mLabel, cfg.unitLabel);
    if cfg.sameCLim, caxis(axM, clim); else, caxis(axM, climM); end
    set(axM,'FontSize',cfg.fontSize,'LineWidth',1);

    outMeniscus = fullfile(cfg.exportDir, cfg.exportMeniscusName);
    exportgraphics(figM, outMeniscus, 'Resolution', cfg.exportDPI);

    % Figure 2: Bottom only
    figB = figure('Color','w','Position',[920 120 760 640]);
    axB = axes(figB); hold(axB,'on'); box(axB,'on');
    trisurf(Fb, Vbz(:,1), Vbz(:,2), Vbz(:,3), ...
        'FaceColor','flat', ...
        'FaceVertexCData', repelem(Ctri_b_disp, 3), ...
        'EdgeColor', cfg.edgeColor, 'EdgeAlpha', cfg.edgeAlpha2D, 'LineWidth', 0.1);
    view(axB, 2);
    axis(axB,'equal'); grid(axB,'on');
    % xlabel(axB, 'x (mm)'); ylabel(axB, 'y (mm)');
    colormap(axB, cfg.cmap);
    cbB = colorbar(axB); cbB.Label.String = sprintf('%s (%s)', bLabel, cfg.unitLabel);
    if cfg.sameCLim, caxis(axB, clim); else, caxis(axB, climB); end
    set(axB,'FontSize',cfg.fontSize,'LineWidth',1);

    outBottom = fullfile(cfg.exportDir, cfg.exportBottomName);
    exportgraphics(figB, outBottom, 'Resolution', cfg.exportDPI);

    fprintf('\nSaved figures:\n  %s\n  %s\n', outMeniscus, outBottom);
end

%% ---------- Global summary from full meniscus + full bottom tables ----------
PinTot  = sum(Tm.Pin_W, 'omitnan');
PoutTot = sum(TbAll.Ponplane_W, 'omitnan');

% Legacy metric (kept for diagnostic comparison): meniscus -> bottom
etaMeniscus = PoutTot / max(PinTot, eps);
lossMeniscus = 1 - etaMeniscus;

% Requested metric: source -> bottom
PsourceTot = compute_source_power(cfg);
if ~isfinite(PsourceTot) || PsourceTot <= 0
    warning('Invalid source power estimate (%.6g W). Falling back to meniscus-based metric.', PsourceTot);
    PsourceTot = PinTot;
end
etaSource = PoutTot / max(PsourceTot, eps);
lossSource = 1 - etaSource;

switch lower(cfg.globalLossMode)
    case 'source_to_bottom'
        etaShown = etaSource;
        lossShown = lossSource;
        refLabel = 'source\rightarrowbottom';
    case 'meniscus_to_bottom'
        etaShown = etaMeniscus;
        lossShown = lossMeniscus;
        refLabel = 'meniscus\rightarrowbottom';
    otherwise
        error('Unknown cfg.globalLossMode: %s', cfg.globalLossMode);
end

fprintf('\n=== Thesis Comparison Summary ===\n');
fprintf('Meniscus triangles plotted: %d\n', height(Tm));
fprintf('Bottom triangles plotted:   %d (valid) / %d (all)\n', height(TbPlot), height(TbAll));
fprintf('Psource,total   = %.6g W\n', PsourceTot);
fprintf('Pin,total       = %.6g W\n', PinTot);
fprintf('Ponplane,total  = %.6g W\n', PoutTot);
fprintf('eta_source      = %.6f\n', etaSource);
fprintf('loss_source     = %.3f %%\n', 100*lossSource);
fprintf('eta_meniscus    = %.6f\n', etaMeniscus);
fprintf('loss_meniscus   = %.3f %%\n', 100*lossMeniscus);

%% ---------------- local functions ----------------
function Psource = compute_source_power(cfg)
if ~isempty(cfg.sourcePowerOverride_W)
    Psource = cfg.sourcePowerOverride_W;
    return;
end

% Plane-wave estimate in medium 1 (air):
% I = |E|^2 / (2*eta0),   P = I * pi*Rf^2
EPS0 = 8.854187817e-12;
MU0 = 1.2566370614359173e-6;
eta0 = sqrt(MU0 / EPS0);

I0 = (abs(cfg.sourceEamp_Vpm)^2) / (2*eta0);
Psource = I0 * pi * cfg.sourceBeamRadius_m^2;
end

function check_columns(T, requiredCols, fileName)
missing = setdiff(requiredCols, T.Properties.VariableNames);
if ~isempty(missing)
    error('File %s is missing required columns: %s', fileName, strjoin(missing, ', '));
end
end

function [V,F] = build_tri_mesh_3d(v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z)
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
