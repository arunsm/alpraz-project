% script to plot surface data in fsaverage5 space; derived from
% fcn_lausannesurf.m 

% INPUT:
% x: 234x1 data vector
% lh_surf_data: surface data as a 10242x1 vector for left hemisphere
% rh_surf_data: surface data as a 10242x1 vector for right hemisphere
% cmap0: user-specified colormap.
% clims, lower/upper limits for color range.

% Outputs:
% f: figure handles (there are two; one for each hemisphere)
% ax: axes handles
% ph: surface handles
% cmap: colormap generated in function

function [f, ax, ph, cmap] = plot_surface_file(x, lh_surf_data, rh_surf_data, cmap0, clims)

% read in vertex and face coordinates in fsaverage5 space

[lh_vertices, lh_faces] = freesurfer_read_surf('/Applications/freesurfer/subjects/fsaverage5/surf/lh.pial');
[rh_vertices, rh_faces] = freesurfer_read_surf('/Applications/freesurfer/subjects/fsaverage5/surf/rh.pial');

if ~exist('cmap0','var')
    cmap0 = jet(256);
end

if ~exist('clims','var')
    clims = [min(x), max(x)];
end

lowval = clims(1) - 0.1*range(clims);
aa = linspace(clims(1),clims(2),size(cmap0,1));
bb = [lowval,clims(1) - 0.01*range(clims), aa];
cmap = interp1(bb',[ones(2,3)*0.25; cmap0],linspace(lowval,clims(2),size(cmap0,1)));

indl = lh_surf_data > 0;
indr = rh_surf_data > 0;

% substituting unknown values (coded as NaN in surface files) with uniform
% color value
lh_surf_data(~indl) = lowval;
rh_surf_data(~indr) = lowval;

f(1) = fcn_rickplot([2,2,4,4]);
ax(1) = axes;
ph(1) = trisurf(lh_faces, lh_vertices(:,1), lh_vertices(:,2), lh_vertices(:,3), lh_surf_data);
axis image
colormap(cmap);

f(2) = fcn_rickplot([2,2,4,4]);
ax(2) = axes;
ph(2) = trisurf(rh_faces, rh_vertices(:,1), rh_vertices(:,2), rh_vertices(:,3), rh_surf_data);
axis image
set(ph,'edgecolor','none');
set(ax,...
    'clim',clims);
colormap(cmap);

function f = fcn_rickplot(pos)
f = figure(...
    'units','inches',...
    'position',pos,...
    'paperpositionmode','auto');

