% script to convert fsaverage5 surface files to a 234x1 vector of Lausanne parcels

% INPUTS:
% lh_surf_data: vertices for left hemisphere
% rh_surf_data: vertices for right hemisphere

% OUTPUT: data - 234x1 vector of Lausanne parcel values

function data = convertLausanneSurface2Data(lh_surf_data, rh_surf_data)

mid = 116; % dividing index for 234-node parcellation
%addpath(genpath('/Applications/freesurfer/matlab'));

% load fsaverage5 files containing parcel IDs for each vertex
load('../../data/Lausanne_surfaceData_fsaverage5/lh_surf_LausanneParcelIDs.mat', 'lh_surf_LausanneParcelIDs');
load('../../data/Lausanne_surfaceData_fsaverage5/rh_surf_LausanneParcelIDs.mat', 'rh_surf_LausanneParcelIDs');

data = nan(234, 1);
for i = 1:234
    % right hemisphere
    if i < mid
        data(i) = mean(rh_surf_data(rh_surf_LausanneParcelIDs == i));
    % left hemisphere
    else
        data(i) = mean(lh_surf_data(lh_surf_LausanneParcelIDs == i));
    end
end