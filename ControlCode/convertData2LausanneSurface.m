% script to convert an 234x1 vector of Lausanne parcels to fsaverage5 surface
% files

% INPUT: data - 234x1 vector of Lausanne parcel values

% OUTPUTS:
% lh_surf_data: vertices for left hemisphere
% rh_surf_data: vertices for right hemisphere

function [lh_surf_data, rh_surf_data] = convertData2LausanneSurface(data)

mid = 116; % dividing index for 234-node parcellation
%addpath(genpath('/Applications/freesurfer/matlab'));

X = importdata('../../lausanne2008/LausanneParcelNames.xlsx');
LausanneParcelNames = X.textdata;

% read annotation of Lausanne parcels in fsaverage5
[Rv, RL, Rct] = read_annotation('../../Lausanne_surfaceData_fsaverage5/rh.myaparc_125.annot');
[Lv, LL, Lct] = read_annotation('../../Lausanne_surfaceData_fsaverage5/lh.myaparc_125.annot');

%% create annotations containing parcel ID for each vertex

lh_surf_LausanneParcelIDs = zeros(size(Lv, 1), 1);
rh_surf_LausanneParcelIDs = zeros(size(Rv, 1), 1);

for i = 1:size(data, 1)
    currentParcelName = LausanneParcelNames{i};
    
    % right hemisphere
    if i < mid
        currentLabel = Rct.table(ismember(Rct.struct_names, currentParcelName), 5);
        if ~isempty(currentLabel)
            rh_surf_LausanneParcelIDs(RL == currentLabel) = i;
        else
            fprintf('parcel %s not found in rh surface annotation\n', currentParcelName);
        end
    % left hemisphere
    else
        currentLabel = Lct.table(ismember(Lct.struct_names, currentParcelName), 5);
        if ~isempty(currentLabel)
            lh_surf_LausanneParcelIDs(LL == currentLabel) = i;
        else
            fprintf('parcel %s not found in lh surface annotation\n', currentParcelName);
        end
    end
end

%save('../../Lausanne_surfaceData_fsaverage5/lh_surf_LausanneParcelIDs.mat', 'lh_surf_LausanneParcelIDs');
%save('../../Lausanne_surfaceData_fsaverage5/rh_surf_LausanneParcelIDs.mat', 'rh_surf_LausanneParcelIDs');

%% extract data from input vector and create corresponding surface file

rh_surf_data = nan(size(rh_surf_LausanneParcelIDs));
lh_surf_data = nan(size(lh_surf_LausanneParcelIDs));

roir = rh_surf_LausanneParcelIDs > 0;
roil = lh_surf_LausanneParcelIDs > 0;

rh_surf_data(roir) = data(rh_surf_LausanneParcelIDs(roir));
lh_surf_data(roil) = data(lh_surf_LausanneParcelIDs(roil));

end