% function to rearrange a matrix with 234 rows (w/brain stem) or 233 (w/o brain stem)
% corresponding to Lausanne parcels by Yeo 7 labels

% expected input: (234 x N) or (233 x N) matrix

function rearrangedMatrix = rearrangeMatrix_Yeo7(inputMatrix)
    load('LindenYeoPurity/yeo7netlabelsLaus125EJC.mat', 'finalLabels');
    uniqueLabels = unique(finalLabels);
    idx_brainstem = 234; % index of brain stem parcel    

    % truncating Lausanne parcel labels if input matrix does not contain
    % brain stem info
    if size(inputMatrix, 1) == 233
        finalLabels(idx_brainstem) = [];
    end
    
    rearrangedMatrix = [];
    for i = 1:numel(uniqueLabels)
        currentLabel = uniqueLabels(i);
        idx_currentLabel = find(finalLabels==currentLabel);
        rearrangedMatrix = [rearrangedMatrix; inputMatrix(idx_currentLabel, :)];
    end
end