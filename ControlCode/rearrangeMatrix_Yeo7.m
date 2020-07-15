% function to rearrange a matrix with 234 rows corresponding to Lausanne parcels by Yeo 7 labels

function rearrangedMatrix = rearrangeMatrix_Yeo7(inputMatrix)
    load('LindenYeoPurity/yeo7netlabelsLaus125EJC.mat', 'finalLabels');
    uniqueLabels = unique(finalLabels);
    rearrangedMatrix = [];
    for i = 1:numel(uniqueLabels)
        currentLabel = uniqueLabels(i);
        idx_currentLabel = find(finalLabels==currentLabel);
        rearrangedMatrix = [rearrangedMatrix; inputMatrix(idx_currentLabel, :)];
    end
end