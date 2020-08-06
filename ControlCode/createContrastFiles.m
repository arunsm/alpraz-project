% script to create separate tables for each contrast (threat, non-threat,
% neutral)

function [] = createContrastFiles(baseName_emotionid, baseName_emotionrec)

contrastNames = {'contrast1_threatcorrectStd', 'contrast3_nonthreatcorrectStd', 'contrast5_neutralcorrectStd'};

readPath = strcat(baseName_emotionid, '.csv');
allControlEnergies_emotionid = readtable(readPath);

for i = 1:numel(contrastNames)
    currentContrast = contrastNames{i};
    currentContrast_idx = strcmp(allControlEnergies_emotionid.contrast, currentContrast);
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(currentContrast_idx, :);
    writetable(allControlEnergies_emotionid_currentContrast, strcat(baseName_emotionid, '_', currentContrast, '.csv'));
end

readPath = strcat(baseName_emotionrec, '.csv');
allControlEnergies_emotionrec = readtable(readPath);

for i = 1:numel(contrastNames)
    currentContrast = contrastNames{i};
    currentContrast_idx = strcmp(allControlEnergies_emotionrec.contrast, currentContrast);
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(currentContrast_idx, :);
    writetable(allControlEnergies_emotionrec_currentContrast, strcat(baseName_emotionrec, '_', currentContrast, '.csv'));
end

end