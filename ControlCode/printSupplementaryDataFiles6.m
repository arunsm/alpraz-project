function [] = printSupplementaryDataFiles6(allControlEnergies_emotionid, allControlEnergies_emotionrec)

parameters
resultsDirCurrentFigure = strcat(resultsDir, filesep, 'SupplementaryDataFiles6', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% mixed models exploring relationship between persistence energy and behavioral efficiency

% emotionid
for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(allControlEnergies_emotionid.contrast, currentContrast), :);
    
    allControlEnergies_emotionid_currentContrast.drug = categorical(allControlEnergies_emotionid_currentContrast.drug);
    allControlEnergies_emotionid_currentContrast.group = categorical(allControlEnergies_emotionid_currentContrast.group);
    allControlEnergies_emotionid_currentContrast.gender = categorical(allControlEnergies_emotionid_currentContrast.gender);
    
    % calculate efficiency = accuracy/reactionTime
    efficiency_corrthreat = allControlEnergies_emotionid_currentContrast.pctcorr_threat./allControlEnergies_emotionid_currentContrast.rtmdn_threatcorr;
    efficiency_corrnonthreat = allControlEnergies_emotionid_currentContrast.pctcorr_nonthreat./allControlEnergies_emotionid_currentContrast.rtmdn_nonthreatcorr;
    efficiency_corrneutral = allControlEnergies_emotionid_currentContrast.pctcorr_neutral./allControlEnergies_emotionid_currentContrast.rtmdn_neutralcorr;
    
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrthreat);
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrnonthreat);
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrneutral);
    
    switch currentContrast
        case 'contrast1_threatcorrectStd'
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast3_nonthreatcorrectStd'
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrnonthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast5_neutralcorrectStd'
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrneutral ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    end
    
    save(strcat(resultsDirCurrentFigure, 'efficiencyModel_emotionid_', currentContrast, '.mat'), 'efficiencyModel_emotionid_currentContrast');
    printFileName = strcat(resultsDirCurrentFigure, 'efficiencyModel_emotionid_', currentContrast, '.xlsx');
    printLinearMixedModel(efficiencyModel_emotionid_currentContrast, printFileName);
end

% emotionrec
for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(allControlEnergies_emotionrec.contrast, currentContrast), :);
    
    allControlEnergies_emotionrec_currentContrast.drug = categorical(allControlEnergies_emotionrec_currentContrast.drug);
    allControlEnergies_emotionrec_currentContrast.group = categorical(allControlEnergies_emotionrec_currentContrast.group);
    allControlEnergies_emotionrec_currentContrast.gender = categorical(allControlEnergies_emotionrec_currentContrast.gender);
    
    % calculating efficiency = accuracy/reactionTime
    efficiency_corrthreat = allControlEnergies_emotionrec_currentContrast.pctcorr_threat./allControlEnergies_emotionrec_currentContrast.rtmdn_threatcorr;
    efficiency_corrnonthreat = allControlEnergies_emotionrec_currentContrast.pctcorr_nonthreat./allControlEnergies_emotionrec_currentContrast.rtmdn_nonthreatcorr;
    efficiency_corrneutral = allControlEnergies_emotionrec_currentContrast.pctcorr_neutral./allControlEnergies_emotionrec_currentContrast.rtmdn_neutralcorr;
    
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrthreat);
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrnonthreat);
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrneutral);

    switch currentContrast
        case 'contrast1_threatcorrectStd'
            efficiencyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'efficiency_corrthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast3_nonthreatcorrectStd'
            efficiencyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'efficiency_corrnonthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast5_neutralcorrectStd'
            efficiencyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'efficiency_corrneutral ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    end
    
    save(strcat(resultsDirCurrentFigure, 'efficiencyModel_emotionrec_', currentContrast, '.mat'), 'efficiencyModel_emotionrec_currentContrast');
    printFileName = strcat(resultsDirCurrentFigure, 'efficiencyModel_emotionrec_', currentContrast, '.xlsx');
    printLinearMixedModel(efficiencyModel_emotionrec_currentContrast, printFileName);
end
end