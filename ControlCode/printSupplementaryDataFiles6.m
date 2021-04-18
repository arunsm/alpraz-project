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
    
    % convert categorical variables to matlab categorical variables
    allControlEnergies_emotionid_currentContrast.drug = categorical(allControlEnergies_emotionid_currentContrast.drug);
    allControlEnergies_emotionid_currentContrast.group = categorical(allControlEnergies_emotionid_currentContrast.group);
    allControlEnergies_emotionid_currentContrast.gender = categorical(allControlEnergies_emotionid_currentContrast.gender);
    
    % normalize continuous variables
    %allControlEnergies_emotionid_currentContrast.age = normalize(allControlEnergies_emotionid_currentContrast.age);
    
    % calculate efficiency = accuracy/reactionTime
    efficiency_corrthreat = allControlEnergies_emotionid_currentContrast.pctcorr_threat./allControlEnergies_emotionid_currentContrast.rtmdn_threatcorr;
    efficiency_corrnonthreat = allControlEnergies_emotionid_currentContrast.pctcorr_nonthreat./allControlEnergies_emotionid_currentContrast.rtmdn_nonthreatcorr;
    efficiency_corrneutral = allControlEnergies_emotionid_currentContrast.pctcorr_neutral./allControlEnergies_emotionid_currentContrast.rtmdn_neutralcorr;
    
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrthreat);
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrnonthreat);
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrneutral);
    
    % split persistence energy into within-person (time variant) and
    % between-person (time invariant) components
%     grand_mean_PE = mean(allControlEnergies_emotionid_currentContrast.persistence_allNodes);
%     grand_mean_centered_PE = allControlEnergies_emotionid_currentContrast.persistence_allNodes - grand_mean_PE;
%     persistence_between_person = zeros(numel(grand_mean_centered_PE), 1);
%     for j = 1:numel(grand_mean_centered_PE)
%         current_subject_ID = allControlEnergies_emotionid_currentContrast.subjectID(j);
%         persistence_between_person(j) = mean(grand_mean_centered_PE(allControlEnergies_emotionid_currentContrast.subjectID == current_subject_ID));
%     end
%     persistence_within_person = grand_mean_centered_PE - persistence_between_person;
%     
%     allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, persistence_between_person);
%     allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, persistence_within_person);
    
    switch currentContrast
        case 'contrast1_threatcorrectStd'
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast3_nonthreatcorrectStd'
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrnonthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast5_neutralcorrectStd'
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrneutral ~  persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    end
    
    save(strcat(resultsDirCurrentFigure, 'efficiencyModel_emotionid_', currentContrast, '.mat'), 'efficiencyModel_emotionid_currentContrast');
    printFileName = strcat(resultsDirCurrentFigure, 'efficiencyModel_emotionid_', currentContrast, '.xlsx');
    printLinearMixedModel(efficiencyModel_emotionid_currentContrast, printFileName);
end

% emotionrec
for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(allControlEnergies_emotionrec.contrast, currentContrast), :);
    
    % convert categorical variables to matlab categorical variables
    allControlEnergies_emotionrec_currentContrast.drug = categorical(allControlEnergies_emotionrec_currentContrast.drug);
    allControlEnergies_emotionrec_currentContrast.group = categorical(allControlEnergies_emotionrec_currentContrast.group);
    allControlEnergies_emotionrec_currentContrast.gender = categorical(allControlEnergies_emotionrec_currentContrast.gender);
        
    % normalize continuous variables
    %allControlEnergies_emotionrec_currentContrast.age = normalize(allControlEnergies_emotionrec_currentContrast.age);
    
    % calculating efficiency = accuracy/reactionTime
    efficiency_corrthreat = allControlEnergies_emotionrec_currentContrast.pctcorr_threat./allControlEnergies_emotionrec_currentContrast.rtmdn_threatcorr;
    efficiency_corrnonthreat = allControlEnergies_emotionrec_currentContrast.pctcorr_nonthreat./allControlEnergies_emotionrec_currentContrast.rtmdn_nonthreatcorr;
    efficiency_corrneutral = allControlEnergies_emotionrec_currentContrast.pctcorr_neutral./allControlEnergies_emotionrec_currentContrast.rtmdn_neutralcorr;
    
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrthreat);
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrnonthreat);
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrneutral);
    
    % split persistence energy into within-person (time variant) and
    % between-person (time invariant) components
%     grand_mean_PE = mean(allControlEnergies_emotionrec_currentContrast.persistence_allNodes);
%     grand_mean_centered_PE = allControlEnergies_emotionrec_currentContrast.persistence_allNodes - grand_mean_PE;
%     persistence_between_person = zeros(numel(grand_mean_centered_PE), 1);
%     for j = 1:numel(grand_mean_centered_PE)
%         current_subject_ID = allControlEnergies_emotionrec_currentContrast.subjectID(j);
%         persistence_between_person(j) = mean(grand_mean_centered_PE(allControlEnergies_emotionrec_currentContrast.subjectID == current_subject_ID));
%     end
%     persistence_within_person = grand_mean_centered_PE - persistence_between_person;
%     
%     allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, persistence_between_person);
%     allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, persistence_within_person);
    
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