function [] = printSupplementaryDataFiles1(allControlEnergies_emotionid, allControlEnergies_emotionrec)

parameters
resultsDirCurrentFigure = strcat(resultsDir, 'SupplementaryDataFiles1', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% create mixed models to examine effects of clinical and demographic variables on persistence_allNodes

contrast = allControlEnergies_emotionid.contrast;
% emotion id
for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(contrast, currentContrast), :);
    
    % converting categorical variables
    allControlEnergies_emotionid_currentContrast.group = categorical(allControlEnergies_emotionid_currentContrast.group);
    allControlEnergies_emotionid_currentContrast.drug = categorical(allControlEnergies_emotionid_currentContrast.drug);
    allControlEnergies_emotionid_currentContrast.gender = categorical(allControlEnergies_emotionid_currentContrast.gender);
    
    % fitting model after checking for normality
    %persistence_allNodes = allControlEnergies_emotionid_currentContrast.persistence_allNodes;
    %h = kstest(zscore(persistence_allNodes));
    %figure, histogram(persistence_allNodes); title(strcat('KS test: ', num2str(h)));
    mixedModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'persistence_allNodes ~ drug + group + drug*group + gender + age + avge_FD + (1|subjectID)', 'FitMethod', 'ML');
    save(strcat(resultsDirCurrentFigure, 'mixedModel_emotionid_', currentContrast, '.mat'), 'mixedModel_emotionid_currentContrast');
    printFileName = strcat(resultsDirCurrentFigure, 'mixedModel_emotionid_', currentContrast, '.xlsx');
    printLinearMixedModel(mixedModel_emotionid_currentContrast, printFileName);
end

contrast = allControlEnergies_emotionrec.contrast;
% emotion memory
for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(contrast, currentContrast), :);
    
    % converting categorical variables
    allControlEnergies_emotionrec_currentContrast.group = categorical(allControlEnergies_emotionrec_currentContrast.group);
    allControlEnergies_emotionrec_currentContrast.drug = categorical(allControlEnergies_emotionrec_currentContrast.drug);
    allControlEnergies_emotionrec_currentContrast.gender = categorical(allControlEnergies_emotionrec_currentContrast.gender);
    
    % fitting model
    %persistence_allNodes = allControlEnergies_emotionrec_currentContrast.persistence_allNodes;
    %h = kstest(zscore(persistence_allNodes));
    %figure, histogram(persistence_allNodes); title(strcat('KS test: ', num2str(h)));
    mixedModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'persistence_allNodes ~ drug + group + drug*group + gender + age + avge_FD + (1|subjectID)', 'FitMethod', 'ML');
    save(strcat(resultsDirCurrentFigure, 'mixedModel_emotionrec_', currentContrast, '.mat'), 'mixedModel_emotionrec_currentContrast');
    printFileName = strcat(resultsDirCurrentFigure, 'mixedModel_emotionrec_', currentContrast, '.xlsx');
    printLinearMixedModel(mixedModel_emotionrec_currentContrast, printFileName);
end

% emotion id
% allControlEnergies_emotionid.group = categorical(allControlEnergies_emotionid.group);
% allControlEnergies_emotionid.drug = categorical(allControlEnergies_emotionid.drug);
% allControlEnergies_emotionid.gender = categorical(allControlEnergies_emotionid.gender);
% allControlEnergies_emotionid.contrast = categorical(allControlEnergies_emotionid.contrast);
% allControlEnergies_emotionid.contrast = reordercats(allControlEnergies_emotionid.contrast, {'contrast5_neutralcorrectStd', 'contrast3_nonthreatcorrectStd', 'contrast1_threatcorrectStd'});
% 
% mixedModel_emotionid = fitlme(allControlEnergies_emotionid, 'persistence_allNodes ~ drug*group*contrast + gender + age + (1|subjectID)', 'FitMethod', 'ML');
% save(strcat(resultsDirCurrentFigure, 'mixedModel_emotionid.mat'), 'mixedModel_emotionid');
% printFileName = strcat(resultsDirCurrentFigure, 'mixedModel_emotionid.xlsx');
% printLinearMixedModel(mixedModel_emotionid, printFileName);
% 
% % emotion rec
% allControlEnergies_emotionrec.group = categorical(allControlEnergies_emotionrec.group);
% allControlEnergies_emotionrec.drug = categorical(allControlEnergies_emotionrec.drug);
% allControlEnergies_emotionrec.gender = categorical(allControlEnergies_emotionrec.gender);
% allControlEnergies_emotionrec.contrast = categorical(allControlEnergies_emotionrec.contrast);
% allControlEnergies_emotionrec.contrast = reordercats(allControlEnergies_emotionrec.contrast, {'contrast5_neutralcorrectStd', 'contrast3_nonthreatcorrectStd', 'contrast1_threatcorrectStd'});
% 
% mixedModel_emotionrec = fitlme(allControlEnergies_emotionrec, 'persistence_allNodes ~ drug*group*contrast + gender + age + (1|subjectID)', 'FitMethod', 'ML');
% save(strcat(resultsDirCurrentFigure, 'mixedModel_emotionrec.mat'), 'mixedModel_emotionrec');
% printFileName = strcat(resultsDirCurrentFigure, 'mixedModel_emotionrec.xlsx');
% printLinearMixedModel(mixedModel_emotionrec, printFileName);
