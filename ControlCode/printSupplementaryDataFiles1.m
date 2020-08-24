function [] = printSupplementaryDataFiles1(allControlEnergies_emotionid, allControlEnergies_emotionrec)

parameters
resultsDirCurrentFigure = strcat(resultsDir, filesep, 'FigureS6', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% create mixed models to examine effects of clinical and demographic variables on persistence_allNodes

% emotion id
contrast = allControlEnergies_emotionid.contrast;
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
    mixedModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'persistence_allNodes ~ drug + group + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    save(strcat(resultsDirCurrentFigure, 'mixedModel_emotionid_', currentContrast, '.mat'), 'mixedModel_emotionid_currentContrast');
    printFileName = strcat(resultsDirCurrentFigure, 'mixedModel_emotionid_', currentContrast, '.xlsx');
    printLinearMixedModel(mixedModel_emotionid_currentContrast, printFileName);
end

% emotion memory
contrast = allControlEnergies_emotionrec.contrast;
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
    mixedModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'persistence_allNodes ~ drug + group + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    save(strcat(resultsDirCurrentFigure, 'mixedModel_emotionrec_', currentContrast, '.mat'), 'mixedModel_emotionrec_currentContrast');
    printFileName = strcat(resultsDirCurrentFigure, 'mixedModel_emotionrec_', currentContrast, '.xlsx');
    printLinearMixedModel(mixedModel_emotionrec_currentContrast, printFileName);
end