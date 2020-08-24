function [] = printSupplementaryDataFiles8()
parameters
nIterations = 500;

resultsDirCurrentFigure = strcat(resultsDir, filesep, 'SupplementaryDataFiles8', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% compare test statistics from mixed model against structural null models

% emotionid

% load results of null model
resultsDir_nullModel = 'Results/minimalControl/avge_FD_thresh_0.5_parcelCoverageThresh_0.5_EuclideanNorm_allNodes_nullModel2_betas/';
mixedModelFolder = strcat(resultsDir, 'SupplementaryDataFiles1');
contrast = allControlEnergies_emotionid.contrast;

for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    fprintf(currentContrast); fprintf('\n');
    
    % load coefficients from actual mixed model
    currentMixedModelPath = strcat(mixedModelFolder, filesep, 'mixedModel_emotionid_', currentContrast, '.mat');
    load(currentMixedModelPath, 'mixedModel_emotionid_currentContrast');
    coefficientEstimates_emotionid_currentContrast = mixedModel_emotionid_currentContrast.Coefficients.Estimate;
    coefficientNames_emotionid_currentContrast = mixedModel_emotionid_currentContrast.Coefficients.Name;
    nCoefficients = numel(coefficientNames_emotionid_currentContrast);
    
    % fit mixed models to null data and compute coefficients
    coefficientEstimates_emotionid_nullModel_currentContrast = zeros(nCoefficients, nIterations);
    pValues_emotionid_nullModel_currentContrast = zeros(nCoefficients, 1);
    for j = 1:nIterations
        fprintf('iteration %d\n', j);
        allControlEnergies_emotionid_nullModel = readtable(strcat(resultsDir_nullModel, 'allControlEnergies_emotionid_iteration', num2str(j), '.csv'));
        allControlEnergies_emotionid_nullModel_currentContrast = allControlEnergies_emotionid_nullModel(strcmp(contrast, currentContrast), :);
        
        % converting categorical variables
        allControlEnergies_emotionid_nullModel_currentContrast.group = categorical(allControlEnergies_emotionid_nullModel_currentContrast.group);
        allControlEnergies_emotionid_nullModel_currentContrast.drug = categorical(allControlEnergies_emotionid_nullModel_currentContrast.drug);
        allControlEnergies_emotionid_nullModel_currentContrast.gender = categorical(allControlEnergies_emotionid_nullModel_currentContrast.gender);
        
        currentModelFormula = strcat('persistence_allNodes ~ drug + group + drug*group + gender + age + (1|subjectID)');
        mixedModel_emotionid_nullModel_currentContrast = fitlme(allControlEnergies_emotionid_nullModel_currentContrast, currentModelFormula, 'FitMethod', 'ML');
        coefficientEstimates_emotionid_nullModel_currentContrast(:, j) = mixedModel_emotionid_nullModel_currentContrast.Coefficients.Estimate;
        pValues_emotionid_nullModel_currentContrast = pValues_emotionid_nullModel_currentContrast + double((abs(coefficientEstimates_emotionid_currentContrast) > abs(coefficientEstimates_emotionid_nullModel_currentContrast(:, j))));
    end
    
    % store coefficients and exact p-values
    pValues_emotionid_nullModel_currentContrast = 1 - pValues_emotionid_nullModel_currentContrast/nIterations;
    T = table(coefficientNames_emotionid_currentContrast, pValues_emotionid_nullModel_currentContrast, coefficientEstimates_emotionid_currentContrast, coefficientEstimates_emotionid_nullModel_currentContrast);
    saveFilePath = strcat(resultsDirCurrentFigure, 'structuralNullResults_emotionid_', currentContrast, '.csv');
    writetable(T, saveFilePath);
end

%% emotionrec

% load results of null model
resultsDir_nullModel = 'Results/minimalControl/avge_FD_thresh_0.5_parcelCoverageThresh_0.5_EuclideanNorm_allNodes_nullModel2_betas/';
mixedModelFolder = strcat(resultsDir, 'SupplementaryDataFiles1');

contrast = allControlEnergies_emotionrec.contrast;
for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    fprintf(currentContrast); fprintf('\n');
    
    % load coefficients from actual mixed model
    currentMixedModelPath = strcat(mixedModelFolder, filesep, 'mixedModel_emotionrec_', currentContrast, '.mat');
    load(currentMixedModelPath, 'mixedModel_emotionrec_currentContrast');
    coefficientEstimates_emotionrec_currentContrast = mixedModel_emotionrec_currentContrast.Coefficients.Estimate;
    coefficientNames_emotionrec_currentContrast = mixedModel_emotionrec_currentContrast.Coefficients.Name;
    nCoefficients = numel(coefficientNames_emotionrec_currentContrast);
    
    % fit mixed models to null data and compute coefficients
    coefficientEstimates_emotionrec_nullModel_currentContrast = zeros(nCoefficients, nIterations);
    pValues_emotionrec_nullModel_currentContrast = zeros(nCoefficients, 1);
    for j = 1:nIterations
        fprintf('iteration %d\n', j);
        allControlEnergies_emotionrec_nullModel = readtable(strcat(resultsDir_nullModel, 'allControlEnergies_emotionrec_iteration', num2str(j), '.csv'));
        allControlEnergies_emotionrec_nullModel_currentContrast = allControlEnergies_emotionrec_nullModel(strcmp(contrast, currentContrast), :);
        
        % converting categorical variables
        allControlEnergies_emotionrec_nullModel_currentContrast.group = categorical(allControlEnergies_emotionrec_nullModel_currentContrast.group);
        allControlEnergies_emotionrec_nullModel_currentContrast.drug = categorical(allControlEnergies_emotionrec_nullModel_currentContrast.drug);
        allControlEnergies_emotionrec_nullModel_currentContrast.gender = categorical(allControlEnergies_emotionrec_nullModel_currentContrast.gender);
        
        currentModelFormula = strcat('persistence_allNodes ~ drug + group + drug*group + gender + age + (1|subjectID)');
        mixedModel_emotionrec_nullModel_currentContrast = fitlme(allControlEnergies_emotionrec_nullModel_currentContrast, currentModelFormula, 'FitMethod', 'ML');
        coefficientEstimates_emotionrec_nullModel_currentContrast(:, j) = mixedModel_emotionrec_nullModel_currentContrast.Coefficients.Estimate;
        pValues_emotionrec_nullModel_currentContrast = pValues_emotionrec_nullModel_currentContrast + double((abs(coefficientEstimates_emotionrec_currentContrast) > abs(coefficientEstimates_emotionrec_nullModel_currentContrast(:, j))));
    end
    
    % store coefficients and exact p-values
    pValues_emotionrec_nullModel_currentContrast = 1 - pValues_emotionrec_nullModel_currentContrast/nIterations;
    T = table(coefficientNames_emotionrec_currentContrast, pValues_emotionrec_nullModel_currentContrast, coefficientEstimates_emotionrec_currentContrast, coefficientEstimates_emotionrec_nullModel_currentContrast);
    saveFilePath = strcat(resultsDirCurrentFigure, 'structuralNullResults_emotionrec_', currentContrast, '.csv');
    writetable(T, saveFilePath);
end
end