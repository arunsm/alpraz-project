function [] = printSupplementaryDataFiles8(allControlEnergies_emotionid, allControlEnergies_emotionrec)
parameters
nIterations = 500;

resultsDirCurrentFigure = strcat(resultsDir, filesep, 'SupplementaryDataFiles8', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

% compare test statistics from mixed model against structural null models

%% emotionid

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

    % getting avgeFD from main results folder (not calculated for null models)
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(contrast, currentContrast), :);
    avge_FD = allControlEnergies_emotionid_currentContrast.avge_FD;
    
    % fit mixed models to null data and compute coefficients
    coefficientEstimates_emotionid_nullModel_currentContrast = zeros(nCoefficients, nIterations);
    pValues_emotionid_nullModel_currentContrast = zeros(nCoefficients, nIterations);
    
    for j = 1:nIterations
        fprintf('iteration %d\n', j);
        
        allControlEnergies_emotionid_nullModel = readtable(strcat(resultsDir_nullModel, 'allControlEnergies_emotionid_iteration', num2str(j), '.csv'));
        allControlEnergies_emotionid_nullModel_currentContrast = allControlEnergies_emotionid_nullModel(strcmp(contrast, currentContrast), :);
        
        % converting categorical variables
        allControlEnergies_emotionid_nullModel_currentContrast.group = categorical(allControlEnergies_emotionid_nullModel_currentContrast.group);
        allControlEnergies_emotionid_nullModel_currentContrast.drug = categorical(allControlEnergies_emotionid_nullModel_currentContrast.drug);
        allControlEnergies_emotionid_nullModel_currentContrast.gender = categorical(allControlEnergies_emotionid_nullModel_currentContrast.gender);
        
        % adding head motion variable
        allControlEnergies_emotionid_nullModel_currentContrast = addvars(allControlEnergies_emotionid_nullModel_currentContrast, avge_FD);
        
        currentModelFormula = strcat('persistence_allNodes ~ drug + group + drug*group + avge_FD + gender + age + (1|subjectID)');
        mixedModel_emotionid_nullModel_currentContrast = fitlme(allControlEnergies_emotionid_nullModel_currentContrast, currentModelFormula, 'FitMethod', 'ML');
        coefficientEstimates_emotionid_nullModel_currentContrast(:, j) = mixedModel_emotionid_nullModel_currentContrast.Coefficients.Estimate;
        pValues_emotionid_nullModel_currentContrast(:, j) = mixedModel_emotionid_nullModel_currentContrast.Coefficients.pValue;
    end
    
    % store coefficients and fraction of significant iterations
    sum_significant_emotionid_nullModel_currentContrast = sum((pValues_emotionid_nullModel_currentContrast < 0.05), 2);
    T = table(coefficientNames_emotionid_currentContrast, sum_significant_emotionid_nullModel_currentContrast, coefficientEstimates_emotionid_currentContrast, coefficientEstimates_emotionid_nullModel_currentContrast);
    saveFilePath = strcat(resultsDirCurrentFigure, 'structuralNullResults_emotionid_', currentContrast, '.csv');
    writetable(T, saveFilePath);
    
    % plot coefficients for group, drug, and groupxdrug interaction
    coefficients2plot_null = [coefficientEstimates_emotionid_nullModel_currentContrast(strcmp(coefficientNames_emotionid_currentContrast, 'group_1'), :); ...
        coefficientEstimates_emotionid_nullModel_currentContrast(strcmp(coefficientNames_emotionid_currentContrast, 'drug_1'), :); ...
        coefficientEstimates_emotionid_nullModel_currentContrast(strcmp(coefficientNames_emotionid_currentContrast, 'group_1:drug_1'), :)];
    
    coefficients2plot_actual = [coefficientEstimates_emotionid_currentContrast(strcmp(coefficientNames_emotionid_currentContrast, 'group_1')); ...
        coefficientEstimates_emotionid_currentContrast(strcmp(coefficientNames_emotionid_currentContrast, 'drug_1')); ...
        coefficientEstimates_emotionid_currentContrast(strcmp(coefficientNames_emotionid_currentContrast, 'group_1:drug_1'))];
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    boxplot(coefficients2plot_null', 'Labels', {'group', 'drug', 'group x drug'}, 'LabelOrientation', 'inline', 'OutlierSize', 8, 'Symbol', 'k.');
    plot(coefficients2plot_actual, 'r*', 'MarkerSize', 10);
    ylim([-0.06, 0.06]);
    
    ylabel('Coefficients');
    saveas(f, strcat(resultsDirCurrentFigure, 'emotionid_', currentContrast, '_structuralNullCoefficients.eps'));
    close(f);
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
    
    % getting avgeFD from main results folder (not calculated for null models)
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(contrast, currentContrast), :);
    avge_FD = allControlEnergies_emotionrec_currentContrast.avge_FD;
    
    % fit mixed models to null data and compute coefficients
    coefficientEstimates_emotionrec_nullModel_currentContrast = zeros(nCoefficients, nIterations);
    pValues_emotionrec_nullModel_currentContrast = zeros(nCoefficients, nIterations);
    for j = 1:nIterations
        fprintf('iteration %d\n', j);
        allControlEnergies_emotionrec_nullModel = readtable(strcat(resultsDir_nullModel, 'allControlEnergies_emotionrec_iteration', num2str(j), '.csv'));
        allControlEnergies_emotionrec_nullModel_currentContrast = allControlEnergies_emotionrec_nullModel(strcmp(contrast, currentContrast), :);
        
        % converting categorical variables
        allControlEnergies_emotionrec_nullModel_currentContrast.group = categorical(allControlEnergies_emotionrec_nullModel_currentContrast.group);
        allControlEnergies_emotionrec_nullModel_currentContrast.drug = categorical(allControlEnergies_emotionrec_nullModel_currentContrast.drug);
        allControlEnergies_emotionrec_nullModel_currentContrast.gender = categorical(allControlEnergies_emotionrec_nullModel_currentContrast.gender);
        
        % adding head motion variable
        allControlEnergies_emotionrec_nullModel_currentContrast = addvars(allControlEnergies_emotionrec_nullModel_currentContrast, avge_FD);
        
        currentModelFormula = strcat('persistence_allNodes ~ drug + group + drug*group + avge_FD + gender + age + (1|subjectID)');
        mixedModel_emotionrec_nullModel_currentContrast = fitlme(allControlEnergies_emotionrec_nullModel_currentContrast, currentModelFormula, 'FitMethod', 'ML');
        coefficientEstimates_emotionrec_nullModel_currentContrast(:, j) = mixedModel_emotionrec_nullModel_currentContrast.Coefficients.Estimate;
        pValues_emotionrec_nullModel_currentContrast(:, j) = mixedModel_emotionrec_nullModel_currentContrast.Coefficients.pValue;
    end
    
    % store coefficients and fraction of significant iterations
    sum_significant_emotionrec_nullModel_currentContrast = sum((pValues_emotionrec_nullModel_currentContrast < 0.05), 2);
    T = table(coefficientNames_emotionrec_currentContrast, sum_significant_emotionrec_nullModel_currentContrast, coefficientEstimates_emotionrec_currentContrast, coefficientEstimates_emotionrec_nullModel_currentContrast);
    saveFilePath = strcat(resultsDirCurrentFigure, 'structuralNullResults_emotionrec_', currentContrast, '.csv');
    writetable(T, saveFilePath);
    
    % plot coefficients for group, drug, and groupxdrug interaction
    coefficients2plot_null = [coefficientEstimates_emotionrec_nullModel_currentContrast(strcmp(coefficientNames_emotionrec_currentContrast, 'group_1'), :); ...
        coefficientEstimates_emotionrec_nullModel_currentContrast(strcmp(coefficientNames_emotionrec_currentContrast, 'drug_1'), :); ...
        coefficientEstimates_emotionrec_nullModel_currentContrast(strcmp(coefficientNames_emotionrec_currentContrast, 'group_1:drug_1'), :)];
    
    coefficients2plot_actual = [coefficientEstimates_emotionrec_currentContrast(strcmp(coefficientNames_emotionrec_currentContrast, 'group_1')); ...
        coefficientEstimates_emotionrec_currentContrast(strcmp(coefficientNames_emotionrec_currentContrast, 'drug_1')); ...
        coefficientEstimates_emotionrec_currentContrast(strcmp(coefficientNames_emotionrec_currentContrast, 'group_1:drug_1'))];
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    boxplot(coefficients2plot_null', 'Labels', {'group', 'drug', 'group x drug'}, 'LabelOrientation', 'inline', 'OutlierSize', 8, 'Symbol', 'k.');
    plot(coefficients2plot_actual, 'r*', 'MarkerSize', 10);
    ylim([-0.06, 0.06]);
    
    saveas(f, strcat(resultsDirCurrentFigure, 'emotionrec_', currentContrast, '_structuralNullCoefficients.eps'));
end
end