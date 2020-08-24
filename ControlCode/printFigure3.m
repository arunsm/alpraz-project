function [] = printFigure3(allControlEnergies_emotionid, allControlEnergies_emotionrec)

parameters
resultsDirCurrentFigure = strcat(resultsDir, filesep, 'Figure3', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% compile coefficients and p-values from mixed models

pValues_persistence = zeros(1, 6); % compile all p-values for 6 tests
coefficients_persistence = zeros(1, 6); % compile all coefficient values for 6 tests
ctr = 1;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    
    % import p-values from mixed models
    mixedModel = strcat(resultsDir, 'SupplementaryDataFiles5', filesep, 'efficiencyModel_emotionid_', currentContrast, '.mat');
    load(mixedModel, 'efficiencyModel_emotionid_currentContrast');
    coefficient_estimates = efficiencyModel_emotionid_currentContrast.Coefficients.Estimate;
    coefficient_names = efficiencyModel_emotionid_currentContrast.Coefficients.Name;
    coefficient_pValues = efficiencyModel_emotionid_currentContrast.Coefficients.pValue;
    
    pValues_persistence(ctr) = coefficient_pValues(strcmp(coefficient_names, 'persistence_allNodes'));
    coefficients_persistence(ctr) = coefficient_estimates(strcmp(coefficient_names, 'persistence_allNodes'));
    ctr = ctr + 1;  
end

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    
    % import p-values from mixed models
    mixedModel = strcat(resultsDir, 'SupplementaryDataFiles5', filesep, 'efficiencyModel_emotionrec_', currentContrast, '.mat');
    load(mixedModel, 'efficiencyModel_emotionrec_currentContrast');
    coefficient_estimates = efficiencyModel_emotionrec_currentContrast.Coefficients.Estimate;
    coefficient_names = efficiencyModel_emotionrec_currentContrast.Coefficients.Name;
    coefficient_pValues = efficiencyModel_emotionrec_currentContrast.Coefficients.pValue;
    
    pValues_persistence(ctr) = coefficient_pValues(strcmp(coefficient_names, 'persistence_allNodes'));
    coefficients_persistence(ctr) = coefficient_estimates(strcmp(coefficient_names, 'persistence_allNodes'));
    ctr = ctr + 1;  
end

pValues_persistence_FDR = mafdr(pValues_persistence, 'BHFDR', true) % correcting for multiple comparisons
coefficients_persistence

%% plot task efficiency versus persistence energy for all contrasts in emotionid

for c = 1:nContrasts
    % import accuracy, reaction time and persistence
    currentContrast = contrastLabels{c};
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(allControlEnergies_emotionid.contrast, currentContrast), :);
    
    contrastString = currentContrast(11:end); contrastString = strrep(contrastString, 'correctStd', ''); 
    accuracy_contrastString = strcat('pctcorr_', contrastString); rtmdn_contrastString = strcat('rtmdn_', contrastString, 'corr');
    accuracy = allControlEnergies_emotionid_currentContrast.(accuracy_contrastString);
    rtmdn = allControlEnergies_emotionid_currentContrast.(rtmdn_contrastString);
    efficiency = accuracy./rtmdn;
    persistence = allControlEnergies_emotionid_currentContrast.persistence_allNodes;
    
    % scatter plots
    f = figure('Visible', 'off'); set(gcf, 'color', 'w'); hold on;
    set(gca, 'FontSize', 20);
    plot(efficiency, persistence, 'k.', 'MarkerSize', 20);
    h = lsline; h.LineWidth = 2; 
    ylabel('persistence energy');
    xlabel('efficiency');
    ylim([0.2, 0.65]);
    xlim([0, 0.6]);
 
    savePath = strcat(resultsDirCurrentFigure, 'emotionid_efficiency_persistenceEnergy_', currentContrast, '.eps'); saveas(f, savePath); close(f);
end

%% plot task efficiency versus persistence energy for all contrasts in emotionrec

for c = 1:nContrasts
    % import accuracy, reaction time and persistence
    currentContrast = contrastLabels{c};
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(allControlEnergies_emotionrec.contrast, currentContrast), :);
    
    contrastString = currentContrast(11:end); contrastString = strrep(contrastString, 'correctStd', ''); 
    accuracy_contrastString = strcat('pctcorr_', contrastString); rtmdn_contrastString = strcat('rtmdn_', contrastString, 'corr');
    accuracy = allControlEnergies_emotionrec_currentContrast.(accuracy_contrastString);
    rtmdn = allControlEnergies_emotionrec_currentContrast.(rtmdn_contrastString);
    efficiency = accuracy./rtmdn;
    persistence = allControlEnergies_emotionrec_currentContrast.persistence_allNodes;
    
    % scatter plots
    f = figure('Visible', 'off'); set(gcf, 'color', 'w'); hold on;
    set(gca, 'FontSize', 20);
    plot(efficiency, persistence, 'k.', 'MarkerSize', 20);
    h = lsline; h.LineWidth = 2; 
    ylabel('persistence energy');
    xlabel('efficiency');
    ylim([0.2, 0.65]);
    xlim([0, 0.4]);
 
    savePath = strcat(resultsDirCurrentFigure, 'emotionrec_efficiency_persistenceEnergy_', currentContrast, '.eps'); saveas(f, savePath); close(f);
end

end