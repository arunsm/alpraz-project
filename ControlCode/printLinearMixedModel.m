% function to print details of LinearMixedModel objects to a spreadsheet
function [] = printLinearMixedModel(inputModel, outputFileName)
    % extract model formula, stats for fixed and random effects
    modelFormula = inputModel.Formula;
    [~, ~, fixedEffectsStats] = inputModel.fixedEffects;
    [~, ~, stats] = inputModel.covarianceParameters; 
    randomEffectsStats = stats{1};
    
    % convert dataset arrays to tables
    fixedEffectsStatsTable = dataset2table(fixedEffectsStats);
    randomEffectsStatsTable = dataset2table(randomEffectsStats);
    
    % print data to file
    writetable(fixedEffectsStatsTable, outputFileName, 'FileType', 'spreadsheet', 'Sheet', 'fixedEffects');
    writetable(randomEffectsStatsTable, outputFileName, 'FileType', 'spreadsheet', 'Sheet', 'randomEffects');
end