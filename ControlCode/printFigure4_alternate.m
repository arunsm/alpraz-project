function [] = printFigure4_alternate(allControlEnergies_emotionid, allControlEnergies_emotionrec)

parameters
importPETdata
importAHBA

resultsDirCurrentFigure = strcat(resultsDir, filesep, 'Figure4', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end
%% plot control input beta coefficients vs PET neurotransmitter maps - emotionid

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlInputs_persistence; % extracting control input
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        current_controlInput_emotionid = trapz(controlInput_emotionid{i}.^2)/nTimeSteps;
        controlInput_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionid;
    end
    
    % evaluate effects of drug and group on control impact of each node
    betas_groupMain = NaN(nNodes, 1); pValues_groupMain = NaN(nNodes, 1);
    betas_drugMain = NaN(nNodes, 1); pValues_drugMain = NaN(nNodes, 1);
    betas_groupDrugInteraction = NaN(nNodes, 1); pValues_groupDrugInteraction = NaN(nNodes, 1);
    
    % looping over nodes, building a model for each
    for i = 1:nNodes
        %fprintf('node%d\n', i);
        controlInput_currentNode = abs(controlInput_emotionid_allSubjects(:, i)); % taking absolute value of total control input for current node
        controlInput_currentNode = zscore(controlInput_currentNode);
        if sum(isnan(controlInput_currentNode)) > 0 % skipping node if one or more subjects does not have parcel coverage in slab
            continue;
        else
            allControlTrajectories_currentNode = addvars(allControlTrajectories_emotionid_currentContrast, controlInput_currentNode); % creating table with additional variable
            allControlTrajectories_currentNode.drug = 1 - allControlTrajectories_currentNode.drug; % NOTE: inverting drug indicator to drug(0, 1) = (placebo, alpraz) for ease of interpretation
            % converting to categorical variables
            allControlTrajectories_currentNode.drug = categorical(allControlTrajectories_currentNode.drug);
            allControlTrajectories_currentNode.group = categorical(allControlTrajectories_currentNode.group);
            controlInputModel = fitlme(allControlTrajectories_currentNode, 'controlInput_currentNode ~ drug + group + drug*group + (1|subjectID)', 'FitMethod', 'ML'); % fitting mixed model
            coefficientNames = controlInputModel.CoefficientNames; coefficients = controlInputModel.Coefficients.Estimate; pValues = controlInputModel.Coefficients.pValue; % extracting coefficients and p-values
            betas_groupMain(i) = coefficients(strcmp(coefficientNames, 'group_1')); betas_drugMain(i) = coefficients(strcmp(coefficientNames, 'drug_1')); betas_groupDrugInteraction(i) = coefficients(strcmp(coefficientNames, 'group_1:drug_1'));
            pValues_groupMain(i) = pValues(strcmp(coefficientNames, 'group_1')); pValues_drugMain(i) = pValues(strcmp(coefficientNames, 'drug_1')); pValues_groupDrugInteraction(i) = pValues(strcmp(coefficientNames, 'group_1:drug_1'));
        end
    end
    
    %betas_groupMain = abs(betas_groupMain); betas_drugMain = abs(betas_drugMain); betas_groupDrugInteraction = abs(betas_groupDrugInteraction); % taking absolute values of coefficients
    
    %%%%% plotting coefficients of drug main effect versus PET neurotransmitter profiles %%%%%
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([10, 80]);
    
    plot(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.8, 20, strcat('\rho=', num2str(rho)), 'Color', 'k', 'FontSize', 12);
    text(0.8, 15, strcat('p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of drug main effect (z-score)');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputDrugEffect_PET_GABAa_FLUMAZENIL_c11.eps'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([10, 80]);
    
    plot(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.8, 20, strcat('\rho=', num2str(rho)), 'Color', 'k', 'FontSize', 12);
    text(0.8, 15, strcat('p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group x drug interaction (z-score)');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputGroupDrugInteraction_PET_GABAa_FLUMAZENIL_c11.eps'));
    
    %%%% plotting coefficients of drug main effect versus GABA receptor expression from Allen atlas %%%%
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_drugMain, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_drugMain, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
   
    lsline;
    xlabel('coefficients of drug main effect (z-score)');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputDrugEffect_AHBA_GABA.eps'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_groupDrugInteraction, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_groupDrugInteraction, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
   
    lsline;
    xlabel('coefficients of group x drug interaction (z-score)');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputGroupDrugInteraction_AHBA_GABA.eps'));
end

%% plot control input beta coefficients vs PET neurotransmitter maps - emotionrec

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlInputs_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = trapz(controlInput_emotionrec{i}.^2)/nTimeSteps;
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    % evaluate effects of drug and group on control impact of each node
    betas_groupMain = NaN(nNodes, 1); pValues_groupMain = NaN(nNodes, 1);
    betas_drugMain = NaN(nNodes, 1); pValues_drugMain = NaN(nNodes, 1);
    betas_groupDrugInteraction = NaN(nNodes, 1); pValues_groupDrugInteraction = NaN(nNodes, 1);
    
    for i = 1:nNodes
        %fprintf('node%d\n', i);
        controlInput_currentNode = abs(controlInput_emotionrec_allSubjects(:, i));
        controlInput_currentNode =  zscore(controlInput_currentNode);
        if sum(isnan(controlInput_currentNode)) > 0 % skipping node if one or more subjects does not have parcel coverage in slab
            continue;
        else
            allControlTrajectories_currentNode = addvars(allControlTrajectories_emotionrec_currentContrast, controlInput_currentNode); % creating table with additional variable
            allControlTrajectories_currentNode.drug = 1 - allControlTrajectories_currentNode.drug; % NOTE: inverting drug indicator to drug(0, 1) = (placebo, alpraz) for ease of interpretation
            % converting to categorical variables
            allControlTrajectories_currentNode.drug = categorical(allControlTrajectories_currentNode.drug);
            allControlTrajectories_currentNode.group = categorical(allControlTrajectories_currentNode.group);
            controlInputModel = fitlme(allControlTrajectories_currentNode, 'controlInput_currentNode ~ drug + group + drug*group + (1|subjectID)', 'FitMethod', 'ML'); % fitting mixed model
            coefficientNames = controlInputModel.CoefficientNames; coefficients = controlInputModel.Coefficients.Estimate; pValues = controlInputModel.Coefficients.pValue; % extracting coefficients and p-values
            betas_groupMain(i) = coefficients(strcmp(coefficientNames, 'group_1')); betas_drugMain(i) = coefficients(strcmp(coefficientNames, 'drug_1')); betas_groupDrugInteraction(i) = coefficients(strcmp(coefficientNames, 'group_1:drug_1'));
            pValues_groupMain(i) = pValues(strcmp(coefficientNames, 'group_1')); pValues_drugMain(i) = pValues(strcmp(coefficientNames, 'drug_1')); pValues_groupDrugInteraction(i) = pValues(strcmp(coefficientNames, 'group_1:drug_1'));
        end
    end
    
    %betas_groupMain = abs(betas_groupMain); betas_drugMain = abs(betas_drugMain); betas_groupDrugInteraction = abs(betas_groupDrugInteraction); % taking absolute values of coefficients
    
%     pValues_groupMain_FDR = mafdr(pValues_groupMain, 'BHFDR', true);
%     pValues_drugMain_FDR = mafdr(pValues_drugMain, 'BHFDR', true);
%     pValues_groupDrugInteraction_FDR = mafdr(pValues_groupDrugInteraction, 'BHFDR', true);
%     
%     sum(pValues_groupMain_FDR < 0.05/234)
%     sum(pValues_drugMain_FDR < 0.05/234)
%     sum(pValues_groupDrugInteraction_FDR < 0.05/234)
    
    % plotting coefficients of group main effect versus PET
    % neurotransmitter profiles
%     f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
%     set(gca, 'FontSize', 20);
%     
%     xlim([0, 1]);
%     ylim([10, 80]);
%     
%     plot(betas_groupMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
%     [rho, pVal_corr] = corr(betas_groupMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
%     text(0.1, 78, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
%     
%     lsline;
%     xlabel('coefficients of group main effect');
%     ylabel('GABA receptor density');
%     saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupEffect_GABA.eps'));
    
    % plotting coefficients of drug main effect versus PET
    % neurotransmitter profiles
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([10, 80]);
    
    plot(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.8, 20, strcat('\rho=', num2str(rho)), 'Color', 'k', 'FontSize', 12);
    text(0.8, 15, strcat('p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of drug main effect (z-score)');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDrugEffect_PET_GABAa_FLUMAZENIL_c11.eps'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([10, 80]);
    
    plot(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.8, 20, strcat('\rho=', num2str(rho)), 'Color', 'k', 'FontSize', 12);
    text(0.8, 15, strcat('p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group x drug interaction (z-scores)');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupDrugInteraction_PET_GABAa_FLUMAZENIL_c11.eps'));

    %%%% plotting coefficients of drug main effect versus GABA receptor expression from Allen atlas %%%%
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_drugMain, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_drugMain, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
   
    lsline;
    xlabel('coefficients of drug main effect (z-score)');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDrugEffect_AHBA_GABA.eps'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_groupDrugInteraction, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_groupDrugInteraction, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
   
    lsline;
    xlabel('coefficients of group x drug interaction (z-score)');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupDrugInteraction_AHBA_GABA.eps'));
end
end