function [] = printFigure4(allControlTrajectories_emotionid, allControlTrajectories_emotionrec)

parameters
importPETdata

resultsDirCurrentFigure = strcat(resultsDir, filesep, 'Figure4', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% plot matrix of PET maps arranged by Yeo7 systems

PET_maps = [PET_5HT1a_WAY_HC36 PET_5HT1b_P943_HC22 PET_5HT2a_ALT_HC19 PET_D1_SCH23390_c11 PET_D2_RACLOPRIDE_c11 PET_DAT_DATSPECT PET_FDOPA_f18 PET_GABAa_FLUMAZENIL_c11 PET_NAT_MRB_c11 PET_SERT_DASB_HC30];
PET_maps(idx_brainstem, :) = [];
PET_maps_Yeo7 = rearrangeMatrix_Yeo7(PET_maps);
nMaps = size(PET_maps_Yeo7, 2);
f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
imagesc(PET_maps_Yeo7); colormap redbluecmap; cb = colorbar; cb.Location = 'westoutside';
set(gca, 'FontSize', 20);
xlim([0.5, nMaps+0.5]); ylim([0, nNodes]);
% refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
% h = refline(0, refLines(1)); h.Color = 'g'; h.LineWidth = 2; %text(nSubjects, mean([0, refLines(1)]), 'Visual');
% h = refline(0, refLines(2)); h.Color = 'g'; h.LineWidth = 2; %text(nSubjects, mean([refLines(1), refLines(2)]), 'Somatomator');
% h = refline(0, refLines(3)); h.Color = 'g'; h.LineWidth = 2; %text(nSubjects, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
% h = refline(0, refLines(4)); h.Color = 'g'; h.LineWidth = 2; %text(nSubjects, mean([refLines(3), refLines(4)]), 'Ventral Attention');
% h = refline(0, refLines(5)); h.Color = 'g'; h.LineWidth = 2; %text(nSubjects, mean([refLines(4), refLines(5)]), 'Limbic');
% h = refline(0, refLines(6)); h.Color = 'g'; h.LineWidth = 2; %text(nSubjects, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
% h = refline(0, refLines(7)); h.Color = 'g'; h.LineWidth = 2; %text(nSubjects, mean([refLines(6), refLines(7)]), 'Default Mode');
% h = refline(0, refLines(8)); h.Color = 'g'; h.LineWidth = 2; %text(nSubjects, mean([refLines(7), refLines(8)]), 'Subcortical');

saveas(f, strcat(resultsDirCurrentFigure, 'PET_allMaps_Yeo7.svg'));
surfacePlots(PET_GABAa_FLUMAZENIL_c11, redbluecmap, [0 100], nifti, subcorticalIndices, resultsDirCurrentFigure, 'PET_GABAa_FLUMAZENIL_c11.svg');

%% creating random maps using spin test
% nRandomizations = 500;
% spinTestSavePath = strcat(resultsDirCurrentFigure, 'spinTestResults');
% fprintf('1\n'); PET_5HT1a_WAY_HC36_random = calculateScrambledBetas(PET_5HT1a_WAY_HC36, strcat(spinTestSavePath, '_PET_5HT1a_WAY_HC36'), nRandomizations);
% fprintf('2\n'); PET_5HT1b_P943_HC22_random = calculateScrambledBetas(PET_5HT1b_P943_HC22, strcat(spinTestSavePath, '_PET_5HT1b_P943_HC22'), nRandomizations);
% fprintf('3\n'); PET_5HT2a_ALT_HC19_random = calculateScrambledBetas(PET_5HT2a_ALT_HC19, strcat(spinTestSavePath, '_PET_5HT2a_ALT_HC19'), nRandomizations);
% fprintf('4\n'); PET_D1_SCH23390_c11_random = calculateScrambledBetas(PET_D1_SCH23390_c11, strcat(spinTestSavePath, '_PET_D1_SCH23390_c11'), nRandomizations);
% fprintf('5\n'); PET_D2_RACLOPRIDE_c11_random = calculateScrambledBetas(PET_D2_RACLOPRIDE_c11, strcat(spinTestSavePath, '_PET_D2_RACLOPRIDE_c11'), nRandomizations);
% fprintf('6\n'); PET_DAT_DATSPECT_random = calculateScrambledBetas(PET_DAT_DATSPECT, strcat(spinTestSavePath, '_PET_DAT_DATSPECT'), nRandomizations);
% fprintf('7\n'); PET_FDOPA_f18_random = calculateScrambledBetas(PET_FDOPA_f18, strcat(spinTestSavePath, '_PET_FDOPA_f18'), nRandomizations);
% fprintf('8\n'); PET_GABAa_FLUMAZENIL_c11_random = calculateScrambledBetas(PET_GABAa_FLUMAZENIL_c11, strcat(spinTestSavePath, '_PET_GABAa_FLUMAZENIL_c11'), nRandomizations);
% fprintf('9\n'); PET_NAT_MRB_c11_random = calculateScrambledBetas(PET_NAT_MRB_c11, strcat(spinTestSavePath, '_PET_NAT_MRB_c11'), nRandomizations);
% fprintf('10\n'); PET_SERT_DASB_HC30_random = calculateScrambledBetas(PET_SERT_DASB_HC30, strcat(spinTestSavePath, '_PET_SERT_DASB_HC30'), nRandomizations);

%% average control input difference w/ and w/o drug vs PET neurotransmitter maps - emotionid

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlInputs_persistence; % extracting control input
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    drug = allControlTrajectories_emotionid_currentContrast.drug;
    group = allControlTrajectories_emotionid_currentContrast.group;
    
    % populate matrix of control input for [nSubjects*2 x nNodes]
    controlInputs_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        current_controlInput_emotionid = trapz(controlInput_emotionid{i}.^2)/nTimeSteps; % calculating total control input for current node as the sum of integral of time-varying control input
        controlInputs_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionid;
    end
    
    % control input w/ alpraz and placebo
    controlInput_emotionid_currentContrast_alpraz = controlInputs_emotionid_allSubjects(drug==0, :);
    controlInput_emotionid_currentContrast_placebo = controlInputs_emotionid_allSubjects(drug==1, :);
    
    controlInputDiff_emotionid_currentContrast = abs(controlInput_emotionid_currentContrast_alpraz - controlInput_emotionid_currentContrast_placebo)'; % calculating absolute value of difference in control input w/ and w/o drug
    avge_controlInputDiff_emotionid_currentContrast = mean(controlInputDiff_emotionid_currentContrast, 2); % average across subjects
    
    % fit linear model
%     T = table(avge_controlInputDiff_emotionid_currentContrast, PET_5HT1a_WAY_HC36, PET_5HT1b_P943_HC22, PET_5HT2a_ALT_HC19, PET_D1_SCH23390_c11, PET_D2_RACLOPRIDE_c11, PET_DAT_DATSPECT, PET_FDOPA_f18, PET_GABAa_FLUMAZENIL_c11, PET_NAT_MRB_c11, PET_SERT_DASB_HC30);
%     linearModel_currentContrast = fitlm(T, 'ResponseVar', 'avge_controlInputDiff_emotionid_currentContrast');
%     spinTestSavePath = strcat(resultsDirCurrentFigure, 'linearModel_emotionid_', currentContrast, '.mat');
%     save(spinTestSavePath, 'linearModel_currentContrast');
    
    % find Lausanne nodes with high difference in control input b/w drug
    % and placebo
    controlInputDiff_emotionid_nodeNames = cell(nNodes, 2);
    for i = 1:nNodes
        controlInputDiff_emotionid_nodeNames{i, 1} = avge_controlInputDiff_emotionid_currentContrast(i);
        controlInputDiff_emotionid_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    controlInputDiff_emotionid_nodeNames = cell2table(controlInputDiff_emotionid_nodeNames, 'VariableNames', {'controlInputDiff_drug_placebo', 'nodeName_Lausanne'});
    controlInputDiff_emotionid_nodeNames = sortrows(controlInputDiff_emotionid_nodeNames, 'controlInputDiff_drug_placebo', 'descend', 'MissingPlacement', 'last'); % sorting nodes by descending value
    writetable(controlInputDiff_emotionid_nodeNames, strcat(resultsDirCurrentFigure, 'controlInputDiff_drug_placebo_emotionid_', currentContrast, '.csv'));
    
    % plot surface maps of differences in control input b/w alpraz and
    % placebo, averaged across subjects
    avge_controlInputDiff_emotionid_currentContrast(isnan(avge_controlInputDiff_emotionid_currentContrast)) = 0;
    surfacePlots([avge_controlInputDiff_emotionid_currentContrast; 0], redbluecmap, [0 0.015], nifti, subcorticalIndices, resultsDirCurrentFigure, strcat('emotionid_', currentContrast, 'controlInputDiff_heatmap_drug_vs_noDrug.svg'));
    
    % compute correlations between control input differences w/ and w/o drug against PET maps
    controlInputDiff_emotionid_currentContrast_PETatlasCorrs = corr(controlInputDiff_emotionid_currentContrast, PET_maps, 'Type', 'Spearman', 'Rows', 'Complete');
    avge_controlInputDiff_emotionid_currentContrast_PETatlasCorrs = mean(controlInputDiff_emotionid_currentContrast_PETatlasCorrs); % average correlation across subjects
    
    % compute correlations between control input difference against PET
    % maps, after randomizing placebo data or PET maps
    %nRandomizations = 500;
    rng(0); % reset random number generator
    randomCorrs = zeros(nRandomizations, numel(PETlabels));
    
%     for i = 1:nRandomizations
%         PET_maps_random = [PET_5HT1a_WAY_HC36_random(:, i) PET_5HT1b_P943_HC22_random(:, i) ...
%             PET_5HT2a_ALT_HC19_random(:, i) PET_D1_SCH23390_c11_random(:, i) PET_D2_RACLOPRIDE_c11_random(:, i) ...
%             PET_DAT_DATSPECT_random(:, i) PET_FDOPA_f18_random(:, i) PET_GABAa_FLUMAZENIL_c11_random(:, i) ...
%             PET_NAT_MRB_c11_random(:, i) PET_SERT_DASB_HC30_random(:, i)]; % compile random maps for i'th iteration
%         PET_maps_random(idx_brainstem, :) = [];
%         controlInputDiff_emotionid_currentContrast_PETatlasCorrs_random = corr(controlInputDiff_emotionid_currentContrast, PET_maps_random, 'Type', 'Spearman', 'Rows', 'Complete');
%         randomCorrs(i, :) = mean(controlInputDiff_emotionid_currentContrast_PETatlasCorrs_random);
%     end
    
        for i = 1:nRandomizations
            PET_maps_random = randomizeMatrix(PET_maps')';
            controlInputDiff_emotionid_currentContrast_PETatlasCorrs_random = corr(controlInputDiff_emotionid_currentContrast, PET_maps_random, 'Type', 'Spearman', 'Rows', 'Complete');
            randomCorrs(i, :) = mean(controlInputDiff_emotionid_currentContrast_PETatlasCorrs_random);
        end
    
    pValues = zeros(1, numel(PETlabels));
    for i = 1:numel(PETlabels)
        pValues(i) = sum(abs(randomCorrs(:, i)) > abs(avge_controlInputDiff_emotionid_currentContrast_PETatlasCorrs(i)));
    end
    pValues = (pValues/nRandomizations);
    pValues_FDR = mafdr(pValues, 'BHFDR', true);
    
    %controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs = atanh(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % Fisher z-transform
    %[h, pValues] = ttest(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % using one-sample t-test to check for significance of subject-wise correlations
    %[pValues, ~, ~, ~, ~] = mult_comp_perm_t1(controlInputDiff_emotionid_currentContrast_PETatlasCorrs);
    h1 = (pValues_FDR < 0.05) & (pValues_FDR > 0.005);
    h2 = (pValues_FDR < 0.005) & (pValues_FDR > 0.0005);
    h3 = (pValues_FDR < 0.0005);
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    ylim([-0.4, 0.4]);
    
    boxplot(controlInputDiff_emotionid_currentContrast_PETatlasCorrs, 'Labels', PETlabels, 'LabelOrientation', 'inline', 'OutlierSize', 8, 'Symbol', 'k.');
    refline(0, 0);
    for i = 1:numel(h1)
        if h1(i)
            text(i, 0, '*', 'Color', 'r', 'FontSize', 20);
        elseif h2(i)
            text(i, 0, '**', 'Color', 'r', 'FontSize', 20);
        elseif h3(i)
            text(i, 0, '***', 'Color', 'r', 'FontSize', 20);
        end
    end
    
    ylabel('Spearman \rho');
    saveas(f, strcat(resultsDirCurrentFigure, 'emotionid_', currentContrast, '_controlInputDiff_drug_vs_noDrug_PETatlases.eps'));
    
end

%% average control input difference w/ and w/o drug vs PET neurotransmitter maps - emotionrec

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlInputs_persistence; % extracting control input
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    drug = allControlTrajectories_emotionrec_currentContrast.drug;
    group = allControlTrajectories_emotionrec_currentContrast.group;
    
    % populate matrix of control input for [nSubjects*2 x nNodes]
    controlInputs_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = trapz(controlInput_emotionrec{i}.^2)/nTimeSteps; % calculating total control input for current node as the sum of integral of time-varying control input
        controlInputs_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    % control input w/ alpraz and placebo
    controlInput_emotionrec_currentContrast_alpraz = controlInputs_emotionrec_allSubjects(drug==0, :);
    controlInput_emotionrec_currentContrast_placebo = controlInputs_emotionrec_allSubjects(drug==1, :);
    
    controlInputDiff_emotionrec_currentContrast = abs(controlInput_emotionrec_currentContrast_alpraz - controlInput_emotionrec_currentContrast_placebo)'; % calculating absolute value of difference in control input w/ and w/o drug
    avge_controlInputDiff_emotionrec_currentContrast = mean(controlInputDiff_emotionrec_currentContrast, 2); % average across subjects
    
    % fit linear model
%     T = table(avge_controlInputDiff_emotionrec_currentContrast, PET_5HT1a_WAY_HC36, PET_5HT1b_P943_HC22, PET_5HT2a_ALT_HC19, PET_D1_SCH23390_c11, PET_D2_RACLOPRIDE_c11, PET_DAT_DATSPECT, PET_FDOPA_f18, PET_GABAa_FLUMAZENIL_c11, PET_NAT_MRB_c11, PET_SERT_DASB_HC30);
%     linearModel_currentContrast = fitlm(T, 'ResponseVar', 'avge_controlInputDiff_emotionrec_currentContrast');
%     spinTestSavePath = strcat(resultsDirCurrentFigure, 'linearModel_emotionrec_', currentContrast, '.mat');
%     save(spinTestSavePath, 'linearModel_currentContrast');
    
    % find Lausanne nodes with high difference in control input b/w drug
    % and placebo
    controlInputDiff_emotionrec_nodeNames = cell(nNodes, 2);
    for i = 1:nNodes
        controlInputDiff_emotionrec_nodeNames{i, 1} = avge_controlInputDiff_emotionrec_currentContrast(i);
        controlInputDiff_emotionrec_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    controlInputDiff_emotionrec_nodeNames = cell2table(controlInputDiff_emotionrec_nodeNames, 'VariableNames', {'controlInputDiff_drug_placebo', 'nodeName_Lausanne'});
    controlInputDiff_emotionrec_nodeNames = sortrows(controlInputDiff_emotionrec_nodeNames, 'controlInputDiff_drug_placebo', 'descend', 'MissingPlacement', 'last'); % sorting nodes by descending value
    writetable(controlInputDiff_emotionrec_nodeNames, strcat(resultsDirCurrentFigure, 'controlInputDiff_drug_placebo_emotionrec_', currentContrast, '.csv'));
    
    % plot surface maps of differences in control input b/w alpraz and
    % placebo, averaged across subjects
    avge_controlInputDiff_emotionrec_currentContrast(isnan(avge_controlInputDiff_emotionrec_currentContrast)) = 0;
    surfacePlots([avge_controlInputDiff_emotionrec_currentContrast; 0], redbluecmap, [0 0.015], nifti, subcorticalIndices, resultsDirCurrentFigure, strcat('emotionrec_', currentContrast, 'controlInputDiff_heatmap_drug_vs_noDrug.svg'));
    
    % compute correlations between control input differences w/ and w/o drug against PET maps
    controlInputDiff_emotionrec_currentContrast_PETatlasCorrs = corr(controlInputDiff_emotionrec_currentContrast, PET_maps, 'Type', 'Spearman', 'Rows', 'Complete');
    avge_controlInputDiff_emotionrec_currentContrast_PETatlasCorrs = mean(controlInputDiff_emotionrec_currentContrast_PETatlasCorrs); % average correlation across subjects
    
    % compute correlations between control input difference against PET
    % maps, after randomizing placebo data or PET maps
    nRandomizations = 500;
    rng(0); % reset random number generator
    randomCorrs = zeros(nRandomizations, numel(PETlabels));
    
%     for i = 1:nRandomizations
%         PET_maps_random = [PET_5HT1a_WAY_HC36_random(:, i) PET_5HT1b_P943_HC22_random(:, i) ...
%             PET_5HT2a_ALT_HC19_random(:, i) PET_D1_SCH23390_c11_random(:, i) PET_D2_RACLOPRIDE_c11_random(:, i) ...
%             PET_DAT_DATSPECT_random(:, i) PET_FDOPA_f18_random(:, i) PET_GABAa_FLUMAZENIL_c11_random(:, i) ...
%             PET_NAT_MRB_c11_random(:, i) PET_SERT_DASB_HC30_random(:, i)]; % compile random maps for i'th iteration
%         PET_maps_random(idx_brainstem, :) = [];
%         controlInputDiff_emotionrec_currentContrast_PETatlasCorrs_random = corr(controlInputDiff_emotionrec_currentContrast, PET_maps_random, 'Type', 'Spearman', 'Rows', 'Complete');
%         randomCorrs(i, :) = mean(controlInputDiff_emotionrec_currentContrast_PETatlasCorrs_random);
%     end
    
    for i = 1:nRandomizations        
        PET_maps_random = randomizeMatrix(PET_maps')';
        controlInputDiff_emotionrec_currentContrast_PETatlasCorrs_random = corr(controlInputDiff_emotionrec_currentContrast, PET_maps_random, 'Type', 'Spearman', 'Rows', 'Complete');
        randomCorrs(i, :) = mean(controlInputDiff_emotionrec_currentContrast_PETatlasCorrs_random);
    end
    
    pValues = zeros(1, numel(PETlabels));
    for i = 1:numel(PETlabels)
        pValues(i) = sum(abs(randomCorrs(:, i)) > abs(avge_controlInputDiff_emotionrec_currentContrast_PETatlasCorrs(i)));
    end
    pValues = (pValues/nRandomizations);
    pValues_FDR = mafdr(pValues, 'BHFDR', true);
    
    %controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs = atanh(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % Fisher z-transform
    %[h, pValues] = ttest(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % using one-sample t-test to check for significance of subject-wise correlations
    %[pValues, ~, ~, ~, ~] = mult_comp_perm_t1(controlInputDiff_emotionid_currentContrast_PETatlasCorrs);
    h1 = (pValues_FDR < 0.05) & (pValues_FDR > 0.005);
    h2 = (pValues_FDR < 0.005) & (pValues_FDR > 0.0005);
    h3 = (pValues_FDR < 0.0005);
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    ylim([-0.4, 0.4]);
    
    boxplot(controlInputDiff_emotionrec_currentContrast_PETatlasCorrs, 'Labels', PETlabels, 'LabelOrientation', 'inline', 'OutlierSize', 8, 'Symbol', 'k.');
    refline(0, 0);
    for i = 1:numel(h1)
        if h1(i)
            text(i, 0, '*', 'Color', 'r', 'FontSize', 20);
        elseif h2(i)
            text(i, 0, '**', 'Color', 'r', 'FontSize', 20);
        elseif h3(i)
            text(i, 0, '***', 'Color', 'r', 'FontSize', 20);
        end
    end
    
    ylabel('Spearman \rho');
    saveas(f, strcat(resultsDirCurrentFigure, 'emotionrec_', currentContrast, '_controlInputDiff_drug_vs_noDrug_PETatlases.eps'));
    
end
end