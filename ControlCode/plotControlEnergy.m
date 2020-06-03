%% loading relevant files

clear all;

addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/beeswarm-master'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/LausanneSurfaceFigures'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/LausanneCoordinates'))

resultsDir = 'Results/minimalControl/avge_FD_thresh_0.5_parcelCoverageThresh_0.5_EuclideanNorm_allNodes_QA/';

allControlEnergies_emotionid = readtable(strcat(resultsDir, 'allControlEnergies_emotionid.csv'));
allControlEnergies_emotionrec = readtable(strcat(resultsDir, 'allControlEnergies_emotionrec.csv'));
load(strcat(resultsDir, 'allControlTrajectories_emotionid.mat'));
load(strcat(resultsDir, 'allControlTrajectories_emotionrec.mat'));
load(strcat(resultsDir, 'structuralAdjacencyMatrix.mat'));

load LindenYeoPurity/yeo7netlabelsLaus125EJC.mat;
subSystemLabels = {'Visual', 'Somatomator', 'DorsalAttention', 'VentralAttention', 'Limbic', 'FrontoparietalControl', 'DefaultMode', 'Subcortical'};
subcorticalIndices = find(finalLabels == 8);
nifti = load_nii('ROIv_scale125_dilated.nii.gz');

X = importdata('../../lausanne2008/LausanneParcelNames.xlsx');
LausanneParcelNames = X.textdata;

nSubSystems = numel(subSystemLabels);
contrastLabels = {'contrast1_threatcorrectStd', 'contrast3_nonthreatcorrectStd', ...
    'contrast5_neutralcorrectStd'};
nContrasts = numel(contrastLabels);

parcelCoverageThresh = 0.5;
nNodes = 234;

%% Table 1 - clinical and demographic variables of cohort

pathToDemographics = '../../Alpraz_subjectDemographics.xlsx';
subjectDemographics = readtable(pathToDemographics);

group = subjectDemographics.Group;

fprintf('percent female controls = %.1f\n', 100*sum(subjectDemographics.sex_M0F1(group==0)/sum(group==0)))
fprintf('proportion female controls: %iF/%iM\n', sum(subjectDemographics.sex_M0F1(group==0)), sum(group==0)-sum(subjectDemographics.sex_M0F1(group==0)))
fprintf('percent female relatives = %.1f\n', 100*sum(subjectDemographics.sex_M0F1(group==1)/sum(group==1)))
fprintf('proportion female controls: %iF/%iM\n', sum(subjectDemographics.sex_M0F1(group==1)), sum(group==1)-sum(subjectDemographics.sex_M0F1(group==1)))
[~, pValue] = ttest2(subjectDemographics.sex_M0F1(group==0), subjectDemographics.sex_M0F1(group==1)) % need to implement different statistical test

fprintf('percent right-handed = %.1f\n', 100-(100*sum(subjectDemographics.hand_R0L1(group==0)/sum(group==0))))
fprintf('proportion right-handed: %iR/%iL\n', sum(group==0)-sum(subjectDemographics.hand_R0L1(group==0)), sum(subjectDemographics.hand_R0L1(group==0)))
fprintf('percent right-handed = %.1f\n', 100-(100*sum(subjectDemographics.hand_R0L1(group==1)/sum(group==1))))
fprintf('proportion right-handed: %iR/%iL\n', sum(group==1)-sum(subjectDemographics.hand_R0L1(group==1)), sum(subjectDemographics.hand_R0L1(group==1)))
[~, pValue] = ttest2(subjectDemographics.hand_R0L1(group==0), subjectDemographics.hand_R0L1(group==1))

fprintf('percent non-smokers = %.1f\n', 100*sum(subjectDemographics.smoke_Y0N1(group==0)/sum(group==0)))
fprintf('proportion smoke: %iN/%iY\n', sum(subjectDemographics.smoke_Y0N1(group==0)), sum(group==0)-sum(subjectDemographics.smoke_Y0N1(group==0)))
fprintf('percent non-smokers = %.1f\n', 100*sum(subjectDemographics.smoke_Y0N1(group==1)/sum(group==1)))
fprintf('proportion smoke: %iN/%iY\n', sum(subjectDemographics.smoke_Y0N1(group==1)), sum(group==1)-sum(subjectDemographics.smoke_Y0N1(group==1)))
[~, pValue] = ttest2(subjectDemographics.smoke_Y0N1(group==0), subjectDemographics.smoke_Y0N1(group==1))

fprintf('mean age (std) = %.1f(%.1f)\n', mean(subjectDemographics.AgeAtFMRI(group==0)), std(subjectDemographics.AgeAtFMRI(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.AgeAtFMRI(group==0)), max(subjectDemographics.AgeAtFMRI(group==0)))
fprintf('mean age (std) = %.1f(%.1f)\n', mean(subjectDemographics.AgeAtFMRI(group==1)), std(subjectDemographics.AgeAtFMRI(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.AgeAtFMRI(group==1)), max(subjectDemographics.AgeAtFMRI(group==1)))
pValue = ranksum(subjectDemographics.AgeAtFMRI(group==0), subjectDemographics.AgeAtFMRI(group==1))

fprintf('mean education (std) = %.1f(%.1f)\n', mean(subjectDemographics.educ(group==0)), std(subjectDemographics.educ(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.educ(group==0)), max(subjectDemographics.educ(group==0)))
fprintf('mean education (std) = %.1f(%.1f)\n', mean(subjectDemographics.educ(group==1)), std(subjectDemographics.educ(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.educ(group==1)), max(subjectDemographics.educ(group==1)))
[~, pValue] = ttest2(subjectDemographics.educ(group==0), subjectDemographics.educ(group==1))

fprintf('mean parental education (std) = %.1f(%.1f)\n', mean(subjectDemographics.par_ed(group==0), 'omitnan'), std(subjectDemographics.par_ed(group==0), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.par_ed(group==0)), max(subjectDemographics.par_ed(group==0)))
fprintf('mean parental education (std) = %.1f(%.1f)\n', mean(subjectDemographics.par_ed(group==1), 'omitnan'), std(subjectDemographics.par_ed(group==1), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.par_ed(group==1)), max(subjectDemographics.par_ed(group==1)))
[~, pValue] = ttest2(subjectDemographics.par_ed(group==0), subjectDemographics.par_ed(group==1))

fprintf('mean height = %.1f(%.1f)\n', mean(subjectDemographics.height(group==0)), std(subjectDemographics.height(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.height(group==0)), max(subjectDemographics.height(group==0)))
fprintf('mean height (std) = %.1f(%.1f)\n', mean(subjectDemographics.height(group==1)), std(subjectDemographics.height(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.height(group==1)), max(subjectDemographics.height(group==1)))
[~, pValue] = ttest2(subjectDemographics.height(group==0), subjectDemographics.height(group==1))

fprintf('mean weight = %.1f(%.1f)\n', mean(subjectDemographics.weight(group==0)), std(subjectDemographics.weight(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.weight(group==0)), max(subjectDemographics.weight(group==0)))
fprintf('mean weight (std) = %.1f(%.1f)\n', mean(subjectDemographics.weight(group==1)), std(subjectDemographics.weight(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.weight(group==1)), max(subjectDemographics.weight(group==1)))
[~, pValue] = ttest2(subjectDemographics.weight(group==0), subjectDemographics.weight(group==1))

fprintf('mean BMI = %.1f(%.1f)\n', mean(subjectDemographics.BMI(group==0)), std(subjectDemographics.BMI(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.BMI(group==0)), max(subjectDemographics.BMI(group==0)))
fprintf('mean BMI (std) = %.1f(%.1f)\n', mean(subjectDemographics.BMI(group==1)), std(subjectDemographics.BMI(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.BMI(group==1)), max(subjectDemographics.BMI(group==1)))
[~, pValue] = ttest2(subjectDemographics.BMI(group==0), subjectDemographics.BMI(group==1))

fprintf('mean SISTOTAL = %.1f(%.1f)\n', mean(subjectDemographics.SISTOTAL(group==0), 'omitnan'), std(subjectDemographics.SISTOTAL(group==0), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.SISTOTAL(group==0)), max(subjectDemographics.SISTOTAL(group==0)))
fprintf('mean SISTOTAL (std) = %.1f(%.1f)\n', mean(subjectDemographics.SISTOTAL(group==1), 'omitnan'), std(subjectDemographics.SISTOTAL(group==1), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.SISTOTAL(group==1)), max(subjectDemographics.SISTOTAL(group==1)))
[~, pValue] = ttest2(subjectDemographics.SISTOTAL(group==0), subjectDemographics.SISTOTAL(group==1))

%% Figure 1

% plot structural matrix

figure; set(gcf, 'color', 'w');
imagesc(A);
colormap('parula');
colorbar;
axis off;
ax = gca;
ax.FontSize = 12;

% Figure depicting Lausanne parcellation

parcelValues = zeros(234, 1);
[~, ax, ph, ~] = fcn_lausannesurf(parcelValues, white, [-20 100]);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');

% plot brain activation maps w/ and w/o alpraz for individual subject

M = importdata('/Users/ArunMahadevan/Documents/BBL/studies/alpraz/derivatives/xcp_output_allTasks2/sub-010707/ses-001989/task-emotionid/norm/sub-010707_ses-001989_task-emotionid_contrast1_threatcorrectStd_lausanne_ROIv_scale125_dilated.txt');
brainStates_alpraz = M.data';
[~, ax, ph, ~] = fcn_lausannesurf(brainStates_alpraz, hot, [-20 100]);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');

%% Supplementary Figure 1 - plot average brain mask coverage
M = importdata('/Users/ArunMahadevan/Documents/BBL/studies/alpraz/slabCoverage_combinedMaskStd_emotionid.txt');
slabCoverage = M.data';
slabCoverage = double(slabCoverage > parcelCoverageThresh);
[~, ax, ph, ~] = fcn_lausannesurf(slabCoverage, parula);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');

%% Figure 2

%% creating mixed models to examine effects of clinical and demographic variables on persistence_allNodes during emotionid

contrast = allControlEnergies_emotionid.contrast;

for i = 1:nContrasts
    currentContrast = contrastLabels{i}
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(contrast, currentContrast), :);
    
    % converting categorical variables
    allControlEnergies_emotionid_currentContrast.group = categorical(allControlEnergies_emotionid_currentContrast.group);
    allControlEnergies_emotionid_currentContrast.drug = categorical(allControlEnergies_emotionid_currentContrast.drug);
    allControlEnergies_emotionid_currentContrast.gender = categorical(allControlEnergies_emotionid_currentContrast.gender);
    
    % fitting model after checking for normality
    persistence_allNodes = allControlEnergies_emotionid_currentContrast.persistence_allNodes; 
    h = kstest(zscore(persistence_allNodes));
    figure, histogram(persistence_allNodes); title(strcat('KS test: ', num2str(h)));
    mixedModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'persistence_allNodes ~ drug + group + drug*group + STAI_TRAIT*drug + SISTOTAL + gender + age + (1|subjectID)', 'FitMethod', 'REML')
    save(strcat(resultsDir, 'mixedModel_emotionid_', currentContrast, '.mat'), 'mixedModel_emotionid_currentContrast');
end

%% creating mixed models to examine effects of clinical and demographic variables on persistence_allNodes during emotionrec

contrast = allControlEnergies_emotionrec.contrast;

for i = 1:nContrasts
    currentContrast = contrastLabels{i}
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(contrast, currentContrast), :);
    
    % converting categorical variables
    allControlEnergies_emotionrec_currentContrast.group = categorical(allControlEnergies_emotionrec_currentContrast.group);
    allControlEnergies_emotionrec_currentContrast.drug = categorical(allControlEnergies_emotionrec_currentContrast.drug);
    allControlEnergies_emotionrec_currentContrast.gender = categorical(allControlEnergies_emotionrec_currentContrast.gender);
    
    % fitting model
    persistence_allNodes = allControlEnergies_emotionrec_currentContrast.persistence_allNodes;
    h = kstest(zscore(persistence_allNodes));
    figure, histogram(persistence_allNodes); title(strcat('KS test: ', num2str(h)));
    mixedModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'persistence_allNodes ~ drug + group + drug*group + STAI_TRAIT*drug + SISTOTAL + gender + age + (1|subjectID)', 'FitMethod', 'REML')
    save(strcat(resultsDir, 'mixedModel_emotionrec_', currentContrast, '.mat'), 'mixedModel_emotionrec_currentContrast');
end

%% plot persistence during emotion id +/- Alpraz for controls and relatives; note that drug(0,1)=(alpraz,placebo)

drug = allControlEnergies_emotionid.drug;
group = allControlEnergies_emotionid.group;
contrast = allControlEnergies_emotionid.contrast;

f = figure; set(gcf, 'color', 'w');
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 6 5];

x = [allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd'))];

groups = [ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    2*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    3*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    4*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    5*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    6*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    7*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    8*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    9*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    10*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    11*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    12*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd'))))];
    
positions = [1, 1.2, 1.4, 1.6, 2, 2.2, 2.4, 2.6, 3, 3.2, 3.4, 3.6];

boxplot(x, groups, 'positions', positions, 'OutlierSize', 8, 'Symbol', 'k.');
set(gca, 'XTickLabel', {'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A'});

%set(gca,'xtick',[])
color = repmat({[0.4660, 0.6740, 0.1880], [0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560], [0.4940, 0.1840, 0.5560]}, 1, 3); % right to left
faceAlpha = repmat([0.75, 0.25], 1, 6); % right to left
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   currentFaceAlpha = faceAlpha(j); 
   patch(get(h(j), 'XData'), get(h(j), 'YData'), color{j}, 'FaceAlpha', currentFaceAlpha);
end

ylabel('persistence energy')
ylim([0.2, 0.65]);
set(gca, 'FontSize', 20);

%% plot persistence during emotion rec +/- Alpraz for controls and relatives; note that drug(0,1)=(alpraz,placebo)

drug = allControlEnergies_emotionrec.drug;
group = allControlEnergies_emotionrec.group;
contrast = allControlEnergies_emotionrec.contrast;

f = figure; set(gcf, 'color', 'w');
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 6 5];

x = [allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd'))];

groups = [ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    2*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    3*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    4*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    5*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    6*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    7*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    8*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    9*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    10*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    11*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    12*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd'))))];
    
positions = [1, 1.2, 1.4, 1.6, 2, 2.2, 2.4, 2.6, 3, 3.2, 3.4, 3.6];

boxplot(x, groups, 'positions', positions, 'OutlierSize', 8, 'Symbol', 'k.');
set(gca, 'XTickLabel', {'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A'});

%set(gca,'xtick',[])
color = repmat({[0.4660, 0.6740, 0.1880], [0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560], [0.4940, 0.1840, 0.5560]}, 1, 3); % right to left
faceAlpha = repmat([0.75, 0.25], 1, 6); % right to left
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   currentFaceAlpha = faceAlpha(j); 
   patch(get(h(j), 'XData'), get(h(j), 'YData'), color{j}, 'FaceAlpha', currentFaceAlpha);
end

ylabel('persistence energy')
ylim([0.2, 0.65]);
set(gca, 'FontSize', 20);

%% calculating and plotting average control impact across all subjects

%% emotionid
avgeControlImpact_emotionid_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlImpact_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionid = controlInput_emotionid{i};
        controlInput_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionid;
    end
    
    avgeControlImpact_emotionid_currentContrast_allSubjects = mean(controlInput_emotionid_allSubjects); % averaging over all subjects
    avgeControlImpact_emotionid_currentContrast_allSubjects(isnan(avgeControlImpact_emotionid_currentContrast_allSubjects)) = 0; % setting NaN values (outside slab) to 0
    avgeControlImpact_emotionid_allSubjects(:, c) = avgeControlImpact_emotionid_currentContrast_allSubjects;
end

avgeControlImpact_emotionid_allSubjects = mean(avgeControlImpact_emotionid_allSubjects, 2); % averaging over all contrasts

emotionid_controlImpact_nodeNames = cell(234, 2);
for i = 1:234
    emotionid_controlImpact_nodeNames{i, 1} = avgeControlImpact_emotionid_allSubjects(i);
    emotionid_controlImpact_nodeNames{i, 2} = LausanneParcelNames{i};
end

emotionid_controlImpact_nodeNames = sortrows(emotionid_controlImpact_nodeNames, 'descend'); % sorting nodes by descending control impact value
emotionid_controlImpact_nodeNames = cell2table(emotionid_controlImpact_nodeNames, 'VariableNames', {'controlImpact_sorted', 'nodeName_Lausanne'});
writetable(emotionid_controlImpact_nodeNames, strcat(resultsDir, 'emotionid_controlImpact_nodeNames.csv'));

[~, ax, ph, ~] = fcn_lausannesurf(avgeControlImpact_emotionid_allSubjects, hot, [0 1.5]);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
figure; [~,f_lat] = plot_subcortvol(avgeControlImpact_emotionid_allSubjects(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, hot, 0, 1.5);
close(f_lat);
[~, ax, ph, ~] = fcn_lausannesurf(avgeControlImpact_emotionid_allSubjects, hot, [0 1.5]);
view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [270, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
figure; [~,f_lat] = plot_subcortvol(avgeControlImpact_emotionid_allSubjects(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, hot, 0, 1.5);
close(f_lat);

%% emotionrec

avgeControlImpact_emotionrec_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlImpact_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = controlInput_emotionrec{i};
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    avgeControlImpact_emotionrec_currentContrast_allSubjects = mean(controlInput_emotionrec_allSubjects);
    avgeControlImpact_emotionrec_currentContrast_allSubjects(isnan(avgeControlImpact_emotionrec_currentContrast_allSubjects)) = 0;
    
    avgeControlImpact_emotionrec_allSubjects(:, c) = avgeControlImpact_emotionrec_currentContrast_allSubjects;
end

avgeControlImpact_emotionrec_allSubjects = mean(avgeControlImpact_emotionrec_allSubjects, 2);

emotionrec_controlImpact_nodeNames = cell(234, 2);
for i = 1:234
    emotionrec_controlImpact_nodeNames{i, 1} = avgeControlImpact_emotionrec_allSubjects(i);
    emotionrec_controlImpact_nodeNames{i, 2} = LausanneParcelNames{i};
end

emotionrec_controlImpact_nodeNames = sortrows(emotionrec_controlImpact_nodeNames, 'descend');
emotionrec_controlImpact_nodeNames = cell2table(emotionrec_controlImpact_nodeNames, 'VariableNames', {'controlImpact_sorted', 'nodeName_Lausanne'});
writetable(emotionrec_controlImpact_nodeNames, strcat(resultsDir, 'emotionrec_controlImpact_nodeNames.csv'));

[~, ax, ph, ~] = fcn_lausannesurf(avgeControlImpact_emotionrec_allSubjects, hot, [0 1.5]);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
figure; [~,f_lat] = plot_subcortvol(avgeControlImpact_emotionrec_allSubjects(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, hot, 0, 1.5);
close(f_lat);
[~, ax, ph, ~] = fcn_lausannesurf(avgeControlImpact_emotionrec_allSubjects, hot, [0 1.5]);
view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [270, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
figure; [~,f_lat] = plot_subcortvol(avgeControlImpact_emotionrec_allSubjects(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, hot, 0, 1.5);
close(f_lat);

%% Figure 3

%% plotting SISTOTAL and STAI_TRAIT against persistence_allNodes for nonthreat emotionid

contrast = allControlEnergies_emotionid.contrast;
allControlEnergies_emotionid_nonthreat = allControlEnergies_emotionid(strcmp(contrast, 'contrast3_nonthreatcorrectStd'), :);
persistence_allNodes_nonthreat = allControlEnergies_emotionid_nonthreat.persistence_allNodes;
allControlEnergies_emotionid_neutral = allControlEnergies_emotionid(strcmp(contrast, 'contrast5_neutralcorrectStd'), :);
persistence_allNodes_neutral = allControlEnergies_emotionid_neutral.persistence_allNodes;
SISTOTAL = allControlEnergies_emotionid_nonthreat.SISTOTAL;
STAI_TRAIT = allControlEnergies_emotionid_nonthreat.STAI_TRAIT;

figure; set(gcf, 'color', 'white');
plot(SISTOTAL, persistence_allNodes_nonthreat, 'k.', 'MarkerSize', 20);
h = lsline; h.LineWidth = 2;
xlabel('SISTOTAL');
ylabel('persistence energy - nonthreat');
set(gca, 'FontSize', 20);
[rho, pValue] = corr(SISTOTAL, persistence_allNodes_nonthreat, 'Rows', 'Complete')
ylim([0.2, 0.65]);

drug = allControlEnergies_emotionid_nonthreat.drug;

figure; set(gcf, 'color', 'w'); set(gca, 'FontSize', 20); hold on;
plot(STAI_TRAIT(drug==1), persistence_allNodes_neutral(drug==1), 'r.', 'MarkerSize', 20);
plot(STAI_TRAIT(drug==0), persistence_allNodes_neutral(drug==0), 'b.', 'MarkerSize', 20);
lsline;
h = lsline; h(1).LineWidth = 2; h(2).LineWidth = 2; 
xlabel('trait anxiety');
ylabel('persistence energy - neutral');
legend('placebo', 'alpraz', 'location', 'northeast'); legend boxoff;
[rho, pValue] = corr(STAI_TRAIT(drug==0), persistence_allNodes_neutral(drug==0), 'Rows', 'Complete')
[rho, pValue] = corr(STAI_TRAIT(drug==1), persistence_allNodes_neutral(drug==1), 'Rows', 'Complete')
ylim([0.2, 0.65]);

%% mixed models for accuracy and reaction time versus persistence for emotionid and emotionrec

for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(allControlEnergies_emotionid.contrast, currentContrast), :);
    
    allControlEnergies_emotionid_currentContrast.drug = categorical(allControlEnergies_emotionid_currentContrast.drug);
    allControlEnergies_emotionid_currentContrast.group = categorical(allControlEnergies_emotionid_currentContrast.group);
    allControlEnergies_emotionid_currentContrast.gender = categorical(allControlEnergies_emotionid_currentContrast.gender);
    
    switch currentContrast
        case 'contrast1_threatcorrectStd'
            figure, histogram(allControlEnergies_emotionid_currentContrast.pctcorr_threat); title('emotionid - pctcorr_threat');
            figure, histogram(allControlEnergies_emotionid_currentContrast.rtmdn_threatcorr); title('emotionid - rtmdn_threatcorr');
            accuracyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'pctcorr_threat ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
            reactionTimeModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'rtmdn_threatcorr ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
        case 'contrast3_nonthreatcorrectStd'
            figure, histogram(allControlEnergies_emotionid_currentContrast.pctcorr_nonthreat); title('emotionid - pctcorr_nonthreat');
            figure, histogram(allControlEnergies_emotionid_currentContrast.rtmdn_nonthreatcorr); title('emotionid - rtmdn_nonthreatcorr');
            accuracyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'pctcorr_nonthreat ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
            reactionTimeModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'rtmdn_nonthreatcorr ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
        case 'contrast5_neutralcorrectStd'
           figure, histogram(allControlEnergies_emotionid_currentContrast.pctcorr_neutral); title('emotionid - pctcorr_neutral');
           figure, histogram(allControlEnergies_emotionid_currentContrast.rtmdn_neutralcorr); title('emotionid - rtmdn_neutralcorr');
           accuracyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'pctcorr_neutral ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
           reactionTimeModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'rtmdn_neutralcorr ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
    end
    
    save(strcat(resultsDir, 'accuracyModel_emotionid_', currentContrast, '.mat'), 'accuracyModel_emotionid_currentContrast');
    save(strcat(resultsDir, 'reactionTimeModel_emotionid_', currentContrast, '.mat'), 'reactionTimeModel_emotionid_currentContrast');
end

for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(allControlEnergies_emotionrec.contrast, currentContrast), :);
    
    allControlEnergies_emotionrec_currentContrast.drug = categorical(allControlEnergies_emotionrec_currentContrast.drug);
    allControlEnergies_emotionrec_currentContrast.group = categorical(allControlEnergies_emotionrec_currentContrast.group);
    allControlEnergies_emotionrec_currentContrast.gender = categorical(allControlEnergies_emotionrec_currentContrast.gender);
    
    switch currentContrast
        case 'contrast1_threatcorrectStd'
            figure, histogram(allControlEnergies_emotionrec_currentContrast.pctcorr_threat); title('emotionrec - pctcorr_threat');
            figure, histogram(allControlEnergies_emotionrec_currentContrast.rtmdn_threatcorr); title('emotionrec - rtmdn_threatcorr');
            accuracyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'pctcorr_threat ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
            reactionTimeModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'rtmdn_threatcorr ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
        case 'contrast3_nonthreatcorrectStd'
            figure, histogram(allControlEnergies_emotionrec_currentContrast.pctcorr_nonthreat); title('emotionrec - pctcorr_nonthreat');
            figure, histogram(allControlEnergies_emotionrec_currentContrast.rtmdn_nonthreatcorr); title('emotionrec - rtmdn_nonthreatcorr');
            accuracyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'pctcorr_nonthreat ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
            reactionTimeModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'rtmdn_nonthreatcorr ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
        case 'contrast5_neutralcorrectStd'
           figure, histogram(allControlEnergies_emotionrec_currentContrast.pctcorr_neutral); title('emotionrec - pctcorr_neutral');
            figure, histogram(allControlEnergies_emotionrec_currentContrast.rtmdn_neutralcorr); title('emotionrec - rtmdn_neutralcorr');
           accuracyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'pctcorr_neutral ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
           reactionTimeModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'rtmdn_neutralcorr ~ persistence_allNodes + drug + group + gender + age + (1|subjectID)', 'FitMethod', 'REML');
    end
    
    save(strcat(resultsDir, 'accuracyModel_emotionrec_', currentContrast, '.mat'), 'accuracyModel_emotionrec_currentContrast');
    save(strcat(resultsDir, 'reactionTimeModel_emotionrec_', currentContrast, '.mat'), 'reactionTimeModel_emotionrec_currentContrast');
end

%% plot task accuracy versus persistence for threat emotionid

allControlEnergies_emotionid_threat = allControlEnergies_emotionid(strcmp(allControlEnergies_emotionid.contrast, 'contrast1_threatcorrectStd'), :);

accuracy = allControlEnergies_emotionid_threat.pctcorr_threat;
persistence = allControlEnergies_emotionid_threat.persistence_allNodes;

drug = allControlEnergies_emotionid_threat.drug;

figure; set(gcf, 'color', 'w'); hold on;
set(gca, 'FontSize', 20);
plot(accuracy(drug==1), persistence(drug==1), 'r.', 'MarkerSize', 20);
plot(accuracy(drug==0), persistence(drug==0), 'b.', 'MarkerSize', 20);
h = lsline; h(1).LineWidth = 2; h(2).LineWidth = 2;
xlabel('accuracy - threat emotion ID');
ylabel('persistence energy');
ylim([0.2, 0.65])
legend('placebo', 'alpraz', 'location', 'northwest'); legend boxoff;

%% plot task accuracy versus persistence for threat emotionrec

allControlEnergies_emotionrec_threat = allControlEnergies_emotionrec(strcmp(allControlEnergies_emotionrec.contrast, 'contrast1_threatcorrectStd'), :);

accuracy = allControlEnergies_emotionrec_threat.pctcorr_threat;
persistence = allControlEnergies_emotionrec_threat.persistence_allNodes;

drug = allControlEnergies_emotionrec_threat.drug;

figure; set(gcf, 'color', 'w'); hold on;
set(gca, 'FontSize', 20);
plot(accuracy(drug==1), persistence(drug==1), 'r.', 'MarkerSize', 20);
plot(accuracy(drug==0), persistence(drug==0), 'b.', 'MarkerSize', 20);
lsline
xlabel('accuracy - threat emotionrec');
ylabel('persistence');
legend('placebo', 'alpraz', 'location', 'northoutside'); legend boxoff;

%% plot PET neurotransmitter profiles

PETlabels = {'5HT1a_WAY_HC36', '5HT1b_P943_HC22', '5HT2a_ALT_HC19', 'D1_SCH23390_c11', 'D2_RACLOPRIDE_c11', 'DAT_DATSPECT.nii', 'FDOPA_f18', 'GABAa_FLUMAZENIL_c11', 'NAT_MRB_c11', 'SERT_DASB_HC30'};

% load PET data
PETdir = '/Users/ArunMahadevan/Documents/BBL/studies/alpraz/PETatlas/';
X = importdata(strcat(PETdir, '5HT1a_WAY_HC36.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT1a_WAY_HC36 = X.data';
X = importdata(strcat(PETdir, '5HT1b_P943_HC22.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT1b_P943_HC22 = X.data';
X = importdata(strcat(PETdir, '5HT2a_ALT_HC19.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT2a_ALT_HC19 = X.data';
X = importdata(strcat(PETdir, 'D1_SCH23390_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_D1_SCH23390_c11 = X.data';
X = importdata(strcat(PETdir, 'D2_RACLOPRIDE_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_D2_RACLOPRIDE_c11 = X.data';
X = importdata(strcat(PETdir, 'DAT_DATSPECT.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_DAT_DATSPECT = X.data';
X = importdata(strcat(PETdir, 'FDOPA_f18.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_FDOPA_f18 = X.data';
X = importdata(strcat(PETdir, 'GABAa_FLUMAZENIL_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_GABAa_FLUMAZENIL_c11 = X.data';
X = importdata(strcat(PETdir, 'NAT_MRB_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_NAT_MRB_c11 = X.data';
X = importdata(strcat(PETdir, 'SERT_DASB_HC30.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_SERT_DASB_HC30 = X.data';

[f, ax, ph, ~] = fcn_lausannesurf(PET_GABAa_FLUMAZENIL_c11, hot, [10 80]); close(f(2)); % plotting only right hemisphere
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
axis(ax(1), 'off'); colorbar;
[g, ax, ph, ~] = fcn_lausannesurf(PET_GABAa_FLUMAZENIL_c11, hot, [10 80]); close(g(2)); % plotting only right hemisphere
view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
axis(ax(1), 'off'); colorbar;

[f, ax, ph, ~] = fcn_lausannesurf(PET_GABAa_FLUMAZENIL_c11, hot, [10 80]); close(f(1)); % plotting only left hemisphere
view(ax(2), [-90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(2), 'off'); colorbar;
[g, ax, ph, ~] = fcn_lausannesurf(PET_GABAa_FLUMAZENIL_c11, hot, [10 80]); close(g(1)); % plotting only left hemisphere
view(ax(2), [-270, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(2), 'off'); colorbar;

figure;
[~,f_lat] = plot_subcortvol(PET_GABAa_FLUMAZENIL_c11(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, hot, 10, 80); % visualizing sub-cortical areas
close(f_lat); colorbar;

%% average control input difference w/ and w/o drug vs PET neurotransmitter maps - emotionid

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlInputs_persistence; % extracting control impact
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    drug = allControlTrajectories_emotionid_currentContrast.drug;
    group = allControlTrajectories_emotionid_currentContrast.group;
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInputs_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        current_controlInput_emotionid = trapz(controlInput_emotionid{i}.^2)/nTimeSteps; % calculating total control input for current node as the sum of integral of time-varying control input
        controlInputs_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionid;
    end
    
    controlInputDiff_emotionid_currentContrast_drug = abs(controlInputs_emotionid_allSubjects(drug==0, :) - controlInputs_emotionid_allSubjects(drug==1, :)); % calculating absolute value of difference in control input w/ and w/o drug
        
    % compute correlations between control impact differences w/ and w/o drug against PET
    % maps
    controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs = [corr(controlInputDiff_emotionid_currentContrast_drug', PET_5HT1a_WAY_HC36, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_5HT1b_P943_HC22, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_5HT2a_ALT_HC19, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_D1_SCH23390_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_D2_RACLOPRIDE_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_DAT_DATSPECT, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_FDOPA_f18, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_NAT_MRB_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug', PET_SERT_DASB_HC30, 'Type', 'Spearman', 'Rows', 'Complete')];
    
    %controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs = atanh(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % Fisher z-transform
    [h, pValues] = ttest(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % using one-sample t-test to check for significance of subject-wise correlations
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    boxplot(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs, 'Labels', PETlabels, 'LabelOrientation', 'inline');
    refline(0, 0);
    for i = 1:numel(h)
        if h(i) == 1
            text(i, 0, '*', 'FontSize', 12);
        end
    end
    
    ylabel('Fisher z (Spearman \rho)');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputDiffDrug_PETatlases.svg'));
    
end

%% average control input w/ and w/o drug vs PET neurotransmitter maps - emotionrec

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlInputs_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    drug = allControlTrajectories_emotionrec_currentContrast.drug;
    group = allControlTrajectories_emotionrec_currentContrast.group;
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = trapz(controlInput_emotionrec{i}.^2)/nTimeSteps; % calculating total control input for current node as the sum of integral of time-varying control input
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    controlInputDiff_emotionrec_currentContrast_drug = abs(controlInput_emotionrec_allSubjects(drug==0, :) - controlInput_emotionrec_allSubjects(drug==1, :));
        
    % compute correlations between control impact differences w/ and w/o drug against PET
    % maps
    controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs = [corr(controlInputDiff_emotionrec_currentContrast_drug', PET_5HT1a_WAY_HC36, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_5HT1b_P943_HC22, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_5HT2a_ALT_HC19, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_D1_SCH23390_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_D2_RACLOPRIDE_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_DAT_DATSPECT, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_FDOPA_f18, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_NAT_MRB_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug', PET_SERT_DASB_HC30, 'Type', 'Spearman', 'Rows', 'Complete')];
    
    controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs = atanh(controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs); % Fisher z-transform
    [h, pValues] = ttest(controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs); % using one-sample t-test to check for significance of subject-wise correlations
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    boxplot(controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs, 'Labels', PETlabels, 'LabelOrientation', 'inline');
    refline(0, 0);
    for i = 1:numel(h)
        if h(i) == 1
            text(i, 0, '*', 'FontSize', 12);
        end
    end
    
    ylabel('Fisher z (Spearman \rho)');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDiffDrug_PETatlases.svg'));
    
end

%% control input beta coefficients vs PET neurotransmitter maps - emotionid

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlInputs_persistence; % extracting control impact
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        current_controlInput_emotionid = trapz(controlInput_emotionid{i}.^2)/nTimeSteps; % calculating total control input for current node as the sum of integral of time-varying control input
        controlInput_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionid;
    end
    
    % evaluate effects of drug and group on control impact of each node
    betas_groupMain = NaN(nNodes, 1); pValues_groupMain = NaN(nNodes, 1);
    betas_drugMain = NaN(nNodes, 1); pValues_drugMain = NaN(nNodes, 1);
    betas_groupDrugInteraction = NaN(nNodes, 1); pValues_groupDrugInteraction = NaN(nNodes, 1);
    
    for i = 1:nNodes
        %fprintf('node%d\n', i);
        controlInput_currentNode = controlInput_emotionid_allSubjects(:, i);
        controlInput_currentNode =  zscore(controlInput_currentNode);
        if sum(isnan(controlInput_currentNode)) > 0 % how much missing data am I willing to tolerate?
            continue;
        else
            allControlTrajectories_currentNode = addvars(allControlTrajectories_emotionid_currentContrast, controlInput_currentNode); % creating table with additional variable
            controlInputModel = fitlme(allControlTrajectories_currentNode, 'controlInput_currentNode ~ drug + group + drug*group + (1|subjectID)', 'FitMethod', 'REML'); % fitting mixed model
            coefficientNames = controlInputModel.CoefficientNames; coefficients = controlInputModel.Coefficients.Estimate; pValues = controlInputModel.Coefficients.pValue; % extracting coefficients and p-values
            betas_groupMain(i) = coefficients(strcmp(coefficientNames, 'group')); betas_drugMain(i) = coefficients(strcmp(coefficientNames, 'drug')); betas_groupDrugInteraction(i) = coefficients(strcmp(coefficientNames, 'group:drug'));
            pValues_groupMain(i) = pValues(strcmp(coefficientNames, 'group')); pValues_drugMain(i) = pValues(strcmp(coefficientNames, 'drug')); pValues_groupDrugInteraction(i) = pValues(strcmp(coefficientNames, 'group:drug'));
        end
    end
    
    betas_groupMain = abs(betas_groupMain); betas_drugMain = abs(betas_drugMain); betas_groupDrugInteraction = abs(betas_groupDrugInteraction); % taking absolute values of coefficients
    
    pValues_groupMain_FDR = mafdr(pValues_groupMain, 'BHFDR', true); 
    pValues_drugMain_FDR = mafdr(pValues_drugMain, 'BHFDR', true);
    pValues_groupDrugInteraction_FDR = mafdr(pValues_groupDrugInteraction, 'BHFDR', true);
    
    sum(pValues_groupMain_FDR < 0.05/234)
    sum(pValues_drugMain_FDR < 0.05/234)
    sum(pValues_groupDrugInteraction_FDR < 0.05/234)
    
    % plotting coefficients of group main effect versus PET
    % neurotransmitter profiles
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);

    xlim([0, 1]);
    ylim([10, 80]);
    
    plot(betas_groupMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.1, 78, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
        
    lsline;
    xlabel('coefficients of group main effect');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputGroupEffect_GABA.svg'));
    
    % plotting coefficients of drug main effect versus PET
    % neurotransmitter profiles
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([0, 1]);
    ylim([10, 80]);
    
    plot(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.1, 78, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of drug main effect');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputDrugEffect_GABA.svg'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([0, 1]);
    ylim([10, 80]);
    
    plot(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.1, 78, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group x drug interaction');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputGroupDrugInteraction_GABA.svg'));
end

%% control impact vs PET neurotransmitter maps - emotionrec

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
        current_controlInput_emotionrec = trapz(controlInput_emotionrec{i}.^2)/nTimeSteps; % calculating total control input for current node as the sum of integral of time-varying control input
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    % evaluate effects of drug and group on control impact of each node
    betas_groupMain = NaN(nNodes, 1); pValues_groupMain = NaN(nNodes, 1);
    betas_drugMain = NaN(nNodes, 1); pValues_drugMain = NaN(nNodes, 1);
    betas_groupDrugInteraction = NaN(nNodes, 1); pValues_groupDrugInteraction = NaN(nNodes, 1);
    
    for i = 1:nNodes
        %fprintf('node%d\n', i);
        controlInput_currentNode = controlInput_emotionrec_allSubjects(:, i);
        controlInput_currentNode =  zscore(controlInput_currentNode);
        if sum(isnan(controlInput_currentNode)) > 0 % how much missing data am I willing to tolerate?
            continue;
        else
            allControlTrajectories_currentNode = addvars(allControlTrajectories_emotionrec_currentContrast, controlInput_currentNode); % creating table with additional variable
            controlInputModel = fitlme(allControlTrajectories_currentNode, 'controlInput_currentNode ~ drug + group + drug*group + (1|subjectID)', 'FitMethod', 'REML'); % fitting mixed model
            coefficientNames = controlInputModel.CoefficientNames; coefficients = controlInputModel.Coefficients.Estimate; pValues = controlInputModel.Coefficients.pValue; % extracting coefficients and p-values
            betas_groupMain(i) = coefficients(strcmp(coefficientNames, 'group')); betas_drugMain(i) = coefficients(strcmp(coefficientNames, 'drug')); betas_groupDrugInteraction(i) = coefficients(strcmp(coefficientNames, 'group:drug'));
            pValues_groupMain(i) = pValues(strcmp(coefficientNames, 'group')); pValues_drugMain(i) = pValues(strcmp(coefficientNames, 'drug')); pValues_groupDrugInteraction(i) = pValues(strcmp(coefficientNames, 'group:drug'));
        end
    end
    
    betas_groupMain = abs(betas_groupMain); betas_drugMain = abs(betas_drugMain); betas_groupDrugInteraction = abs(betas_groupDrugInteraction); % taking absolute values of coefficients
    
    pValues_groupMain_FDR = mafdr(pValues_groupMain, 'BHFDR', true); 
    pValues_drugMain_FDR = mafdr(pValues_drugMain, 'BHFDR', true);
    pValues_groupDrugInteraction_FDR = mafdr(pValues_groupDrugInteraction, 'BHFDR', true);
    
    sum(pValues_groupMain_FDR < 0.05/234)
    sum(pValues_drugMain_FDR < 0.05/234)
    sum(pValues_groupDrugInteraction_FDR < 0.05/234)
    
    % plotting coefficients of group main effect versus PET
    % neurotransmitter profiles
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);

    xlim([0, 1]);
    ylim([10, 80]);
    
    plot(betas_groupMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.1, 78, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
        
    lsline;
    xlabel('coefficients of group main effect');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupEffect_GABA.svg'));
    
    % plotting coefficients of drug main effect versus PET
    % neurotransmitter profiles
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([0, 1]);
    ylim([10, 80]);
    
    plot(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.1, 78, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of drug main effect');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDrugEffect_GABA.svg'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([0, 1]);
    ylim([10, 80]);
    
    plot(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.1, 78, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group x drug interaction');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupDrugInteraction_GABA.svg'));
end

%% plot GABA(A) receptor expression from Allen atlas

% calculate GABA gene expression from Allen human brain atlas data
load('/Users/ArunMahadevan/Documents/BBL/studies/alpraz/GeneExpression/ParcellatedGeneExpressionLausanne125.mat');
GABRA1 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRA1'));
GABRA2 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRA2'));
GABRA3 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRA3'));
GABRA5 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRA5'));
GABRB1 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRB1'));
GABRB2 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRB2'));
GABRB3 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRB3'));
GABRG2 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRG2'));

% % calculating expression of type 1 and type 2 receptors as minimum of
% % associated subunit expression
% GABAA_type1 = min([GABRA1, GABRB1, GABRB2, GABRB3, GABRG2], [], 2); %GABAA_type1(idx_brainstem) = [];
% GABAA_type2 = min([GABRA2, GABRA3, GABRA5, GABRB1, GABRB2, GABRB3, GABRG2], [], 2); %GABAA_type2(idx_brainstem) = [];
% 
% %GABAA_type1(GABAA_type1 < prctile(GABAA_type1, 80) | isnan(GABAA_type1)) = 0; % keeping only top 20th percentile and non-NaN values
% %GABAA_type2(GABAA_type1 < prctile(GABAA_type2, 80) | isnan(GABAA_type2)) = 0; % keeping only top 20th percentile and non-NaN values
% 
% [f, ax, ph, ~] = fcn_lausannesurf(GABAA_type1, hot, [0 0.58]); close(f(2)); % plotting only right hemisphere
% view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
% axis(ax(1), 'off'); colorbar;
% [g, ax, ph, ~] = fcn_lausannesurf(GABAA_type1, hot, [0 0.58]); close(g(2)); % plotting only right hemisphere
% view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
% axis(ax(1), 'off'); colorbar;
% figure;
% [f_ant,f_lat] = plot_subcortvol(GABAA_type1(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, hot, 0, 0.58); % visualizing sub-cortical areas
% close(f_lat); colorbar;
% 
% [f, ax, ph, ~] = fcn_lausannesurf(GABAA_type2, hot, [0 0.58]); close(f(2)); % plotting only right hemisphere
% view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
% axis(ax(1), 'off'); colorbar;
% [g, ax, ph, ~] = fcn_lausannesurf(GABAA_type2, hot, [0 0.58]); close(g(2)); % plotting only right hemisphere
% view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
% axis(ax(1), 'off'); colorbar;
% figure;
% [f_ant,f_lat] = plot_subcortvol(GABAA_type2(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, hot, 0, 0.58); % visualizing sub-cortical areas
% close(f_lat); colorbar;

%% control impact vs GABA expression calculations

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlImpact_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = controlInput_emotionrec{i};
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    % evaluate effects of drug and group on control impact of each node
    betas_groupMain = NaN(nNodes, 1); pValues_groupMain = NaN(nNodes, 1);
    betas_drugMain = NaN(nNodes, 1); pValues_drugMain = NaN(nNodes, 1);
    betas_groupDrugInteraction = NaN(nNodes, 1); pValues_groupDrugInteraction = NaN(nNodes, 1);
    
    for i = 1:nNodes
        %fprintf('node%d\n', i);
        controlInput_currentNode = controlInput_emotionrec_allSubjects(:, i);
        controlInput_currentNode =  zscore(controlInput_currentNode);
        if sum(isnan(controlInput_currentNode)) > 0 % how much missing data am I willing to tolerate?
            continue;
        else
            allControlTrajectories_currentNode = addvars(allControlTrajectories_emotionrec_currentContrast, controlInput_currentNode); % creating table with additional variable
            controlInputModel = fitlme(allControlTrajectories_currentNode, 'controlImpact_currentNode ~ drug + group + drug*group + (1|subjectID)', 'FitMethod', 'REML'); % fitting mixed model
            coefficientNames = controlInputModel.CoefficientNames; coefficients = controlInputModel.Coefficients.Estimate; pValues = controlInputModel.Coefficients.pValue; % extracting coefficients and p-values
            betas_groupMain(i) = coefficients(strcmp(coefficientNames, 'group')); betas_drugMain(i) = coefficients(strcmp(coefficientNames, 'drug')); betas_groupDrugInteraction(i) = coefficients(strcmp(coefficientNames, 'group:drug'));
            pValues_groupMain(i) = pValues(strcmp(coefficientNames, 'group')); pValues_drugMain(i) = pValues(strcmp(coefficientNames, 'drug')); pValues_groupDrugInteraction(i) = pValues(strcmp(coefficientNames, 'group:drug'));
        end
    end
    
    pValues_groupMain_FDR = mafdr(pValues_groupMain, 'BHFDR', true); 
    pValues_drugMain_FDR = mafdr(pValues_drugMain, 'BHFDR', true);
    pValues_groupDrugInteraction_FDR = mafdr(pValues_groupDrugInteraction, 'BHFDR', true);
    
    sum(pValues_groupMain_FDR < 0.05/234)
    sum(pValues_drugMain_FDR < 0.05/234)
    sum(pValues_groupDrugInteraction_FDR < 0.05/234)
    
    % plotting coefficients of group main effect versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([-1, 1]);
    ylim([0, 1]);
    
    plot(betas_groupMain, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupMain, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_groupMain, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupMain, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
    
    plot(betas_groupMain, GABRA3, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupMain, GABRA3, 'Rows', 'complete');
    text(-0.9, 0.92, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group main effect');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', '\alpha3', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlImpactGroupEffect_GABA.svg'));
    
    % plotting coefficients of drug main effect versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([-1, 1]);
    ylim([0, 1]);
    
    plot(betas_drugMain, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_drugMain, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
    
    plot(betas_drugMain, GABRA3, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA3, 'Rows', 'complete');
    text(-0.9, 0.92, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of drug main effect');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', '\alpha3', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlImpactDrugEffect_GABA.svg'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([-1, 1]);
    ylim([0, 1]);
    
    plot(betas_groupDrugInteraction, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_groupDrugInteraction, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
    
    plot(betas_groupDrugInteraction, GABRA3, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA3, 'Rows', 'complete');
    text(-0.9, 0.92, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group x drug interaction');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', '\alpha3', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlImpactGroupDrugInteraction_GABA.svg'));
end

%% control input vs GABA calculations

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    fprintf(currentContrast); fprintf('\n');
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlInputs_persistence; % extracting control input
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    
    % populate matrix of control inout for [nSubjects*2 x nNodes]
    controlInput_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = trapz(controlInput_emotionrec{i}.^2)/nTimeSteps; % calculating control input for each node as the sum of integral of squared energy trajectories, divided by number of time steps
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    % evaluate effects of drug and group on control impact of each node
    betas_groupMain = NaN(nNodes, 1); pValues_groupMain = NaN(nNodes, 1);
    betas_drugMain = NaN(nNodes, 1); pValues_drugMain = NaN(nNodes, 1);
    betas_groupDrugInteraction = NaN(nNodes, 1); pValues_groupDrugInteraction = NaN(nNodes, 1);
    
    for i = 1:nNodes
        %fprintf('node%d\n', i);
        controlInput_currentNode = controlInput_emotionrec_allSubjects(:, i);
        controlInput_currentNode =  zscore(controlInput_currentNode);
        if sum(isnan(controlInput_currentNode)) > 0 % how much missing data am I willing to tolerate?
            continue;
        else
            allControlTrajectories_currentNode = addvars(allControlTrajectories_emotionrec_currentContrast, controlInput_currentNode); % creating table with additional variable
            controlInputModel = fitlme(allControlTrajectories_currentNode, 'controlInput_currentNode ~ drug + group + drug*group + (1|subjectID)', 'FitMethod', 'REML'); % fitting mixed model
            coefficientNames = controlInputModel.CoefficientNames; coefficients = controlInputModel.Coefficients.Estimate; pValues = controlInputModel.Coefficients.pValue; % extracting coefficients and p-values
            betas_groupMain(i) = coefficients(strcmp(coefficientNames, 'group')); betas_drugMain(i) = coefficients(strcmp(coefficientNames, 'drug')); betas_groupDrugInteraction(i) = coefficients(strcmp(coefficientNames, 'group:drug'));
            pValues_groupMain(i) = pValues(strcmp(coefficientNames, 'group')); pValues_drugMain(i) = pValues(strcmp(coefficientNames, 'drug')); pValues_groupDrugInteraction(i) = pValues(strcmp(coefficientNames, 'group:drug'));
        end
    end
    
    pValues_groupMain_FDR = mafdr(pValues_groupMain, 'BHFDR', true); 
    pValues_drugMain_FDR = mafdr(pValues_drugMain, 'BHFDR', true);
    pValues_groupDrugInteraction_FDR = mafdr(pValues_groupDrugInteraction, 'BHFDR', true);
    
    sum(pValues_groupMain_FDR < 0.05/234)
    sum(pValues_drugMain_FDR < 0.05/234)
    sum(pValues_groupDrugInteraction_FDR < 0.05/234)
    
    % plotting coefficients of drug main effect versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_drugMain, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_drugMain, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
    
    plot(betas_drugMain, GABRA3, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA3, 'Rows', 'complete');
    text(-0.9, 0.92, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of drug main effect');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', '\alpha3', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDrugEffect_GABA.svg'));
    
    % plotting coefficients of group main effect versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_groupMain, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupMain, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_groupMain, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupMain, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
    
    plot(betas_groupMain, GABRA3, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupMain, GABRA3, 'Rows', 'complete');
    text(-0.9, 0.92, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group main effect');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', '\alpha3', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupEffect_GABA.svg'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    
    xlim([-1, 1]);
    ylim([0, 1]);
    
    plot(betas_groupDrugInteraction, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_groupDrugInteraction, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
    
    plot(betas_groupDrugInteraction, GABRA3, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA3, 'Rows', 'complete');
    text(-0.9, 0.92, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group x drug interaction');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', '\alpha3', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupDrugInteraction_GABA.svg'));
end

%% structural null model analysis

nullResultsDir = 'Results/minimalControl/avge_FD_thresh_0.5_parcelCoverageThresh_0.5_EuclideanNorm_allNodes_nullModel/';

allControlEnergies_emotionid_null = readtable(strcat(nullResultsDir, 'allControlEnergies_emotionid.csv'));
allControlEnergies_emotionrec_null = readtable(strcat(nullResultsDir, 'allControlEnergies_emotionrec.csv'));

contrast = allControlEnergies_emotionid.contrast;

for i = 1:nContrasts
    currentContrast = contrastLabels{i}
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(contrast, currentContrast), :);
    allControlEnergies_emotionid_currentContrast_null = allControlEnergies_emotionid_null(strcmp(contrast, currentContrast), :);
    
    drug = allControlEnergies_emotionid_currentContrast.drug;
    group = allControlEnergies_emotionid_currentContrast.group;
    
    [~, p, ci, stats] = ttest(allControlEnergies_emotionid_currentContrast.persistence_allNodes(drug==1), allControlEnergies_emotionid_currentContrast_null.persistence_allNodes(drug==1))
    figure, boxplot([allControlEnergies_emotionid_currentContrast.persistence_allNodes(drug==1), allControlEnergies_emotionid_currentContrast_null.persistence_allNodes(drug==1)]);
end

%% plot control trajectories

%%%%%% emotionid %%%%%%
% re-ordering nodewise control trajectories by Yeo7 sub-network
idx = [find(finalLabels==1); find(finalLabels==2); find(finalLabels==3); find(finalLabels==4); find(finalLabels==5); find(finalLabels==6); find(finalLabels==7); find(finalLabels==8)];
controlInputs_stability_emotionid_drug = controlInputs_stability_emotionid_drug(:, idx, :); controlInputs_stability_emotionid_noDrug = controlInputs_stability_emotionid_noDrug(:, idx, :);
controlInputs_stability_emotionrec_drug = controlInputs_stability_emotionrec_drug(:, idx, :); controlInputs_stability_emotionrec_noDrug = controlInputs_stability_emotionrec_noDrug(:, idx, :);

% plotting correlation between energy cost with drug and GABA(A) receptor expression
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
plot(GABAA_type1(idx_type1), avge_energy_cost_emotionid_drug(idx_type1), 'r.', 'MarkerSize', 12);
[rho, pVal] = corr(GABAA_type1, avge_energy_cost_emotionid_drug, 'Rows', 'complete');
text(0.6, max(avge_energy_cost_emotionid_drug), strcat('type 1, \rho=', num2str(rho), ', p=', num2str(pVal)), 'Color', 'r');
plot(GABAA_type2(idx_type2), avge_energy_cost_emotionid_drug(idx_type2), 'b.', 'MarkerSize', 12);
[rho, pVal] = corr(GABAA_type2, avge_energy_cost_emotionid_drug, 'Rows', 'complete');
text(0.6, max(avge_energy_cost_emotionid_drug)-0.1, strcat('type 2, \rho=', num2str(rho), ', p=', num2str(pVal)), 'Color', 'b');
xlabel('GABA(A) receptor expression');
ylabel('control energy');
ax = gca;
ax.FontSize = 20;
xlim([0 1]);
currentSaveFileName = strcat(resultsDir, 'controlEnergy_GABA_emotionid_', currentContrastLabel, '_drug.svg');
saveas(f, currentSaveFileName);
close(f);

% plotting control trajectories with drug
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
imagesc(mean(U_opt_stability_emotionid_drug, 3)'); cb = colorbar; cb.Location = 'westoutside';
xlabel('Time (a.u.)');
ylabel('Control input by region');
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
xlim([0, nTimeSteps]); ylim([0, nNodes]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(7), refLines(8)]), 'Subcortical');
ax = gca;
ax.FontSize = 20;
currentSaveFileName = strcat(resultsDir, 'stateTrajectories_emotionid_drug_', currentContrastLabel, '.svg');
saveas(f, currentSaveFileName);
close(f);

% plotting control trajectories without drug
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
imagesc(mean(U_opt_stability_emotionid_noDrug, 3)'); cb = colorbar; cb.Location = 'westoutside';
xlabel('Time (a.u.)');
ylabel('Control input by region');
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(1001, mean([0, refLines(1)]), 'Visual');
xlim([0, nTimeSteps]); ylim([0, nNodes]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(7), refLines(8)]), 'Subcortical');
ax = gca;
ax.FontSize = 20;
currentSaveFileName = strcat(resultsDir, 'stateTrajectories_emotionid_noDrug_', currentContrastLabel, '.svg');
saveas(f, currentSaveFileName);
close(f);

%%%%%% emotionrec %%%%%%

% plotting correlation between energy cost and GABA(A) receptor expression
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
plot(GABAA_type1(idx_type1), avge_energy_cost_emotionrec_drug(idx_type1), 'r.', 'MarkerSize', 12);
[rho, pVal] = corr(GABAA_type1, avge_energy_cost_emotionrec_drug, 'Rows', 'complete');
text(0.6, max(avge_energy_cost_emotionrec_drug), strcat('type 1, \rho=', num2str(rho), ', p=', num2str(pVal)), 'Color', 'r');
plot(GABAA_type2(idx_type2), avge_energy_cost_emotionrec_drug(idx_type2), 'b.', 'MarkerSize', 12);
[rho, pVal] = corr(GABAA_type2, avge_energy_cost_emotionrec_drug, 'Rows', 'complete');
text(0.6, max(avge_energy_cost_emotionrec_drug)-0.1, strcat('type 2, \rho=', num2str(rho), ', p=', num2str(pVal)), 'Color', 'b');
xlabel('GABA(A) receptor expression');
ylabel('control energy');
ax = gca;
ax.FontSize = 20;
xlim([0 1]);
currentSaveFileName = strcat(resultsDir, 'controlEnergy_GABA_emotionrec_', currentContrastLabel, '_drug.svg');
saveas(f, currentSaveFileName);
close(f);

% plotting control trajectories with drug
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
imagesc(mean(U_opt_stability_emotionrec_drug, 3)'); cb = colorbar; cb.Location = 'westoutside';
xlabel('Time (a.u.)');
ylabel('Control input by region');
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
xlim([0, nTimeSteps]); ylim([0, nNodes]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(7), refLines(8)]), 'Subcortical');
ax = gca;
ax.FontSize = 20;
currentSaveFileName = strcat(resultsDir, 'stateTrajectories_emotionrec_drug_', currentContrastLabel, '.svg');
saveas(f, currentSaveFileName);
close(f);

% plotting control trajectories without drug
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
imagesc(mean(U_opt_stability_emotionrec_noDrug, 3)'); cb = colorbar; cb.Location = 'westoutside';
xlabel('Time (a.u.)');
ylabel('Control input by region');
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
xlim([0, nTimeSteps]); ylim([0, nNodes]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(7), refLines(8)]), 'Subcortical');
ax = gca;
ax.FontSize = 20;
currentSaveFileName = strcat(resultsDir, 'stateTrajectories_emotionrec_noDrug_', currentContrastLabel, '.svg');
saveas(f, currentSaveFileName);
close(f);