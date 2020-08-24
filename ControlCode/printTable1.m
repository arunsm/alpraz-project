function [] = printTable1()
parameters
group = subjectDemographics.Group;

fprintf('percent female controls = %.1f\n', 100*sum(subjectDemographics.sex_M0F1(group==0)/sum(group==0)))
fprintf('proportion female controls: %iF/%iM\n', sum(subjectDemographics.sex_M0F1(group==0)), sum(group==0)-sum(subjectDemographics.sex_M0F1(group==0)))
fprintf('percent female relatives = %.1f\n', 100*sum(subjectDemographics.sex_M0F1(group==1)/sum(group==1)))
fprintf('proportion female controls: %iF/%iM\n', sum(subjectDemographics.sex_M0F1(group==1)), sum(group==1)-sum(subjectDemographics.sex_M0F1(group==1)))
x = table([sum(subjectDemographics.sex_M0F1(group==0)) sum(subjectDemographics.sex_M0F1(group==1)); ...
     sum(group==0)-sum(subjectDemographics.sex_M0F1(group==0)) sum(group==1)-sum(subjectDemographics.sex_M0F1(group==1))]);
[~, pValue] = fishertest(x)

fprintf('percent right-handed = %.1f\n', 100-(100*sum(subjectDemographics.hand_R0L1(group==0)/sum(group==0))))
fprintf('proportion right-handed: %iR/%iL\n', sum(group==0)-sum(subjectDemographics.hand_R0L1(group==0)), sum(subjectDemographics.hand_R0L1(group==0)))
fprintf('percent right-handed = %.1f\n', 100-(100*sum(subjectDemographics.hand_R0L1(group==1)/sum(group==1))))
fprintf('proportion right-handed: %iR/%iL\n', sum(group==1)-sum(subjectDemographics.hand_R0L1(group==1)), sum(subjectDemographics.hand_R0L1(group==1)))
x = table([sum(subjectDemographics.hand_R0L1(group==0)) sum(subjectDemographics.hand_R0L1(group==1)); ...
     sum(group==0)-sum(subjectDemographics.hand_R0L1(group==0)) sum(group==1)-sum(subjectDemographics.hand_R0L1(group==1))]);
[~, pValue] = fishertest(x)

fprintf('percent non-smokers = %.1f\n', 100*sum(subjectDemographics.smoke_Y0N1(group==0)/sum(group==0)))
fprintf('proportion smoke: %iN/%iY\n', sum(subjectDemographics.smoke_Y0N1(group==0)), sum(group==0)-sum(subjectDemographics.smoke_Y0N1(group==0)))
fprintf('percent non-smokers = %.1f\n', 100*sum(subjectDemographics.smoke_Y0N1(group==1)/sum(group==1)))
fprintf('proportion smoke: %iN/%iY\n', sum(subjectDemographics.smoke_Y0N1(group==1)), sum(group==1)-sum(subjectDemographics.smoke_Y0N1(group==1)))
x = table([sum(subjectDemographics.smoke_Y0N1(group==0)) sum(subjectDemographics.smoke_Y0N1(group==1)); ...
     sum(group==0)-sum(subjectDemographics.smoke_Y0N1(group==0)) sum(group==1)-sum(subjectDemographics.smoke_Y0N1(group==1))]);
[~, pValue] = fishertest(x)

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

% using different table for STAI_TRAIT and alpraz_levels
allSubjectInfo_emotionID = readtable('../../data/Alpraz_emotionid.xlsx');
subjectIDs = subjectDemographics.bblid;

STAI_TRAIT = zeros(numel(subjectIDs), 1);
for i = 1:numel(subjectIDs)
    currentSubjectID = subjectIDs(i)
    STAI_TRAIT(i) = mean(allSubjectInfo_emotionID.STAI_TRAIT(allSubjectInfo_emotionID.bblid == currentSubjectID)); % averaging anxiety values over sessions
end
fprintf('mean STAI_TRAIT = %.1f(%.1f)\n', mean(STAI_TRAIT(group==0)), std(STAI_TRAIT(group==0)))
fprintf('range: %.1f-%.1f\n', min(STAI_TRAIT(group==0)), max(STAI_TRAIT(group==0)))
fprintf('mean STAI_TRAIT (std) = %.1f(%.1f)\n', mean(STAI_TRAIT(group==1)), std(STAI_TRAIT(group==1)))
fprintf('range: %.1f-%.1f\n', min(STAI_TRAIT(group==1)), max(STAI_TRAIT(group==1)))
pValue = ranksum(STAI_TRAIT(group==0), STAI_TRAIT(group==1))

drug = allSubjectInfo_emotionID.drug;
alpraz_levels = allSubjectInfo_emotionID.alpraz_levels(drug == 0); % using alpraz_levels in session where drug was administered
fprintf('mean alpraz_levels = %.1f(%.1f)\n', mean(alpraz_levels(group==0), 'omitnan'), std(alpraz_levels(group==0), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(alpraz_levels(group==0)), max(alpraz_levels(group==0)))
fprintf('mean alpraz_levels (std) = %.1f(%.1f)\n', mean(alpraz_levels(group==1), 'omitnan'), std(alpraz_levels(group==1), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(alpraz_levels(group==1)), max(alpraz_levels(group==1)))
[~, pValue] = ttest2(alpraz_levels(group==0), alpraz_levels(group==1))
end