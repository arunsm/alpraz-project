function [] = printFigure2(allControlEnergies_emotionid, allControlEnergies_emotionrec, ...
    allControlTrajectories_emotionid, allControlTrajectories_emotionrec)

parameters
resultsDirCurrentFigure = strcat(resultsDir, filesep, 'Figure2', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% plot persistence energy +/- Alpraz for controls and relatives; note that drug(0,1)=(alpraz,placebo)

% emotion id
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

savePath = strcat(resultsDirCurrentFigure, 'boxPlot_persistenceEnergy_emotionID.svg');
saveas(f, savePath);
close(f);

% emotion memory
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

savePath = strcat(resultsDirCurrentFigure, 'boxPlot_persistenceEnergy_emotionrec.svg');
saveas(f, savePath);
close(f);

%% subject-average control impact maps

% emotionid
avgeControlImpact_emotionid_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    group = allControlTrajectories_emotionid_currentContrast.group;
    controlImpact_emotionid = allControlTrajectories_emotionid_currentContrast.controlImpact_persistence; % extracting control impact
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlImpact_emotionid);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlImpact_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        current_controlImpact_emotionid = controlImpact_emotionid{i};
        controlImpact_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlImpact_emotionid;
    end
    
    avgeControlImpact_emotionid_currentContrast_allSubjects = mean(controlImpact_emotionid_allSubjects); % averaging over all subjects
    avgeControlImpact_emotionid_allSubjects(:, c) = avgeControlImpact_emotionid_currentContrast_allSubjects;
    avgeControlImpact_emotionid_allSubjects = round(avgeControlImpact_emotionid_allSubjects, 2); % rounding to 2 decimal places
    
    emotionid_controlImpact_nodeNames = cell(nNodes, 2);
    for i = 1:nNodes
        emotionid_controlImpact_nodeNames{i, 1} = avgeControlImpact_emotionid_allSubjects(i, c);
        emotionid_controlImpact_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    emotionid_controlImpact_nodeNames = cell2table(emotionid_controlImpact_nodeNames, 'VariableNames', {'controlImpact', 'nodeName_Lausanne'});
    emotionid_controlImpact_nodeNames = sortrows(emotionid_controlImpact_nodeNames, 'controlImpact', 'descend', 'MissingPlacement', 'last'); % sorting nodes by descending value
    writetable(emotionid_controlImpact_nodeNames, strcat(resultsDirCurrentFigure, 'emotionid_controlImpact_nodeNames_', currentContrast, '.csv'));
end

avgeControlImpact_emotionid_allSubjects(isnan(avgeControlImpact_emotionid_allSubjects)) = 0; % setting NaN values (outside slab) to 0 for visualization
cmap = redbluecmap;
%cmap = cmap(6:end, :);
C = avgeControlImpact_emotionid_allSubjects(:, 1);
C(234, 1) = 0;
surfacePlots(C, cmap, [min(C) - 0.1*range(C), max(C) + 0.1*range(C)], nifti, ...
    subcorticalIndices, resultsDirCurrentFigure, 'avgeControlImpact_emotionid_threat');

% emotionrec
avgeControlImpact_emotionrec_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlImpact_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlImpact_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlImpact_emotionrec);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlImpact_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlImpact_emotionrec = controlImpact_emotionrec{i};
        controlImpact_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlImpact_emotionrec;
    end
    
    avgeControlImpact_emotionrec_currentContrast_allSubjects = mean(controlImpact_emotionrec_allSubjects); % averaging over all subjects
    avgeControlImpact_emotionrec_allSubjects(:, c) = avgeControlImpact_emotionrec_currentContrast_allSubjects;
    avgeControlImpact_emotionrec_allSubjects = round(avgeControlImpact_emotionrec_allSubjects, 2); % rounding to 2 decimal places
    
    emotionrec_controlImpact_nodeNames = cell(nNodes, 2);
    for i = 1:nNodes
        emotionrec_controlImpact_nodeNames{i, 1} = avgeControlImpact_emotionrec_allSubjects(i, c);
        emotionrec_controlImpact_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    emotionrec_controlImpact_nodeNames = cell2table(emotionrec_controlImpact_nodeNames, 'VariableNames', {'controlImpact', 'nodeName_Lausanne'});
    emotionrec_controlImpact_nodeNames = sortrows(emotionrec_controlImpact_nodeNames, 'controlImpact', 'descend', 'MissingPlacement', 'last'); % sorting nodes by descending value
    writetable(emotionrec_controlImpact_nodeNames, strcat(resultsDirCurrentFigure, 'emotionrec_controlImpact_nodeNames_', currentContrast, '.csv'));
end

avgeControlImpact_emotionrec_allSubjects(isnan(avgeControlImpact_emotionrec_allSubjects)) = 0; % setting NaN values (outside slab) to 0
cmap = redbluecmap;
%cmap = cmap(6:end, :);
C = avgeControlImpact_emotionrec_allSubjects(:, 1);
C(234, 1) = 0;
surfacePlots(C, cmap, [min(C) - 0.1*range(C), max(C) + 0.1*range(C)], nifti, ...
    subcorticalIndices, resultsDirCurrentFigure, 'avgeControlImpact_emotionrec_threat');
end
