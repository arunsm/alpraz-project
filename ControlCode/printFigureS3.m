function [] = printFigureS3(allControlTrajectories_emotionid, allControlTrajectories_emotionrec)

parameters
resultsDirCurrentFigure = strcat(resultsDir, filesep, 'FigureS3', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% subject-average beta maps

% emotionid
avgeBetas_emotionid_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    group = allControlTrajectories_emotionid_currentContrast.group;
    betas_emotionid = allControlTrajectories_emotionid_currentContrast.xf; % extracting betas
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(betas_emotionid);
    
    % populate matrix of betas for [nSubjects*2 x nNodes]
    betas_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        currentBetas_emotionid = betas_emotionid{i};
        betas_emotionid_allSubjects(i, current_parcelsToInclude_idx) = currentBetas_emotionid;
    end
    
    avgeBetas_emotionid_currentContrast_allSubjects = mean(betas_emotionid_allSubjects); % averaging over all subjects
    avgeBetas_emotionid_allSubjects(:, c) = avgeBetas_emotionid_currentContrast_allSubjects;
    avgeBetas_emotionid_allSubjects = round(avgeBetas_emotionid_allSubjects, 2); % rounding to 2 decimal places
    
    emotionid_betas_nodeNames = cell(nNodes, 2);
    for i = 1:nNodes
        emotionid_betas_nodeNames{i, 1} = avgeBetas_emotionid_allSubjects(i, c);
        emotionid_betas_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    emotionid_betas_nodeNames = cell2table(emotionid_betas_nodeNames, 'VariableNames', {'betas', 'nodeName_Lausanne'});
    emotionid_betas_nodeNames = sortrows(emotionid_betas_nodeNames, 'betas', 'descend', 'MissingPlacement',  'last'); % sorting nodes by descending value
    writetable(emotionid_betas_nodeNames, strcat(resultsDirCurrentFigure, 'emotionid_betas_nodeNames_', currentContrast, '.csv'));
end

avgeBetas_emotionid_allSubjects(isnan(avgeBetas_emotionid_allSubjects)) = -1000; % setting NaN values (outside slab) to large negative value for plotting
surfacePlots(avgeBetas_emotionid_allSubjects(:, 1), redbluecmap, [-0.2 0.2], nifti, ...
    subcorticalIndices, resultsDirCurrentFigure, 'avge_betas_emotionid_threat');

% emotionrec
avgeBetas_emotionrec_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    betas_emotionrec = allControlTrajectories_emotionrec_currentContrast.xf; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(betas_emotionrec);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    betas_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_betas_emotionrec = betas_emotionrec{i};
        betas_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_betas_emotionrec;
    end
    
    avgeBetas_emotionrec_currentContrast_allSubjects = mean(betas_emotionrec_allSubjects); % averaging over all subjects
    avgeBetas_emotionrec_allSubjects(:, c) = avgeBetas_emotionrec_currentContrast_allSubjects;
    avgeBetas_emotionrec_allSubjects = round(avgeBetas_emotionrec_allSubjects, 2); % rounding to 2 decimal places
    
    emotionrec_betas_nodeNames = cell(nNodes, 2);
    for i = 1:nNodes
        emotionrec_betas_nodeNames{i, 1} = avgeBetas_emotionrec_allSubjects(i, c);
        emotionrec_betas_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    emotionrec_betas_nodeNames = cell2table(emotionrec_betas_nodeNames, 'VariableNames', {'betas', 'nodeName_Lausanne'});
    emotionrec_betas_nodeNames = sortrows(emotionrec_betas_nodeNames, 'betas', 'descend', 'MissingPlacement',  'last'); % sorting nodes by descending value
    writetable(emotionrec_betas_nodeNames, strcat(resultsDirCurrentFigure, 'emotionrec_betas_nodeNames_', currentContrast, '.csv'));
end

avgeBetas_emotionrec_allSubjects(isnan(avgeBetas_emotionrec_allSubjects)) = -1000; % setting NaN values (outside slab) to 0
surfacePlots(avgeBetas_emotionrec_allSubjects(:, 1), redbluecmap, [-0.2 0.2], nifti, ...
    subcorticalIndices, resultsDirCurrentFigure, 'avgeBetas_emotionrec_threat');
end