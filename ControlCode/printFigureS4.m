function [] = printFigureS4(allControlTrajectories_emotionid, allControlTrajectories_emotionrec)

parameters
resultsDirCurrentFigure = strcat(resultsDir, filesep, 'FigureS4', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% subject-average control input maps

% emotionid
avgeControlInput_emotionid_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    group = allControlTrajectories_emotionid_currentContrast.group;
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlInputs_persistence; % extracting betas
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    
    % populate matrix of betas for [nSubjects*2 x nNodes]
    controlInput_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        currentControlInput_emotionid = trapz(controlInput_emotionid{i}.^2)/nTimeSteps;
        controlInput_emotionid_allSubjects(i, current_parcelsToInclude_idx) = currentControlInput_emotionid;
    end
    
    avgeControlInput_emotionid_currentContrast_allSubjects = mean(controlInput_emotionid_allSubjects); % averaging over all subjects
    avgeControlInput_emotionid_allSubjects(:, c) = avgeControlInput_emotionid_currentContrast_allSubjects;
    %avgeControlInput_emotionid_allSubjects = round(avgeControlInput_emotionid_allSubjects, 2); % rounding to 2 decimal places
    
    emotionid_controlInput_nodeNames = cell(nNodes, 2);
    for i = 1:nNodes
        emotionid_controlInput_nodeNames{i, 1} = avgeControlInput_emotionid_allSubjects(i, c);
        emotionid_controlInput_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    emotionid_controlInput_nodeNames = cell2table(emotionid_controlInput_nodeNames, 'VariableNames', {'avgeControlInput', 'nodeName_Lausanne'});
    emotionid_controlInput_nodeNames = sortrows(emotionid_controlInput_nodeNames, 'avgeControlInput', 'descend', 'MissingPlacement',  'last'); % sorting nodes by descending value
    writetable(emotionid_controlInput_nodeNames, strcat(resultsDirCurrentFigure, 'emotionid_controlInput_nodeNames_', currentContrast, '.csv'));
end

avgeControlInput_emotionid_allSubjects(isnan(avgeControlInput_emotionid_allSubjects)) = -2000; % setting NaN values (outside slab) to large negative value for plotting
surfacePlots(avgeControlInput_emotionid_allSubjects(:, 1), redbluecmap, [0 0.02], nifti, ...
    subcorticalIndices, resultsDirCurrentFigure, 'avgeControlInput_emotionid_threat');

% emotionrec
avgeControlInput_emotionrec_allSubjects = zeros(nNodes, nContrasts);

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
        currentControlInput_emotionrec = trapz(controlInput_emotionrec{i}.^2)/nTimeSteps;
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = currentControlInput_emotionrec;
    end
    
    avgeControlInput_emotionrec_currentContrast_allSubjects = mean(controlInput_emotionrec_allSubjects); % averaging over all subjects
    avgeControlInput_emotionrec_allSubjects(:, c) = avgeControlInput_emotionrec_currentContrast_allSubjects;
    %avgeControlInput_emotionrec_allSubjects = round(avgeControlInput_emotionrec_allSubjects, 2); % rounding to 2 decimal places
    
    emotionrec_controlInput_nodeNames = cell(nNodes, 2);
    for i = 1:nNodes
        emotionrec_controlInput_nodeNames{i, 1} = avgeControlInput_emotionrec_allSubjects(i, c);
        emotionrec_controlInput_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    emotionrec_controlInput_nodeNames = cell2table(emotionrec_controlInput_nodeNames, 'VariableNames', {'avgeControlInput', 'nodeName_Lausanne'});
    emotionrec_controlInput_nodeNames = sortrows(emotionrec_controlInput_nodeNames, 'avgeControlInput', 'descend', 'MissingPlacement',  'last'); % sorting nodes by descending value
    writetable(emotionrec_controlInput_nodeNames, strcat(resultsDirCurrentFigure, 'emotionrec_controlInput_nodeNames_', currentContrast, '.csv'));
end

avgeControlInput_emotionrec_allSubjects(isnan(avgeControlInput_emotionrec_allSubjects)) = -2000; % setting NaN values (outside slab) to 0
surfacePlots(avgeControlInput_emotionrec_allSubjects(:, 1), redbluecmap, [0 0.02], nifti, ...
    subcorticalIndices, resultsDirCurrentFigure, 'avgeControlInput_emotionrec_threat');
end