%% Script to perform sensitivity analysis on control parameters

%% set parameters

clear all;
currentNormalizationMethod = 'EuclideanNorm';
current_FD_thresh = 0.5;
currentParcelCoverageThresh = 0.5;
currentStructuralEdges = 'QA';
currentControlInputs = 'allNodes';
controlAlgorithm = 'optimalControl';

parameters = 0.1:0.01:1;
nParameters = numel(parameters);

switch currentStructuralEdges
    case 'StreamlineCounts'
        A = getAverageStreamlineCounts(); % calculating average streamline counts from Q7 data
    case 'QA'
        A = getQA(); % calculating quantitative anisotropy from Q7 data
    case 'zero'
        A = zeros(234, 234); % creating a null adjacency matrix
end

nNodes = size(A, 1); % number of nodes in parcellation
nTimeSteps = 1000; % number of time steps in simulation
T = 3;
%rho = 100; % parameter that weights energy and distance constraints, set as in Betzel et al, Sci. Rep (2016)

% loading parcel Yeo sub-network label for Lausanne-234 parcels
% 1=Visual, 2=Somatomator, 3=Dorsal Attention, 4=Ventral Attention, 5=Limbic, 6=Frontoparietal Control, 7=Default Mode, 8=Subcortical
load('LindenYeoPurity/yeo7netlabelsLaus125EJC.mat', 'finalLabels');
%finalLabels(idx_brainstem) = []; % removing brain stem parcel
subSystemLabels = {'Visual', 'Somatomator', 'DorsalAttention', 'VentralAttention', 'Limbic', 'FrontoparietalControl', 'DefaultMode', 'Subcortical'};
nSubSystems = numel(subSystemLabels);

% setting control inputs
switch currentControlInputs
    case 'allNodes'
        B = eye(size(A)); % all nodes as control inputs
    case 'cognitiveControlRegions'
        d = double(finalLabels==3 | finalLabels==4 | finalLabels==6); % dorsal and ventral attention networks + frontoparietal control network are control inputs
        d(d==0) = 0.00001; % adding jitter to non-control nodes to reduce numerical error
        B = diag(d); % creating diagonal matrix
    case 'visual_subcortex'
        d = double(finalLabels==1 | finalLabels==8); % visual and sub-cortical networks
        d(d==0) = 0.00001; % adding jitter to non-control nodes to reduce numerical error
        B = diag(d); % creating diagonal matrix
end

% creating results folder and writing parameters as text file
resultsDir = strcat('/data/joy/BBL/studies/alpraz/scripts/ControlCode/Results/', ...
    date, '_parameterAnalysis_avge_FD_thresh_', num2str(current_FD_thresh), ...
    '_parcelCoverageThresh_', num2str(currentParcelCoverageThresh), ...
    '_', currentNormalizationMethod, '_', currentControlInputs, '_', currentStructuralEdges, '/');

if ~exist(resultsDir)
    mkdir(resultsDir);
end
fid = fopen(strcat(resultsDir, 'parameterAnalysis_parameters.csv'), 'w');
fprintf(fid, 'controlAlgorithm,normalizationMethod,avge_FD_thresh,parcelCoverageThresh,T,controlInputs,structuralEdges\n%s,%s,%d,%d,%d,%s,%s', controlAlgorithm, currentNormalizationMethod, current_FD_thresh, currentParcelCoverageThresh, T, currentControlInputs, currentStructuralEdges);
fclose(fid);

%% set paths and read in subject data

rootFolder_contrasts = '/data/joy/BBL/studies/alpraz/rawData/derivatives/xcp_output_allTasks2/';
contrastLabels = {'contrast1_threatcorrectStd', 'contrast3_nonthreatcorrectStd', ...
    'contrast5_neutralcorrectStd'};
nContrasts = numel(contrastLabels);
atlasLabel = 'lausanne_ROIv_scale125_dilated';

pathToMotionFile = '/data/joy/BBL/studies/alpraz/rawData/derivatives/fmriprep/subjectMotion.csv';
pathToDemographics = '/data/joy/BBL/studies/alpraz/rawData/derivatives/stimDesigns/Alpraz_subjectDemographics.xlsx';
pathToID = '/data/joy/BBL/studies/alpraz/rawData/derivatives/stimDesigns/Alpraz_emotionid.xlsx';
pathToRec = '/data/joy/BBL/studies/alpraz/rawData/derivatives/stimDesigns/Alpraz_emotionrec.xlsx';

allSubjectInfo_emotionID = readtable(pathToID);
allSubjectInfo_emotionrec = readtable(pathToRec);
subjectDemographics = readtable(pathToDemographics);
subjectMotion = readtable(pathToMotionFile);

allSubjectIDs = subjectDemographics.bblid;
nSubjects = numel(allSubjectIDs);

%% optimal control energy calculations

taskIDs = {'task-emotionid', 'task-emotionrec'};
nTaskIDs = numel(taskIDs);

for k = 1:nTaskIDs
    currentTaskID = taskIDs{k};
    fprintf('%s\n', currentTaskID);
    
    totalEnergyCost_allParameters_threat = [];
    totalEnergyCost_allParameters_nonthreat = [];
    totalEnergyCost_allParameters_neutral = [];
    
    % iterating over all time horizons for current task
    for p = 1:nParameters
        rho = parameters(p);
        fprintf('\trho=%d\n', rho);
        
        totalEnergyCost_threat = [];
        totalEnergyCost_nonthreat = [];
        totalEnergyCost_neutral = [];
        
        for i = 1:nSubjects
            currentSubjectID = allSubjectIDs(i, 1); currentSubjectID_str = strcat('sub-0', num2str(currentSubjectID));
            
            if strcmp(currentTaskID, 'task-emotionid')
                currentSubjectInfo = allSubjectInfo_emotionID(allSubjectInfo_emotionID.bblid == currentSubjectID, :);
            else
                currentSubjectInfo = allSubjectInfo_emotionrec(allSubjectInfo_emotionrec.bblid == currentSubjectID, :);
            end
            
            allSessionIDs_currentSubject = unique(currentSubjectInfo.fmriid);
            %fprintf('\t%s\n', currentSubjectID_str);
            
            for contrast_idx = 1:nContrasts
                currentContrastLabel = contrastLabels{contrast_idx};
                %fprintf('%s\n', currentContrastLabel);
                
                for j = 1:numel(allSessionIDs_currentSubject)
                    currentSessionID = allSessionIDs_currentSubject(j); currentSessionID_str = strcat('ses-00', num2str(currentSessionID));
                    %fprintf('%s\n', currentSessionID_str);
                    
                    % check motion, skip if exceeds thresh
                    current_avge_FD = subjectMotion.avge_FD(strcmp(subjectMotion.id0, currentSubjectID_str) & strcmp(subjectMotion.id1, currentSessionID_str) & strcmp(subjectMotion.id2, currentTaskID));
                    if current_avge_FD > current_FD_thresh
                        %fprintf('Skipping subject %s, task %s: avge FD exceeds threshold in session %s\n', currentSubjectID_str, currentTaskID, currentSessionID_str);
                        break;
                    else
                        
                        pathToMask = strcat(rootFolder_contrasts, currentSubjectID_str, filesep, currentSessionID_str, filesep, currentTaskID, '/norm/', currentSubjectID_str, '_', currentSessionID_str, '_', currentTaskID, '_maskStd_', atlasLabel, '.txt');
                        M = importdata(pathToMask);
                        parcelCoverage = M.data';
                        parcelsToInclude = parcelCoverage > currentParcelCoverageThresh; % boolean list of parcels whose coverage exceeds the threshold
                        parcelsToInclude_idx = find(parcelCoverage > currentParcelCoverageThresh);
                        nNodes_currentSubject = numel(parcelsToInclude_idx);
                        
                        % truncating A and B matrices to be within imaging slab
                        A_slab = A(parcelsToInclude, parcelsToInclude);
                        B_slab = B(parcelsToInclude, parcelsToInclude);
                        finalLabels_slab = finalLabels(parcelsToInclude_idx);
                        
                        pathToContrast = strcat(rootFolder_contrasts, currentSubjectID_str, filesep, currentSessionID_str, filesep, currentTaskID, '/norm/', currentSubjectID_str, '_', currentSessionID_str, '_', currentTaskID, '_', currentContrastLabel, '_', atlasLabel, '.txt');
                        X = importdata(pathToContrast);
                        if isstruct(X)
                            xf = X.data'; % final state = contrasts for task; note that this also contains values for nodes outside imaging slab
                            %xf(idx_brainstem) = []; % removing brain stem parcel
                            xf = xf(parcelsToInclude); % only keeping nodes within imaging slab
                            
                            switch currentNormalizationMethod
                                case 'EuclideanNorm'
                                    xf = xf/norm(xf); % normalizing brain states by Euclidean norm
                                case 'none'
                                    xf = xf; % no normalization applied
                            end
                        else
                            %fprintf('data for %s, %s, %s not available\n', currentSubjectID_str, currentSessionID_str, currentTaskID);
                            break;
                        end
                        
                        S = eye(numel(xf), numel(xf));
                        
                        %fprintf('calculating persistence...\n')
                        % calculating persistence of brain state - control energy needed
                        % to maintain state xf
                        [~, controlInputs_persistence, numericalError] = opt_eng_cont(A_slab, T, B_slab, xf, xf, rho, S, true);
                        energyCost_persistence = sum(trapz(controlInputs_persistence.^2)/nTimeSteps); % calculating the total energetic cost as the sum of integral of squared energy trajectories, divided by number of time steps
                        
                        switch currentContrastLabel
                            case 'contrast1_threatcorrectStd'
                                totalEnergyCost_threat = [totalEnergyCost_threat; energyCost_persistence];
                            case 'contrast3_nonthreatcorrectStd'
                                totalEnergyCost_nonthreat = [totalEnergyCost_nonthreat; energyCost_persistence];
                            case 'contrast5_neutralcorrectStd'
                                totalEnergyCost_neutral = [totalEnergyCost_neutral; energyCost_persistence];
                        end
                    end
                end
            end
        end
        
        totalEnergyCost_allParameters_threat = [totalEnergyCost_allParameters_threat, totalEnergyCost_threat];
        totalEnergyCost_allParameters_nonthreat = [totalEnergyCost_allParameters_nonthreat, totalEnergyCost_nonthreat];
        totalEnergyCost_allParameters_neutral = [totalEnergyCost_allParameters_neutral, totalEnergyCost_neutral];
    end
    
    f = figure('visible', 'off'); set(gcf, 'color', 'w'); set(gca, 'FontSize', 20);
    imagesc(flip(corr(totalEnergyCost_allParameters_threat))); colorbar;
    xticks([1, nParameters]); yticks([1, nParameters]); yticklabels({'1', '0.1'}); xticklabels({'0.1', '1'});
    xlabel('\rho'); ylabel('\rho');
    currentSaveFileName = strcat(resultsDir, 'energyCost_rho_', currentTaskID, '_contrast1_threatcorrectStd', '.svg');
    saveas(f, currentSaveFileName);
    close(f);
    
    f = figure('visible', 'off'); set(gcf, 'color', 'w'); set(gca, 'FontSize', 20);
    imagesc(flip(corr(totalEnergyCost_allParameters_nonthreat))); colorbar;
    xticks([1, nParameters]); yticks([1, nParameters]); yticklabels({'1', '0.1'}); xticklabels({'0.1', '1'});
    xlabel('\rho'); ylabel('\rho');
    currentSaveFileName = strcat(resultsDir, 'energyCost_rho_', currentTaskID, '_contrast3_nonthreatcorrectStd', '.svg');
    saveas(f, currentSaveFileName);
    close(f);
    
    f = figure('visible', 'off'); set(gcf, 'color', 'w'); set(gca, 'FontSize', 20);
    imagesc(flip(corr(totalEnergyCost_allParameters_neutral))); colorbar;
    xticks([1, nParameters]); yticks([1, nParameters]); yticklabels({'1', '0.1'}); xticklabels({'0.1', '1'});
    xlabel('\rho'); ylabel('\rho');
    currentSaveFileName = strcat(resultsDir, 'energyCost_rho_', currentTaskID, '_contrast5_neutralcorrectStd', '.svg');
    saveas(f, currentSaveFileName);
    close(f);
end
