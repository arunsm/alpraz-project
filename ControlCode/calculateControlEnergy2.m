%%%%% script to calculate control energy for brain state transitions in the
%%%%% alpraz dataset

%% set parameters

clear all;

addpath(genpath('~/matlab/BCT'));

brainStates = {'spinTest'};
normalizationMethods = {'EuclideanNorm'};
controlAlgorithm = 'minimalControl';
avge_FD_thresh = [0.5]; % threshold for eliminating subjects based on motion (in mm)
parcelCoverageThresh = [0.5]; % threshold for keeping a parcel based on coverage in subject mask
structuralEdges = {'QA'};
controlInputs = {'allNodes'};

%% iterating through parameter choices
for n = 1:numel(normalizationMethods)
    currentNormalizationMethod = normalizationMethods{n};
    for b = 1:numel(brainStates)
        currentBrainStates = brainStates{b};
        for f = 1:numel(avge_FD_thresh)
            current_FD_thresh = avge_FD_thresh(f);
            for p = 1:numel(parcelCoverageThresh)
                currentParcelCoverageThresh = parcelCoverageThresh(p);
                for se = 1:numel(structuralEdges)
                    currentStructuralEdges = structuralEdges{se};
                    for c = 1:numel(controlInputs)
                        currentControlInputs = controlInputs{c};
                        switch currentStructuralEdges
                            case 'StreamlineCounts'
                                A = getAverageStreamlineCounts(); % calculating average streamline counts from Q7 data
                            case 'QA'
                                A = getQA(); % calculating quantitative anisotropy from Q7 data
                            case 'nullModel1'
                                A = getQA();
                                A = null_model_und_sign(A); % calling BCT function to randomize matrix while preserving degree and strength distributions
                            case 'nullModel2'
                                A = getQA();
                                A = randmio_und(A, 100); % calling BCT function to randomize matrix while preserving degree distribution but not strength
                        end
                        idx_brainstem = 234; % index of brain stem parcel
                        %A(idx_brainstem, :) = []; A(:, idx_brainstem) = []; % removing brain stem parcel
                        nNodes = size(A, 1);
                        T = 3; % control horizon
                        nTimeSteps = 1000; % number of time steps in simulation
                        rho = 100; % parameter that weights energy and distance constraints, set as in Betzel et al, Sci. Rep (2016)
                        
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
                            date, '_avge_FD_thresh_', num2str(current_FD_thresh), ...
                            '_parcelCoverageThresh_', num2str(currentParcelCoverageThresh), ...
                            '_', currentNormalizationMethod, '_', currentControlInputs, '_', ...
                            currentStructuralEdges, '_', currentBrainStates, '/');
                        
                        if ~exist(resultsDir)
                            mkdir(resultsDir);
                        end

                        fid = fopen(strcat(resultsDir, 'calculateOptEnergy_parameters.csv'), 'w');
                        fprintf(fid, 'controlAlgorithm,normalizationMethod,avge_FD_thresh,parcelCoverageThresh,T,rho,controlInputs,structuralEdges,brainStates\n%s,%s,%d,%d,%d,%d,%s,%s,%s', ...
                            controlAlgorithm, currentNormalizationMethod, current_FD_thresh, currentParcelCoverageThresh, T, rho, currentControlInputs, currentStructuralEdges, currentBrainStates);
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
                        
                        allControlEnergies_emotionid = {};
                        allControlEnergies_emotionrec = {};
                        allControlTrajectories_emotionid = {};
                        allControlTrajectories_emotionrec = {};
                        
                        taskIDs = {'task-emotionid', 'task-emotionrec'};
                        nTaskIDs = numel(taskIDs);
                        
                        for k = 1:nTaskIDs
                            currentTaskID = taskIDs{k};
                            fprintf('%s\n', currentTaskID);
                            
                            ctr = 1;
                            for i = 1:nSubjects
                                currentSubjectID = allSubjectIDs(i, 1); currentSubjectID_str = strcat('sub-0', num2str(currentSubjectID));
                                
                                if strcmp(currentTaskID, 'task-emotionid')
                                    currentSubjectInfo = allSubjectInfo_emotionID(allSubjectInfo_emotionID.bblid == currentSubjectID, :);
                                else
                                    currentSubjectInfo = allSubjectInfo_emotionrec(allSubjectInfo_emotionrec.bblid == currentSubjectID, :);
                                end
                                
                                allSessionIDs_currentSubject = unique(currentSubjectInfo.fmriid);
                                fprintf('%s\n', currentSubjectID_str);
                                
                                for contrast_idx = 1:nContrasts
                                    currentContrastLabel = contrastLabels{contrast_idx};
                                    fprintf('%s\n', currentContrastLabel);
                                    
                                    for j = 1:numel(allSessionIDs_currentSubject)
                                        currentSessionID = allSessionIDs_currentSubject(j); currentSessionID_str = strcat('ses-00', num2str(currentSessionID));
                                        fprintf('%s\n', currentSessionID_str);
                                        
                                        % check motion, skip if exceeds thresh
                                        current_avge_FD = subjectMotion.avge_FD(strcmp(subjectMotion.id0, currentSubjectID_str) & strcmp(subjectMotion.id1, currentSessionID_str) & strcmp(subjectMotion.id2, currentTaskID));
                                        if current_avge_FD > current_FD_thresh
                                            fprintf('Skipping subject %s, task %s: avge FD exceeds threshold in session %s\n', currentSubjectID_str, currentTaskID, currentSessionID_str);
                                            break;
                                        else
                                            
                                            % looking up group and drug info
                                            currentGroup = currentSubjectInfo.group(currentSubjectInfo.fmriid == currentSessionID); % looking up group info (0=control, 1=relative)
                                            currentDrug = currentSubjectInfo.drug(currentSubjectInfo.fmriid == currentSessionID); % looking up drug info (0=alpraz, 1=placebo)
                                            % looking up accuracy
                                            current_pctcorr_threat = currentSubjectInfo.pctcorr_threat(currentSubjectInfo.fmriid == currentSessionID);
                                            current_pctcorr_nonthreat = currentSubjectInfo.pctcorr_nonthreat(currentSubjectInfo.fmriid == currentSessionID);
                                            current_pctcorr_neutral = currentSubjectInfo.pctcorrneutral(currentSubjectInfo.fmriid == currentSessionID);
                                            % looking up median reaction times
                                            current_rtmdn_threatcorr = currentSubjectInfo.rtmdn_threatcorr(currentSubjectInfo.fmriid == currentSessionID);
                                            current_rtmdn_nonthreatcorr = currentSubjectInfo.rtmdn_nonthreatcorr(currentSubjectInfo.fmriid == currentSessionID);
                                            current_rtmdn_neutralcorr = currentSubjectInfo.rtmdneutralcorr(currentSubjectInfo.fmriid == currentSessionID);
                                            % looking up STAI-TRAIT scores and
                                            % alpraz levels
                                            current_STAI_TRAIT = currentSubjectInfo.STAI_TRAIT(currentSubjectInfo.fmriid == currentSessionID);
                                            current_alpraz_levels = currentSubjectInfo.alpraz_levels(currentSubjectInfo.fmriid == currentSessionID);
                                            % looking up demographic variables
                                            currentAge = subjectDemographics.AgeAtFMRI(subjectDemographics.bblid==currentSubjectID); % looking up age
                                            currentGender = subjectDemographics.sex_M0F1(subjectDemographics.bblid==currentSubjectID); % looking up gender
                                            currentSISTOTAL = subjectDemographics.SISTOTAL(subjectDemographics.bblid==currentSubjectID); % looking up SISTOTAL score
                                            
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
                                                x = X.data'; % final state = contrasts for task; note that this also contains values for nodes outside imaging slab
                                                %xf(idx_brainstem) = []; % removing brain stem parcel
                                                x = x(parcelsToInclude); % only keeping nodes within imaging slab
                                                
                                                switch currentBrainStates
                                                    case 'betas'
                                                        xf = x; % using beta values from first-level GLM as brain states
                                                    case 'spinTest'
                                                        scrambled_xf_savePath = strcat('spinTestResults', filesep, currentSubjectID_str, '_', currentSessionID_str, '_', currentTaskID, '_', currentContrastLabel, '_', atlasLabel);
                                                        xf = calculateScrambledBetas(x, scrambled_xf_savePath); % using scrambled activation patterns calculated using spin-test
                                                        % replace subcortical data with betas; spin-test only scrambles cortical data
                                                        subcorticalIndices = isnan(xf);
                                                        xf(subcorticalIndices) = x(subcorticalIndices);
                                                end
                                                
                                                switch currentNormalizationMethod
                                                    case 'EuclideanNorm'
                                                        xf = x/norm(x); % normalizing brain states by Euclidean norm
                                                    case 'none'
                                                        xf = x; % no normalization applied
                                                end
                                            else
                                                fprintf('data for %s, %s, %s not available\n', currentSubjectID_str, currentSessionID_str, currentTaskID);
                                                break;
                                            end
                                            
                                            % setting initial state to 0's
                                            %x0 = zeros(size(xf));
                                            
                                            S = eye(numel(xf), numel(xf));
                                            
                                            % setting diagonal elements of S to include only nodes within slab
                                            %S = diag(parcelsToInclude);
                                            %S(idx_brainstem, :) = []; S(:, idx_brainstem) = []; % removing brain stem parcel
                                            
                                            fprintf('calculating persistence...\n')
                                            
                                            switch controlAlgorithm
                                                case 'optimalControl'
                                                    % calculating persistence of brain state - control energy needed
                                                    % to maintain state xf
                                                    [stateTrajectories_persistence, controlInputs_persistence, numericalError] = opt_eng_cont(A_slab, T, B_slab, xf, xf, rho, S, true);
                                                    stateTrajectories_persistence = stateTrajectories_persistence(:, 1:nNodes_currentSubject);
                                                    energyCost_persistence = trapz(controlInputs_persistence.^2)/nTimeSteps; % calculating the total energetic cost as the sum of integral of squared energy trajectories, divided by number of time steps
                                                    
                                                case 'minimalControl'
                                                    % calculating persistence of brain state - control energy needed
                                                    % to maintain state xf
                                                    [stateTrajectories_persistence, controlInputs_persistence, numericalError] = min_eng_cont(A_slab, T, B_slab, xf, xf, true);
                                                    stateTrajectories_persistence = stateTrajectories_persistence(:, 1:nNodes_currentSubject);
                                                    energyCost_persistence = trapz(controlInputs_persistence.^2)/nTimeSteps; % calculating the total energetic cost as the sum of integral of squared energy trajectories, divided by number of time steps
                                            end
                                            
                                            persistence_allNodes = sum(energyCost_persistence);
                                            persistence_visual = sum(energyCost_persistence(finalLabels_slab==1));
                                            persistence_ventralAttention = sum(energyCost_persistence(finalLabels_slab==4));
                                            persistence_limbic = sum(energyCost_persistence(finalLabels_slab==5));
                                            persistence_frontoparietal = sum(energyCost_persistence(finalLabels_slab==6));
                                            persistence_DMN = sum(energyCost_persistence(finalLabels_slab==7));
                                            persistence_subcortex = sum(energyCost_persistence(finalLabels_slab==8));
                                            
                                            fprintf('calculating control impact...\n')
                                            
                                            controlImpact_currentSubject = zeros(1, nNodes_currentSubject);
                                            switch controlAlgorithm
                                                case 'optimalControl'
                                                    parfor a = 1:nNodes_currentSubject
                                                        B_controlImpact = B_slab;
                                                        if B_controlImpact(a, a) == 1
                                                            B_controlImpact(a, a) = 0;
                                                            [~, U_opt_stability_controlImpact, ~] = opt_eng_cont(A_slab, T, B_controlImpact, xf, xf, rho, S, true);
                                                            energyCost_persistence_controlImpact = trapz(U_opt_stability_controlImpact.^2)/nTimeSteps;
                                                            controlImpact_currentSubject(a) = log10(sum(energyCost_persistence_controlImpact)/sum(energyCost_persistence));
                                                        end
                                                    end
                                                    
                                                case 'minimalControl'
                                                    parfor a = 1:nNodes_currentSubject
                                                        B_controlImpact = B_slab;
                                                        if B_controlImpact(a, a) == 1
                                                            B_controlImpact(a, a) = 0;
                                                            [~, U_opt_stability_controlImpact, ~] = min_eng_cont(A_slab, T, B_controlImpact, xf, xf, true);
                                                            energyCost_persistence_controlImpact = trapz(U_opt_stability_controlImpact.^2)/nTimeSteps;
                                                            controlImpact_currentSubject(a) = log10(sum(energyCost_persistence_controlImpact)/sum(energyCost_persistence));
                                                        end
                                                    end
                                            end
                                            
                                            switch currentTaskID
                                                case 'task-emotionid'
                                                    allControlEnergies_emotionid(ctr, :) = {currentSubjectID, currentSessionID, currentContrastLabel, ...
                                                        persistence_allNodes, persistence_visual, persistence_ventralAttention, ...
                                                        persistence_limbic, persistence_frontoparietal, persistence_DMN, persistence_subcortex, ...
                                                        currentGroup, currentDrug, currentAge, currentGender, currentSISTOTAL, ...
                                                        current_pctcorr_threat, current_pctcorr_nonthreat, current_pctcorr_neutral, ...
                                                        current_rtmdn_threatcorr, current_rtmdn_nonthreatcorr, current_rtmdn_neutralcorr, ...
                                                        current_STAI_TRAIT, current_alpraz_levels};
                                                    allControlTrajectories_emotionid(ctr, :) = {currentSubjectID, currentSessionID, currentContrastLabel, ...
                                                        parcelsToInclude_idx, xf, controlInputs_persistence, stateTrajectories_persistence, controlImpact_currentSubject, ...
                                                        numericalError, currentGroup, currentDrug, currentAge, currentGender};
                                                    
                                                case 'task-emotionrec'
                                                    allControlEnergies_emotionrec(ctr, :) = {currentSubjectID, currentSessionID, currentContrastLabel, ...
                                                        persistence_allNodes, persistence_visual, persistence_ventralAttention, ...
                                                        persistence_limbic, persistence_frontoparietal, persistence_DMN, persistence_subcortex, ...
                                                        currentGroup, currentDrug, currentAge, currentGender, currentSISTOTAL, ...
                                                        current_pctcorr_threat, current_pctcorr_nonthreat, current_pctcorr_neutral, ...
                                                        current_rtmdn_threatcorr, current_rtmdn_nonthreatcorr, current_rtmdn_neutralcorr, ...
                                                        current_STAI_TRAIT, current_alpraz_levels};
                                                    allControlTrajectories_emotionrec(ctr, :) = {currentSubjectID, currentSessionID, currentContrastLabel, ...
                                                        parcelsToInclude_idx, xf, controlInputs_persistence, stateTrajectories_persistence, controlImpact_currentSubject, ...
                                                        numericalError, currentGroup, currentDrug, currentAge, currentGender};
                                            end
                                            ctr = ctr + 1;
                                        end
                                    end
                                end
                            end
                        end
                        
                        allControlEnergies_emotionid = cell2table(allControlEnergies_emotionid, 'VariableNames', {'subjectID', 'sessionID', 'contrast', ...
                            'persistence_allNodes', 'persistence_visual', 'persistence_ventralAttention', ...
                            'persistence_limbic', 'persistence_frontoparietal', 'persistence_DMN', 'persistence_subcortex', ...
                            'group', 'drug', 'age', 'gender', 'SISTOTAL', 'pctcorr_threat', 'pctcorr_nonthreat', 'pctcorr_neutral', ...
                            'rtmdn_threatcorr', 'rtmdn_nonthreatcorr', 'rtmdn_neutralcorr', 'STAI_TRAIT', 'alpraz_levels'});
                        
                        allControlEnergies_emotionrec = cell2table(allControlEnergies_emotionrec, 'VariableNames', {'subjectID', 'sessionID', 'contrast', ...
                            'persistence_allNodes', 'persistence_visual', 'persistence_ventralAttention', ...
                            'persistence_limbic', 'persistence_frontoparietal', 'persistence_DMN', 'persistence_subcortex', ...
                            'group', 'drug', 'age', 'gender', 'SISTOTAL', 'pctcorr_threat', 'pctcorr_nonthreat', 'pctcorr_neutral', ...
                            'rtmdn_threatcorr', 'rtmdn_nonthreatcorr', 'rtmdn_neutralcorr', 'STAI_TRAIT', 'alpraz_levels'});
                        
                        allControlTrajectories_emotionid = cell2table(allControlTrajectories_emotionid, 'VariableNames', {'subjectID', 'sessionID', 'contrast', ...
                            'parcelsToInclude_idx', 'xf', 'controlInputs_persistence', 'stateTrajectories_persistence', 'controlImpact_persistence', 'numericalError', ...
                            'group', 'drug', 'age', 'gender'});
                        
                        allControlTrajectories_emotionrec = cell2table(allControlTrajectories_emotionrec, 'VariableNames', {'subjectID', 'sessionID', 'contrast', ...
                            'parcelsToInclude_idx', 'xf', 'controlInputs_persistence', 'stateTrajectories_persistence', 'controlImpact_persistence', 'numericalError', ...
                            'group', 'drug', 'age', 'gender'});
                        
                        % removing subjects with incomplete data (fewer than 2 sessions)
                        for i = 1:nSubjects
                            currentSubjectID = allSubjectIDs(i);
                            X = allControlEnergies_emotionid(allControlEnergies_emotionid.subjectID == currentSubjectID, :);
                            allSessions_currentSubject = unique(X.sessionID);
                            if numel(allSessions_currentSubject) < 2
                                allControlEnergies_emotionid(allControlEnergies_emotionid.subjectID == currentSubjectID, :) = [];
                                allControlTrajectories_emotionid(allControlTrajectories_emotionid.subjectID == currentSubjectID, :) = [];
                            end
                        end
                        
                        for i = 1:nSubjects
                            currentSubjectID = allSubjectIDs(i);
                            X = allControlEnergies_emotionrec(allControlEnergies_emotionrec.subjectID == currentSubjectID, :);
                            allSessions_currentSubject = unique(X.sessionID);
                            if numel(allSessions_currentSubject) < 2
                                allControlEnergies_emotionrec(allControlEnergies_emotionrec.subjectID == currentSubjectID, :) = [];
                                allControlTrajectories_emotionrec(allControlTrajectories_emotionrec.subjectID == currentSubjectID, :) = [];
                            end
                        end
                        
                        writetable(allControlEnergies_emotionid, strcat(resultsDir, 'allControlEnergies_emotionid.csv'));
                        writetable(allControlEnergies_emotionrec, strcat(resultsDir, 'allControlEnergies_emotionrec.csv'));
                        
                        createContrastFiles(resultsDir);
                        
                        save(strcat(resultsDir, 'allControlTrajectories_emotionid.mat'), 'allControlTrajectories_emotionid');
                        save(strcat(resultsDir, 'allControlTrajectories_emotionrec.mat'), 'allControlTrajectories_emotionrec');
                        save(strcat(resultsDir, 'structuralAdjacencyMatrix.mat'), 'A');
                    end
                end
            end
        end
    end
end