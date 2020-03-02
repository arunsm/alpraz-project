%%%%% script to calculate control energy for brain state transitions in the
%%%%% alpraz dataset

%% set parameters

clear all;
controlAlgorithm = 'optimalControl';
avge_FD_thresh = [0.25, 0.5]; % threshold for eliminating subjects based on motion (in mm)
parcelCoverageThresh = [0.25, 0.5, 0.75]; % threshold for keeping a parcel based on coverage in subject mask
structuralEdges = {'QA', 'StreamlineCounts', 'zero'};

%% iterating through parameter choices
for f = 1:numel(avge_FD_thresh)
    current_FD_thresh = avge_FD_thresh(f);
    for p = 1:numel(parcelCoverageThresh)
        currentParcelCoverageThresh = parcelCoverageThresh(p);
        for se = 1:numel(structuralEdges)
            currentStructuralEdges = structuralEdges{se};
            
            switch currentStructuralEdges
                case 'StreamlineCounts'
                    A = getAverageStreamlineCounts(); % calculating average streamline counts from Q7 data
                case 'QA'
                    A = getQA(); % calculating quantitative anisotropy from Q7 data
                case 'zero'
                    A = zeros(234, 234); % creating a null adjacency matrix
            end
            idx_brainstem = 234; % index of brain stem parcel
            %A(idx_brainstem, :) = []; A(:, idx_brainstem) = []; % removing brain stem parcel
            nNodes = size(A, 1);
            T = 3; % control horizon
            nTimeSteps = 1000; % number of time steps in simulation
            rho = 100; % parameter that weights energy and distance constraints, set as in Betzel et al, Sci. Rep (2016)
            controlInputs = 'allNodes'; % parameter to set control inputs; possible values are 'allNodes', 'cognitiveControlRegions'
            
            % loading parcel Yeo sub-network label for Lausanne-234 parcels
            % 1=Visual, 2=Somatomator, 3=Dorsal Attention, 4=Ventral Attention, 5=Limbic, 6=Frontoparietal Control, 7=Default Mode, 8=Subcortical
            load('LindenYeoPurity/yeo7netlabelsLaus125EJC.mat', 'finalLabels');
            %finalLabels(idx_brainstem) = []; % removing brain stem parcel
            subSystemLabels = {'Visual', 'Somatomator', 'DorsalAttention', 'VentralAttention', 'Limbic', 'FrontoparietalControl', 'DefaultMode', 'Subcortical'};
            nSubSystems = numel(subSystemLabels);
            
            % setting control inputs
            switch controlInputs
                case 'allNodes'
                    B = eye(size(A)); % all nodes as control inputs
                case 'cognitiveControlRegions'
                    B = diag(finalLabels==3 | finalLabels==4 | finalLabels==6); % dorsal and ventral attention networks + frontoparietal control network are control inputs
            end
            
            % creating results folder and writing parameters as text file
            resultsDir = strcat('/data/joy/BBL/studies/alpraz/scripts/ControlCode/Results/', ...
                date, '_avge_FD_thresh_', num2str(current_FD_thresh), '_parcelCoverageThresh_', num2str(currentParcelCoverageThresh), '_', currentStructuralEdges, '/');
            
            if ~exist(resultsDir)
                mkdir(resultsDir);
            end
            fid = fopen(strcat(resultsDir, 'calculateOptEnergy_parameters.csv'), 'w');
            fprintf(fid, 'controlAlgorithm,avge_FD_thresh,parcelCoverageThresh,T,rho,controlInputs,structuralEdges\n%s,%d,%d,%d,%d,%s,%s', controlAlgorithm, current_FD_thresh, currentParcelCoverageThresh, T, rho, controlInputs, currentStructuralEdges);
            fclose(fid);
            
            %% set paths and read in subject data
            
            pathToQualityFile = '/data/joy/BBL/studies/alpraz/rawData/derivatives/xcp_output_allTasks2/group/n191_quality.csv';
            pathToMotionFile = '/data/joy/BBL/studies/alpraz/rawData/derivatives/fmriprep/subjectMotion.csv';
            rootFolder_contrasts = '/data/joy/BBL/studies/alpraz/rawData/derivatives/xcp_output_allTasks2/';
            
            contrastLabels = {'contrast1_threatcorrectStd', 'contrast3_nonthreatcorrectStd', ...
                'contrast5_neutralcorrectStd'};
            
            nContrasts = numel(contrastLabels);
            atlasLabel = 'lausanne_ROIv_scale125_dilated';
            
            QA = xlsread('/data/joy/BBL/studies/alpraz/rawData/derivatives/stimDesigns/AlprazQA_forArun.xlsx');
            drug_group_info = [QA(:, 2) QA(:, 3) QA(:, 6), QA(:, 8)]; % look-up table with [subjectID, sessionID, group(0=control, 1=relative), drug(0=no,1=yes)]
            
            motion = readtable(pathToMotionFile);
            avge_FD = motion.avge_FD;
            
            Q = readtable(pathToQualityFile);
            totalSubjects_and_Sessions = size(Q, 1);
            allSubjectIDs = unique(Q.id0);
            nSubjects = numel(allSubjectIDs);
            
            %% calculate optimal control energy
            
            % initializing cell arrays containing control energy and stability values for
            % all contrasts, averaged across each individual sub-system and overall
            
            allControlEnergies_emotionid = cell(8, nContrasts+1);
            allControlEnergies_emotionrec = cell(8, nContrasts+1);
            allControlTrajectories_emotionid = cell(5, nContrasts+1);
            allControlTrajectories_emotionrec = cell(5, nContrasts+1);
            
            allControlEnergies_emotionid{2, 1} = 'Total Control Energy (Alpraz+)'; allControlEnergies_emotionrec{2, 1} = 'Total Control Energy (Alpraz+)';
            allControlEnergies_emotionid{3, 1} = 'Total Control Energy (Alpraz-)'; allControlEnergies_emotionrec{3, 1} = 'Total Control Energy (Alpraz-)';
            allControlEnergies_emotionid{4, 1} = 'Persistence (Alpraz+)'; allControlEnergies_emotionrec{4, 1} = 'Persistence (Alpraz+)';
            allControlEnergies_emotionid{5, 1} = 'Persistence (Alpraz-)'; allControlEnergies_emotionrec{5, 1} = 'Persistence (Alpraz-)';
            allControlEnergies_emotionid{6, 1} = 'Control Impact (Alpraz+)'; allControlEnergies_emotionrec{6, 1} = 'Control Impact (Alpraz+)';
            allControlEnergies_emotionid{7, 1} = 'Control Impact (Alpraz-)'; allControlEnergies_emotionrec{7, 1} = 'Control Impact (Alpraz-)';
            allControlEnergies_emotionid{8, 1} = 'Group (0=control, 1=relative)'; allControlEnergies_emotionrec{8, 1} = 'Group (0=control, 1=relative)';
            
            allControlTrajectories_emotionid{2, 1} = 'Control Input (Alpraz+)'; allControlTrajectories_emotionrec{2, 1} = 'Control Input (Alpraz+)';
            allControlTrajectories_emotionid{3, 1} = 'Control Input (Alpraz-)'; allControlTrajectories_emotionrec{3, 1} = 'Control Input (Alpraz-)';
            allControlTrajectories_emotionid{4, 1} = 'State Trajectories (Alpraz+)'; allControlTrajectories_emotionrec{4, 1} = 'State Trajectories (Alpraz+)';
            allControlTrajectories_emotionid{5, 1} = 'State Trajectories (Alpraz-)'; allControlTrajectories_emotionrec{5, 1} = 'State Trajectories (Alpraz-)';
            
            for contrast_idx = 1:nContrasts
                currentContrastLabel = contrastLabels{contrast_idx};
                allControlEnergies_emotionid{1, contrast_idx+1} = currentContrastLabel; allControlEnergies_emotionrec{1, contrast_idx+1} = currentContrastLabel;
                allControlTrajectories_emotionid{1, contrast_idx+1} = currentContrastLabel; allControlTrajectories_emotionrec{1, contrast_idx+1} = currentContrastLabel;
                
                fprintf('%s\n', currentContrastLabel)
                
                group_idx_emotionid = NaN(nSubjects, 1); group_idx_emotionrec = NaN(nSubjects, 1);
                
                logsum_energy_cost_emotionid_drug = NaN(nSubjects, nSubSystems+1); logsum_energy_cost_emotionid_noDrug = NaN(nSubjects, nSubSystems+1);
                logsum_energy_cost_emotionrec_drug = NaN(nSubjects, nSubSystems+1); logsum_energy_cost_emotionrec_noDrug = NaN(nSubjects, nSubSystems+1);
                
                persistence_emotionid_drug = NaN(nSubjects, nSubSystems+1); persistence_emotionid_noDrug = NaN(nSubjects, nSubSystems+1);
                persistence_emotionrec_drug = NaN(nSubjects, nSubSystems+1); persistence_emotionrec_noDrug = NaN(nSubjects, nSubSystems+1);
                
                controlImpact_persistence_emotionid_drug = NaN(nSubjects, nNodes); controlImpact_persistence_emotionid_noDrug = NaN(nSubjects, nNodes);
                controlImpact_persistence_emotionrec_drug = NaN(nSubjects, nNodes); controlImpact_persistence_emotionrec_noDrug = NaN(nSubjects, nNodes);
                
                controlInputs_persistence_emotionid_drug = NaN(nTimeSteps+1, nNodes, nSubjects); controlInputs_persistence_emotionid_noDrug = NaN(nTimeSteps+1, nNodes, nSubjects);
                controlInputs_persistence_emotionrec_drug = NaN(nTimeSteps+1, nNodes, nSubjects); controlInputs_persistence_emotionrec_noDrug = NaN(nTimeSteps+1, nNodes, nSubjects);
                
                stateTrajectories_persistence_emotionid_drug = NaN(nTimeSteps+1, nNodes, nSubjects); stateTrajectories_persistence_emotionid_noDrug = NaN(nTimeSteps+1, nNodes, nSubjects);
                stateTrajectories_persistence_emotionrec_drug = NaN(nTimeSteps+1, nNodes, nSubjects); stateTrajectories_persistence_emotionrec_noDrug = NaN(nTimeSteps+1, nNodes, nSubjects);
                
                numericalError_emotionid = []; numericalError_emotionrec = [];
                
                for i = 1:totalSubjects_and_Sessions
                    currentSubjectID = Q.id0{i}; currentSubjectID_num = str2double(extractAfter(currentSubjectID, 'sub-'));
                    currentSessionID = Q.id1{i}; currentSessionID_num = str2double(extractAfter(currentSessionID, 'ses-'));
                    currentTaskID = Q.id2{i};
                    fprintf('%s\t%s\t%s\n', currentSubjectID, currentSessionID, currentTaskID);
                    
                    subject_idx = find(strcmp(currentSubjectID, allSubjectIDs));
                    group = drug_group_info(drug_group_info(:, 1)==currentSubjectID_num & drug_group_info(:, 2)==currentSessionID_num , 3); % looking up group info
                    drug = drug_group_info(drug_group_info(:, 1)==currentSubjectID_num & drug_group_info(:, 2)==currentSessionID_num , 4); % looking up drug info
                    
                    current_avge_FD = avge_FD(i);
                    if current_avge_FD > current_FD_thresh
                        fprintf('Skipping subject %s, session %s, task %s: avge FD exceeds threshold\n', currentSubjectID, currentSessionID, currentTaskID);
                        continue;
                    else
                        
                        pathToMask = strcat(rootFolder_contrasts, currentSubjectID, filesep, currentSessionID, filesep, currentTaskID, '/norm/', currentSubjectID, '_', currentSessionID, '_', currentTaskID, '_maskStd_', atlasLabel, '.txt');
                        M = importdata(pathToMask);
                        parcelCoverage = M.data';
                        parcelsToInclude = parcelCoverage > currentParcelCoverageThresh; % boolean list of parcels whose coverage exceeds the threshold
                        parcelsToInclude_idx = find(parcelCoverage > currentParcelCoverageThresh);
                        nNodes_currentSubject = numel(parcelsToInclude_idx);
                        
                        % truncating A and B matrices to be within imaging slab
                        A_slab = A(parcelsToInclude, parcelsToInclude);
                        B_slab = B(parcelsToInclude, parcelsToInclude);
                        finalLabels_slab = finalLabels(parcelsToInclude_idx);
                        
                        pathToContrast = strcat(rootFolder_contrasts, currentSubjectID, filesep, currentSessionID, filesep, currentTaskID, '/norm/', currentSubjectID, '_', currentSessionID, '_', currentTaskID, '_', currentContrastLabel, '_', atlasLabel, '.txt');
                        X = importdata(pathToContrast);
                        if isstruct(X)
                            xf = X.data'; % final state = contrasts for task; note that this also contains values for nodes outside imaging slab
                            %xf(idx_brainstem) = []; % removing brain stem parcel
                            xf = xf(parcelsToInclude); % only keeping nodes within imaging slab
                            xf = xf/norm(xf); % IMPORTANT - normalizing brain states by Euclidean norm
                        else
                            fprintf('data for %s, %s, %s not available\n', currentSubjectID, currentSessionID, currentTaskID);
                            continue;
                        end
                        
                        % setting initial state to 0's
                        x0 = zeros(size(xf));
                        
                        S = eye(numel(xf), numel(xf));
                        
                        % setting diagonal elements of S to include only nodes within slab
                        %S = diag(parcelsToInclude);
                        %S(idx_brainstem, :) = []; S(:, idx_brainstem) = []; % removing brain stem parcel
                        
                        switch controlAlgorithm
                            case 'optimalControl'
                                fprintf('calculating control energy for transition from 0-state...\n')
                                
                                % calculating optimal control
                                [stateTrajectories_x0_xf, controlInputs_x0_xf, n_err] = opt_eng_cont(A_slab, T, B_slab, x0, xf, rho, S, true);
                                stateTrajectories_x0_xf = stateTrajectories_x0_xf(:, 1:nNodes_currentSubject);
                                energyCost = trapz(controlInputs_x0_xf.^2)/nTimeSteps; % calculating the total energetic cost as the sum of integral of squared energy trajectories
                                sumEnergyCost = zeros(1, nSubSystems+1);
                                for s = 1:nSubSystems
                                    sumEnergyCost(s) = sum(energyCost(finalLabels_slab==s)); % summing over all nodes within sub-system and included in slab
                                end
                                sumEnergyCost(end) = sum(energyCost); % last entry stores energy cost summed over all control nodes in imaging slab
                                
                                fprintf('calculating persistence...\n')
                                % calculating stability of brain state - control energy needed
                                % to maintain state xf
                                [stateTrajectories_persistence, controlInputs_persistence, ~] = opt_eng_cont(A_slab, T, B_slab, xf, xf, rho, S, true);
                                stateTrajectories_persistence = stateTrajectories_persistence(:, 1:nNodes_currentSubject);
                                energyCost_persistence = trapz(controlInputs_persistence.^2)/nTimeSteps; % calculating the total energetic cost as the sum of integral of squared energy trajectories
                                
                                persistence = zeros(1, nSubSystems+1);
                                for s = 1:nSubSystems
                                    persistence(s) = 1/log10(sum(energyCost_persistence(finalLabels_slab==s)));
                                end
                                persistence(end) = 1/log10(sum(energyCost_persistence)); % last entry stores stability summed over all control inputs
                                
                                fprintf('calculating control impact...\n')
                                controlImpact_currentSubject = zeros(1, nNodes_currentSubject);
                                parfor j = 1:nNodes_currentSubject
                                    B_controlImpact = B_slab; B_controlImpact(j, j) = 0;
                                    [~, U_opt_stability_controlImpact, ~] = opt_eng_cont(A_slab, T, B_controlImpact, xf, xf, rho, S, true);
                                    energyCost_persistence_controlImpact = trapz(U_opt_stability_controlImpact.^2)/nTimeSteps;
                                    controlImpact_currentSubject(j) = log10(sum(energyCost_persistence_controlImpact)/sum(energyCost_persistence));
                                end
                                
                            case 'minimalControl'
                                fprintf('calculating control energy for transition from 0-state...\n')
                                
                                % calculating minimal control
                                [stateTrajectories_x0_xf, controlInputs_x0_xf, n_err] = min_eng_cont(A_slab, T, B_slab, x0, xf, true);
                                stateTrajectories_x0_xf = stateTrajectories_x0_xf(:, 1:nNodes_currentSubject);
                                energyCost = trapz(controlInputs_x0_xf.^2)/nTimeSteps; % calculating the total energetic cost as the sum of integral of squared energy trajectories
                                sumEnergyCost = zeros(1, nSubSystems+1);
                                for s = 1:nSubSystems
                                    sumEnergyCost(s) = sum(energyCost(finalLabels_slab==s));
                                end
                                sumEnergyCost(end) = sum(energyCost); % last entry stores energy cost summed over all control inputs
                                
                                fprintf('calculating stability...\n')
                                % calculating stability of brain state - control energy needed
                                % to maintain state xf
                                [stateTrajectories_persistence, controlInputs_persistence, ~] = min_eng_cont(A_slab, T, B_slab, xf, xf, true);
                                stateTrajectories_persistence = stateTrajectories_persistence(:, 1:nNodes_currentSubject);
                                energyCost_persistence = trapz(controlInputs_persistence.^2)/nTimeSteps; % calculating the total energetic cost as the sum of integral of squared energy trajectories
                                
                                persistence = zeros(1, nSubSystems+1);
                                for s = 1:nSubSystems
                                    persistence(s) = 1/log10(sum(energyCost_persistence(finalLabels_slab==s)));
                                end
                                persistence(end) = 1/log10(sum(energyCost_persistence)); % last entry stores stability summed over all control inputs
                                
                                %                                 fprintf('calculating control impact...\n')
                                %                                 controlImpact_currentSubject = zeros(1, nNodes);
                                %                                 parfor j = 1:nNodes_currentSubject
                                %                                     B_controlImpact = B_slab; B_controlImpact(j, j) = 0;
                                %                                     [~, U_opt_stability_controlImpact, ~] = opt_eng_cont(A_slab, T, B_controlImpact, xf, xf, rho, S, true);
                                %                                     energyCost_stability_controlImpact = trapz(U_opt_stability_controlImpact.^2)/nTimeSteps;
                                %                                     controlImpact_currentSubject(j) = log10(sum(energyCost_stability_controlImpact)/sum(energyCost_stability));
                                %                                 end
                        end
                        
                        switch currentTaskID
                            case 'task-emotionid'
                                if drug
                                    logsum_energy_cost_emotionid_drug(subject_idx, :) = log10(sumEnergyCost); % storing log of total energy cost
                                    persistence_emotionid_drug(subject_idx, :) = persistence; % storing persistence
                                    controlInputs_persistence_emotionid_drug(:, parcelsToInclude_idx, subject_idx) = controlInputs_persistence; % storing control input trajectories for persistence calculation
                                    stateTrajectories_persistence_emotionid_drug(:, parcelsToInclude_idx, subject_idx) = stateTrajectories_persistence; % storing state trajectories for persistence calculation
                                    controlImpact_persistence_emotionid_drug(subject_idx, parcelsToInclude_idx) = controlImpact_currentSubject;
                                else
                                    logsum_energy_cost_emotionid_noDrug(subject_idx, :) = log10(sumEnergyCost); % storing log of total energy cost
                                    persistence_emotionid_noDrug(subject_idx, :) = persistence; % storing stability
                                    controlInputs_persistence_emotionid_noDrug(:, parcelsToInclude_idx, subject_idx) = controlInputs_persistence; % storing control input trajectories for persistence calculation
                                    stateTrajectories_persistence_emotionid_noDrug(:, parcelsToInclude_idx, subject_idx) = stateTrajectories_persistence; % storing state trajectories for persistence calculation
                                    controlImpact_persistence_emotionid_noDrug(subject_idx, parcelsToInclude_idx) = controlImpact_currentSubject;
                                end
                                numericalError_emotionid = [numericalError_emotionid, n_err];
                                group_idx_emotionid(subject_idx) = group;
                                
                            case 'task-emotionrec'
                                if drug
                                    logsum_energy_cost_emotionrec_drug(subject_idx, :) = log10(sumEnergyCost); % storing log of total energy cost
                                    persistence_emotionrec_drug(subject_idx, :) = persistence; % storing stability
                                    controlInputs_persistence_emotionrec_drug(:, parcelsToInclude_idx, subject_idx) = controlInputs_persistence; % storing control input trajectories for persistence calculation
                                    stateTrajectories_persistence_emotionrec_drug(:, parcelsToInclude_idx, subject_idx) = stateTrajectories_persistence; % storing state trajectories for persistence calculation
                                    controlImpact_persistence_emotionrec_drug(subject_idx, parcelsToInclude_idx) = controlImpact_currentSubject;
                                else
                                    logsum_energy_cost_emotionrec_noDrug(subject_idx, :) = log10(sumEnergyCost); % storing log of total energy cost
                                    persistence_emotionrec_noDrug(subject_idx, :) = persistence; % storing stability
                                    controlInputs_persistence_emotionrec_noDrug(:, parcelsToInclude_idx, subject_idx) = controlInputs_persistence; % storing control input trajectories for persistence calculation
                                    stateTrajectories_persistence_emotionrec_noDrug(:, parcelsToInclude_idx, subject_idx) = stateTrajectories_persistence; % storing state trajectories for persistence calculation
                                    controlImpact_persistence_emotionrec_noDrug(subject_idx, parcelsToInclude_idx) = controlImpact_currentSubject;
                                end
                                numericalError_emotionrec = [numericalError_emotionrec, n_err];
                                group_idx_emotionrec(subject_idx) = group;
                        end
                    end
                end
                
                % removing NaN rows
                idx = [find(isnan(logsum_energy_cost_emotionid_drug(:, end))); find(isnan(logsum_energy_cost_emotionid_noDrug(:, end)))];
                logsum_energy_cost_emotionid_drug(idx, :) = []; logsum_energy_cost_emotionid_noDrug(idx, :) = [];
                controlInputs_persistence_emotionid_drug(:, :, idx) = []; controlInputs_persistence_emotionid_noDrug(:, :, idx) = [];
                stateTrajectories_persistence_emotionid_drug(:, :, idx) = []; stateTrajectories_persistence_emotionid_noDrug(:, :, idx) = [];
                controlImpact_persistence_emotionid_drug(idx, :) = []; controlImpact_persistence_emotionid_noDrug(idx, :) = [];
                group_idx_emotionid(idx) = [];
                
                idx = [find(isnan(logsum_energy_cost_emotionrec_drug(:, end))); find(isnan(logsum_energy_cost_emotionrec_noDrug(:, end)))];
                logsum_energy_cost_emotionrec_drug(idx, :) = []; logsum_energy_cost_emotionrec_noDrug(idx, :) = [];
                controlInputs_persistence_emotionrec_drug(:, :, idx) = []; controlInputs_persistence_emotionrec_noDrug(:, :, idx) = [];
                stateTrajectories_persistence_emotionrec_drug(:, :, idx) = []; stateTrajectories_persistence_emotionrec_noDrug(:, :, idx) = [];
                controlImpact_persistence_emotionrec_drug(idx, :) = []; controlImpact_persistence_emotionrec_noDrug(idx, :) = [];
                group_idx_emotionrec(idx) = [];
                
                idx = [find(isnan(persistence_emotionid_drug(:, end))); find(isnan(persistence_emotionid_noDrug(:, end)))];
                persistence_emotionid_drug(idx, :) = []; persistence_emotionid_noDrug(idx, :) = [];
                idx = [find(isnan(persistence_emotionrec_drug(:, end))); find(isnan(persistence_emotionrec_noDrug(:, end)))];
                persistence_emotionrec_drug(idx, :) = []; persistence_emotionrec_noDrug(idx, :) = [];
                
                allControlEnergies_emotionid{2, contrast_idx+1} = logsum_energy_cost_emotionid_drug; allControlEnergies_emotionrec{2, contrast_idx+1} = logsum_energy_cost_emotionrec_drug;
                allControlEnergies_emotionid{3, contrast_idx+1} = logsum_energy_cost_emotionid_noDrug; allControlEnergies_emotionrec{3, contrast_idx+1} = logsum_energy_cost_emotionrec_noDrug;
                allControlEnergies_emotionid{4, contrast_idx+1} = persistence_emotionid_drug; allControlEnergies_emotionrec{4, contrast_idx+1} = persistence_emotionrec_drug;
                allControlEnergies_emotionid{5, contrast_idx+1} = persistence_emotionid_noDrug; allControlEnergies_emotionrec{5, contrast_idx+1} = persistence_emotionrec_noDrug;
                allControlEnergies_emotionid{6, contrast_idx+1} = controlImpact_persistence_emotionid_drug; allControlEnergies_emotionrec{6, contrast_idx+1} = controlImpact_persistence_emotionrec_drug;
                allControlEnergies_emotionid{7, contrast_idx+1} = controlImpact_persistence_emotionid_noDrug; allControlEnergies_emotionrec{7, contrast_idx+1} = controlImpact_persistence_emotionrec_noDrug;
                allControlTrajectories_emotionid{2, contrast_idx+1} = controlInputs_persistence_emotionid_drug; allControlTrajectories_emotionrec{2, contrast_idx+1} = controlInputs_persistence_emotionrec_drug;
                allControlTrajectories_emotionid{3, contrast_idx+1} = controlInputs_persistence_emotionid_noDrug; allControlTrajectories_emotionrec{3, contrast_idx+1} = controlInputs_persistence_emotionrec_noDrug;
                allControlTrajectories_emotionid{4, contrast_idx+1} = stateTrajectories_persistence_emotionid_drug; allControlTrajectories_emotionrec{4, contrast_idx+1} = stateTrajectories_persistence_emotionrec_drug;
                allControlTrajectories_emotionid{5, contrast_idx+1} = stateTrajectories_persistence_emotionid_noDrug; allControlTrajectories_emotionrec{5, contrast_idx+1} = stateTrajectories_persistence_emotionrec_noDrug;
                
                currentSaveFileName = strcat(resultsDir, 'numericalError_emotionid_', currentContrastLabel, '.mat');
                save(currentSaveFileName, 'numericalError_emotionid');
                
                currentSaveFileName = strcat(resultsDir, 'numericalError_emotionrec_', currentContrastLabel, '.mat');
                save(currentSaveFileName, 'numericalError_emotionrec');
            end
            
            allControlEnergies_emotionid{8, 2} = group_idx_emotionid; allControlEnergies_emotionrec{8, 2} = group_idx_emotionrec;
            
            save(strcat(resultsDir, 'allControlEnergies_emotionid.mat'), 'allControlEnergies_emotionid');
            save(strcat(resultsDir, 'allControlEnergies_emotionrec.mat'), 'allControlEnergies_emotionrec');
            save(strcat(resultsDir, 'allControlTrajectories_emotionid.mat'), 'allControlTrajectories_emotionid');
            save(strcat(resultsDir, 'allControlTrajectories_emotionrec.mat'), 'allControlTrajectories_emotionrec');
            save(strcat(resultsDir, 'structuralAdjacencyMatrix.mat'), 'A');
        end
    end
end