% heuristic to create a cohort file to feed into xcpEngine
% format:
% id0, id1, id2, img, task_design
% subjectID, sessionID, taskID, path_to_preprocessed_BOLD_image, path_to_design_file

baseDir_stimDesigns = '/data/joy/BBL/studies/alpraz/rawData/derivatives/stimDesigns/';

% collecting all fmriprep output
baseDir_fmriprep = '/data/joy/BBL/studies/alpraz/rawData/derivatives/fmriprep/';
d_sub = dir(strcat(baseDir_fmriprep, filesep, 'sub*.html'));
fnme_sub = {d_sub.name};
nSubjects = numel(fnme_sub);

QA = xlsread('/data/joy/BBL/studies/alpraz/rawData/derivatives/stimDesigns/AlprazQA_forArun.xlsx');
form_info = [QA(:, 2) QA(:, 3) QA(:, 7)]; % look-up table with [subjectID, sessionID, form(A0,B1)]

fid = fopen("cohortFile.csv", 'a');
fprintf(fid, '%s,%s,%s,%s,%s\n', 'id0', 'id1', 'id2', 'img', 'task_design');

for i = 1:nSubjects
    currentSubject = strrep(fnme_sub{i}, '.html', '');
    fprintf(currentSubject); fprintf('\n');
    currentSubjectID = str2num(extractAfter(currentSubject, 'sub-'));
    d_ses = dir(strcat(baseDir_fmriprep, filesep, currentSubject, filesep, 'ses*'));
    fnme_ses = {d_ses.name};
    nSessions = numel(fnme_ses);
    
    for j = 1:nSessions
        currentSession = fnme_ses{j};
        fprintf(currentSession); fprintf('\n');
        currentSessionID = str2num(extractAfter(currentSession, 'ses-'));
        d_BOLD_files = dir(strcat(baseDir_fmriprep, filesep, currentSubject, filesep, currentSession, filesep, 'func', filesep, '*preproc_bold.nii.gz'));
        fnme_BOLD_files = {d_BOLD_files.name};
        nBOLD_files = numel(fnme_BOLD_files);
        
        for k = 1:nBOLD_files
            currentFile = fnme_BOLD_files{k};
            if contains(currentFile, 'task-emotionid')
                currentTaskType = 'task-emotionid';
                path_to_design_template = strcat(baseDir_stimDesigns, 'emotionid_design_template.fsf');
            elseif contains(currentFile, 'task-emotionrec')
                currentTaskType = 'task-emotionrec';
                path_to_design_template = strcat(baseDir_stimDesigns, 'emotionrec_design_template.fsf');
            end
            
            % extracting motion parameters and storing as text file
            path_to_confounds = strcat(baseDir_fmriprep, filesep, currentSubject, filesep, currentSession, ...
                    filesep, 'func', filesep, currentSubject, '_', currentSession, '_', currentTaskType, '_desc-confounds_regressors.tsv');
            allConfounds = tdfread(path_to_confounds);
            motionParameters = [allConfounds.trans_x allConfounds.trans_y allConfounds.trans_z allConfounds.rot_x allConfounds.rot_y allConfounds.rot_z];
            path_to_motionParameters = strcat(baseDir_fmriprep, filesep, currentSubject, filesep, currentSession, ...
                    filesep, 'func', filesep, currentSubject, '_', currentSession, '_', currentTaskType, '_motionParameters.txt');
            dlmwrite(path_to_motionParameters, motionParameters);
                
            % creating subject-specific design file
            fid_template = fopen(path_to_design_template, 'r');
            f = fread(fid_template,'*char')';
            fclose(fid_template);
            f = strrep(f, 'sessionID', strcat('0', num2str(currentSessionID))); % replacing the placeholder 'sessionID' with actual session ID padded with one zero
            f = strrep(f, 'subjectID_char', currentSubject); % replacing the placeholder 'subjectID_char' with actual subjectID
            f = strrep(f, 'sesID_char', currentSession); % replacing the placeholder 'sessionID_char' with actual sessionID
            f = strrep(f, 'taskID', currentTaskType); % replacing the placeholder 'taskID' with actual taskID
            currentDesignFilePath = strrep(path_to_design_template, 'template', currentSession);
            fid_design = fopen(currentDesignFilePath, 'w');
            fprintf(fid_design, '%s', f);
            fclose(fid_design);
            
            % assigning appropriate design form - use this section if assigning uniform design files (form A/B)
            %if sum(form_info(:, 1) == currentSubjectID) == 0
            %continue; % skipping if subject ID not found in QA list
            %elseif form_info(form_info(:, 1) == currentSubjectID & form_info(:, 2) == currentSessionID, 3) == 0
            %form = 'formA';
            %elseif form_info(form_info(:, 1) == currentSubjectID & form_info(:, 2) == currentSessionID, 3) == 1
            %form = 'formB';
            %end
            
            currentFilePath = strcat(baseDir_fmriprep, currentSubject, filesep, currentSession, filesep, 'func', filesep, currentFile);
            fprintf(fid, '%s,%s,%s,%s,%s\n', currentSubject, currentSession, currentTaskType, currentFilePath, currentDesignFilePath);
        end
    end
end

fclose(fid);
