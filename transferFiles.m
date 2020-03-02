baseDir_stimDesigns = '/data/joy/BBL/studies/alpraz/rawData/derivatives/stimDesigns/';

% collecting all fmriprep output
baseDir_fmriprep = '/data/joy/BBL/studies/alpraz/rawData/derivatives/fmriprep/';
d_sub = dir(strcat(baseDir_fmriprep, filesep, 'sub*.html'));
fnme_sub = {d_sub.name};
nSubjects = numel(fnme_sub);

QA = xlsread('/data/joy/BBL/studies/alpraz/rawData/derivatives/stimDesigns/AlprazQA_forArun.xlsx');

for i = 1:nSubjects
    currentSubject = strrep(fnme_sub{i}, '.html', '');
    fprintf(currentSubject); fprintf('\n');
    d_ses = dir(strcat(baseDir_fmriprep, filesep, currentSubject, filesep, 'ses*'));
    fnme_ses = {d_ses.name};
    nSessions = numel(fnme_ses);
    
    for j = 1:nSessions
        currentSession = fnme_ses{j};
        fprintf(currentSession); fprintf('\n');
        d_BOLD_files = dir(strcat(baseDir_fmriprep, filesep, currentSubject, filesep, currentSession, filesep, 'func', filesep, '*preproc_bold.nii.gz'));
        fnme_BOLD_files = {d_BOLD_files.name};
        nBOLD_files = numel(fnme_BOLD_files);
        
        for k = 1:nBOLD_files
            currentFile = fnme_BOLD_files{k};
            if contains(currentFile, 'task-emotionid')
                currentTaskType = 'task-emotionid';
            elseif contains(currentFile, 'task-emotionrec')
                currentTaskType = 'task-emotionrec';
            end
            
            currentSource = strcat(currentSubject, '_', currentSession, '_', currentTaskType, '_desc-confounds_regressors.tsv');
            currentDestination = strcat(baseDir_fmriprep, currentSubject, filesep, currentSession, filesep, 'func/');
            
            fprintf('moving file %s\n', currentSource);
            
            movefile currentCopyPath currentDestination;
        end
    end
end