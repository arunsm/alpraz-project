pathToQualityFile = '/Users/ArunMahadevan/Documents/BBL/studies/alpraz/derivatives/xcp_output_allTasks2/group/n191_quality.csv';
Q = readtable(pathToQualityFile);

fid = fopen("subjectMotion.csv", 'w');
fprintf(fid, '%s,%s,%s,%s\n', 'id0', 'id1', 'id2', 'avge_FD');

nVols2Skip = 6; % specify number of non-steady state volumes to skip

totalSubjects_and_Sessions = size(Q, 1);

for i = 1:totalSubjects_and_Sessions
    currentSubjectID = Q.id0{i}; fprintf('%s\n', currentSubjectID);
    currentSessionID = Q.id1{i}; fprintf('%s\n', currentSessionID);
    currentTaskID = Q.id2{i}; fprintf('%s\n', currentTaskID);
    
    currentPath = strcat('/Users/ArunMahadevan/Documents/BBL/studies/alpraz/derivatives/fmriprep/', ...
        currentSubjectID, filesep, currentSessionID, filesep, 'func/', currentSubjectID, '_', currentSessionID, '_', currentTaskID, '_desc-confounds_regressors.tsv');
    
    confounds = readmatrix(currentPath, 'FileType', 'text');
    FD = confounds(nVols2Skip+1:end, 6);
    avge_FD = mean(FD);
    
    fprintf(fid, '%s,%s,%s,%d\n', currentSubjectID, currentSessionID, currentTaskID, avge_FD);
end