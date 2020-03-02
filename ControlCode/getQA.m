% script to calculate average streamline counts from multiple subject DTI data

function QA = getQA()
	subjectIDs = {'1342338080', '1904002815', '2411119678', '3258774044', '3788244813', '4484682335', '7715135642', '8201504234', '8567397981', '9325872444'};
	nSubjects = numel(subjectIDs);
	rootFolder = '../../DTI_data/connectivity_QA/connectivity/';
    %rootFolder = '/data/joy/BBL/studies/alpraz/DTI_data/connectivity_QA/connectivity/';
	nNodes = 234;
	QA = zeros(nNodes, nNodes);
	
	for i = 1:nSubjects
		currentSubjectID = subjectIDs{i};
		currentFilePath = strcat(rootFolder, currentSubjectID, filesep, currentSubjectID, '.src.gz.odf8.f5.gqi.1.25.fib.gz.ROIv_scale125_dilated.qa.end.connectivity.mat');
		load(currentFilePath, 'connectivity');
        %connectivity = connectivity/max(connectivity(:)); % scaling each subject's connectivity to [0 1]
		QA = QA + connectivity;
	end
	
	QA = QA/nSubjects;
end