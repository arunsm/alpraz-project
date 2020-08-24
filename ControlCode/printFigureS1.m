function [] = printFigureS1()

parameters
resultsDirCurrentFigure = strcat(resultsDir, filesep, 'FigureS1', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

M = importdata('../../data/slabCoverage_combinedMaskStd_emotionid.txt');
slabCoverage = M.data';
slabCoverage = double(slabCoverage > parcelCoverageThresh);

surfacePlots(slabCoverage, parula, [0 1], nifti, subcorticalIndices, resultsDirCurrentFigure, strcat('slabCoverage.svg'));

end