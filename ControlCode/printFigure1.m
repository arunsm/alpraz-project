function [] = printFigure1(A, allControlTrajectories_emotionid)

parameters
resultsDirCurrentFigure = strcat(resultsDir, filesep, 'Figure1', filesep);
if ~exist(resultsDirCurrentFigure)
    mkdir(resultsDirCurrentFigure)
end

%% plot heatmap of Q7 structural matrix
f = figure; set(gcf, 'color', 'w'); hold on;
imagesc(A);
colormap('parula');
colorbar;
axis off;
ax = gca;
ax.FontSize = 12;
savePath = strcat(resultsDirCurrentFigure, 'Q7_averageStructuralMatrix_QA.svg');
saveas(f, savePath); close(f);

%% figure depicting Lausanne parcellation
parcelValues = zeros(234, 1);
[f, ax, ph, ~] = fcn_lausannesurf(parcelValues, white, [-20 100]);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDirCurrentFigure, 'LausanneParcellation1.svg');
saveas(f(1), savePath); close(f(1)); close(f(2)); 

%% plot brain activation maps w/ and w/o alpraz for individual subject
brainStates_alpraz = zeros(234, 1);
xf = allControlTrajectories_emotionid.xf{2};
parcelsToInclude_idx = allControlTrajectories_emotionid.parcelsToInclude_idx{1};
brainStates_alpraz(parcelsToInclude_idx) = xf;

[f, ax, ph, ~] = fcn_lausannesurf(brainStates_alpraz, redbluecmap);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDirCurrentFigure, 'sampleBrainStates_Placebo.svg');
saveas(f(1), savePath); close(f(1)); close(f(2));

brainStates_alpraz = zeros(234, 1);
xf = allControlTrajectories_emotionid.xf{1};
parcelsToInclude_idx = allControlTrajectories_emotionid.parcelsToInclude_idx{1};
brainStates_alpraz(parcelsToInclude_idx) = xf;

[f, ax, ph, ~] = fcn_lausannesurf(brainStates_alpraz, redbluecmap);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDirCurrentFigure, 'sampleBrainStates_Alpraz.svg');
saveas(f(1), savePath); close(f(1)); close(f(2));
end