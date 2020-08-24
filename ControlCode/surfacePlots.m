% function to generate and save 4 cortical and 2 subcortical plots for scalar
% metrics

function [] = surfacePlots(metricToPlot, cmap, clims, nifti, subcorticalIndices, resultsDir, savePathBaseName)

% cortex - left hemisphere (2 plots)
[f, ax, ph, ~] = fcn_lausannesurf(metricToPlot, cmap, clims);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDir, savePathBaseName, '_cortex1.svg'); saveas(f(1), savePath); close(f(1));
savePath = strcat(resultsDir, savePathBaseName, '_cortex2.svg'); saveas(f(2), savePath); close(f(2)); 

% cortex - right hemisphere (2 plots)
[f, ax, ph, ~] = fcn_lausannesurf(metricToPlot, cmap, clims);
view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [270, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDir, savePathBaseName, '_cortex3.svg'); saveas(f(1), savePath); close(f(1));
savePath = strcat(resultsDir, savePathBaseName, '_cortex4.svg'); saveas(f(2), savePath); close(f(2));

% subcortex (2 plots)
figure; [f1, f2] = plot_subcortvol(metricToPlot(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, cmap, clims(1), clims(2));
savePath = strcat(resultsDir, savePathBaseName, '_subcortex1.svg'); saveas(f1, savePath); close(f1);
savePath = strcat(resultsDir, savePathBaseName, '_subcortex2.svg'); saveas(f2, savePath); close(f2);

end