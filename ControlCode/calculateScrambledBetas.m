function xf = calculateScrambledBetas(x, savePath)

%% add spin-test toolbox to path
addpath(genpath('~/matlab/spin-test-master'))
%addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/spin-test-master'));
nPermutations = 1;

%% call function to convert node values to surface data
[lh_surf_data, rh_surf_data] = convertData2LausanneSurface(x);

%% confirm that surface plots look as expected
% [f, ax, ph, ~] = plot_surface_file(x, lh_surf_data, rh_surf_data, jet, [0 8]); close(f(2));
% view(ax(1), [270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
% axis(ax(1), 'off'); 

%% perform spin test
SpinPermuFS(lh_surf_data, rh_surf_data, nPermutations, savePath);
load(savePath, 'bigrotl', 'bigrotr');

%% confirm that surface plots for scrambled data look as expected
% load(savePath, 'bigrotl', 'bigrotr');
% [f, ax, ph, ~] = plot_surface_file(x, bigrotl(nPermutations, :)', bigrotr(nPermutations, :)', jet, [0 8]); close(f(2));
% view(ax(1), [270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
% axis(ax(1), 'off'); 

%% return scrambled node values from last permutation
xf = convertLausanneSurface2Data(bigrotl(nPermutations, :)', bigrotr(nPermutations, :)');

end