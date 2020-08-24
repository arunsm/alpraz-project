%% loading relevant files

clear all; clc;
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/LausanneSurfaceFigures'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/LausanneCoordinates'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/mult_comp_perm_t1'));

parameters % run script to set parameters

allControlEnergies_emotionid = readtable(strcat(resultsDir, 'allControlEnergies_emotionid.csv'));
allControlEnergies_emotionrec = readtable(strcat(resultsDir, 'allControlEnergies_emotionrec.csv'));
load(strcat(resultsDir, 'allControlTrajectories_emotionid.mat'));
load(strcat(resultsDir, 'allControlTrajectories_emotionrec.mat'));
load(strcat(resultsDir, 'structuralAdjacencyMatrix.mat'));

%% table 1
printTable1()

%% figure 1 - structural matrix heatmap, brain maps for sample subject, Lausanne parcellation for schematic
printFigure1(A, allControlTrajectories_emotionid)

%% figure 2 - box plots of persistence energy, control impact brain maps, sorted Lausanne parcels of control impact
printFigure2(allControlEnergies_emotionid, allControlEnergies_emotionrec, ...
    allControlTrajectories_emotionid, allControlTrajectories_emotionrec)

%% figure 3 - persistence energy vs task efficiency
printFigure3(allControlEnergies_emotionid, allControlEnergies_emotionrec)

%% figure 4 - spatial maps of drug differences in control input, PET maps, correlations between PET maps and drug differences
printFigure4(allControlTrajectories_emotionid, allControlTrajectories_emotionrec)

%% Supplementary Data Files 1 -  mixed models with categorical variables
printSupplementaryDataFiles1(allControlEnergies_emotionid, allControlEnergies_emotionrec)

%% Supplementary Data Files 2 - mixed models with alpraz_levels
printSupplementaryDataFiles2(allControlEnergies_emotionid, allControlEnergies_emotionrec)

%% Supplementary Data Files 3 - mixed models with SISTOTAL
printSupplementaryDataFiles3(allControlEnergies_emotionid, allControlEnergies_emotionrec)

%% Supplementary Data Files 5 - mixed models for efficiency vs PE
printSupplementaryDataFiles5(allControlEnergies_emotionid, allControlEnergies_emotionrec)

%% figure S1 - brain map of slab coverage
printFigureS1()

%% figure S2
printFigureS2()

%% figure S3 - brain maps of GLM betas, sorted Lausanne parcels of GLM beta values
printFigureS3(allControlTrajectories_emotionid, allControlTrajectories_emotionrec)

%% figure S4 - brain maps of control inputs, sorted Lausanne parcels of control inputs
printFigureS4(allControlTrajectories_emotionid, allControlTrajectories_emotionrec)

%% figure S5 - 
printFigureS5()
