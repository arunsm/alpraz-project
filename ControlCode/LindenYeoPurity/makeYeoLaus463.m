%%
clear all; close all; clc
addpath(genpath('/Applications/freesurfer/matlab'));
savedir = '~/Dropbox/Cornblath_Bassett_Projects/code/shared_code/LindenYeoPurity/';

[YLv,YLL,YLct] = read_annotation('lh.Yeo2011_7Networks_N1000.annot');
[YRv,YRL,YRct] = read_annotation('rh.Yeo2011_7Networks_N1000.annot');

lausanneScale = 250;

[Lv,LL,Lct] = read_annotation(['lh.myaparc_',num2str(lausanneScale),'.annot']);
[Rv,RL,Rct] = read_annotation(['rh.myaparc_',num2str(lausanneScale),'.annot']);

%%

RnodeAssignments = zeros((Rct.numEntries),7);

for N = 1:(Rct.numEntries)
    for Y = 1:7
        RnodeAssignments(N,Y) = sum(YRL(RL == Rct.table(N,5)) == YRct.table(Y+1,5));
    end
end

LnodeAssignments = zeros((Lct.numEntries),7);

for N = 1:(Lct.numEntries)
    for Y = 1:7
        LnodeAssignments(N,Y) = sum(YLL(LL == Lct.table(N,5)) == YLct.table(Y+1,5));
    end
end

%% calculate purity of each label

[~,LYeoLabels] = max(LnodeAssignments,[],2);
[~,RYeoLabels] = max(RnodeAssignments,[],2);

LnodeProportion = LnodeAssignments ./ sum(LnodeAssignments,2);
LnodePurity = zeros(length(LYeoLabels),1);
for i = 1:length(LYeoLabels)
    LnodePurity(i) = LnodeProportion(i,LYeoLabels(i));
end

RnodeProportion = RnodeAssignments ./ sum(RnodeAssignments,2);
RnodePurity = zeros(length(RYeoLabels),1);

for i = 1:length(RYeoLabels)
    RnodePurity(i) = RnodeProportion(i,RYeoLabels(i));
end
%% get labels by node index

load human_regionNames.mat
if lausanneScale == 125
    lausnames = roinames{3};
elseif lausanneScale == 250
    lausnames = roinames{4};
elseif lausanneScale == 60
    lausnames = roinames{2};
elseif lausanneScale == 36
    lausnames = roinames{1};
end

Ls = size(Lct.struct_names,1); Rs = size(Rct.struct_names,1);

for i = 1:Ls
    Lct.struct_names{i} = ['lh_',Lct.struct_names{i}];
end

for i = 1:Rs
    Rct.struct_names{i} = ['rh_',Rct.struct_names{i}];
end

[~,lind] = ismember(Lct.struct_names,lausnames);
[~,rind] = ismember(Rct.struct_names,lausnames);

finalLabels = zeros((Ls + Rs - 4 + 15),1);

for i = 1:Ls
    if lind(i) ~= 0
        finalLabels(lind(i)) = LYeoLabels(i);
    end
end

for i = 1:Rs
    if rind(i) ~= 0
        finalLabels(rind(i)) = RYeoLabels(i);
    end
end

finalLabels(finalLabels == 0) = 8;

load([savedir,'yeo7netlabelsLaus125.mat']);
if lausanneScale == 125
    f=figure; imagesc([finalLabels,network7labels]); 
    xticks([1:2]);xticklabels({'Eli','Graham'});
    saveas(f,fullfile(savedir,'EliVsGrahamLaus125.pdf'));
end

cd(savedir);
save(['yeo7netlabelsLaus',num2str(lausanneScale),'EJC.mat'],'finalLabels','names');

%% 
%
figure; 
subplot(1,2,1); imagesc(LnodeAssignments);
xlabel('Cognitive System'); ylabel('Node'); title('L');
colorbar('southoutside');
subplot(1,2,2); imagesc(RnodeAssignments);
xlabel('Cognitive System'); ylabel('Node'); title('R');
colorbar('southoutside');
figure;
subplot(1,2,1); imagesc(LnodePurity);
h = colorbar('southoutside'); ylabel(h,'Purity');
subplot(1,2,2); imagesc(RnodePurity);
h = colorbar('southoutside'); ylabel(h,'Purity');

%}