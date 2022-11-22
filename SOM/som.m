close all;

% Create a Self-Organizing Map
load('Workspace_Script_UtrechtListmat.mat')
Name = ["ROInr","Cell_Area","Nucl_Area","X80ArAr_Mean","sSMA113_Mean","aSMA_Mean","X124Xe_Mean","X127I_Mean","X129Xe_Mean","X132Xe_Mean","RORgt_Mean","Ecad_Mean","CD20_Mean","GATA3_Mean","GATA3_Nucl_Mean","Tbet_Mean","Tbet_Nucl_Mean","CD16_Mean","BetaCatenin_Mean","PanKeratin_Mean","CTLA4_Mean","PDL1_Mean","IFNg_Mean","CD45RO_Mean","AKT_Mean","HLA_DR-DP-DQ_Mean","FOXP3_Mean","FOXP3_Nucl_Mean","CD4_Mean","CD103_Mean","pSTAT3_Mean","CD68_Mean","IL10_Mean","CD45_Mean","CD8a_Mean","ICOS_Mean","pS6_Mean","PD1_Mean","NFkB_Mean","NFkB_Nucl_Mean","IL17a_Mean","Ki67_Mean","Ki67_Nucl_Mean","GranzymeB_Mean","CD3_Mean","pERK_Mean","Cleaved Caspase 3_Mean","ERK_Mean","TCRgd_Mean","pAKT_Mean","H3_Mean","H3_Nucl_Mean","Ir193_Mean","Ir193_Nucl_Mean"];
r = [];
for i = 1: 54
    if contains(Name(i),"CD45")
        r = [r, i];
    end
end
%remove CD45 DATA
for i = 1 : size(r,2)
    data(:, r(i)-(i-1)) = [];
end

dimension1 = 5;
dimension2 = 5;
net = selforgmap([dimension1 dimension2], 500, 0);

% Train the Network
batch = data(randperm(size(data, 1)), :);
[net,tr] = train(net,batch(1:round(size(data,1)/5), :)');

% Test the Network

% View the Network

%figure, plotsomplanes(net)
center = net.IW{1};
id = 1;
for i = 1 : size(center, 1)
    for y = 1 : size(center, 1)
        dis(id) = pdist([center(i, :); center(y, :)], 'minkowski');
        id = id + 1;
    end
end
m = max(dis);
bestScore = 0;
for i =1:length(fval(:,1))
    gates = x(i, :);

    cell_types(:,1)=(data_top(:,1)>gates(1));
    cell_types(:,2)=(data_top(:,2)>gates(2)| data_top(:,5)>gates(5));
    cell_types(:,3)=(data_top(:,3)>gates(3));
    cell_types(:,4)=(data_top(:,4)>gates(4));
    cell_types(:,5)=(data_top(:,6)>gates(6));
    cell_types(:,6)=(data_top(:,7)>gates(7) & data_top(:,1)<gates(1));
    cell_types(:,7)=(data_top(:,8)>gates(8));

    s = (sum(sum(cell_types,2)==1))./size(cell_types,1);
    if bestScore < s
        bestScore = s;
        bestGates = x(i, :);
    end
end
ald = 0;

gates = bestGates;
cell_types(:,1)=(data_top(:,1)>gates(1));
cell_types(:,2)=(data_top(:,2)>gates(2)| data_top(:,5)>gates(5));
cell_types(:,3)=(data_top(:,3)>gates(3));
cell_types(:,4)=(data_top(:,4)>gates(4));
cell_types(:,5)=(data_top(:,6)>gates(6));
cell_types(:,6)=(data_top(:,7)>gates(7) & data_top(:,1)<gates(1));
cell_types(:,7)=(data_top(:,8)>gates(8));

for i = 1 : 7
    %size(dataORGINAL(cell_types(:,i), :))
    yp = net(data(cell_types(:,i)==0, :)');
    [argvalue, argmax] = max(yp);
    c = center(unique(argmax), :);
    ald = ald + mean(pdist(c', 'minkowski'))/m;
end
t = ["Fibroblasts", "Epithelium", "Bcells", "Monocytes", "Macrophages", "IL17", "T cells", "Others"];
for i=1: 7
    figure, plotsomhits(net, data(cell_types(:,i), :)')
    title(t(i));
end