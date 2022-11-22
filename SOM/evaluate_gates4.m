function [scores,cell_types] = evaluate_gates4(gates, data)
%%aSMA, Ecad, CD20, CD16, pankeratin, CD68, IL17a, CD3
%evaluate_gates - this is a fitness function for multiobjective
%optimisation.  It evaluates a single individual - so given a dataset and
%some gate thresholds it returns the percentage of cells unlabelled and
%percentage in multiple classes.
%Inputs
%   10 inputs one per marker (cd3, cd4, cd8, gdTCR, cd68, cd14, cd20,
%   ecadherin, pankeratin, asma) giving the threshold value for that
%   marker.
%Outputs
%   scores(1) - the percentage of cells unlabelled
%   scores(2) - the percentage of cells labelled as multiple classes
%   cell_types - the matrix of which cells are labelled as which type using
%   these gates.

%TO DO
% * Only do top level classification on first pass - then next level on a
% secon pass - done!
% * Calculate upper and lower bounds from data - 50%-99%? - median and
% prctile will calculate these from data.
% * Adjust data so it only has columns needed for this task and columns are
% in order corresponding to gates - done

%aSMA, Ecad, CD20, CD16, pankeratin, CD68, IL17a, CD3
%Fibroblasts, Epithelium, Bcells, Monocytes, Epithelium, Macrophages, IL17, T cells

%Calculate the cell types using the gates
cell_types(:,1)=(data(:,1)>gates(1));
%cell_types(:,2)=(data(:,1)>gates(1) & data(:,2)>gates(2));
%cell_types(:,3)=(data(:,1)>gates(1) & data(:,3)>gates(3));
%cell_types(:,4)=(data(:,1)>gates(1) & data(:,4)>gates(4));
cell_types(:,2)=(data(:,2)>gates(2)| data(:,5)>gates(5));
cell_types(:,3)=(data(:,3)>gates(3));
cell_types(:,4)=(data(:,4)>gates(4));
cell_types(:,5)=(data(:,6)>gates(6));
cell_types(:,6)=(data(:,7)>gates(7) & data(:,1)<gates(1));
cell_types(:,7)=(data(:,8)>gates(8));

%look at numbers of cells unlabelled and in overlaps.
scores(1)=(sum(sum(cell_types,2)==0))./size(cell_types,1);
%top_level_cell_types=[1,5,6,7,8,9];
scores(2)=1 - (sum(sum(cell_types,2)==1))./size(cell_types,1);
%{
%modif arthur
scores(1)= 1 - ((sum(sum(cell_types,2)==1))./size(cell_types,1));
if scores >.4
    scores(2)=(sum(sum(cell_types,2)>1))./size(cell_types,1);
else
    center = net.IW{1};
    % cell that are label only ones
    % distance between same cell type
    ald = 0;
    for i = 1 : 7
        %size(dataORGINAL(cell_types(:,i), :))
        yp = net(dataORGINAL(cell_types(:,i), :)');
        [argvalue, argmax] = max(yp);
        c = center(argmax);
        ald = ald + mean(pdist(c', 'minkowski'))/m;
    end
    scores(2) = ald/7;
    scores(2)
end
%}