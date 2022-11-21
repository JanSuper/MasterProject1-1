%Setup the options for the GA – tell it to plot the graph of the pareto front when the GA is run.
% Read the data for one stage.
%data = csvread('data.Uninflamed.csv',1,0);
% aSMA 6, Ecad 12, Cd20 13, CD16 18, pan keratin 20, cd68 32, Il17a 41, CD3 45
% Limit to the column you would like to see.
data_top = data(:,[6,12,13,18,20,32,41,45]);
f = @(gates)evaluate_gates4(gates,data_top);


options = optimoptions('gamultiobj','PlotFcn',@gaplotpareto);
lb=median(data_top);
ub=prctile(data_top,99);

%Run the GA.

%aSMA, Ecad, CD20, CD16, pankeratin, CD68, IL17a, CD3
%Fibroblasts, Epithelium, Bcells, Monocytes, Epithelium, Macrophages, IL17, T cells

[x,fval,e,o]=gamultiobj(f,8,[],[],[],[],lb,ub,options)
x

%To run the pareto search you use;
[xp,fvalp,ep,op]=paretosearch(f,8,[],[],[],[],lb,ub);

%In both cases the 7 indicates that we are optimizing 7 parameter

%obj1 = fval(:,1)

for i =1:length(fval(:,1))
    gates = x(i, :);
    cell_types(:,1)=(data_top(:,1)>gates(1));
    cell_types(:,2)=(data_top(:,2)>gates(2)| data_top(:,5)>gates(5));
    cell_types(:,3)=(data_top(:,3)>gates(3));
    cell_types(:,4)=(data_top(:,4)>gates(4));
    cell_types(:,5)=(data_top(:,6)>gates(6));
    cell_types(:,6)=(data_top(:,7)>gates(7) & data_top(:,1)<gates(1));
    cell_types(:,7)=(data_top(:,8)>gates(8));

    ct = zeros(1,10);

    for ii=1:size(cell_types,1)
        ct(8) = fval(i,1);
        ct(9) = fval(i,2);
        ct(10) = abs(fval(i,1)-fval(i,2));
        if sum(cell_types(ii,:)) == 1
            for c=1:size(cell_types,2)
                if cell_types(ii,c)
                    ct(c)=ct(c)+1;
                end
            end
        end
    end
    %ct_prop_all = 100*ct(1:8)/sum(ct(1:8));

    ct_prop = [100*ct(1:7)/sum(ct(1:7)) ct(10)];
    if i == 1
        resg = ct;
        resg_prop = [ct_prop gates];
        %resg_prop_all = [ct_prop_all gates];
    else
        resg = [resg; ct];
        resg_prop = [resg_prop; ct_prop gates];
        %resg_prop_all = [resg_prop_all; ct_prop_all gates];
    end
    resg = unique(resg, 'rows');
    resg_prop = unique(resg_prop, 'rows');
    %resg_prop_all = unique(resg_prop_all, 'rows');
end

resg_prop = sortrows(resg_prop,8);
tmp = resg_prop(1,8);
difff = resg_prop(:,8);

figure(1)
tiledlayout(1,2)
nexttile
bar(resg_prop(1, 1:7))
nexttile
plot(fvalp(:,1), fvalp(:,2),'m*')

%to_plot = res_prop(res_prop(:,7)>0.19 & res_prop(:,7)<0.2, :);
figure(2)
bar(resg_prop(1, 1:7));

%%%%%% USING THE LIST FROM UTRECHT 

for i =1:length(fval(:,1))
    gates = x(i, :);
    cell_types(:,1)=(data_top(:,1)>gates(1));
    cell_types(:,2)=(data_top(:,2)>gates(2)| data_top(:,5)>gates(5));
    cell_types(:,3)=(data_top(:,3)>gates(3));
    cell_types(:,4)=(data_top(:,4)>gates(4));
    cell_types(:,5)=(data_top(:,6)>gates(6));
    cell_types(:,6)=(data_top(:,7)>gates(7) & data_top(:,1)<gates(1));
    cell_types(:,7)=(data_top(:,8)>gates(8));

    ct = zeros(1,10);

    for ii=1:size(cell_types,1)
        %obj 1 - unassigned
        ct(8) = fval(i,1);
        %obj 2 - multiple
        ct(9) = fval(i,2);
        ct(10) = abs(fval(i,1)-fval(i,2));
        if sum(cell_types(ii,:)) == 1
            for c=1:size(cell_types,2)
                if cell_types(ii,c)
                    ct(c)=ct(c)+1;
                end
            end
        elseif sum(cell_types(ii,:)) > 1
            if cell_types(ii,6) ==1
                ct(6)=ct(6)+1;
            elseif cell_types(ii,7) ==1
                ct(7)=ct(7)+1;
            elseif cell_types(ii,3) ==1
                ct(1)=ct(3)+1;
            elseif cell_types(ii,4) ==1
                ct(1)=ct(4)+1;
            elseif cell_types(ii,5) ==1
                ct(1)=ct(5)+1;
            elseif cell_types(ii,1) ==1
                ct(1)=ct(1)+1;
            else
                ct(2)=ct(2)+1;
            end
        end
    end
    %ct_prop_all = 100*ct(1:8)/sum(ct(1:8));

    ct_prop = [100*ct(1:8)/sum(ct(1:8)) ct(8) ct(10)];
    if i == 1
        resg = ct;
        resg_prop = [ct_prop gates];
        %resg_prop_all = [ct_prop_all gates];
    else
        resg = [resg; ct];
        resg_prop = [resg_prop; ct_prop gates];
        %resg_prop_all = [resg_prop_all; ct_prop_all gates];
    end
    resg = unique(resg, 'rows');
    resg_prop = unique(resg_prop, 'rows');
    %resg_prop_all = unique(resg_prop_all, 'rows');
end

resg_prop = sortrows(resg_prop,10);
tmp = resg_prop(1,10);
difff = resg_prop(:,10);

figure(1)
tiledlayout(1,2)
nexttile
bar(resg_prop(1, 1:8))
nexttile
plot(fvalp(:,1), fvalp(:,2),'m*')

%to_plot = res_prop(res_prop(:,7)>0.19 & res_prop(:,7)<0.2, :);
%Point Rachel has selected for her presentation
figure(2)
bar(resg_prop(1, 1:8));
