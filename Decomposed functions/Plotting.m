%% Plot



grey = [0.4,0.4,0.4];


%Plot all demand / production profiles
for i=1:Noutputs
figure()
title(strcat(Outputs{i},' Demand / Production profiles'))
xlabel('Timestep')
ylabel('Power [kW]')
hold on
stdnam=[{'Storage Discharge'};{'Purchased from network'};{'Storage Charge'}; {'Sold on network'}; {'Dissipated'}];
prodnames=[Pmat{i,2}' Pmat{i,10}']';
stdval={Pmat{i,4}' Pmat{i,6}'};
k = bar([Pmat{i,3}' Pmat{i,11}' stdval{L{i}(1:2)}],0.6,'stacked');   
stdval={Pmat{i,5}' -Pmat{i,7}' -Pmat{i,8}' -Pmat{i,13}'};
h = bar([stdval{L{i}(3:5)} -Pmat{i,13}'],0.6,'stacked');
plot(Dall{2}(i,:),'k','LineWidth',2)
names=[prodnames ;stdnam(L{i}(1:5)); Pmat{i,12}];
color=parula(max(size(prodnames,1),8));
%color=copper(size(prodnames,1)+10);
for j=1:size(prodnames,1)
    set(k(j),'facecolor',color(j,:));
end
if L{i}(2)~=0
    set(k(end),'facecolor','y')
end
if L{i}(1)~=0
    set(k(end-L{i}(2)),'facecolor','g')
end
if L{i}(3)~=0
    set(h(1),'facecolor','r')
end
if L{i}(4)~=0
    set(h(1+L{i}(3)),'facecolor',[1 .5 0])
end
if L{i}(5)~=0
    set(h(L{i}(3)+L{i}(4)+L{i}(5)),'facecolor',grey)
end

color=cool(max(size(Pmat{i,12},1),7));
for j=1:size(Pmat{i,12},1)
    set(h(L{i}(3)+L{i}(4)+L{i}(5)+j),'facecolor',color(j,:))
end
legend(gca,names)
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
hold off
end

fig=i;
nst=0;

%Plot all storages charge / discharge profiles
for i=1:Noutputs
if ~isempty(Pmat{i,9})
fig=fig+1;
nst=nst+1;
nst=find(ismember({Storages{:,1}},Outputs(i)));
figure()
hold on
yyaxis left
h=bar([Pmat{i,9}' -Pmat{i,14}'],1,'stacked');               %ESS energy level 
set(h(1),'facecolor','b')
if L{i}(7)~=0
    set(h(2),'facecolor','y')
end
axis([0,ntimes,-500,Storages{nst,2}*Storages{nst,8}/100*1.1])
title(strcat(Outputs{i},' Storage level and charge/discharge profiles'))
xlabel('Timestep') % x-axis label
ylabel('Storage energy content') % y-axis label
set(gca,'XTick', 0:4:24*days);
yyaxis right
ylabel('Power injected (green) / withdrawn (red) from storage')
h=bar([-Pmat{i,5}' Pmat{i,15}'],0.6,'stacked');               %ESS charge / discharge 
set(h(1),'facecolor','g')
if L{i}(8)~=0
    set(h(2),'facecolor',[0 .5 0])
end
h=bar([-Pmat{i,4}' -Pmat{i,16}'],0.6,'stacked');               %ESS charge / discharge 
set(h(1),'facecolor','r')
if L{i}(9)~=0
    set(h(2),'facecolor',[0.5 0 0])
end
axis([0,length(Pmat{1,4}),-500,Storages{nst,2}*Storages{nst,8}/100*1.1])
yyaxis left
legnames=[{'Storage Charge Level'};{'Storage self-discharge'};{'Storage charge'};{'Storage charge loss'};{'Storage discharge'};{'Storage discharge loss'}];
logictag=logical([1 L{i}(7) 1 L{i}(8) 1 L{i}(9)]);
legend(gca,legnames(logictag))

% ax(1)=gca                   %da sistemare poich� non funzica
% yyaxis left
% ax(2)=gca
% linkaxes(ax,'y')

plot(get(gca,'xlim'), [Storages{nst,2}*Storages{nst,9}/100 Storages{nst,2}*Storages{nst,9}/100],'--k','Linewidth',2)
plot(get(gca,'xlim'), [Storages{nst,2}*Storages{nst,8}/100 Storages{nst,2}*Storages{nst,8}/100],'--k','Linewidth',2)
% linkaxes([h k p],'y')
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
hold off
end
end


%Objective=sum(sum([Fuels{:,2}]'.*(fuelusage.*timestep)))+sum(sum([Networks{:,4}]'.*(NETWORKbought.*timestep)))-sum(sum([Networks{:,5}]'.*(NETWORKsold.*timestep)));



figure()
hold on
h=bar([costvals{:}],1,'stacked');              %ESS energy level 
color=parula(max(size(costvals,2),8));
%color=copper(size(prodnames,1)+10);
for j=1:size([costvals{:}],2)
    set(h(j),'facecolor',color(j,:));
end
if can~=0
    k=bar(-[gainvals{:}],1,'stacked');
    color=hot(max(size(gainvals,2),8));
    for j=1:size([gainvals{:}],2)
        set(k(j),'facecolor',color(j,:));
    end
end
title('Cost function explosion')
xlabel('Timestep') % x-axis label
ylabel('Timestep cost') % y-axis label
set(gca,'XTick', 0:4:24*days);
plot(Obj,'k','LineWidth',2)
tags{end+1}='Overall cost function';
legend(gca,[tags{:}])
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
hold off
