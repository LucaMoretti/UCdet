%% Plot

%Prepare data for each produced good
%Pmat: 
%column1 --> Output name
%column2 --> Tags of dispatchable machines producing that good
%column3 --> Output profile of each dispatchable machine
%column4 --> Storage discharging power profile
%column5 --> Storage charging power profile
%column6 --> Purchased from network
%column7 --> Sold on network
%column8 --> Dissipated
%column9 --> Storage content
%column10 --> Tags of dispatchable machines producing that good
%column11 --> Output profile of each dispatchable machine
%column12 --> Utility consumers tag
%column13 --> Utility consumers consumption profile
nstor=0;
nnet=0;
Pmat=cell(Noutputs,16);
for i=1:Noutputs
    nprod=0;
    Pmat{i,1}=Outputs{i};
    f=0;
    l=0;
    for j=1:Nmachines
        %Production by machine i of good j
        for h=1:numel(Machines{j,3})
            if isequal(Outputs(i),Machines{j,3}(h)) && sum(value(OUTPUT{j}(h,:)))>1e-6
                if f==0
                    Pmat{i,2}=Machines{j,1};
                    Pmat{i,3}=[value(OUTPUT{j}(h,:))];
                else
                    Pmat{i,2}=[Pmat{i,2};Machines{j,1}];
                    Pmat{i,3}=[Pmat{i,3};value(OUTPUT{j}(h,:))];
                end
                f=f+1;
            end
        end
        %Consumption by machine i of good j
        for h=1:numel(Machines{j,2})
            if isequal(Outputs(i),Machines{j,2}(h)) && sum(value(INPUT{j}(h,:)))~=0
                if l==0
                    Pmat{i,12}=Machines{j,1};
                    Pmat{i,13}=[value(INPUT{j}(h,:))];
                else
                    Pmat{i,12}=[Pmat{i,12};Machines{j,1}];
                    Pmat{i,13}=[Pmat{i,13};value(INPUT{j}(h,:))];
                end
                l=l+1;
            end
        end
    end
    f=0;
    for j=1:Nundisp
        %Production by undispatchable machine i of good j
        for h=1:numel(UndProd{j,2})
            if isequal(Outputs(i),UndProd{j,2}(h)) && sum(UndProd{j,3}(h,:))~=0
                if f==0
                    Pmat{i,10}=UndProd{j,1};
                    Pmat{i,11}=UndProd{j,3}(h,:);
                else
                    Pmat{i,10}=[Pmat{i,10};UndProd{j,1}];
                    Pmat{i,11}=[Pmat{i,11};UndProd{j,3}(h,:)];
                end
                f=f+1;
            end
        end
    end
    C=logical(zeros(9,1));
    Pmat{i,4}=zeros(1,ntimes);
    Pmat{i,5}=zeros(1,ntimes);
    Pmat{i,6}=zeros(1,ntimes);
    Pmat{i,7}=zeros(1,ntimes);
    Pmat{i,8}=zeros(1,ntimes);
    if ismember(Outputs{i},{Storages{:,1}}) %strcmp({Storages{:,1}},Outputs{i}))>0 
        nstor=ismember({Storages{:,1}},Outputs{i});
        C(1)=1;                                 %C: 1=Sdisch 2=NetPurch 3=Stcharge 4=NetSold 5=Dissipation 6=Storagelevel
        C(3)=1;                                 % 7=self disch 8=charge loss 9=discharge loss
        C(6)=1;
        if Storages{nstor,7}~=0
            C(7)=1;
            Pmat{i,14}=value(STORAGEcharge(nstor,:))*Storages{nstor,7};
        end
        Storprof=value(STORAGEpower(nstor,:));
        Storprof(Storprof>=0)=0;
        Pmat{i,4}=value(STORAGEpower(nstor,:))-Storprof;
        Pmat{i,5}=Storprof;
        Pmat{i,9}=value(STORAGEcharge(nstor,:));
        if Storages{nstor,4}~=1
            C(8)=1;
            Pmat{i,15}=Pmat{i,5}*(1-Storages{nstor,4});
        end
        if Storages{nstor,6}~=1
            C(9)=1;
            Pmat{i,16}=Pmat{i,4}*(1-Storages{nstor,6});
        end
    end
    if sum(strcmp({Networks{:,1}},Outputs{i}))>0
        nnet=nnet+1;
        if sum(value(NETWORKbought(nnet,:)))>1e-6
            C(2)=1;
            Pmat{i,6}=value(NETWORKbought(nnet,:));
        end
        if sum(value(NETWORKsold(nnet,:)))>1e-6
            C(4)=1;
            Pmat{i,7}=value(NETWORKsold(nnet,:));
        end
end
    

    Pmat{i,8}=value(Diss(i,:));
    if abs(sum(value(Diss(i,:))))>0.001
        C(5)=1;
    end
    L{i}=C;

end

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
plot(D{2}(i,:),'k','LineWidth',2)
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
axis([0,ntimes,-500,Storages{nst,2}*Storages{nst,8}/100*1.1])
yyaxis left
legnames=[{'Storage Charge Level'};{'Storage self-discharge'};{'Storage charge'};{'Storage charge loss'};{'Storage discharge'};{'Storage discharge loss'}];
logictag=logical([1 L{i}(7) 1 L{i}(8) 1 L{i}(9)]);
legend(gca,legnames(logictag))

% ax(1)=gca                   %da sistemare poichè non funzica
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


Objective=sum(sum([Fuels{:,2}]'.*(fuelusage.*timestep)))+sum(sum([Networks{:,4}]'.*(NETWORKbought.*timestep)))-sum(sum([Networks{:,5}]'.*(NETWORKsold.*timestep)));

cont=0;
can=0;
tags={};
costvals={};
gainvals={};

for i=1:Nmachines
    if ismember(Machines{i,2},[Fuels{:,1}])&&sum(value(INPUT{i}))>1e-6
        cont=cont+1;
        tags{cont}=strcat(Machines{i,1},' fuel');
        costvals{cont}=value(INPUT{i})'.*Fuels{ismember([Fuels{:,1}], Machines{i,2}),2};
    end
end
for i=1:Nnetworks
    if sum(value(NETWORKbought(i,:)))>0
        cont=cont+1;
        tags{cont}=strcat(Networks{costorder(i),1}, 'bought from network');
        costvals{cont}=value(NETWORKbought(i,:)').*Networks{costorder(i),4};
    end
end
for i=1:Nnetworks
    if sum(value(NETWORKsold(i,:)))>0
        can=can+1;
        cont=cont+1;
        tags{cont}=strcat(Networks{costorder(i),1}, 'sold on network');
        gainvals{can}=value(NETWORKsold(i,:)').*Networks{costorder(i),5};
    end
end
for i=1:Noutputs
    if sum(value(slacks(i,:))) > 0
        cont=cont+1;
        tags{cont}=strcat(Outputs{i}, ' unmet demand');
        costvals{cont}=value(slacks(i,:))'.*slackcost(i);
    end
end

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
Obj=value((sum([Fuels{:,2}]'.*(fuelusage.*timestep),1))+(sum([Networks{costorder,4}]'.*(NETWORKbought.*timestep),1))-(sum([Networks{costorder,5}]'.*(NETWORKsold.*timestep),1)));
plot(Obj,'k','LineWidth',2)
tags{end+1}='Overall cost function';
legend(gca,[tags{:}])
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
hold off
