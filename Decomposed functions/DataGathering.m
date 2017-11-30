%% CREATION OF SOLUTION CELL ARRAY (Pmat)

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


if symtype == 1 

Pmat=cell(Noutputs,16);
for i=1:Noutputs
    nprod=0;
    Pmat{i,1}=Outputs{i};
    f=0;
    l=0;
    for j=1:Nmachines
        %Production by machine i of good j
        for h=1:numel(Machines{j,3})
            if isequal(Outputs(i),Machines{j,3}(h)) && sum((OUTPUT{j}(h,:)))>1e-6
                if f==0
                    Pmat{i,2}=Machines{j,1};
                    Pmat{i,3}=[(OUTPUT{j}(h,:))];
                else
                    Pmat{i,2}=[Pmat{i,2};Machines{j,1}];
                    Pmat{i,3}=[Pmat{i,3};(OUTPUT{j}(h,:))];
                end
                f=f+1;
            end
        end
        %Consumption by machine i of good j
        for h=1:numel(Machines{j,2})
            if isequal(Outputs(i),Machines{j,2}(h)) && sum((INPUT{j}(h,:)))~=0
                if l==0
                    Pmat{i,12}=Machines{j,1};
                    Pmat{i,13}=[(INPUT{j}(h,:))];
                else
                    Pmat{i,12}=[Pmat{i,12};Machines{j,1}];
                    Pmat{i,13}=[Pmat{i,13};(INPUT{j}(h,:))];
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
    Pmat{i,4}=zeros(1,ntimestot);
    Pmat{i,5}=zeros(1,ntimestot);
    Pmat{i,6}=zeros(1,ntimestot);
    Pmat{i,7}=zeros(1,ntimestot);
    Pmat{i,8}=zeros(1,ntimestot);
    
    if ismember(Outputs{i},{Storages{:,1}}) %strcmp({Storages{:,1}},Outputs{i}))>0 
        nstor=ismember({Storages{:,1}},Outputs{i});
        C(1)=1;                                 %C: 1=Sdisch 2=NetPurch 3=Stcharge 4=NetSold 5=Dissipation 6=Storagelevel
        C(3)=1;                                 % 7=self disch 8=charge loss 9=discharge loss
        C(6)=1;
        if Storages{nstor,7}~=0
            C(7)=1;
            Pmat{i,14}=(STORAGEcharge(nstor,:))*Storages{nstor,7};
        end
        Storprof=(STORAGEpower(nstor,:));
        Storprof(Storprof>=0)=0;
        Pmat{i,4}(1:ntimes)=(STORAGEpower(nstor,:))-Storprof;
        Pmat{i,5}(1:ntimes)=Storprof;
        Pmat{i,9}=(STORAGEcharge(nstor,:));
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
        [~,~,nnet]=intersect(Outputs(i),Networks(:,1));
        if sum((NETWORKbought(nnet,:)))>1e-6
            C(2)=1;
            Pmat{i,6}(1:ntimes)=(NETWORKbought(nnet,:));
        end
        if sum((NETWORKsold(nnet,:)))>1e-6
            C(4)=1;
            Pmat{i,7}(1:ntimes)=(NETWORKsold(nnet,:));
        end
    end
    

    Pmat{i,8}(1:ntimes)=(Diss(i,:));
    if abs(sum((Diss(i,:))))>0.001
        C(5)=1;
    end
    L{i}=C;

end

cont=0;
can=0;
tags={};
costvals={};
gainvals={};

for i=1:Nmachines
    if ismember(Machines{i,2},[Fuels{:,1}])&&sum((INPUT{i}))>1e-6
        cont=cont+1;
        tags{cont}=strcat(Machines{i,1},' fuel');
        costvals{cont}=(INPUT{i})'.*Fuels{ismember([Fuels{:,1}], Machines{i,2}),2}*timestep(1);
        cont=cont+1;
        tags{cont}=strcat(Machines{i,1},' start-up cost');
        costvals{cont}=startflags(i,:)'.*SUcosts(i).*Fuels{ismember([Fuels{:,1}], Machines{i,2}),2}*timestep(1);
    end
end
for i=1:Nnetworks
    if sum((NETWORKbought(i,:)))>0
        cont=cont+1;
        tags{cont}=strcat(Networks{i,1}, 'bought from network');
        costvals{cont}=(NETWORKbought(i,:)').*Networks{i,4}*timestep(1);
    end
end
for i=1:Noutputs
    if sum((slacks(i,:))) > 0
        cont=cont+1;
        tags{cont}=strcat(Outputs{i}, ' unmet demand');
        costvals{cont}=(slacks(i,:))'.*slackcost(i)*timestep(1);
    end
end
for i=1:Nnetworks
    if sum((NETWORKsold(i,:)))>0
        can=can+1;
        cont=cont+1;
        tags{cont}=strcat(Networks{i,1}, 'sold on network');
        gainvals{can}=(NETWORKsold(i,:)').*Networks{i,5}*timestep(1),;
    end
end


Obj=((sum([Fuels{:,2}]'.*(fuelusage.*repmat(timestep',[Nfuels 1])),1))+(sum([Networks{:,4}]'.*(NETWORKbought.*repmat(timestep',[Nnetworks 1])),1))-(sum([Networks{:,5}]'.*(NETWORKsold.*repmat(timestep',[Nnetworks 1])),1)));

elseif symtype == 2 || symtype == 3
    
    if symtype == 2 
        numsaves=ntimes;
    else
        numsaves=roladvance;
    end
    
if runcount == 1
    
Pmat=cell(Noutputs,16);
for i=1:Noutputs
    nprod=0;
    Pmat{i,1}=Outputs{i};
    f=0;
    l=0;
    for j=1:Nmachines
        %Production by machine i of good j
        for h=1:numel(Machines{j,3})
            if isequal(Outputs(i),Machines{j,3}(h)) %&& sum((OUTPUT{j}(h,:)))>1e-6
                if f==0
                    Pmat{i,2}=Machines{j,1};
                    Pmat{i,3}=[(OUTPUT{j}(h,:))];
                else
                    Pmat{i,2}=[Pmat{i,2};Machines{j,1}];
                    Pmat{i,3}=[Pmat{i,3};(OUTPUT{j}(h,:))];
                end
                f=f+1;
            end
        end
        %Consumption by machine i of good j
        for h=1:numel(Machines{j,2})
            if isequal(Outputs(i),Machines{j,2}(h)) %&& sum((INPUT{j}(h,:)))~=0
                if l==0
                    Pmat{i,12}=Machines{j,1};
                    Pmat{i,13}=[(INPUT{j}(h,:))];
                else
                    Pmat{i,12}=[Pmat{i,12};Machines{j,1}];
                    Pmat{i,13}=[Pmat{i,13};(INPUT{j}(h,:))];
                end
                l=l+1;
            end
        end
    end
    f=0;
    for j=1:Nundisp
        %Production by undispatchable machine i of good j
        for h=1:numel(UndProd{j,2})
            if isequal(Outputs(i),UndProd{j,2}(h)) %&& sum(UndProd{j,3}(h,:))~=0
                if f==0
                    Pmat{i,10}=UndProd{j,1};
                    Pmat{i,11}=UndProd{j,3}(h,1:numsaves);
                else
                    Pmat{i,10}=[Pmat{i,10};UndProd{j,1}];
                    Pmat{i,11}=[Pmat{i,11};UndProd{j,3}(h,1:numsaves)];
                end
                f=f+1;
            end
        end
    end
    C=logical(zeros(9,1));
    Pmat{i,4}=zeros(1,ntimestot);
    Pmat{i,5}=zeros(1,ntimestot);
    Pmat{i,6}=zeros(1,ntimestot);
    Pmat{i,7}=zeros(1,ntimestot);
    Pmat{i,8}=zeros(1,ntimestot);
    
    if Nstorages~=0&&ismember(Outputs{i},{Storages{:,1}}) %strcmp({Storages{:,1}},Outputs{i}))>0 
        nstor=ismember({Storages{:,1}},Outputs{i});
        C(1)=1;                                 %C: 1=Sdisch 2=NetPurch 3=Stcharge 4=NetSold 5=Dissipation 6=Storagelevel
        C(3)=1;                                 % 7=self disch 8=charge loss 9=discharge loss
        C(6)=1;
        if Storages{nstor,7}~=0
            C(7)=1;
            Pmat{i,14}=(STORAGEcharge(nstor,:))*Storages{nstor,7};
        end
        Storprof=(STORAGEpower(nstor,:));
        Storprof(Storprof>=0)=0;
        Pmat{i,4}=(STORAGEpower(nstor,:))-Storprof;
        Pmat{i,5}=Storprof;
        Pmat{i,9}=(STORAGEcharge(nstor,:));
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
        [~,~,nnet]=intersect(Outputs(i),Networks(:,1));
        %if sum((NETWORKbought(nnet,:)))>1e-6
            %C(2)=1;
            Pmat{i,6}=(NETWORKbought(nnet,:));
        %end
        %if sum((NETWORKsold(nnet,:)))>1e-6
            %C(4)=1;
            Pmat{i,7}=(NETWORKsold(nnet,:));
        %end
    end
    

    Pmat{i,8}=(Diss(i,:));
    %if abs(sum((Diss(i,:))))>0.001
        %C(5)=1;
    %end
    L{i}=C;

end

cont=0;
can=0;
tags={};
costvals={};
gainvals={};

for i=1:Nmachines
    if ismember(Machines{i,2},[Fuels{:,1}]) %&&sum((INPUT{i}))>1e-6
        cont=cont+1;
        tags{cont}=strcat(Machines{i,1},' fuel');
        pr=Fuels{ismember([Fuels{:,1}], Machines{i,2}),2};
        costvals{cont}=(INPUT{i})'.*pr(1:numsaves);
        cont=cont+1;
        tags{cont}=strcat(Machines{i,1},' start-up cost');
        costvals{cont}=startflags(i,:)'.*SUcosts(i).*pr(1:numsaves);
    end
end
for i=1:Nnetworks
    %if sum((NETWORKbought(i,:)))>0
        cont=cont+1;
        tags{cont}=strcat(Networks{i,1}, 'bought from network');
        pr=Networks{i,4};
        costvals{cont}=(NETWORKbought(i,:)').*pr(1:numsaves);
    %end
end
for i=1:Noutputs
    %if sum((slacks(i,:))) > 0
        cont=cont+1;
        tags{cont}=strcat(Outputs{i}, ' unmet demand');
        costvals{cont}=(slacks(i,:))'.*slackcost(i);
    %end
end
for i=1:Nnetworks
    %if sum((NETWORKsold(i,:)))>0
        can=can+1;
        cont=cont+1;
        tags{cont}=strcat(Networks{i,1}, 'sold on network');
        pr=Networks{i,5};
        gainvals{can}=(NETWORKsold(i,:)').*pr(1:numsaves);
    %end
end

fuelprice=cell2mat(cellfun(@(x) x(1:numsaves,:)',Fuels(:,2),'UniformOutput',false));
if Nnetworks~=0
    netprice=cell2mat(cellfun(@(x) x(1:numsaves,:)',Networks(:,4),'UniformOutput',false));
    netval=cell2mat(cellfun(@(x) x(1:numsaves,:)',Networks(:,5),'UniformOutput',false));
else
    netprice=0;
    netval=0;
end

%Obj=[Obj ((sum([Fuels{:,2}]'.*(fuelusage.*timestep),1))+(sum([Networks{:,4}]'.*(NETWORKbought.*timestep),1))-(sum([Networks{:,5}]'.*(NETWORKsold.*timestep),1)))];
Obj=(sum(fuelprice.*(fuelusage.*repmat(timestep(1:numsaves)',[Nfuels 1])),1))+(sum(netprice.*(NETWORKbought.*repmat(timestep(1:numsaves)',[Nnetworks 1])),1))-(sum(netval.*(NETWORKsold.*repmat(timestep(1:numsaves)',[Nnetworks 1])),1));

else
   
pos=tstart:(tstart+numsaves-1);
for i=1:Noutputs
    f=0;
    l=0;
    for j=1:Nmachines
        %Production by machine i of good j
        for h=1:numel(Machines{j,3})
            if isequal(Outputs(i),Machines{j,3}(h)) %&& sum((OUTPUT{j}(h,:)))>1e-6  %condizione da rivedere anche in fase creazione matrice
                f=f+1;
                Pmat{i,3}(f,pos)=[(OUTPUT{j}(h,:))];
            end
        end
        %Consumption by machine i of good j
        for h=1:numel(Machines{j,2})
            if isequal(Outputs(i),Machines{j,2}(h)) %&& sum((INPUT{j}(h,:)))~=0     %condizione da rivedere
                l=l+1;
                Pmat{i,13}(l,pos)=[(INPUT{j}(h,:))];
            end
        end
    end
    f=0;
    for j=1:Nundisp
        %Production by undispatchable machine i of good j
        for h=1:numel(UndProd{j,2})
            if isequal(Outputs(i),UndProd{j,2}(h)) %&& sum(UndProd{j,3}(h,:))~=0
                f=f+1;
                Pmat{i,11}(f,pos)=UndProd{j,3}(h,1:numsaves);
            end
        end
    end
    if Nstorages~=0&&ismember(Outputs{i},{Storages{:,1}}) %strcmp({Storages{:,1}},Outputs{i}))>0 
        nstor=ismember({Storages{:,1}},Outputs{i});
        if Storages{nstor,7}~=0
            Pmat{i,14}(pos)=(STORAGEcharge(nstor,:))*Storages{nstor,7};
        end
        Storprof=(STORAGEpower(nstor,:));
        Storprof(Storprof>=0)=0;
        Pmat{i,4}(pos)=(STORAGEpower(nstor,:))-Storprof;
        Pmat{i,5}(pos)=Storprof;
        Pmat{i,9}(pos)=(STORAGEcharge(nstor,:));
        if Storages{nstor,4}~=1
            Pmat{i,15}(pos)=Pmat{i,5}(pos)*(1-Storages{nstor,4});
        end
        if Storages{nstor,6}~=1
            Pmat{i,16}(pos)=Pmat{i,4}(pos)*(1-Storages{nstor,6});
        end
    end
    if sum(strcmp({Networks{:,1}},Outputs{i}))>0
        [~,~,nnet]=intersect(Outputs(i),Networks(:,1));
        %if sum((NETWORKbought(nnet,:)))>1e-6       %condizione da rivedere
            Pmat{i,6}(pos)=(NETWORKbought(nnet,:));
        %end
        %if sum((NETWORKsold(nnet,:)))>1e-6         %condizione da rivedere
            Pmat{i,7}(pos)=(NETWORKsold(nnet,:));
        %end
end
    

    Pmat{i,8}(pos)=(Diss(i,:));
    %if abs(sum(Pmat{i,8}))>0.001
    %    C(5)=1;
    %end
    %L{i}(5)=C(5);

end

cont=0;
can=0;

for i=1:Nmachines
    if ismember(Machines{i,2},[Fuels{:,1}])
        cont=cont+1;
        pr=Fuels{ismember([Fuels{:,1}], Machines{i,2}),2};
        costvals{cont}(pos,1)=(INPUT{i})'.*pr(1:numsaves);
        cont=cont+1;
        tags{cont}=strcat(Machines{i,1},' start-up cost');
        costvals{cont}(pos,1)=startflags(i,:)'.*SUcosts(i).*pr(1:numsaves);        
    end
end
for i=1:Nnetworks
    %if sum((NETWORKbought(i,:)))>0
        cont=cont+1;
        pr=Networks{i,4};
        costvals{cont}(pos,1)=(NETWORKbought(i,:)').*pr(1:numsaves);
    %end
end
for i=1:Noutputs
    %if sum((slacks(i,:))) > 0
        cont=cont+1;
        costvals{cont}(pos,1)=((slacks(i,:))'.*slackcost(i))';
    %end
end
for i=1:Nnetworks
    %if sum((NETWORKsold(i,:)))>0
        can=can+1;
        cont=cont+1;
        pr=Networks{i,5};
        gainvals{can}(pos,1)=(NETWORKsold(i,:)').*pr(1:numsaves);
    %end
end

fuelprice=cell2mat(cellfun(@(x) x(1:numsaves,:)',Fuels(:,2),'UniformOutput',false));
if Nnetworks~=0
    netprice=cell2mat(cellfun(@(x) x(1:numsaves,:)',Networks(:,4),'UniformOutput',false));
    netval=cell2mat(cellfun(@(x) x(1:numsaves,:)',Networks(:,5),'UniformOutput',false));
else
    netprice=0;
    netval=0;
end

%Obj=[Obj ((sum([Fuels{:,2}]'.*(fuelusage.*timestep),1))+(sum([Networks{:,4}]'.*(NETWORKbought.*timestep),1))-(sum([Networks{:,5}]'.*(NETWORKsold.*timestep),1)))];
Obj=[Obj ((sum(fuelprice.*(fuelusage.*repmat(timestep(1:numsaves)',[Nfuels 1])),1))+(sum(netprice.*(NETWORKbought.*repmat(timestep(1:numsaves)',[Nnetworks 1])),1))-(sum(netval.*(NETWORKsold.*repmat(timestep(1:numsaves)',[Nnetworks 1])),1)))];

end
end