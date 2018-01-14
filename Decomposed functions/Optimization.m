
%% CREATION OF OPTIMIZATION MODEL

if suppresswarning==1
    warning('off','all')
end

%Variable time resolution
%timestep=sdpvar(ntimes-1,1,'full');

%Binary variables for machines
Z=binvar(Nmachines,ntimes,'full');
indicator=sdpvar(Nmachines,ntimes+histdepth-1);
%On/Off hystory
OnOffHist=binvar(Nmachines,histdepth,'full');

%for continuous ramp constraints considering only 1d.o.f. machines (1 value for each machine)
LastProd=sdpvar(1,Nmachines,'full');

D{2}=sdpvar(size(Dall{1},1),ntimes,'full');
for h=1:size(Fuels,1)
    Fuels{h,2}=sdpvar(ntimes,1,'full');
end
if Nnetworks~=0
    for h=1:Nnetworks
        Networks{h,4}=sdpvar(ntimes,1,'full');
        Networks{h,5}=sdpvar(ntimes,1,'full');
    end
else
    Networks{1,4}=sdpvar(1,1,'full');
    Networks{1,5}=sdpvar(1,1,'full');
end
if Nundisp~=0
    for h=1:size(UndProd,1)
        UndProd{h,3}=sdpvar(size(UndProd{h,3},1),ntimes,'full');    
    end
else
    UndProd=cell(1,3);
    UndProd{1,3}=sdpvar(1,1,'full');
end
    


%Storages variables

if Nstorages~=0
    STORstart=sdpvar(Nstorages,1,'full');
    STORAGEpower=sdpvar(Nstorages,ntimes,'full'); %create two variables for each storage and each timestep: storage content and charge/discharge for that timestep
    STORAGEpin=sdpvar(Nstorages,ntimes,'full');
    STORAGEpout=sdpvar(Nstorages,ntimes,'full');
    zSTOR=binvar(Nstorages,ntimes,'full');
    STORAGEcharge=sdpvar(Nstorages,ntimes,'full');
else
    STORstart=sdpvar(1,1);
    STORAGEpower=sdpvar(1,1);
    STORAGEcharge=sdpvar(1,1);
    STORAGEpin=sdpvar(1,1);
    STORAGEpout=sdpvar(1,1);
end

%Networks variables
if Nnetworks~=0
    NETWORKbought=sdpvar(Nnetworks,ntimes,'full');
    NETWORKsold=sdpvar(Nnetworks,ntimes,'full');
else
    NETWORKbought=sdpvar(1,1,'full');
    NETWORKsold=sdpvar(1,1,'full');
end
zNET=binvar(Nnetworks,ntimes,'full');


%Slack variables
slacks=sdpvar(size(D{2},1),ntimes,'full');                       
slackcost=ones(1,size(D{2},1))*1e2;
Constr = slacks >= 0;


%%%%%%%%%%%%%%%%%%%%%%%%
% MACHINES CONSTRAINTS %
%%%%%%%%%%%%%%%%%%%%%%%%

%Machines characteristic curves
for i=1:Nmachines
    INPUT{i}=sdpvar(numel(Machines{i,2}),ntimes,'full');           %Input variables for each machine
    OUTPUT{i}=sdpvar(numel(Machines{i,3}),ntimes,'full');           %Output variables for each machine
    coeffs{i}=cell(ntimes,1);
    slopev{i}=sdpvar(J(i)-1,numel(Machines{i,3}),ntimes,'full');
    interceptv{i}=sdpvar(J(i)-1,numel(Machines{i,3}),ntimes,'full'); 
    for j=1:ntimes
        coeffs{i}{j} = sdpvar(J(i),numel(Machines{i,2})+numel(Machines{i,3}),'full'); %coefficients of operating points, variables since they will depend on T
    end
    

    
    if convexflag(i)&&convcheck==true     %characteristic function IS CONVEX in all outputs
        %convexflag{i}=true;
        %only if the convexcheck has been passed (I need to make references
        %to variable coefficients, so I can't make the check here)
        for k=1:numel(Machines{i,3})
            for f=1:(J(i)-1)
                Constr=[Constr INPUT{i}(:)'>=OUTPUT{i}(k,:).*permute(slopev{i}(f,k,:),[1,3,2])+permute(interceptv{i}(f,k,:),[1,3,2]).*Z(i,:)];
                Constr=[Constr (OUTPUT{i}>=0):'Output Positivity'];
            end
        end
        
    else   %characteristic function IS NOT CONVEX in all outputs 
            
        alfas{i}=sdpvar(J(i),ntimes,'full');           %alfa variables for each machine and each timestep
        Constr=[Constr (0<=alfas{i}<=1):'Alfas limitation'];                             %alfa values always in between 0 and 1
        Constr=[Constr (sum(alfas{i},1)==Z(i,:)):'Alfas only if machine on'];          %sum of alfas in one timestep equal to one if machine on
        betas{i}=binvar(J(i)-1,ntimes);          %beta variables for each machine and each timestep
        Constr=[Constr sum(betas{i},1)== 1];                       %only one segment active in each timestep
        for j=2:(J(i)-1)
        Constr=[Constr alfas{i}(j,:)<=betas{i}(j-1,:)+betas{i}(j,:)];        %only alfas which are edge of selected beta segment are greater than zero
        end
        Constr=[Constr alfas{i}(1,:)<=betas{i}(1,:)];                       %first alfa constraint
        Constr=[Constr alfas{i}(J(i),:)<=betas{i}((J(i)-1),:)];             %last alfa constraint

        %Relation between alfas, betas and input/output(s)
            colP=1;     %starting column for reading sampled points matrix

        for j=1:numel(Machines{i,2})
            temp=cellfun(@(x) x(:,colP),coeffs{i},'uni',false);
            temp=[temp{:}];
            Constr=[Constr INPUT{i}(j,:) == sum(temp.*alfas{i},1)];
            %Constr=[Constr (INPUT{i}(j,h)==coeffs{i}{h}(:,colP)'*alfas{i}(:,h))];
            colP=colP+1;
        end
        for j=1:numel(Machines{i,3})
            temp=cellfun(@(x) x(:,colP),coeffs{i},'uni',false);
            temp=[temp{:}];                
            Constr=[Constr OUTPUT{i}(j,:) == sum(temp.*alfas{i},1)];
            %Constr=[Constr OUTPUT{i}(j,h)==coeffs{i}{h}(:,colP)'*alfas{i}(:,h)];
            colP=colP+1;
        end
 
    end

end



%Operation limits, ramp limits(1) and on/off variables
%(1) ramp limits are assumed to be in x/h --> to obtain constraint coherent
%with current timestep we just multiply that limit for timestep duration in
%hours
for i=1:Nmachines
    
    if onoffflag(i)
    Constr=[Constr Z(i)==1];
    end
    
    Constr=[Constr Z(i,:)*Machines{i,5}(1)<= INPUT{i}(:)' <= Z(i,:)*Machines{i,5}(2)];    %operation limit always referred to INPUT
    %Big M setted coherently with primary output of machine i (it is used
    %in ramp limit, which is referred to that output only
    if all(~(strcmp(D{1},Machines{i,3}(1))))~=1
        Mbig=M(strcmp(D{1},Machines{i,3}(1)));
    else
        Mbig=Mbigspare;
    end
    %Vincolo in salita con accensione macchina a carico arbitrario
    Constr=[Constr OUTPUT{i}(1,2:end) - OUTPUT{i}(1,1:(end-1)) <= Machines{i,6}(2)*timestep(1:end-1)' + (1-Z(i,1:(end-1)))*Mbig]; %again, ramp limit on main machine output
    Constr=[Constr OUTPUT{i}(1,1) - LastProd(i) <= Machines{i,6}(2)*basetimestep + (1-OnOffHist(i,end))*Mbig]; %Initial ramping constraint 
    %Vincolo in salita con accensione macchina a minimo carico --> il
    %minimo carico può essere sostituito con un nuovo parametro ad hoc che
    %consenta alla macchina in fase di avviamento di lavorare anche sotto al minimo carico normale 
%     Constr=[Constr OUTPUT{i}(1,2:end) - OUTPUT{i}(1,1:(end-1)) <= Machines{i,6}(2)*timestep.*Z(i,1:(end-1)) + (Z(i,2:end)-Z(i,1:(end-1)))*Machines{i,5}(1)];
    
    %Vincolo in discesa 
    Constr=[Constr OUTPUT{i}(1,1:(end-1)) - OUTPUT{i}(1,2:end) <= Machines{i,6}(1)*timestep(1:end-1)' + (1-Z(i,2:end))*Mbig]; %again, ramp limit on main machine output
    Constr=[Constr LastProd(i) - OUTPUT{i}(1,1) <= Machines{i,6}(1)*timestep + (1-Z(i,1))*Mbig]; %Initial ramping constraint
    %Vincolo in discesa con spegnimento macchina a minimo carico (da
    %completare)
%     Constr=[Constr OUTPUT{i}(1,1:(end-1)) - OUTPUT{i}(1,2:end) <= Machines{i,6}(1)*timestep.*Z(i,2:end) + (Z(i,1:(end-1))-Z(i,2:end)).*(OUTPUT{i}(1,1:(end-1))]; 

end



%Exclusive operation: among machines tagged with the same flag only one
%must be allowed to work
logicflags=logical(zeros(Nmachines,exclusivegroups));
for i=1:exclusivegroups
    for j=1:Nmachines
        if Machines{j,7}(1)==exclusivetags(i)
            logicflags(j,i)=true;
        end
    end
end
            
for i=1:exclusivegroups
    Constr=[Constr sum(Z(logicflags(:,i),:))<=1];
end



% Min uptime/downtime constraints

%Extended binary matrix
Zext=binvar(Nmachines,ntimes+histdepth);
Constr=[Constr Zext==[OnOffHist Z]];
%On/Off indicator
Constr=[Constr indicator(:,:) == Zext(:,2:end)-Zext(:,1:(end-1))];
extendedtimestep=[ones(histdepth,1)'*basetimestep timestep'];

%Uptime
for unit = 1:Nmachines
    baseminupsteps=ceil(Machines{unit,7}(2)/basetimestep);
    for k = (histdepth-baseminupsteps+2):(histdepth+ntimes)
        % indicator will be 1 only when switched on
        l=0;
        while (k+l+1)<(histdepth+ntimes)&&(sum(extendedtimestep(k:(k+l)))<=Machines{unit,7}(2))
            l=l+1;
        end
        range = k:(k+l);
        % range = k:min(histdepth+ntimes,k+minupsteps-1);
        % Constraints will be redundant unless indicator = 1
        Constr = [Constr, Zext(unit,range) >= indicator(unit,k-1)];
    end
end

%Downtime
for unit = 1:Nmachines
    basemindownsteps=ceil(Machines{unit,7}(3)/basetimestep);
    for k = (histdepth-basemindownsteps+2):(histdepth+ntimes)
        % indicator will be 1 only when switched on
        l=0;
        while (k+l+1)<(histdepth+ntimes)&&(sum(extendedtimestep(k:(k+l)))<=Machines{unit,7}(3))
            l=l+1;
        end
        range = k:(k+l);
        % range = k:min(histdepth+ntimes,k+mindownsteps-1);
        % Constraints will be redundant unless indicator = 1
        Constr = [Constr, Zext(unit,range) <= indicator(unit,k-1) + 1];
    end
end

% Startup flag
delta=sdpvar(Nmachines, ntimes, 'full');
Constr=[Constr 0<= delta<= 1];
Constr=[Constr Zext(:,(histdepth+1):end)-Zext(:,histdepth:(end-1))<=delta];

%%%%%%%%%%%%%%%%%%%%%%%
% STORAGE CONSTRAINTS %
%%%%%%%%%%%%%%%%%%%%%%%

Constr=[Constr STORAGEpin>=0];
Constr=[Constr STORAGEpout>=0];
  
FinStorCharge = sdpvar(1,Nstorages);

for i=1:Nstorages
    Constr=[Constr STORAGEpower(i,:)==STORAGEpout(i,:)*Storages{i,6}-STORAGEpin(i,:)/Storages{i,4}];
    Constr=[Constr STORAGEpout(i,:)<=Storages{i,5}*zSTOR(i,:)];    %limit in charge/discharge
    Constr=[Constr STORAGEpin(i,:)<=Storages{i,3}*(1-zSTOR(i,:))];    %limit in charge/discharge
    Constr=[Constr STORAGEcharge(i,1)==STORstart(i)];    %energy content initial condition
    Constr=[Constr STORAGEcharge(i,2:end)==STORAGEcharge(i,1:(end-1)).*(1-Storages{i,7}.*timestep(1:end-1)')-(STORAGEpout(i,1:(end-1))-STORAGEpin(i,1:(end-1))).*timestep(1:end-1)']; %energy content evolution in time
    % Energy(k+1)=Energy(k)[kWh]-Power(k)[kW]*?t[h] <--- assicurati che
    % unità di misura siano coerenti!!!
%     if symtype~=3
%         Constr=[Constr STORAGEcharge(i,end).*(1-Storages{i,7}.*timestep(end))-(STORAGEpout(i,end)-STORAGEpin(i,end)).*timestep(end)==STORstart(i)]; %cyclic storage charge condition  
        Constr=[Constr STORAGEcharge(i,end).*(1-Storages{i,7}.*timestep(end))-(STORAGEpout(i,end)-STORAGEpin(i,end)).*timestep(end) == FinStorCharge(i)]; %final storage valorization
%     end
    Constr=[Constr Storages{i,9}./100<=STORAGEcharge(i,:)./Storages{i,2}<=Storages{i,8}./100];    %SOC constraints
    Constr=[Constr Storages{i,9}./100<=FinStorCharge(i)./Storages{i,2}<=Storages{i,8}./100];    %SOC constraints    
end



%%%%%%%%%%%%%%%%%%%%
% GRID CONSTRAINTS %
%%%%%%%%%%%%%%%%%%%%

if Nnetworks~=0
    Constr=[Constr NETWORKbought >= 0];
    Constr=[Constr NETWORKsold >= 0];
    for i=1:Nnetworks
        %Setting big M in accordance with good of network (if possible)
        if all(~(strcmp(D{1},Networks{i,1})))==0
            maxsold=min(M(strcmp(D{1},Networks{i,1})),Networks{i,3}); %check if there are more restrictive limits on exchange with network
            maxbought=min(M(strcmp(D{1},Networks{i,1})),Networks{i,2});
        else
            maxsold=min(Mbigspare,Networks{i,3});         
            maxbought=min(Mbigspare,Networks{i,2});
        end  
        Constr=[Constr NETWORKsold(i,:) <= zNET(i,:) .* maxsold.*timestep'];
        Constr=[Constr NETWORKbought(i,:) <= (1-zNET(i,:)) .* maxbought.*timestep'];
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMAND / PRODUCTION BALANCES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netprod=sdpvar(Noutputs,ntimes,'full');     %Net units production for each Output (already includes internal consumption of each Output)

Diss=sdpvar(Noutputs,ntimes,'full');                                %NB: CURRENTLY NO UPPER LIMIT SET ON DISSIPATION!!!
Constr=[Constr Diss>=0];
%Constr=[Constr Diss<=repmat([0 inf inf]',1,ntimes)];

for i=1:Noutputs
    nprod=0;
    ncons=0;
    for j=1:Nmachines
        %for each good and each machine, it checks whether that machine
        %produces/consumes that good
        nprod=nprod+ismember(Outputs(i),Machines{j,3});
        ncons=ncons+ismember(Outputs(i),Machines{j,2});
    end
    
    %Here it creates the variables consumption and production of good i for
    %each producer/consumer and each timestep
    machcons{i}=sdpvar(ncons,ntimes,'full');
    machprod{i}=sdpvar(nprod,ntimes,'full');
    if Nundisp~=0; undmachprod{i}=sdpvar(sum(ismember([UndProd{:,2}],Outputs(i))) ,ntimes,'full'); end
end


%THERE MIGHT BE SMARTER WAY OF DOING THIS, WITHOUT AUX VARIABLES --> WOULD
%IT IMPROVE COMPUTATION TIME?
for j=1:Noutputs
    f=0;
    l=0;
    for i=1:Nmachines
        %Production by dispatchable machine i of good j
        if ismember(Outputs(j),Machines{i,3})
            f=f+1;
            Constr=[Constr machprod{j}(f,:)==OUTPUT{i}(ismember(Machines{i,3},Outputs(j)),:)];
        end
        %Consumption by dispatchable machine i of good j
        if ismember(Outputs(j),Machines{i,2})
            l=l+1;
            Constr=[Constr machcons{j}(l,:)==INPUT{i}(ismember(Machines{i,2},Outputs(j)),:)+delta(i,:)*SUcosts(i)];
        end

    end
    f=0;
    for i=1:Nundisp
        %Production by undispatchable machine i of good j
        if ismember(Outputs(j),UndProd{i,2})
            f=f+1;
            Constr=[Constr undmachprod{j}(f,:)==UndProd{i,3}(ismember(UndProd{i,2},Outputs(j)),:)];
        end
    end
    
    if Nundisp~=0
        Constr=[Constr netprod(j,:)==sum(machprod{j},1)+sum(undmachprod{j},1)-sum(machcons{j},1)];
    else
        Constr=[Constr netprod(j,:)==sum(machprod{j},1)-sum(machcons{j},1)];
    end
    
end


%Consider or not storage and network in balance for each good
store=zeros(1,Noutputs);
net=zeros(1,Noutputs);
if Nstorages~=0
    store=ismember(Outputs(:),Storages(:,1));
end
if Nnetworks~=0
    net=ismember(Outputs(:),Networks(:,1));
end

for j=1:Noutputs
    %balances in each scenario
    if store(j)==1&&net(j)==1
        [~,~,netind]=intersect(Outputs(j),Networks(:,1));
        [~,~,storeind]=intersect(Outputs(j),Storages(:,1));
        Constr=[Constr [netprod(j,:) + STORAGEpower(storeind,:) + NETWORKbought(netind,:) - NETWORKsold(netind,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
    elseif store(j)==1&&net(j)==0
        [~,~,storeind]=intersect(Outputs(j),Storages(:,1));
        Constr=[Constr [netprod(j,:) + + STORAGEpower(storeind,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
    elseif store(j)==0&&net(j)==1
        [~,~,netind]=intersect(Outputs(j),Networks(:,1));
        Constr=[Constr [netprod(j,:) + NETWORKbought(netind,:) - NETWORKsold(netind,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
    else
        Constr=[Constr [netprod(j,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
    end
end

for i=1:Nfuels
    ncons=0;
    for j=1:Nmachines
        %for each fuel and each machine, it checks whether that machine
        %consumes that good
        cont=cell2mat(strfind(Machines{j,2}(:),Fuels{i}));
        cont(isempty(cont))=0;
        ncons=ncons+cont;
    end
    %Here it creates the variables consumption and production of good i for
    %each producer/consumer and each timestep
    fuelcons{i}=sdpvar(ncons,ntimes,'full');
end

fuelusage=sdpvar(Nfuels,ntimes,'full');

for j=1:Nfuels
    l=0;
    for i=1:Nmachines
        %usage of fuel j by machine i
        for h=1:numel(Machines{i,2})
            if isequal(Fuels{j,1},Machines{i,2}(h))
                l=l+1;
                Constr=[Constr fuelcons{j}(l,:)==INPUT{i}(h,:)+delta(i,:)*SUcosts(i)];
            end
        end
    end
    
    Constr=[Constr fuelusage(j,:)==sum(fuelcons{j},1)];

end


if suppresswarning==1
    warning('on','all')
end

%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE FUNCTION %
%%%%%%%%%%%%%%%%%%%%%%
%[~,~,costorder]=intersect(Outputs,{Networks{:,1}},'stable'); %cost and networks might be listed differently
if Nnetworks~=0
    Objective=sum(sum([Fuels{:,2}]'.*(fuelusage.*repmat(timestep',[Nfuels 1]))))+sum(sum([Networks{:,4}]'.*(NETWORKbought.*repmat(timestep',[Nnetworks 1]))))-sum(sum([Networks{:,5}]'.*(NETWORKsold.*repmat(timestep',[Nnetworks 1]))))+sum(slackcost*(slacks.*repmat(timestep',[size(D{2},1) 1])));
else
    Objective=sum(sum([Fuels{:,2}]'.*(fuelusage.*repmat(timestep',[Nfuels 1]))))+sum(sum(repmat(slackcost',1,ntimes).*slacks.*repmat(timestep',[size(D{2},1) 1])))-sum(FinStorCharge)*0.0001;
end

coefcontainer=[coeffs{:}];
Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart coefcontainer{:} slopev{:} interceptv{:}};
Param_opt=Param;
% Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart};
if symtype==2||symtype==1
    advance=ntimes+1;
elseif symtype==3
    advance=roladvance+1;
end


Outs={INPUT{:} OUTPUT{:} STORAGEcharge STORAGEpower NETWORKbought NETWORKsold Diss slacks Zext(:,advance:(advance+histdepth-1)) fuelusage delta netprod FinStorCharge};

ops = sdpsettings('solver','gurobi','gurobi.MIPGap',0.005,'gurobi.MIPGapAbs',5e-2,'verbose',3);

Model=optimizer(Constr,Objective,ops,Param,Outs);

%sol=optimize(Constr,Objective,ops)

