%% CREATION OF OPTIMIZATION PROBLEM

if suppresswarning==1
    warning('off','all')
end

%%ATTIVA CONTROLLO CONVESSITà 
convcheck=true;


%Binary variables for machines
Z=binvar(Nmachines,ntimes); 


%Storages variables
STORAGEpower=sdpvar(Nstorages,ntimes); %create two variables for each storage and each timestep: storage content and charge/discharge for that timestep
STORAGEcharge=sdpvar(Nstorages,ntimes);
finstore=sdpvar(Nstorages);         %final storage state
%Networks variables
NETWORKbought=sdpvar(Nnetworks,ntimes);
NETWORKsold=sdpvar(Nnetworks,ntimes);
zNET=binvar(Nnetworks,ntimes);


%Slack variables
slacks=sdpvar(size(D{2},1),ntimes);                       
slackcost=ones(1,size(D{2},1))*1e7;
Constr = slacks >= 0;

%%%%%%%%%%%%%%%%%%%%%%%%
% MACHINES CONSTRAINTS %
%%%%%%%%%%%%%%%%%%%%%%%%

%Machines characteristic curves
for i=1:Nmachines
    INPUT{i}=sdpvar(numel(Machines{i,2}),ntimes);           %Input variables for each machine
    OUTPUT{i}=sdpvar(numel(Machines{i,3}),ntimes);           %Output variables for each machine
    J(i)=size(Machines{i,4},1);                                 %Number of sampled points

    %Convexity check
    slope=0;
    intercept=0;
    for f=1:(J(i)-1)                    %for each segment
        for k=1:numel(Machines{i,3})    %and for each output (NB:we assume we have only one input)
            slope(f,k)=(Machines{i,4}(f+1,1)-Machines{i,4}(f,1))/(Machines{i,4}(f+1,1+k)-Machines{i,4}(f,1+k));
            intercept(f,k)=Machines{i,4}(f+1,1)-slope(f,k)*Machines{i,4}(f+1,1+k);
        end
    end
    

    
    if all(all(slope(2:end,:)>=slope(1:(end-1),:)))&&convcheck==true     %characteristic function IS CONVEX in all outputs
        convexflag{i}=true;
        for k=1:numel(Machines{i,3})
            for f=1:(J(i)-1)
                Constr=[Constr INPUT{i}(:)'>=OUTPUT{i}(k,:).*slope(f,k)+intercept(f,k).*Z(i,:)];
                Constr=[Constr OUTPUT{i}>=0];
            end
        end
        
    else   %characteristic function IS NOT CONVEX in all outputs 
            
        alfas{i}=sdpvar(J(i),ntimes);           %alfa variables for each machine and each timestep
        Constr=[Constr 0<=alfas{i}<=1];                             %alfa values always in between 0 and 1
        Constr=[Constr sum(alfas{i},1)==Z(i,:)];          %sum of alfas in one timestep equal to one if machine on
        betas{i}=binvar(J(i)-1,ntimes);          %beta variables for each machine and each timestep
        Constr=[Constr sum(betas{i},1)== 1];                       %only one segment active in each timestep
        for j=2:(J(i)-1)
        Constr=[Constr alfas{i}(j,:)<=betas{i}(j-1,:)+betas{i}(j,:)];        %only alfas which are edge of selected beta segment are greater than zero
        end
        Constr=[Constr alfas{i}(1,:)<=betas{i}(1,:)];                       %first alfa constraint
        Constr=[Constr alfas{i}(J(i),:)<=betas{i}((J(i)-1),:)];             %last alfa constraint

        %Relation between alfas, betas and input/output(s)
        colP=1;       %starting column for reading sampled points matrix
        for j=1:numel(Machines{i,2})
            Constr=[Constr INPUT{i}(j,:)==Machines{i,4}(:,colP)'*alfas{i}(:,:)];
            colP=colP+1;
        end
        for j=1:numel(Machines{i,3})
            Constr=[Constr OUTPUT{i}(j,:)==Machines{i,4}(:,colP)'*alfas{i}(:,:)];
            colP=colP+1;
        end
    end

end



%Operation limits, ramp limits(1) and on/off variables
%(1) ramp limits are assumed to be in x/h --> to obtain constraint coherent
%with current timestep we just multiply that limit for timestep duration in
%hours
for i=1:Nmachines
    Constr=[Constr Z(i,:)*Machines{i,5}(1)<= OUTPUT{i}(1,:) <= Z(i,:)*Machines{i,5}(2)];    %operation limit always referred to main machine output
    %Big M setted coherently with primary output of machine i (it is used
    %in ramp limit, which is referred to that output only
    if all(~(strcmp(D{1},Machines{i,3}(1))))~=1
        Mbig=M(strcmp(D{1},Machines{i,3}(1)));
    else
        Mbig=Mbigspare;
    end
    %Vincolo in salita con accensione macchina a carico arbitrario
    Constr=[Constr OUTPUT{i}(1,2:end) - OUTPUT{i}(1,1:(end-1)) <= Machines{i,6}(2)*timestep + (1-Z(i,1:(end-1)))*Mbig]; %again, ramp limit on main machine output
    %Vincolo in salita con accensione macchina a minimo carico --> il
    %minimo carico può essere sostituito con un nuovo parametro ad hoc che
    %consenta alla macchina in fase di avviamento di lavorare anche sotto al minimo carico normale 
%     Constr=[Constr OUTPUT{i}(1,2:end) - OUTPUT{i}(1,1:(end-1)) <= Machines{i,6}(2)*timestep.*Z(i,1:(end-1)) + (Z(i,2:end)-Z(i,1:(end-1)))*Machines{i,5}(1)];
    
    %Vincolo in discesa 
    Constr=[Constr OUTPUT{i}(1,1:(end-1)) - OUTPUT{i}(1,2:end) <= Machines{i,6}(1)*timestep + (1-Z(i,2:end))*Mbig]; %again, ramp limit on main machine output
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
for k = 2:ntimes
 for unit = 1:Nmachines
  % indicator will be 1 only when switched on
  indicator = Z(unit,k)-Z(unit,k-1);
  range = k:min(ntimes,k+ceil(Machines{unit,7}(3)/timestep)-1);
  % Constraints will be redundant unless indicator = 1
  Constr = [Constr, Z(unit,range) >= indicator];
 end
end

for k = 2:ntimes
 for unit = 1:Nmachines
  % indicator will be 1 only when switched on
  indicator = Z(unit,k-1)-Z(unit,k);
  range = k:min(ntimes,k+ceil(Machines{unit,7}(4)/timestep)-1);
  % Constraints will be redundant unless indicator = 1
  Constr = [Constr, Z(unit,range) <= 1-indicator];
 end
end


%%%%%%%%%%%%%%%%%%%%%%%
% STORAGE CONSTRAINTS %
%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nstorages
    Constr=[Constr -Storages{i,3}<=STORAGEpower(i,:)<=Storages{i,5}];    %limit in charge/discharge
    Constr=[Constr STORAGEcharge(i,1)==Storages{i,10}];    %energy content initial condition
    Constr=[Constr STORAGEcharge(i,2:end)==STORAGEcharge(i,1:(end-1)).*(1-Storages{i,7}.*timestep)-STORAGEpower(i,1:(end-1)).*timestep]; %energy content evolution in time
    % Energy(k+1)=Energy(k)[kWh]-Power(k)[kW]*?t[h] <--- assicurati che
    % unità di misura siano coerenti!!!
    Constr=[Constr finstore(i)==STORAGEcharge(i,end)-STORAGEpower(i,end)*timestep];    %final charge condition (initial condition of next segment
    Constr=[Constr finstore(i)==STORAGEcharge(i,1)];    %cyclic storage charge condition
    Constr=[Constr Storages{i,9}./100<=STORAGEcharge(i,:)./Storages{i,2}<=Storages{i,8}./100];    %SOC constraints
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
        Constr=[Constr NETWORKsold(i,:) <= zNET(i,:) .* maxsold];
        Constr=[Constr NETWORKbought(i,:) <= (1-zNET(i,:)) .* maxbought];
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMAND / PRODUCTION BALANCES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netprod=sdpvar(Noutputs,ntimes);     %Net units production for each Output (already includes internal consumption of each Output)

Diss=sdpvar(Noutputs,ntimes);                                %NB: CURRENTLY NO UPPER LIMIT SET ON DISSIPATION!!!
Constr=[Constr Diss>=0];

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
    machcons{i}=sdpvar(ncons,ntimes);
    machprod{i}=sdpvar(nprod,ntimes);
    undmachprod{i}=zeros(1,ntimes);
end

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
            Constr=[Constr machcons{j}(l,:)==INPUT{i}(ismember(Machines{i,2},Outputs(j)),:)];
        end

    end
    f=0;
    for i=1:Nundisp
        %Production by undispatchable machine i of good j
        if ismember(Outputs(j),UndProd{i,2})
            f=f+1;
            undmachprod{j}(f,:)=UndProd{i,3}(ismember(UndProd{i,2},Outputs(j)),:);
        end
    end
    
    Constr=[Constr netprod(j,:)==sum(machprod{j},1)+sum(undmachprod{j},1)-sum(machcons{j},1)];

end


%Consider or not storage and network in balance for each good
store=ismember(Outputs(:),Storages(:,1));
net=ismember(Outputs(:),Networks(:,1));

for j=1:Noutputs
    %balances in each scenario
    if store(j)==1&&net(j)==1
        [~,~,netind]=intersect(Outputs(j),Networks(:,1));
        [~,~,storeind]=intersect(Outputs(j),Storages(:,1));
        Constr=[Constr [netprod(j,:) + STORAGEpower(storeind,:) + NETWORKbought(netind,:) - NETWORKsold(netind,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
    elseif store(j)==1&&net(j)==0
        [~,~,storeind]=intersect(Outputs(j),Storages(:,1));
        Constr=[Constr [netprod(j,:) + STORAGEpower(storeind,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
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
    fuelcons{i}=sdpvar(ncons,ntimes);
end

fuelusage=sdpvar(Nfuels,ntimes);

for j=1:Nfuels
    l=0;
    for i=1:Nmachines
        %usage of fuel j by machine i
        for h=1:numel(Machines{i,2})
            if isequal(Fuels{j,1},Machines{i,2}(h))
                l=l+1;
                Constr=[Constr fuelcons{j}(l,:)==INPUT{i}(h,:)];
            end
        end
    end
    
    Constr=[Constr fuelusage(j,:)==sum(fuelcons{j},1)];

end


if suppresswarning==1
    warning('on','all')
end

waitbar(0.6,han,'Problem Solution')

%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE FUNCTION %
%%%%%%%%%%%%%%%%%%%%%%
%[~,~,costorder]=intersect(Outputs,{Networks{:,1}},'stable'); %cost and networks might be listed differently
Objective=sum(sum([Fuels{:,2}]'.*(fuelusage.*timestep)))+sum(sum([Networks{:,4}]'.*(NETWORKbought.*timestep)))-sum(sum([Networks{:,5}]'.*(NETWORKsold.*timestep)))+sum(slackcost*slacks);



