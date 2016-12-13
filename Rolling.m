%% DATA READING

clear
clc

addpath(genpath('D:\Dottorato\Ottimizzazione\Filone moretti!'))

Filename='Input.xlsm';

[~,range]=xlsread(Filename,'Demand','datarange');
[Dtot,Dnames]=xlsread(Filename,'Demand',range{1});
[~,range]=xlsread(Filename,'Undispatch','ddatarange');
[Undisp,Undispnames]=xlsread(Filename,'Undispatch',range{1});
[~,range]=xlsread(Filename,'Prices','pricerange');
[Prices,Pricetags]=xlsread(Filename,'Prices',range{1});

x=find(~cellfun(@isempty,Undispnames(1,:)));
Nundisp=length(x);
x(end+1)=length(Undispnames(2,:))+1;

for i=1:Nundisp
    UndProd{i,1}=Undispnames(1,x(i));
    UndProd{i,2}=Undispnames(2,x(i):(x(i+1)-1));
    UndProd{i,3}=Undisp(:,x(i):(x(i+1)-1))';
end
    
Objective=0;
Constr=0;

    
%Simulation horizon and timestep settings
timestep = xlsread(Filename,'Demand','tdur');   % simulation timestep [h]
ntimes=10;                                      % rolling horizon 
ntimestot=size(Dtot,1);                         % total number of timesteps
days=ntimes*timestep/24;                        % simulation days

updateT=1;                                      % period for calling again solution of problem

D=Dtot(ntimes,:);                               % Demand limited to first ntimes timesteps
D={{Dnames{:}}' D'};                            %Demand for each good: row=good ; columns=timestep

rolcont=0;

% FROM HERE WE START TO USE D (problem sized on ntimes)

slacks=sdpvar(size(D,2),ntimes);                       % slack variables
slackcost=ones(1,size(D,2))*1e7;
Constr = slacks >= 0;                                  % slacks are positive

% Create object.
ExcelApp = actxserver('Excel.Application');
ExcelApp.Visible = 1;

% Open file located in the current folder.
ExcelApp.Workbooks.Open(fullfile(pwd,'\',Filename));
% Run Macro1, defined in "ThisWorkBook" with one parameter. A return value cannot be retrieved.
Nmachines=ExcelApp.Run('machinesref');  %Total number of machines
% ExcelApp.workbooks.Saved = 1;
ExcelApp.Workbooks.Item(Filename).Save 
ExcelApp.Quit;
ExcelApp.release;

Nmachines=double(Nmachines);
[~,range]=xlsread(Filename,'DATA','I2');
[~,refranges]=xlsread(Filename,'DATA',range{1});

%INPUT DATA: 
%column1 --> machine name; 
%column2 --> input type; 
%column3 --> output/s type/s; 
%column4 --> sampled operating points; 
%column5 --> operation limits (on first output)
%column6 --> ramp limits (on first output)
%column7 --> SU/SD times, penalty and min up times                          NB: SU time substituted with exclusive on flag (column 1 of group 7)

%Each row corresponds to a different machine
for i=1:Nmachines
    for j=1:3
    [~, Machines{i,j}]=xlsread(Filename,'Machines',refranges{i,j});
    end
    for j=4:7
    Machines{i,j}=xlsread(Filename,'Machines',refranges{i,j});
    end
    flagsvector(i)=Machines{i,7}(1);
end

exclusivetags=unique(flagsvector(flagsvector~=0));
exclusivegroups=length(exclusivetags);


%Inputs lists all possible machine inputs. 
Inputs=[Machines{:,2}];
Inputs=unique(Inputs, 'Stable');    %Gets rid of duplicates


Outputs=([D{:,1}])';                                    % outputs (goods for which a balance is performed) are only taken from goods which have a demand profile; for 
Noutputs=size(Outputs,2);                               % good which only have internal consumption it is necessary to define a fictitious demand equal to zero.
                                                        

%All possible goods (without duplicates)
Goods=unique([Inputs Outputs],'stable');
Ngoods=size(Goods,2);


%STORAGE DATA: 
%column1 --> stored good; 
%column2 --> storage capacity;
%column3 --> max charge rate; 
%column4 --> charge efficiency
%column5 --> max discharge rate;
%column6 --> discharge efficiency;
%column7 --> self discharge rate;
%column8 --> max SOC
%column9 --> min SOC
%column10 --> storage initial condition (== final condition)

%sarebbe bello avere reference variabile su inizio nomi storage ma vabbè
i=1;
[~,storname]=xlsread(Filename,'Stor&Net',strcat('A',num2str(2+i)));
while ~isempty(storname)
    Storages(i,1)=storname;
    storvals=xlsread(Filename,'Stor&Net',strcat('B',num2str(2+i),':J',num2str(2+i)));
    for j=2:10
        Storages{i,j}=storvals(j-1);
    end
    i=i+1;
    [~,storname]=xlsread(Filename,'Stor&Net',strcat('A',num2str(2+i)));
end    

 
Nstorages=size(Storages,1);
STORAGEpower=sdpvar(Nstorages,ntimes); %create two variables for each storage and each timestep: storage content and charge/discharge for that timestep
STORAGEcharge=sdpvar(Nstorages,ntimes);
finstore=sdpvar(Nstorages);         %final storage state

%Networks data: 
%column1 --> network's good
%column2 --> max withdrawal rate
%column3 --> max injection rate
i=1;
[k,netname]=xlsread(Filename,'Stor&Net',strcat('L',num2str(2+i)));
while ~isempty(netname)
    Networks(i,1)=netname;
    netvals{1}=xlsread(Filename,'Stor&Net',strcat('M',num2str(2+i),':M',num2str(2+i)));
    netvals{2}=xlsread(Filename,'Stor&Net',strcat('N',num2str(2+i),':N',num2str(2+i)));
    for j=2:3
        if ~isempty(netvals{j-1})
            Networks{i,j}=netvals{j-1};
        else
            Networks{i,j}=Inf;
        end        
    end
    i=i+1;
    [k,netname]=xlsread(Filename,'Stor&Net',strcat('L',num2str(2+i)));
end    


for i=1:Noutputs
    M(i)=max(D{2}(i,:))*10;  %M1 related to maximum good demand (an order of magnitude higher
end
Mbigspare=max(M);

%controllo qualità su customizzazione big M
% M(:)=Mbigspare;
%Risultato del test: miglioramento dell'1.4% nella velocità di calcolo! Immagino diventi davvero rilevante solo quando gli ordini di grandezza nelle domande dei vari beni sono molto diverse 

% strcmp(D{1},Networks{1,1})    %Comando fondamentale confronto stringa / vettore
% ismember(D{1},Machines{i,3}) nel caso in cui siano due vettori

%Setting exchange limits on various networks: limit might be set by user or
%will be set to big M of network good
for i=1:size(Networks,1)      
    if ~isnumeric(Networks{i,2})||~isfinite(Networks{i,2})  %if limit for exchange is infinite or not a number
        if all(~(strcmp(D{1},Networks{i,1})))~=0
            Networks{i,2}=M(strcmp(D{1},Networks{i,1}));
        else
            Networks{i,2}=Mbigspare;
        end
    end
    if ~isnumeric(Networks{i,3})||~isfinite(Networks{i,3})  %if limit for exchange is infinite or not a number
        if all(~(strcmp(D{1},Networks{i,1})))~=1
            Networks{i,3}=M(strcmp(D{1},Networks{i,1}));
        else
            Networks{i,3}=Mbigspare;
        end
    end
end

Nnetworks=size(Networks,1);
NETWORKbought=sdpvar(Nnetworks,ntimes);
NETWORKsold=sdpvar(Nnetworks,ntimes);
zNET=binvar(Nnetworks,ntimes);

x=find(~cellfun(@isempty,Pricetags(1,:)));  %

if length(x)==2
    for i=1:x(2)-1      %these are the fuels and the prices
        Fuels{i,1}=Pricetags(2,i);
        Fuels{i,2}=Prices(:,i);
    end
    for i=1:Nnetworks
        a=find(strcmp(Networks{i,1},Pricetags(2,x(2):end)))-1; %column for pricing of network i
        Networks{i,4}=Prices(:,x(2)+a);
        Networks{i,5}=Prices(:,x(2)+a+1);
    end
end
       
Nfuels=size(Fuels,2);

system('taskkill /F /IM EXCEL.EXE');


Z=binvar(Nmachines,ntimes); %Binary variables for machines


%% MACHINES CONSTRAINTS
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
    
    %%ATTIVA CONTROLLO CONVESSITà 
    convcheck=false;
    
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
    Constr=[Constr OUTPUT{i}(1,2:end) - OUTPUT{i}(1,1:(end-1)) <= Machines{i,6}(2)*timestep + (1-Z(i,1:(end-1)))*Mbig]; %again, ramp limit on main machine output
    Constr=[Constr OUTPUT{i}(1,1:(end-1)) - OUTPUT{i}(1,2:end) <= Machines{i,6}(1)*timestep + (1-Z(i,2:end))*Mbig]; %again, ramp limit on main machine output
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



%% STORAGE CONSTRAINTS
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

%% GRID CONSTRAINTS


if Nnetworks~=0
    Constr=[Constr NETWORKbought >= 0];
    Constr=[Constr NETWORKsold >= 0];
    for i=1:Nnetworks
        %Setting big M in accordance with good of network (if possible)
        if all(~(strcmp(D{1},Networks{i,1})))~=0
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

netprod=sdpvar(Noutputs,ntimes);     %Net units production for each Output (already includes internal consumption of each Output)

Diss=sdpvar(Noutputs,ntimes);                                %NB: CURRENTLY NO UPPER LIMIT SET ON DISSIPATION!!!
Constr=[Constr Diss>=0];

for i=1:Noutputs
    nprod=0;
    ncons=0;
    for j=1:Nmachines
        %for each good and each machine, it checks whether that machine
        %produces/consumes that good
        cont=cell2mat(strfind(Machines{j,3}(:),Outputs{i}));
        cont(isempty(cont))=0;
        nprod=nprod+cont;
        cont=cell2mat(strfind(Machines{j,2}(:),Outputs{i}));
        cont(isempty(cont))=0;
        ncons=ncons+cont;
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
        for h=1:numel(Machines{i,3})
            if isequal(Outputs(j),Machines{i,3}(h))
                f=f+1;
                Constr=[Constr machprod{j}(f,:)==OUTPUT{i}(h,:)];
            end
        end
        %Consumption by dispatchable machine i of good j
        for h=1:numel(Machines{i,2})
            if isequal(Outputs(j),Machines{i,2}(h))
                l=l+1;
                Constr=[Constr machcons{j}(l,:)==INPUT{i}(h,:)];
            end
        end
    end
    f=0;
    for i=1:Nundisp
        %Production by undispatchable machine i of good j
        for h=1:numel(UndProd{i,2})
            if isequal(Outputs(j),UndProd{i,2}(h))
                f=f+1;
                undmachprod{j}(f,:)=UndProd{i,3}(h,:);
            end
        end
    end
    
    Constr=[Constr netprod(j,:)==sum(machprod{j},1)+sum(undmachprod{j},1)-sum(machcons{j},1)];

end


%Consider or not storage and network in balance for each good
store=zeros(Noutputs,1);
net=zeros(Noutputs,1);
for j=1:Noutputs
    for h=1:Nstorages
        if isequal(Storages{h,1},Outputs{j})
            store(j)=1;
            break           %this way h keeps track of which row in Storages we should consider
        end
    end
    for k=1:Nnetworks
        if isequal(Networks{k,1},Outputs{j})
            net(j)=1;
            break           %this way k keeps track of which row in Networks we should consider
        end
    end
    %balances in each scenario
    if store(j)==1&&net(j)==1
        Constr=[Constr [netprod(j,:) + STORAGEpower(h,:) + NETWORKbought(k,:) - NETWORKsold(k,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
    elseif store(j)==1&&net(j)==0
        Constr=[Constr [netprod(j,:) + STORAGEpower(h,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
    elseif store(j)==0&&net(j)==1
        Constr=[Constr [netprod(j,:) + NETWORKbought(k,:) - NETWORKsold(k,:) + slacks(j,:) == D{2}(j,:)+Diss(j,:)]];
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

Objective=sum(sum([Fuels{:,2}]'.*(fuelusage.*timestep)))+sum(sum([Networks{:,4}]'.*(NETWORKbought.*timestep)))-sum(sum([Networks{:,5}]'.*(NETWORKsold.*timestep)))+sum(slackcost*slacks);

ops = sdpsettings('solver','gurobi','gurobi.MIPGap',0.005);


%% Rolling

sol=optimize(Constr,Objective,ops)

for i=updateT:updateT:ntimestot
    rolcont=rolcont+i;
    D=Dtot(rolcont:(rolcont+ntimes),:);                               % Demand limited to first ntimes timesteps
    D={{Dnames{:}}' D'};                            %Demand for each good: row=good ; columns=timestep
    sol=optimize(Constr,Objective,ops)


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
            if isequal(Outputs(i),Machines{j,3}(h)) && sum(value(OUTPUT{j}(h,:)))~=0
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
    if sum(ismember((Outputs{i}),{Storages{:,1}}))>0 
        nstor=nstor+1;
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
    if sum(ismember((Outputs{i}),[Networks{:,1}]))==length(Outputs{i})
        nnet=nnet+1;
        C(2)=1;
        C(4)=1;
        Pmat{i,6}=value(NETWORKbought(nnet,:));
        Pmat{i,7}=value(NETWORKsold(nnet,:));
    end
    

    Pmat{i,8}=value(Diss(i,:));
    if abs(sum(value(Diss(i,:))))>0.001
        C(5)=1;
    end
    L{i}=C;

end

grey = [0.4,0.4,0.4];

%column1 --> Output name
%column2 --> Tags of dispatchable machines producing that good
%column3 --> Output profile of each dispatchable machine
%column4 --> Storage discharging power profile
%column5 --> Storage charging power profile
%column6 --> Purchased from network
%column7 --> Sold on network
%column8 --> Dissipated
%column9 --> Storage content
%column10 --> Tags of undispatchable machines producing that good
%column11 --> Output profile of each dispatchable machine
%column12 --> Utility consumers tag
%column13 --> Utility consumers consumption profile
%column14 --> Storage self discharge profile
%column15 --> Storage charge losses
%column16 --> Storage discharge losses

%Plot all demand / production profiles
for i=1:Noutputs
figure()
title(strcat(Outputs{i},' Demand / Production profiles'))
xlabel('Timestep')
ylabel('Power [kW]')
hold on
stdnam=[{'Storage Discharge'};{'Purchased'};{'Storage Charge'}; {'Sold'}; {'Dissipated'}];
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
figure()
hold on
yyaxis left
h=bar([Pmat{i,9}' -Pmat{i,14}'],1,'stacked')               %ESS energy level 
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
    if ismember(Machines{i,2},[Fuels{:,1}])&&sum(value(INPUT{i}))>0
        cont=cont+1;
        tags{cont}=strcat(Machines{i,1},' fuel');
        costvals{cont}=value(INPUT{i})'.*Fuels{ismember([Fuels{:,1}], Machines{i,2}),2};
    end
end
for i=1:Nnetworks
    if sum(value(NETWORKbought(i,:)))>0
        cont=cont+1;
        tags{cont}=strcat(Networks{i,1}, 'bought from network');
        costvals{cont}=value(NETWORKbought(i,:)').*Networks{i,4};
    end
    if sum(value(NETWORKsold(i,:)))>0
        can=can+1;
        cont=cont+1;
        tags{cont}=strcat(Networks{i,1}, 'sold on network');
        gainvals{can}=value(NETWORKsold(i,:)').*Networks{i,5};
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
for j=1:size(costvals,2)
    set(h(j),'facecolor',color(j,:));
end
if can~=0
    k=bar(-[gainvals{:}],1,'stacked');
    color=hot(max(size(gainvals,2),8));
    for j=1:size(gainvals,2)
        set(k(j),'facecolor',color(j,:));
    end
end
title('Cost function explosion')
xlabel('Timestep') % x-axis label
ylabel('Timestep cost') % y-axis label
set(gca,'XTick', 0:4:24*days);
Obj=value((sum([Fuels{:,2}]'.*(fuelusage.*timestep),1))+(sum([Networks{:,4}]'.*(NETWORKbought.*timestep),1))-(sum([Networks{:,5}]'.*(NETWORKsold.*timestep),1)));
plot(Obj,'k','LineWidth',2)
tags{end+1}='Overall cost function';
legend(gca,[tags{:}])
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
hold off







