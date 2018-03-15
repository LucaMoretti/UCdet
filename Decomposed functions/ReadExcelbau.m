%% DATA READING FROM EXCEL FILE
global convcheck
%%IMPORTANTE!! Bisogna fare un unico accesso all'excel, e smistare poi in
%%Matlab i profili alle varie variabili. In questo modo velocizziamo di
%%brutto questa fase della simulazione

%addpath(genpath('D:\Dottorato\Ottimizzazione\Filone moretti!'))

excelpath=fileparts(pwd());

Filename='Input - V_PVT_inter.xlsm';
Filepath=strcat(excelpath,'\',Filename);

% Create object.
ExcelApp = actxserver('Excel.Application');
ExcelApp.Visible = 0;
ExcelApp.Workbooks.Open(fullfile(Filepath));
Workbook = ExcelApp.ActiveWorkbook;
Worksheets = Workbook.sheets;


% Open file located in the current folder.
Nmachines=ExcelApp.Run('machinesref');  %Total number of machines

%READING DEMAND PROFILES
Worksheets.Item('Demand').Activate;
temp=get(ExcelApp.ActiveSheet,'Range',get(ExcelApp.ActiveSheet,'Range','datarange'));
temp = temp.value;
Dall = {temp(1,:)',cell2mat(temp(2:end,:))'};
% [~,range]=xlsread(Filepath,'Demand','datarange');
% [Dall,Dnames]=xlsread(Filepath,'Demand',range{1});

%READING UNDISPATCHABLE PRODUCTION PROFILES
Worksheets.Item('Undispatch').Activate;
range=get(ExcelApp.ActiveSheet,'Range','ddatarange');
% [~,range]=xlsread(Filepath,'Undispatch','ddatarange');
if range.value~=0
    temp=get(ExcelApp.ActiveSheet,'Range',range);
    temp=temp.value;
    Undisp=temp(3:end,:);
    Undispnames=temp(1:2,:);
    
%     [Undisp,Undispnames]=xlsread(Filepath,'Undispatch',range{1});
    x=find(~cellfun(@isempty,Undispnames(1,:)));
    Nundisp=length(x);
    x(end+1)=length(Undispnames(2,:))+1;
else
    Nundisp=0;
    UndProdall=cell(1,1);
end

for i=1:Nundisp
    UndProdall{i,1}=Undispnames(1,x(i));
    UndProdall{i,2}=Undispnames(2,x(i):(x(i+1)-1));
    UndProdall{i,3}=Undisp(:,x(i):(x(i+1)-1))';
end


%READING PRICES PROFILES
Worksheets.Item('Prices').Activate;
temp=get(ExcelApp.ActiveSheet,'Range',get(ExcelApp.ActiveSheet,'Range','pricerange'));
temp=temp.value;
Prices=cell2mat(temp(4:end,:));
Pricetags=temp(1:3,:);
Pricetags(cell2mat(cellfun(@(x) ~ischar(x),Pricetags,'Uniformoutput',0)))={''};

% [~,range]=xlsread(Filepath,''s','pricerange');
% [Prices,Pricetags]=xlsread(Filepath,'Prices',range{1});

    
%Simulation horizon and timestep settings
Worksheets.Item('Demand').Activate;
basetimestep=get(ExcelApp.ActiveSheet,'Range','tdur');
basetimestep=basetimestep.value;
% basetimestep = xlsread(Filepath,'Demand','tdur');   % simulation timestep [h]
ntimestot=size(Dall{2},2);                               % total number of timesteps
days=ntimestot*basetimestep/24;                        % simulation days

%Temperature Profiles %THIS WILL BECOME GENERAL NON-CONTROLLABLE DEGREES OF
%FREEDOM
ndof=get(ExcelApp.ActiveSheet,'Range','ambdof');
ndof=ndof.value;
% ndof=xlsread(Filepath,'Demand','ambdof');


%READING AMBIENT VARIABLES
temp=get(ExcelApp.ActiveSheet,'Range',get(ExcelApp.ActiveSheet,'Range','ambrange'));
temp=temp.value;
ambvar=cell2mat(temp(2:end,:));
ambnames=temp(1,:);
% [~,range]=xlsread(Filepath,'Demand','ambrange');
% [ambvar,ambnames]=xlsread(Filepath,'Demand',range{1});

%killing measure unit (if present)
ambnames=cellfun(@(x) strtok(x),ambnames,'uniformoutput',false);


 

% %Demand for each good: row=good ; columns=timestep
% Dall={{Dnames{:}}' Dall'}; 

Nmachines=double(Nmachines);

Worksheets.Item('DATA').Activate;
refranges=get(ExcelApp.ActiveSheet,'Range',get(ExcelApp.ActiveSheet,'Range','rangeread'));
refranges=refranges.value;
% [~,range]=xlsread(Filepath,'DATA','rangeread');
% [~,refranges]=xlsread(Filepath,'DATA',range{1});

%INPUT DATA: 
%column1 --> machine name;
%column2 --> input type;
%column3 --> output/s type/s;
%column4 --> sampled operating points;
%column5 --> operation limits (on first output)
%column6 --> ramp limits (on first output)
%column7 --> SU/SD times, penalty and min up times                          NB: SU time substituted with exclusive on flag (column 1 of group 7)
%column9 --> operating coefficients
%column9 --> Ambient parameters

%Each row corresponds to a different machine
Worksheets.Item('Machines').Activate;
for i=1:Nmachines
    for j=1:3
        temp=get(ExcelApp.ActiveSheet,'Range',refranges{i,j});
        if ischar(temp.value)
            Machines{i,j}={temp.value};
        else
            Machines{i,j}=temp.value;
        end
%     [~, Machines{i,j}]=xlsread(Filepath,'Machines',refranges{i,j});
    end
    %ambient tags 
    if refranges{i,8} ~= "N/A"
        temp=get(ExcelApp.ActiveSheet,'Range',refranges{i,8});
        Machines{i,9}=(temp.value);
%         [~, Machines{i,9}]=xlsread(Filepath,'Machines',refranges{i,8});
    end
    for j=4:7
        temp=get(ExcelApp.ActiveSheet,'Range',refranges{i,j});
        Machines{i,j}=cell2mat(temp.value);
%     Machines{i,j}=xlsread(Filepath,'Machines',refranges{i,j});
    end
    flagsvector(i)=Machines{i,7}(1);
end


exclusivetags=unique(flagsvector(flagsvector~=0));
exclusivegroups=length(exclusivetags);


a=cellfun(@(x) x(:,2:3),Machines(:,7),'UniformOutput',false);
histdepth=max([ceil(max([a{:}])/basetimestep) 1]);                                          %already in number of timesteps



%Rearranging operating points
for i=1:Nmachines
    nambs=length(Machines{i,9});
    matrix=Machines{i,4};
    if nambs ~= 0
        ntemps=sum(~isnan(matrix(:,1:nambs)),1); %temperatures at which characteristic curves is available
        J(i)=size(matrix,1)/ntemps(end);      %number of sampled points for each machine
    else
        J(i)=size(matrix,1);
    end
    Machines{i,4}=cell(1,2);
    t=0;
    for h=1:nambs
        for l=2:size(matrix,1)
            if isnan(matrix(l,h))
                matrix(l,h)=matrix(l-1,h);
            end
        end
    end
    if nambs~=0
        for h=1:ntemps(end)
            Machines{i,4}{1}(h,:)=matrix(1+(h-1)*J(i),1:nambs);
            Machines{i,4}{2}{h,:}=matrix(1+(h-1)*J(i):h*J(i),nambs+1:end);
        end    
    else
        Machines{i,4}{2}=matrix;
    end
end

 C = cell(1,1,ntimestot);

%Operating Coefficients for each machine            N.B. It could be possible to extrapolate as well, and to change the interpolation method
                                                    %from linear
for i=1:Nmachines
    t=Machines{i,4}{1};
    if size(t,1)>1
        influencer=ismember(ambnames,Machines{i,9});
        ambprof=ambvar(:,influencer);
        x=cat(3,Machines{i,4}{2}{:});
        if size(Machines{i,9},2)>1
            x=num2cell(x,3);
            F=cellfun(@(x) scatteredInterpolant(t,squeeze(x),'linear','none'),x,'UniformOutput',false);
            Machines{i,8}=cell2mat(cellfun(@(x) permute(x(ambprof),[2 3 1]),F,'UniformOutput',false));
            tempvar=Machines{i,8}(:,1,:);
            tempvar(isnan(tempvar))=0;
            Machines{i,8}(:,1,:)=tempvar;
            tempvar=Machines{i,8}(:,2:end,:);
            tempvar(isnan(tempvar))=0;
            Machines{i,8}(:,2:end,:)=tempvar;
            [x1 x2 x3] = size(Machines{i,8});
            Machines{i,8}=squeeze(mat2cell(Machines{i,8}, x1, x2, ones(1,x3)));
%             Machines{i,8}=cellfun(@(x) mat2cell(x),Machines{i,8},'UniformOutput',false);
        elseif size(Machines{i,9},2)==1
            x = permute(x,[3 1 2]);
            C=interp1(t,x,ambprof,'linear',0);
            C = permute(C,[2 3 1]);
            C=mat2cell(C,size(C,1),size(C,2),ones(size(C,3),1));
            Machines{i,8}=squeeze(C);
        end
    else
        [C(:)]=deal(Machines{i,4}(2));
        Machines{i,8}=squeeze(C);
    end
end

%StartUp cost vector
SUcosts=cell2mat(cellfun(@(x) x(4),Machines(:,7),'UniformOutput',false));

%Inputs lists all possible machine inputs. 
Inputs=[Machines{:,2}];
Inputs=unique(Inputs, 'Stable');    %Gets rid of duplicates

%All possible outputs (without duplicates)
                                                        %Now defined on the basis of demand profiles, but what if there is an
                                                        %internal consumption that you don't see in the demand profiles?
Outputs=([Dall{:,1}])';
% Outputs=unique([Machines{:,3}]);                                                          

Noutputs=size(Outputs,2);

%All possible goods (without duplicates)
% Goods=unique([Inputs Outputs],'stable');
slope=cell(Nmachines,1);
intercept=cell(Nmachines,1);

% Convexity check
for i=1:Nmachines
    for h=1:ntimestot                            %for each time instant
        for f=1:(J(i)-1)                    %for each segment
            for k=1:numel(Machines{i,3})    %and for each output (NB:we assume we have only one input)
                if convcheck
                    slope{i}(f,k,h)=(Machines{i,8}{h}(f+1,1)-Machines{i,8}{h}(f,1))/(Machines{i,8}{h}(f+1,1+k)-Machines{i,8}{h}(f,1+k));
                    intercept{i}(f,k,h)=Machines{i,8}{h}(f+1,1)-slope{i}(f,k,h)*Machines{i,8}{h}(f+1,1+k);
                else
                    slope{i}(f,k,h)=0;
                    intercept{i}(f,k,h)=0;
                end
            end
        end
    end
    convexflag(i)=all(all(all((slope{i}(2:end,:,:)>=slope{i}(1:(end-1),:,:)))));
end



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

%sarebbe bello avere reference variabile su inizio nomi storage ma vabb�
i=0;
Storages={[]};
Worksheets.Item('Stor&Net').Activate;
storname=get(ExcelApp.ActiveSheet,'Range',strcat('A',num2str(3+i)));
storname=storname.value;
% [~,storname]=xlsread(Filepath,'Stor&Net',strcat('A',num2str(3+i)));
while ~isnan(storname)
    i=i+1;
    Storages(i,1)={storname};
    temp=get(ExcelApp.ActiveSheet,'Range',strcat('B',num2str(2+i),':J',num2str(2+i)));
    storvals=cell2mat(temp.value);
%     storvals=xlsread(Filepath,'Stor&Net',strcat('B',num2str(2+i),':J',num2str(2+i)));
    for j=2:10
        Storages{i,j}=storvals(j-1);
    end
    storname=get(ExcelApp.ActiveSheet,'Range',strcat('A',num2str(3+i)));
    storname=storname.value;
%     [~,storname]=xlsread(Filepath,'Stor&Net',strcat('A',num2str(3+i)));
end    

Nstorages=i;

Networksall={[] [] [] 0 0};
%Networks data: 
%column1 --> network's good
%column2 --> max withdrawal rate
%column3 --> max injection rate
i=0;
netname=get(ExcelApp.ActiveSheet,'Range',strcat('L',num2str(3+i)));
netname=netname.value;
% [k,netname]=xlsread(Filepath,'Stor&Net',strcat('L',num2str(3+i)));
while ~isnan(netname)
    i=i+1;
    Networksall(i,1)={netname};
    netvals{1}=get(ExcelApp.ActiveSheet,'Range',strcat('M',num2str(2+i),':M',num2str(2+i)));
    netvals{1}=netvals{1}.value;
%     netvals{1}=xlsread(Filepath,'Stor&Net',strcat('M',num2str(2+i),':M',num2str(2+i)));
    netvals{2}=get(ExcelApp.ActiveSheet,'Range',strcat('N',num2str(2+i),':N',num2str(2+i)));
    netvals{2}=netvals{2}.value;
%     netvals{2}=xlsread(Filepath,'Stor&Net',strcat('N',num2str(2+i),':N',num2str(2+i)));
    for j=2:3
        if ~ischar(netvals{j-1})
            Networksall{i,j}=netvals{j-1};
        else
            Networksall{i,j}=Inf;
        end        
    end
    
    netname=get(ExcelApp.ActiveSheet,'Range',strcat('L',num2str(3+i)));
    netname=netname.value;
%     [k,netname]=xlsread(Filepath,'Stor&Net',strcat('L',num2str(3+i)));
end    

Nnetworks=i;

for i=1:Noutputs
M(i)=max(Dall{2}(i,:))*100;  %M1 related to maximum good demand (an order of magnitude higher
end
Mbigspare=max(M);

%controllo qualit� su customizzazione big M
% M(:)=Mbigspare;
%Risultato del test: miglioramento dell'1.4% nella velocit� di calcolo! Immagino diventi davvero rilevante solo quando gli ordini di grandezza nelle domande dei vari beni sono molto diverse 

% strcmp(D{1},Networks{1,1})    %Comando fondamentale confronto stringa / vettore
% ismember(D{1},Machines{i,3}) nel caso in cui siano due vettori

%Setting exchange limits on various networks: limit might be set by user or
%will be set to big M of network good
for i=1:Nnetworks      
    if ~isnumeric(Networksall{i,2})||~isfinite(Networksall{i,2})  %if limit for exchange is infinite or not a number
        if all(~(strcmp(Dall{1},Networksall{i,1})))~=0
            Networksall{i,2}=M(strcmp(Dall{1},Networksall{i,1}));
        else
            Networksall{i,2}=Mbigspare;
        end
    end
    if ~isnumeric(Networksall{i,3})||~isfinite(Networksall{i,3})  %if limit for exchange is infinite or not a number
        if all(~(strcmp(Dall{1},Networksall{i,1})))~=1
            Networksall{i,3}=M(strcmp(Dall{1},Networksall{i,1}));
        else
            Networksall{i,3}=Mbigspare;
        end
    end
end


ExcelApp.Workbooks.Item(Filename).Save 
ExcelApp.Quit;
ExcelApp.release;


x=find(~cellfun(@isempty,Pricetags(1,:)));  %check whether we have both networks and fuels

if length(x)==2         %both fuels and networks
    for i=1:x(2)-1      %these are the fuels and the prices
        Fuelsall{i,1}=Pricetags(2,i);
        Fuelsall{i,2}=Prices(:,i);
    end
    for i=1:Nnetworks
        a=find(strcmp(Networksall{i,1},Pricetags(2,x(2):end)))-1; %column for pricing of network i
        Networksall{i,4}=Prices(:,x(2)+a);
        Networksall{i,5}=Prices(:,x(2)+a+1);
    end
else                                    %only fuel prices
    for i=1:length(Pricetags(2,:))      %these are the fuels and the prices
        Fuelsall{i,1}=Pricetags(2,i);
        Fuelsall{i,2}=Prices(:,i);
    end
end   
Nfuels=size(Fuelsall,1);
% 
system('taskkill /F /IM EXCEL.EXE');            %alterantives to brutally murder all excel tasks?