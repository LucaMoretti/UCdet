%% DATA READING FROM EXCEL FILE

clear
clc

%addpath(genpath('D:\Dottorato\Ottimizzazione\Filone moretti!'))

excelpath=fileparts(pwd());

Filename=strcat(excelpath,'\Input.xlsm');

[~,range]=xlsread(Filename,'Demand','datarange');
[D,Dnames]=xlsread(Filename,'Demand',range{1});
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
    
%Simulation horizon and timestep settings
timestep = xlsread(Filename,'Demand','tdur');   % simulation timestep [h]
ntimes=size(D,1);                               % total number of timesteps
days=ntimes*timestep/24;                        % simulation days

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
 

%Demand for each good: row=good ; columns=timestep
D={{Dnames{:}}' D'}; 

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

%All possible outputs (without duplicates)
                                                        %Now defined on the basis of demand profiles, but what if there is an
                                                        %internal consumption?
Outputs=([D{:,1}])';
% Outputs=unique([Machines{:,3}]);                                                          

Noutputs=size(Outputs,2);

%All possible goods (without duplicates)
% Goods=unique([Inputs Outputs],'stable');
% Ngoods=size(Goods,2);

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


x=find(~cellfun(@isempty,Pricetags(1,:)));  %check whether we have both networks and fuels

if length(x)==2         %both fuels and networks
    for i=1:x(2)-1      %these are the fuels and the prices
        Fuels{i,1}=Pricetags(2,i);
        Fuels{i,2}=Prices(:,i);
    end
    for i=1:Nnetworks
        a=find(strcmp(Networks{i,1},Pricetags(2,x(2):end)))-1; %column for pricing of network i
        Networks{i,4}=Prices(:,x(2)+a);
        Networks{i,5}=Prices(:,x(2)+a+1);
    end
else                                    %only fuel prices
    for i=1:length(Pricetags(2,:))      %these are the fuels and the prices
        Fuels{i,1}=Pricetags(2,i);
        Fuels{i,2}=Prices(:,i);
    end
end   
Nfuels=size(Fuels,2);

system('taskkill /F /IM EXCEL.EXE');
