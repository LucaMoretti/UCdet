clear
clc

%SELECTION OF SIMULATION TYPE
% 1 --> single batch
% 2 --> contiguous batches
% 3 --> rolling horizon

symtype = 2;
nbatches = 52;
roltsteps = 10;

%%Convexity check on/off 
convcheck=true;

han = waitbar(0.05,'Simulation initialization');
ReadExcel
suppresswarning=1;

%Variables initialization (required to understand variables structure)
D=Dall;
Fuels=Fuelsall;
UndProd=UndProdall;
Networks=Networksall;

% You might need to act on: ntimes, D, Fuels, Networks, UndProd

waitbar(0.6,han,'Problem solution')

%SINGLE BATCH
if symtype==1
    
    %Simulation time horizon
    ntimes=ntimestot;
    %Creation of problem structure
    Optimization
    %Implementation of parameters values for current simulation instance
    D=Dall;
    Fuels=Fuelsall;
    if Nundisp~=0
        UndProd=UndProdall;
    else
        UndProd=cell(1,3);
        UndProd{1,3}=0;
    end
    OnOffHist=zeros(Nmachines,histdepth);
    Networks=Networksall;
    %Creation of parameters input vector
    Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist};
    %Problem solution and data gathering
    Solution
    DataGathering
    
%CONTIGUOUS BATCHES
elseif symtype==2
    %Calculation of time horizon for each batch (last one is longer)
    tdur(1:(nbatches-1))= floor(ntimestot/nbatches);
    tdur(nbatches)= ntimestot - floor(ntimestot/nbatches)*(nbatches-1);
    tstart=1;
    %Time horizon for all simulation apart from last one
    ntimes=tdur(1);
    
    %Preliminary variables initialization (required to understand variables structure)
    D(2) = cellfun(@(x) x(:,tstart:(tstart+tdur(1)-1)),Dall(2),'UniformOutput',false);
    Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+tdur(1)-1),:),Fuelsall(:,2),'UniformOutput',false);
    if Nundisp~=0
        UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+tdur(1)-1)),UndProdall(:,3),'UniformOutput',false);
    else
        UndProd=cell(1,3);
        UndProd{1,3}=0;
    end
    OnOffHist=zeros(Nmachines,histdepth);
    Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(1)-1)),Networksall(:,4:5),'UniformOutput',false);
    %Creation of problem structure
    Optimization
    
    
    for runcount=1:(nbatches-1)
        
        waitbar(0.6,han,strcat('Solution batch n°',num2str(runcount)))
        
        %Implementation of parameters values for current simulation instance
        D(2) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),Dall(2),'UniformOutput',false);
        Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1),:),Fuelsall(:,2),'UniformOutput',false);
        if Nundisp~=0
            UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),UndProdall(:,3),'UniformOutput',false);
        else
            UndProd=cell(1,3);
            UndProd{1,3}=0;
        end
        OnOffHist=zeros(Nmachines,histdepth);
        Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),Networksall(:,4:5),'UniformOutput',false);
        %Creation of parameters input vector
        Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist};
        %Problem solution and data gathering
        Solution
        DataGathering        
        %starting time update update
        tstart=tstart+tdur(runcount);
        
    end
    
    %final run --> we need to recreate the problem since the size of the
    %simulation has changed!!
    runcount=runcount+1;
    ntimes=tdur(runcount);
    %Preliminary variables initialization (required to understand variables structure)
    D(2) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),Dall(2),'UniformOutput',false);
    Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1),:),Fuelsall(:,2),'UniformOutput',false);
    if Nundisp~=0
        UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),UndProdall(:,3),'UniformOutput',false);
    else
        UndProd=cell(1,3);
        UndProd{1,3}=0;
    end
    OnOffHist=zeros(Nmachines,histdepth);
    Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),Networksall(:,4:5),'UniformOutput',false);
    %Creation of problem structure
    Optimization    
    %Implementation of parameters values for current simulation instance
    D(2) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),Dall(2),'UniformOutput',false);
    Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1),:),Fuelsall(:,2),'UniformOutput',false);
    if Nundisp~=0
        UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),UndProdall(:,3),'UniformOutput',false);
    else
        UndProd=cell(1,3);
        UndProd{1,3}=0;
    end
    OnOffHist=zeros(Nmachines,histdepth);
    Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),Networksall(:,4:5),'UniformOutput',false);
    %Creation of parameters input vector
    Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist};
    %Problem solution and data gathering
    Solution
    DataGathering        
    
elseif symtype==3
    
end

waitbar(0.9,han,'Solution plotting')
Plotting
waitbar(1,han,'Simulation completed!')
close(han) 