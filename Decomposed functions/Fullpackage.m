%% MAIN ROUTINE

clear
clc

%SELECTION OF SIMULATION TYPE
% 1 --> single batch
% 2 --> contiguous batches
% 3 --> rolling horizon

symtype = 3;

%DATA FOR SYMTYPE #2
nbatches = 2;

%DATA FOR SYMTYPE #3
roltsteps = 24;
roladvance = 1;

%%Convexity check on/off 
convcheck=false;

han = waitbar(0.05,'Simulation initialization');
ReadExcel
suppresswarning=1;

%Variables initialization (required to understand variables structure)
D=Dall;
Fuels=Fuelsall;
UndProd=UndProdall;
Networks=Networksall;
STORstart=0;


% You might need to act on: ntimes, D, Fuels, Networks, UndProd

waitbar(0.4,han,'Problem solution')

                        %%%%%%%%%%%%%%
                        %SINGLE BATCH%
                        %%%%%%%%%%%%%%
if symtype==1 || (symtype==3&&roltsteps>=ntimestot)
    
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
    LastProd=zeros(1,Nmachines);
    Networks=Networksall;
    if Nstorages~=0
        STORstart=cell2mat(Storages(:,10));
    end
    %Creation of parameters input vector
    Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart};
    %Problem solution and data gathering
    Solution
    DataGathering
  
                        %%%%%%%%%%%%%%%%%%%%
                        %CONTIGUOUS BATCHES%
                        %%%%%%%%%%%%%%%%%%%%
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
    Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(1)-1)),Networksall(:,4:5),'UniformOutput',false);
    %Creation of problem structure
    Optimization
    OnOffHist=zeros(Nmachines,histdepth);
    LastProd=zeros(1,Nmachines);    
    
    for runcount=1:(nbatches-1)
        
        waitbar(0.4+0.5/nbatches*runcount,han,strcat('Solution batch n°',num2str(runcount)))
        
        %Implementation of parameters values for current simulation instance
        D(2) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),Dall(2),'UniformOutput',false);
        Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1),:),Fuelsall(:,2),'UniformOutput',false);
        if Nundisp~=0
            UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),UndProdall(:,3),'UniformOutput',false);
        else
            UndProd=cell(1,3);
            UndProd{1,3}=0;
        end
        Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),Networksall(:,4:5),'UniformOutput',false);
        %Creation of parameters input vector
        Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd cell2mat(Storages(:,10))};
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
    Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),Networksall(:,4:5),'UniformOutput',false);
    %Creation of problem structure
    OnOffHisttemp=OnOffHist;
    LastProdtemp=LastProd;
    Optimization    
    OnOffHist=OnOffHisttemp;
    LastProd=LastProdtemp;
    %Implementation of parameters values for current simulation instance
    D(2) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),Dall(2),'UniformOutput',false);
    Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1),:),Fuelsall(:,2),'UniformOutput',false);
    if Nundisp~=0
        UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),UndProdall(:,3),'UniformOutput',false);
    else
        UndProd=cell(1,3);
        UndProd{1,3}=0;
    end
    Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),Networksall(:,4:5),'UniformOutput',false);
    %Creation of parameters input vector
    Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd cell2mat(Storages(:,10))};
    %Problem solution and data gathering
    Solution
    DataGathering        
    
                        %%%%%%%%%%%%%%%%%
                        %ROLLING HORIZON%
                        %%%%%%%%%%%%%%%%%
elseif symtype==3
    
    %Calculation of simulation numbers and duration
    nsims=ceil(ntimestot/roladvance);
    tstart=1;    
    ntimes=roltsteps;
    
    %Preliminary variables initialization (required to understand variables structure)
    D(2) = cellfun(@(x) x(:,tstart:(tstart+roltsteps-1)),Dall(2),'UniformOutput',false);
    Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+roltsteps-1),:),Fuelsall(:,2),'UniformOutput',false);
    if Nundisp~=0
        UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+roltsteps-1)),UndProdall(:,3),'UniformOutput',false);
    else
        UndProd=cell(1,3);
        UndProd{1,3}=0;
    end
    if Nnetworks~=0
        Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+roltsteps-1)),Networksall(:,4:5),'UniformOutput',false);
    end
    %Creation of problem structure
    convexflag=zeros(Nmachines,1);
    coeffs=cellfun(@(x) x(tstart:(tstart+roltsteps-1)),{Machines{:,8}},'UniformOutput',false);
    Optimization
    OnOffHist=zeros(Nmachines,histdepth);
    LastProd=zeros(1,Nmachines);
    %STORAGEcharge=bsxfun(@times,ones(size(STORAGEcharge)),cell2mat(Storages(:,10)));
    if Nstorages~=0
        STORstart=cell2mat(Storages(:,10));
    else
        STORstart=0;
    end
    runcount=0;
    
    while (ntimestot-tstart)>roltsteps;
   
        runcount=runcount+1;
        %Implementation of parameters values for current simulation instance
        D(2) = cellfun(@(x) x(:,tstart:(tstart+roltsteps-1)),Dall(2),'UniformOutput',false);
        Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+roltsteps-1),:),Fuelsall(:,2),'UniformOutput',false);
        if Nundisp~=0
            UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+roltsteps-1)),UndProdall(:,3),'UniformOutput',false);
        else
            UndProd=cell(1,3);
            UndProd{1,3}=0;
        end
        if Nnetworks~=0
            Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+roltsteps-1)),Networksall(:,4:5),'UniformOutput',false);
        else
            Networks{:,4}=0;
            Networks{:,5}=0;
        end
        %Creation of parameters input vector
%         Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart actualcoeffs{:}};
        Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart};
        %Problem solution and data gathering
        Solution
        DataGathering        
        %starting time update update
        tstart=tstart+roladvance;   
        
        OnOffHisttemp=OnOffHist;
        LastProdtemp=LastProd;
        temp=STORstart;
        coeffs=cellfun(@(x) x(tstart:(tstart+roltsteps-1)),{Machines{:,8}},'UniformOutput',false);
        Optimization    
        STORstart=temp;
        OnOffHist=OnOffHisttemp;
        LastProd=LastProdtemp;
        
    end
    
    %while tstart<=ntimestot
        runcount=runcount+1;
        %final runs --> we need to recreate the problem since the size of the
        %rolling horizon exceeds the simulation horizon!!
        minhorizon=max(histdepth+2,6);
        %if ntimestot-roladvance-tstart+1<=minhorizon
            ntimes=ntimestot-tstart+1;
            roladvance=ntimestot-tstart+1;
        %else
        %    ntimes=ntimestot-tstart+1;
        %end
        %Preliminary variables initialization (required to understand variables structure)
        D(2) = cellfun(@(x) x(:,tstart:(tstart+ntimes-1)),Dall(2),'UniformOutput',false);
        Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+ntimes-1),:),Fuelsall(:,2),'UniformOutput',false);
        if Nundisp~=0
            UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+ntimes-1)),UndProdall(:,3),'UniformOutput',false);
        else
            UndProd=cell(1,3);
            UndProd{1,3}=0;
        end
        if Nnetworks~=0
            Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+ntimes-1)),Networksall(:,4:5),'UniformOutput',false);
        else
            Networks{:,4}=0;
            Networks{:,5}=0;
        end
        %Creation of problem structure
        OnOffHisttemp=OnOffHist;
        LastProdtemp=LastProd;
        temp=STORstart;
        coeffs=cellfun(@(x) x(tstart:(tstart+ntimes-1)),{Machines{:,8}},'UniformOutput',false);
        Optimization    
        STORstart=temp;
        OnOffHist=OnOffHisttemp;
        LastProd=LastProdtemp;
        %Implementation of parameters values for current simulation instance
        D(2) = cellfun(@(x) x(:,tstart:(tstart+ntimes-1)),Dall(2),'UniformOutput',false);
        Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+ntimes-1),:),Fuelsall(:,2),'UniformOutput',false);
        if Nundisp~=0
            UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+ntimes-1)),UndProdall(:,3),'UniformOutput',false);
        else
            UndProd=cell(1,3);
            UndProd{1,3}=0;
        end
        if Nnetworks~=0
            Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+ntimes-1)),Networksall(:,4:5),'UniformOutput',false);
        else
            Networks{:,4}=0;
            Networks{:,5}=0;
        end
        %Creation of parameters input vector
        Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart};
        %Problem solution and data gathering
        Solution
        DataGathering   
        tstart=tstart+roladvance;
    %end
    
end

waitbar(0.9,han,'Solution plotting')
Plotting
waitbar(1,han,'Simulation completed!')
close(han) 