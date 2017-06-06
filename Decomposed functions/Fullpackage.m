%% MAIN ROUTINE
clc
read=true;

tic

if read
    clear
    han = waitbar(0.05,'Simulation initialization');
    ReadExcel
else
    han = waitbar(0.05,'Data already read');
end

suppresswarning=1;

varstep=true;

%SELECTION OF SIMULATION TYPE
% 1 --> single batch
% 2 --> contiguous batches
% 3 --> rolling horizon

symtype = 2;

%DATA FOR SYMTYPE #2
nbatches = 52;

%DATA FOR SYMTYPE #3
roltsteps = 24;
roladvance = 5;

%%Convexity check on/off 
convcheck=false;


%Variables initialization (required to understand variables structure)
D=Dall;
Fuels=Fuelsall;
if Nundisp~=0
    UndProd=UndProdall;
else
    UndProd=cell(1,3);
    UndProd{1,3}=0;
end

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
    timestep=ones(ntimes,1)*basetimestep;
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
    actualcoeffs={Machines{:,8}};
    actualslope=slope(:);
    actualintercept=intercept(:);
    actualcoeffs=[actualcoeffs{:}];
    %Creation of parameters input vector
    Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart actualcoeffs{:} actualslope{:} actualintercept{:}};
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
    timestep=ones(ntimes,1)*basetimestep;
    
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
    if Nstorages~=0
        STORstart=cell2mat(Storages(:,10));
    else
        STORstart=0;
    end
    
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
        actualcoeffs=cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),{Machines{:,8}},'UniformOutput',false);
        actualslope=cellfun(@(x) x(:,:,tstart:(tstart+tdur(runcount)-1)),slope(:),'UniformOutput',false);
        actualintercept=cellfun(@(x) x(:,:,tstart:(tstart+tdur(runcount)-1)),intercept(:),'UniformOutput',false);
        actualcoeffs=[actualcoeffs{:}];
        %Creation of parameters input vector
        Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart actualcoeffs{:} actualslope{:} actualintercept{:}};
        %Problem solution and data gathering
        Solution
        DataGathering        
        %starting time update update
        tstart=tstart+tdur(runcount);
        
    end
    
    %final run --> we need to recreate the problem if the size of the
    %simulation has changed!!
    runcount=runcount+1;
    ntimes=tdur(runcount);
    timestep=ones(ntimes,1)*basetimestep;
    
    if tdur(end)~=tdur(end-1)
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
        STORstarttemp=STORstart;
        LastProdtemp=LastProd;
        Optimization    
        STORstart=STORstarttemp;
        OnOffHist=OnOffHisttemp;
        LastProd=LastProdtemp;
    end
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
    actualcoeffs=cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),{Machines{:,8}},'UniformOutput',false);
    actualslope=cellfun(@(x) x(:,:,tstart:(tstart+tdur(runcount)-1)),slope(:),'UniformOutput',false);
    actualintercept=cellfun(@(x) x(:,:,tstart:(tstart+tdur(runcount)-1)),intercept(:),'UniformOutput',false);
    actualcoeffs=[actualcoeffs{:}];
    %Creation of parameters input vector
    Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart actualcoeffs{:} actualslope{:} actualintercept{:}};
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
    if varstep
        T=roltsteps*basetimestep;   %Time horizon [h]
        N=[9 2 2 1];
        deltatimes=[basetimestep basetimestep*2 basetimestep*3 basetimestep*5];
        if sum(N.*deltatimes) ~= T
            error('Selected timestep mesh is not consistent with selected time horizon')
        end 
        nsteps=deltatimes/basetimestep; 
        n=cell2mat(arrayfun(@(i) ones(1,N(i))*deltatimes(i),1:length(N),'UniformOutput',false));
        %n=[n{:}];
        ntimes=sum(N);
        timestep=n'*basetimestep;
    else
        ntimes=roltsteps;
        timestep=ones(ntimes,1)*basetimestep;
    end
    

    
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
    end
    
    
    %Creation of problem structure
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
    
    while (ntimestot-tstart)>roltsteps
   
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

        
        
        %Reshape profiles if variable timestep mesh is on                                 METTI TUTTO IN UNA FUNZIONE A PARTE
        if varstep    
            
            temp=arrayfun(@(i) mean(D{2}(:,(sum(n(1:(i-1)))+1:sum(n(1:i)))),2),1:length(n),'UniformOutput',false);
            D{2}=[temp{:}];
            Fuels(:,2)=cellfun(@(k) cell2mat(k)', cellfun(@(x) arrayfun(@(i) mean(x(sum(n(1:(i-1)))+1:sum(n(1:i)))),1:length(n),'UniformOutput',false),Fuels(:,2),'UniformOutput',false),'UniformOutput',false);
            if Nundisp~=0
                UndProd(:,3)=cellfun(@(k) cell2mat(k), cellfun(@(x) arrayfun(@(i) mean(x(:,sum(n(1:(i-1)))+1:sum(n(1:i))),2),1:length(n),'UniformOutput',false),UndProd(:,3),'UniformOutput',false),'UniformOutput',false);
            end
            if Nnetworks~=0
                Networks(:,4:5)=cellfun(@(k) cell2mat(k)', cellfun(@(x) arrayfun(@(i) mean(x(sum(n(1:(i-1)))+1:sum(n(1:i))),1),1:length(n),'UniformOutput',false),Networks(:,4:5),'UniformOutput',false),'UniformOutput',false);
            end
            
            T4run=Tprof(tstart:(tstart+roltsteps-1));
            T4run=arrayfun(@(i) mean(T4run(sum(n(1:(i-1)))+1:sum(n(1:i)))),1:length(n),'UniformOutput',false);
            T4run=[T4run{:}];
            
            C = cell(1,1,ntimes);
            actualcoeffs=cell(ntimes,Nmachines);
            for i=1:Nmachines
            t=Machines{i,4}{1};
            if size(t,2)>1
                x=cat(3,Machines{i,4}{2}{:});
                x = permute(x,[3 1 2]);
                C=interp1(t,x,T4run);
                C = permute(C,[2 3 1]);
                C=mat2cell(C,size(C,1),size(C,2),ones(size(C,3),1));
            else
                [C(:)]=deal(Machines{i,4}{2});
            end
            actualcoeffs(:,i)=(squeeze(C));          
            end
            
            for i=1:Nmachines
                for h=1:ntimes                            %for each time instant
                    for f=1:(J(i)-1)                    %for each segment
                        for k=1:numel(Machines{i,3})    %and for each output (NB:we assume we have only one input)
                            actualslope{i,1}(f,k,h)=(actualcoeffs{h,i}(f+1,1)-actualcoeffs{h,i}(f,1))/(actualcoeffs{h,i}(f+1,1+k)-actualcoeffs{h,i}(f,1+k));
                            actualintercept{i,1}(f,k,h)=actualcoeffs{h,i}(f+1,1)-slope{i}(f,k,h)*actualcoeffs{h,i}(f+1,1+k);
                        end
                    end
                end
            end
            
        else
            
            actualcoeffs=(cellfun(@(x) x(tstart:(tstart+roltsteps-1)),{Machines{:,8}},'UniformOutput',false));
            actualcoeffs=[actualcoeffs{:}];
            actualslope=cellfun(@(x) x(:,:,tstart:(tstart+roltsteps-1)),slope(:),'UniformOutput',false);
            actualintercept=cellfun(@(x) x(:,:,tstart:(tstart+roltsteps-1)),intercept(:),'UniformOutput',false);
       
        end   

        %Creation of parameters input vector
        Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart actualcoeffs{:} actualslope{:} actualintercept{:}};
        %Problem solution and data gathering
        Solution
        DataGathering        
        %starting time update update
        tstart=tstart+roladvance;   
    end
    
    %while tstart<=ntimestot                    %if data changes it might be better to shrink the horizon until you reach the end of the time horizon
        runcount=runcount+1;
            
%         if varstep
%             N=[9 2 2 1];
%             deltatimes=[basetimestep basetimestep*2 basetimestep*3 basetimestep*5];
%             nsteps=deltatimes/basetimestep; 
%             n=cell2mat(arrayfun(@(i) ones(1,N(i))*deltatimes(i),1:length(N),'UniformOutput',false));
%             ntimes=sum(N);
%             timestep=n'*basetimestep;  
%         else
%         
        %final runs --> we need to recreate the problem since the size of the
        %rolling horizon exceeds the simulation horizon!!
        minhorizon=max(histdepth+2,6);
        %if ntimestot-roladvance-tstart+1<=minhorizon
            ntimes=ntimestot-tstart+1;
            roladvance=ntimestot-tstart+1;
            timestep=ones(ntimes,1)*basetimestep;
%         end
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
        
        
%         if varstep    
%             
%             temp=arrayfun(@(i) mean(D{2}(:,(sum(n(1:(i-1)))+1:sum(n(1:i)))),2),1:length(n),'UniformOutput',false);
%             D{2}=[temp{:}];
%             Fuels(:,2)=cellfun(@(k) cell2mat(k)', cellfun(@(x) arrayfun(@(i) mean(x(sum(n(1:(i-1)))+1:sum(n(1:i)))),1:length(n),'UniformOutput',false),Fuels(:,2),'UniformOutput',false),'UniformOutput',false);
%             if Nundisp~=0
%                 UndProd(:,3)=cellfun(@(k) cell2mat(k)', cellfun(@(x) arrayfun(@(i) mean(x(sum(n(1:(i-1)))+1:sum(n(1:i))),2),1:length(n),'UniformOutput',false),UndProd(:,3),'UniformOutput',false),'UniformOutput',false);
%             end
%             if Nnetworks~=0
%                 Networks(:,4:5)=cellfun(@(k) cell2mat(k)', cellfun(@(x) arrayfun(@(i) mean(x(sum(n(1:(i-1)))+1:sum(n(1:i))),1),1:length(n),'UniformOutput',false),Networks(:,4:5),'UniformOutput',false),'UniformOutput',false);
%             end
%             
%             T4run=Tprof(tstart:(tstart+roltsteps-1));
%             T4run=arrayfun(@(i) mean(T4run(sum(n(1:(i-1)))+1:sum(n(1:i)))),1:length(n),'UniformOutput',false);
%             T4run=[T4run{:}];
%             
%             C = cell(1,1,ntimes);
%             actualcoeffs=cell(ntimes,Nmachines);
%             for i=1:Nmachines
%             t=Machines{i,4}{1};
%             if size(t,2)>1
%                 x=cat(3,Machines{i,4}{2}{:});
%                 x = permute(x,[3 1 2]);
%                 C=interp1(t,x,T4run);
%                 C = permute(C,[2 3 1]);
%                 C=mat2cell(C,size(C,1),size(C,2),ones(size(C,3),1));
%             else
%                 [C(:)]=deal(Machines{i,4}{2});
%             end
%             actualcoeffs(:,i)=(squeeze(C));          
%             end
%             
%             for i=1:Nmachines
%                 for h=1:ntimes                            %for each time instant
%                     for f=1:(J(i)-1)                    %for each segment
%                         for k=1:numel(Machines{i,3})    %and for each output (NB:we assume we have only one input)
%                             actualslope{i,1}(f,k,h)=(actualcoeffs{h,i}(f+1,1)-actualcoeffs{h,i}(f,1))/(actualcoeffs{h,i}(f+1,1+k)-actualcoeffs{h,i}(f,1+k));
%                             actualintercept{i,1}(f,k,h)=actualcoeffs{h,i}(f+1,1)-slope{i}(f,k,h)*actualcoeffs{h,i}(f+1,1+k);
%                         end
%                     end
%                 end
%             end
%             
%         else
            
            actualcoeffs=(cellfun(@(x) x(tstart:(tstart+ntimes-1)),{Machines{:,8}},'UniformOutput',false));
            actualcoeffs=[actualcoeffs{:}];
            actualslope=cellfun(@(x) x(:,:,tstart:(tstart+ntimes-1)),slope(:),'UniformOutput',false);
            actualintercept=cellfun(@(x) x(:,:,tstart:(tstart+ntimes-1)),intercept(:),'UniformOutput',false);
       
%         end
        
%         actualcoeffs=cellfun(@(x) x(tstart:(tstart+ntimes-1)),{Machines{:,8}},'UniformOutput',false);
%         actualcoeffs=[actualcoeffs{:}];
%         actualslope=cellfun(@(x) x(:,:,tstart:(tstart+ntimes-1)),slope(:),'UniformOutput',false);
%         actualintercept=cellfun(@(x) x(:,:,tstart:(tstart+ntimes-1)),intercept(:),'UniformOutput',false);
        %Creation of parameters input vector
        Param={D{2} Fuels{:,2} Networks{:,4:5} UndProd{:,3} OnOffHist LastProd STORstart actualcoeffs{:} actualslope{:} actualintercept{:}};
        %Problem solution and data gathering
        Solution
        DataGathering   
        tstart=tstart+roladvance;
    %end
    
end

totalfo=sum(Obj)

waitbar(0.9,han,'Solution plotting')
Plotting
waitbar(1,han,'Simulation completed!')
close(han) 
toc
