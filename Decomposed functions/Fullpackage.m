clear
clc

%SELECTION OF SIMULATION TYPE
% 1 --> single batch
% 2 --> contiguous batches
% 3 --> rolling horizon

symtype = 2;
nbatches = 52;
roltsteps = 10;
UndProdall=cell(1,3);

han = waitbar(0.05,'Simulation initialization');
ReadExcel
D=Dall;
Fuels=Fuelsall;
UndProd=UndProdall;
Networks=Networksall;
suppresswarning=1;

waitbar(0.4,han,'Problem formulation')
% You might need to act on: ntimes, D, Fuels, days (not even used I think),
% Networks, UndProd

waitbar(0.6,han,'Problem solution')
if symtype==1
    ntimes=ntimestot;
    Optimization
    Solution
    DataGathering
elseif symtype==2
    tdur(1:(nbatches-1))= floor(ntimestot/nbatches);
    tdur(nbatches)= ntimestot - floor(ntimestot/nbatches)*(nbatches-1);
    tstart=1;
    for runcount=1:nbatches
        
        waitbar(0.6,han,strcat('Solution batch n°',num2str(runcount)))
        
        ntimes=tdur(runcount);
        D(2) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),Dall(2),'UniformOutput',false);
        Fuels(:,2) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1),:),Fuelsall(:,2),'UniformOutput',false);
        if Nundisp~=0
            UndProd(:,3) = cellfun(@(x) x(:,tstart:(tstart+tdur(runcount)-1)),UndProdall(:,3),'UniformOutput',false);
        else
            UndProd=UndProdall;
        end
        Networks(:,4:5) = cellfun(@(x) x(tstart:(tstart+tdur(runcount)-1)),Networksall(:,4:5),'UniformOutput',false);
        tic
        Optimization
        toc
        Solution
        DataGathering
        tstart=tstart+tdur(runcount);
    end
        
    
elseif symtype==3
    
end

waitbar(0.9,han,'Solution plotting')
Plotting
waitbar(1,han,'Simulation completed!')
close(han) 