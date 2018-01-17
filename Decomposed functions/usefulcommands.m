x=cat(3,Machines{i,4}{2}{:});
x = permute(x,[3 1 2]);
t=Machines{i,4}{1};
C=interp1(t,x,Tprof);
C = permute(C,[2 3 1]);
C=mat2cell(C,size(C,1),size(C,2),ones(size(C,3),1));
Machines(i,8)=squeeze(C);



    for h=1:ntimes                            %for each time instant
        for f=1:(J(i)-1)                    %for each segment
            for k=1:numel(Machines{i,3})    %and for each output (NB:we assume we have only one input)
                slope(f,k,h)=(Machines{i,8}{h}(f+1,1)-Machines{i,8}{h}(f,1))/(Machines{i,8}{h}(f+1,1+k)-Machines{i,8}{h}(f,1+k));
                intercept(f,k,h)=Machines{i,8}{h}(f+1,1)-slope(f,k)*Machines{i,8}{h}(f+1,1+k);
            end
        end
    end
    
    all(all(all((slope(2:end,:,:)>=slope(1:(end-1),:,:)))))
    
T=24;
delt=1/2;
N=1;
n=10;
while n(end)>1
    N=N+1;
    n=0;
    Tsp=[];
    for i=1:N
        n(i)=T/N/delt/2^(i-1);
        deli(i)=delt*2^(i-1);
        Tsp=[Tsp ones(1,round(n(i)))*deli(i)];
    end
end

N=0;
Ttot=0;
while Ttot<T
    N=N+1;
    for i=1:N
        n(i)=2^(i-1);
    end
    Ttot=n(end)*delt*N;
end
    
% n(i+1)=floor(T/N/delt/2^(i));
% deli(i+1)=(T-sum(Tsp))/n(end);
% Tsp=[Tsp ones(1,round(n(end)))*deli(end)];

data=ones(48/delt);

Time=datetime(1900,1,1,0,0,0):hours(1):datetime(1900,12,31,23,59,0);
TT=timetable(Time',D{2}');
Time2=datetime(1900,1,1,0,0,0):minutes(1):datetime(1900,12,31,23,59,0);
TT=retime(TT,Time2,'linear');
D{2}=TT.Var1(:)';

TT=timetable(Time',Undisp);
TT=retime(TT,Time2,'linear');
Undisp=TT.Undisp;

%% Dispatchment analysis
Period_in=1440;
Period_fin=1440*2;
Tot_cons=sum(cellfun(@(x)sum(x(Period_in:Period_fin)),costvals,'UniformOutput',true))/100;
DG_op=[mean(Pmat{1,3}(1,Pmat{1,3}(1,Period_in:Period_fin)>0)/60);mean(Pmat{1,3}(2,Pmat{1,3}(2,Period_in:Period_fin)>0)/90);mean(Pmat{1,3}(3,Pmat{1,3}(3,Period_in:Period_fin)>0)/150);0];
DG_temp=[sum(Pmat{1,3}(1,Period_in:Period_fin)>0);sum(Pmat{1,3}(2,Period_in:Period_fin)>0);sum(Pmat{1,3}(3,Period_in:Period_fin)>0)]/60;%;sum(Pmat{1,3}(3,Pmat{1,3}(3,:)>0));sum(Pmat{1,3}(4,Pmat{1,3}(4,:)>0))];
Demand=sum(Dall{1,2}(1,Period_in:Period_fin))/60;
DG_prod=sum(sum(Pmat{1,3}(1:2,Period_in:Period_fin)))/60;
PV_prod=(Demand-DG_prod)/(0.97);
PV_tot=sum(sum([Pmat{2,11}(Period_in:Period_fin) Pmat{3,11}(Period_in:Period_fin)]))/60;
Diss=sum(sum([sum(Pmat{1,8}(Period_in:Period_fin)),sum(Pmat{1,8}(Period_in:Period_fin)),sum(Pmat{1,8}(Period_in:Period_fin))]))/60;
Penetr=(Demand-DG_prod)/PV_tot;
Stor_fin=Pmat{2,9}(1,Period_in:Period_fin)+Pmat{3,9}(1,Period_in:Period_fin);
Outs=[Demand;DG_prod;PV_tot;Penetr;DG_op;Tot_cons];

%Diss=sum(Pmat{1,8}(:))+sum(Pmat{2,8}(:))+sum(Pmat{3,8}(:))/60;

Powerbalance=(sum(sum(Pmat{1,3}(Period_in:Period_fin)))-sum(sum(Pmat{1,13}(Period_in:Period_fin)))-sum(sum((Dall{1,2}(Period_in:Period_fin)))))/60;%-sum(sum([Pmat{1:3,8}]))

%% SOC
figure(7)
plot([1:10081],Pmat{2,9}+Pmat{3,9},'r',[1:10081],Pmat{2,9},'b-',[1:10081],Pmat{3,9},'b-.')

Disch=sum(sum([Pmat{2,4}(Period_in:Period_fin), Pmat{3,4}(Period_in:Period_fin)]))/60;
Charg=sum(sum([Pmat{2,5}(Period_in:Period_fin), Pmat{3,5}(Period_in:Period_fin)]))/60;
title('Battery level summer Batches')
xlabel('time [min]')
ylabel('kWh')
%%
FinalBESS=[Pmat{2,9}(Period_fin),Pmat{3,9}(Period_fin)];
TotBESS=sum(FinalBESS);

%%
plot(Period_in:Period_fin,Pmat{2,9}(Period_in:Period_fin)+Pmat{3,9}(Period_in:Period_fin))