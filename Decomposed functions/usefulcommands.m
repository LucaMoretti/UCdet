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
