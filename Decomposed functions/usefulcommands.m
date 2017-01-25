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