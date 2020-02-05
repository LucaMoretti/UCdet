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
            Machines{i,8}(:,1,:)=changem(Machines{i,8}(:,1,:),100,NaN);
            Machines{i,8}(:,2:end,:)=changem(Machines{i,8}(:,2:end,:),0,NaN);
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