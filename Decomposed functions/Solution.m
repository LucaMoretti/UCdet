%% PROBLEM SOLUTION AND VARIABLES VALUE EXTRACTION

sol=Model{Param};

if symtype==1||symtype==2
    st=1;
    en=Nmachines;
    INPUT=sol(st:en);
    st=Nmachines+1;
    en=2*Nmachines;
    OUTPUT=sol(st:en);
    en=en+1;
    STORAGEcharge=sol{en};
    en=en+1;
    STORAGEpower=sol{en};
    en=en+1;
    NETWORKbought=sol{en};
    en=en+1;
    NETWORKsold=sol{en};
    en=en+1;
    Diss=sol{en};
    en=en+1;
    slacks=sol{en};
    en=en+1;
    OnOffHist=sol{en};
    en=en+1;
    fuelusage=sol{en};
    en=en+1;
    startflags=sol{en};
    STORstart=sol{en+2}';
elseif  symtype==3
    st=1;
    en=Nmachines;
    INPUT=cellfun(@(x) x(1:roladvance),sol(st:en),'UniformOutput',false);
    st=Nmachines+1;
    en=2*Nmachines;
    OUTPUT=cellfun(@(x) x(:,1:roladvance),sol(st:en),'UniformOutput',false);
    en=en+1;
    if Nstorages~=0
        STORAGEcharge=sol{en}(:,1:roladvance);
        STORstart=sol{en}(:,min(roladvance+1,end));
        STORstart=fix(STORstart*1e5)/1e5;
        en=en+1;
        STORAGEpower=sol{en}(:,1:roladvance);
    else
        en=en+1;
    end
    en=en+1;
    if Nnetworks~=0  
        NETWORKbought=sol{en}(:,1:roladvance);
        en=en+1;
        NETWORKsold=sol{en}(:,1:roladvance);
        en=en+1;
    else
        en=en+2;
    end
    Diss=sol{en}(:,1:roladvance);
    en=en+1;
    slacks=sol{en}(:,1:roladvance);
    en=en+1;
    OnOffHist=sol{en};
    en=en+1;
    fuelusage=sol{en}(:,1:roladvance);
    en=en+1;
    startflags=sol{en}(:,1:roladvance);
end

LastProd=cellfun(@(x) x(1,end),OUTPUT(:),'UniformOutput',false);
LastProd=[LastProd{:}];