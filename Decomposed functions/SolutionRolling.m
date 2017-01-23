%%%%%%%%%%%%%%%%%%%%
% PROBLEM SOLUTION %
%%%%%%%%%%%%%%%%%%%%
% ops = sdpsettings('solver','gurobi','gurobi.MIPGap',0.005);
% sol=optimize(Constr,Objective,ops)
sol=Model{Param};
st=1;
en=Nmachines;
INPUT=cellfun(@(x) x(1:roladvance),sol(st:en),'UniformOutput',false);
st=Nmachines+1;
en=2*Nmachines;
OUTPUT=cellfun(@(x) x(:,1:roladvance),sol(st:en),'UniformOutput',false);
en=en+1;
STORAGEcharge=sol{en}(1:roladvance);
en=en+1;
STORAGEpower=sol{en}(1:roladvance);
en=en+1;
NETWORKbought=sol{en}(1:roladvance);
en=en+1;
NETWORKsold=sol{en}(1:roladvance);
en=en+1;
Diss=sol{en}(1:roladvance);
en=en+1;
slacks=sol{en}(1:roladvance);
en=en+1;
OnOffHist=sol(en);
en=en+1;
fuelusage=sol{en}(1:roladvance);

LastProd=cellfun(@(x) x(1,1),OUTPUT(:),'UniformOutput',false);
LastProd=[LastProd{:}];