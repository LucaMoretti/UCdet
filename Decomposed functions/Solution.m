%%%%%%%%%%%%%%%%%%%%
% PROBLEM SOLUTION %
%%%%%%%%%%%%%%%%%%%%
% ops = sdpsettings('solver','gurobi','gurobi.MIPGap',0.005);
% sol=optimize(Constr,Objective,ops)
sol=Model{Param};
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
OnOffHist=sol(en);
en=en+1;
fuelusage=sol{en};
