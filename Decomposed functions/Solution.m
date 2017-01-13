%%%%%%%%%%%%%%%%%%%%
% PROBLEM SOLUTION %
%%%%%%%%%%%%%%%%%%%%
tic
ops = sdpsettings('solver','gurobi','gurobi.MIPGap',0.005);
sol=optimize(Constr,Objective,ops)
toc