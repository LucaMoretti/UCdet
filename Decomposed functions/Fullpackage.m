clear
clc
han = waitbar(0.05,'Simulation initialization');
ReadExcel
suppresswarning=1;
waitbar(0.4,han,'Problem formulation')
% You might need to act on: ntimes, D, Fuels, days (not even used I think),
% Networks, Prices, timestep, UndProd
Optimization
waitbar(0.9,han,'Solution plotting')
Plotting
waitbar(1,han,'Simulation completed!')
close(han) 