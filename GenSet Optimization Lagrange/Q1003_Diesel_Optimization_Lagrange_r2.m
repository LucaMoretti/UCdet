%---------------------------------%
%---DO NOT CHANGE THIS PROGRAM!---%
%---------------------------------%
%If you have any issues please do ask system engineers!

%------|File name: Q1003_Diesel_Optimization_Lagrange|--------%
%----------|--Revision--|----Date--|----Author--|-------------%
%----------|     0      |2018/24/01|  L.Meraldi |-------------%
%----------|     1      |2018/24/01|  L.Meraldi |-------------%


%This program test standard optimization using Lagrange multipliers
%in MATLAB environment for solving
%minimum fuel consumption of a group of diesel generators

%Revision 1: Spinning Reserve constraints added, number of turned on
%generator automatically chosen

clc
clear variables

%The variables number must be consistent no-checks are performed
%The variables number must be consistent no-checks are performed
DG_set = 100;
DG_SR = 0;
DG_Pn = [60 100 40 60];
DG_PMIN = [20 20 30 10];
DG_Eff_out_perc1 = 0.2*[0 10 21 33 46 60];                                 %consumed fuel @  DG_Eff_in_perc(Percentage with respect to the nominal value)
DG_Eff_out_perc2 = 1.75*[0 5 12 20 30 45];
DG_Eff_out_perc3 = 2*[0 5 12 20 30 45];
DG_Eff_out_perc4 = 0.2*[0 10 21 33 46 60];
DG_Eff_out_perc = [DG_Eff_out_perc1;
    DG_Eff_out_perc2;
    DG_Eff_out_perc3;
    DG_Eff_out_perc4
    ];

DG_Eff_in_perc = [0 0.2 0.4 0.6 0.8 1;
    0 0.2 0.4 0.6 0.8 1;
    0 0.2 0.4 0.6 0.8 1;
    0 0.2 0.4 0.6 0.8 1];

%Input management
DG_num = length(DG_Pn);                                                    %Count number of diesel generator
DG_Eff_seg_num = length(DG_Eff_in_perc)-1;                                 %Count number of linear segments in DG efficiencies curves
DG_Eff_in =  DG_Eff_in_perc.*DG_Pn';
DG_Eff_out = DG_Eff_out_perc;

%Choose configuration to be tried
MaxDG_P = DG_set+DG_SR;
Config = zeros(2^DG_num-1,DG_num);
Config(1,1)=1;
for i=2:2^DG_num-1
    Config(i,:)= Config(i-1,:);
    Config(i,1)= Config(i-1,1) + 1;
    for t=1:DG_num-1
        if (Config(i,t)>1)
            Config(i,t)=0;
            Config(i,t+1)=Config(i,t+1)+1;
        end
    end
end

t=0;
if (MaxDG_P<sum(DG_Pn) && MaxDG_P>min(DG_PMIN))
    for i=1:2^DG_num-1
        DG_Pachiv = Config(i,:)*DG_Pn';
        DG_PMINSet = Config(i,:)*DG_PMIN';
        if (DG_Pachiv>MaxDG_P)
            t=t+1;
            Configuration(t,:) = Config(i,:);
            if (DG_set < DG_PMINSet)
                DG_set_Rew(t) = DG_PMINSet;
            else
                DG_set_Rew(t) = DG_set;
            end
        end
    end
else
    disp('Problem unfeseable');
end
Trials = size(Configuration);

for j=1:Trials(1)
    clear DG_PMIN_Conf DG_Pn_Conf DG_Eff_in_conf DG_Eff_out_conf;
    
    pos = 0;
    for z=1:DG_num
        if Configuration(j,z)==1
            pos = pos+1;
            DG_PMIN_Conf(pos) = DG_PMIN(z);
            DG_Pn_Conf(pos) = DG_Pn(z);
            DG_Eff_in_Conf(pos,:) = DG_Eff_in(z,:);
            DG_Eff_out_Conf(pos,:) = DG_Eff_out(z,:);
        end
    end
    
    DG_num_Conf = length(DG_Pn_Conf);
    
    %Create linear system to be solved
    Coeff = zeros(DG_num_Conf,3);
    
    for i=1:DG_num_Conf
        Coeff(i,:) = polyfit(DG_Eff_in_Conf(i,:),DG_Eff_out_Conf(i,:),2);
    end
    
    %Create constraints & inequalities matrix
    CountConstraints = DG_num_Conf + 1;
    A = zeros (CountConstraints);
    b = zeros (CountConstraints,1);
    
    for i=1:DG_num_Conf
        Coeff(i,:) = polyfit(DG_Eff_in_Conf(i,:),DG_Eff_out_Conf(i,:),2);
        %         x = linspace(DG_Eff_in_Conf(1),DG_Eff_in_Conf(end),100);
        %         y = polyval(Coeff(i,:),x);
        %         plot(x,y,DG_Eff_in_Conf(i,:),DG_Eff_out_Conf(i,:));
        A(i,i)=2*Coeff(i,1);
        b(i)=-Coeff(i,2);
        A(i,DG_num_Conf+1)=-1;
    end
    
    A(DG_num_Conf+1,1:DG_num_Conf)=1;
    b(DG_num_Conf+1)=DG_set_Rew(j);
    
    X = linsolve(A,b);
    Pset = zeros(DG_num_Conf,1);
    Pset = (X(1:DG_num_Conf))';
    Pset_out_max = zeros(1,DG_num_Conf);
    Pset_out_min = zeros(1,DG_num_Conf);
    
    %Check validity of solutions
    for i=1:DG_num_Conf
        if Pset(i)>DG_Pn_Conf(i)
            b(DG_num_Conf+1)=b(DG_num_Conf+1)-DG_Pn_Conf(i);
            A(DG_num_Conf+1,i)=0;
            A(i,i)=1;
            b(i) = 0;
            X = linsolve(A,b);
            Pset(1:DG_num_Conf) = X(1:DG_num_Conf)';
            Pset_out_max(i) = DG_Pn_Conf(i);
        end
    end
    
    for i=1:DG_num_Conf
        if Pset_out_max(i)~=0
            Pset(i) = Pset_out_max(i);
        end
    end
    
    while (sum(Pset<DG_PMIN_Conf)>=1)
        
        for i=1:DG_num_Conf
            if Pset(i)<DG_PMIN_Conf(i)
                b(DG_num_Conf+1)=b(DG_num_Conf+1)-DG_PMIN_Conf(i);
                A(DG_num_Conf+1,i)=0;
                A(i,i)=1;
                b(i) = 0;
                X = linsolve(A,b);
                Pset(1:DG_num_Conf) = X(1:DG_num_Conf)';
                Pset_out_min(i) = DG_PMIN_Conf(i);
            end
        end
        
        for i=1:DG_num_Conf
            if (Pset_out_min(i)~=0)
                Pset(i) = Pset_out_min(i);
                %Pset(i) = Pset_out_max(i);
            end
        end
        
        for i=1:DG_num_Conf
            if Pset_out_max(i)~=0
                Pset(i) = Pset_out_max(i);
            end
        end
        
    end
    
    Flag_Conf = ones(1,DG_num_Conf);
    if (b(DG_num_Conf+1)<0)
        Flag_Conf = Flag_Conf*0;
        for i=1:DG_num_Conf
            A(i,i)=2*Coeff(i,1);
            b(i)=-Coeff(i,2);
            A(i,DG_num_Conf+1)=-1;
            A(DG_num_Conf+1,1:DG_num_Conf)=1;
            b(DG_num_Conf+1)=DG_set_Rew(j);
        end
        for i=1:DG_num_Conf
            if  Pset_out_min(i)~=0
                b(DG_num_Conf+1)=b(DG_num_Conf+1)-DG_PMIN_Conf(i);
                A(DG_num_Conf+1,i)=0;
                A(i,i)=1;
                Flag_Conf(i) = 1;
            end
        end
    end
    X = linsolve(A,b);
    
    for i=1:DG_num_Conf
        if Flag_Conf(i) ==0
            Pset(i) = X(i);
        end
    end
    
    clear Y
    for i=1:DG_num_Conf
        Y(i) = polyval(Coeff(i,:),Pset(i));
    end
    Fuel_Lagrange(j) = sum(Y);
    
    pos=0;
    for z=1:DG_num
        if Configuration(j,z)==1
            pos = pos+1;
            Pset_Matrix_out(j,z) = Pset(pos);
        end
    end
    
end

%Find best solution
clc
Fuel_Lagrange'
Pset_Matrix_out

[Fuel_Lagrange_min,Position] = min(Fuel_Lagrange)
Pset_min = Pset_Matrix_out(Position,:)

% Pset_Lagrange = Pset
%
% Fuel_Lagrange = sum(Y)
%
% Pset_MILP = [60.0000   80.0000  100.0000]
%
% for i=1:DG_num
%     Y(i) = polyval(Coeff(i,:),Pset_MILP(i));
% end
%
% Fuel_MILP = sum(Y)

