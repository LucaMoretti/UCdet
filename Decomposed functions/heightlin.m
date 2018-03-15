% creation of Taylor expansions in different points
clear

syms h x0
flux = 1/h;
Texp=taylor(flux,h,'ExpansionPoint',x0,'Order',2);

exprack=[0.000001:0.000001:1];

hdisc=[0.000001:0.000001:1];

    
x0=exprack(1);
funcs=subs(Texp);
epcoef{1} = double(coeffs(funcs));
j=1;
for i=2:length(exprack)
    i
    fluxreal = 1/hdisc(i);
    fluxlin=epcoef{j}(1) + epcoef{j}(2)*hdisc(i);
    if abs(fluxreal-fluxlin)>500
        j=j+1;
        x0=exprack(i);
        funcs=subs(Texp);
        epcoef{j} = double(coeffs(funcs));
    end
end


h=[0.000001:0.000001:1];
plot(h,1./h);
hold on
for i=1:length(epcoef)
    plot(h,epcoef{i}(1)+epcoef{i}(2).*h);
end
 ylim([0 10000])
 xlim([0.3 0.4])