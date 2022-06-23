%run <x>reallign first

figure(2)
hold on
rates = [];
for i = 1:length(chargeCycles)
    thisChg = chargeCycles{i};
    hsi = find(thisChg(:,6)<.013);%find the row corresponding to the start of the bottom hold
    crate = round(abs(thisChg(2,7)/C_rate),2);
    t = thisChg(1,1);
    times = (thisChg(:,1)-t).*crate;
    voltage = thisChg(:,6);
    CR = round(crate);
    rates(i) = round(crate,1);
    CR = CR+1;
    clr = {'r-','g-','b-','g-'};
    plot(times(1:hsi),voltage(1:hsi),clr{CR}) 
end
title('Voltage vs. SOC for various C-rates')
'C-rates are:'
string(unique(rates))
