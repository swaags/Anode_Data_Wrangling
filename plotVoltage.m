%run <x>reallign first

figure(2)
hold on
rates = [];
for i = 1:length(chargeCycles)
    thisChg = chargeCycles{i};
    hsi = find(thisChg(:,6)<.013);%find the row corresponding to the start of the bottom hold
    crate = round(abs(thisChg(2,7)/C_rate),2)
    t = thisChg(1,1);
    times = (thisChg(:,1)-t).*crate;
    times = times(1:hsi);
    voltage = thisChg(:,6);
    voltage = voltage(1:hsi);
    if hsi(1) == 1 || sum(isnan(voltage)) == length(voltage)
        continue
    end
    timeStep = times(2) - times(1);
    if timeStep < 9
        times = round(downsample(times,round(10/timeStep)),-1);
        voltage = downsample(voltage,round(10/timeStep));
    end
    CR = round(crate);
    rates(i) = round(crate,1);
    CR = CR+1;
    clr = {'r-','g-','b-','g-'};
    plot(times,voltage,clr{CR}) 
    names = ["times",strcat("volt ",string(i)," C: ",string(crate))];
    NormVoltCell{i} = table(times,voltage,'VariableNames',names);
end
NormVoltCell=NormVoltCell(~cellfun('isempty',NormVoltCell));
NormVolt = NormVoltCell{1};

for i = 2:length(NormVoltCell)
    NormVolt = outerjoin(NormVolt,NormVoltCell{i},'MergeKeys',true);
end
%change the pink stuff in quotes to whatever you want the filename to be
%writetable(NormVolt,fullfile(projdir,'Pattern6NormalizedVoltages.csv'))
title('Voltage vs. SOC for various C-rates')
'C-rates are:'
string(unique(rates))
