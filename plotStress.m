%run <x>reallign first

figure(2)
hold on
rates = [];
for i = 1:length(chargeCycles)
    thisChg = chargeCycles{i};
    if min(thisChg(:,6)) < 0
        %hsi = find(thisChg(:,6)< -0.097) %find the row corresponding to the start of the bottom hold (plating)
        hsi = length(thisChg);
    else
        hsi = find(thisChg(:,6)<0.013) %find the row corresponding to the start of the bottom hold
    end
    crate = round(abs(thisChg(2,7)/C_rate),2);
    t = thisChg(1,1);
    times = (thisChg(:,1)-t).*crate;
    times = times(1:hsi);
    startStrs = thisChg(1,2);%build stress differences for each section of cycle
    stress = thisChg(:,2)-startStrs;
    stress = stress(1:hsi);
    if hsi(1) == 1 || sum(isnan(stress)) == length(stress)
        continue
    end
    timeStep = times(2) - times(1);
    if timeStep < 9
        times = round(downsample(times,round(10/timeStep)),-1);
        stress = downsample(stress,round(10/timeStep));
    end
    CR = round(crate*2);
    rates(i) = round(crate,1);
    CR = CR+1
    clr = {'r-','g-','b-','c-','m-','y-','k-'};
    plot(times,stress,clr{CR})
    names = ["times",strcat("stress ",string(i)," C: ",string(crate))];
    NormStressCell{i} = table(times,stress,'VariableNames',names);
end
NormStressCell=NormStressCell(~cellfun('isempty',NormStressCell));
NormStress = NormStressCell{1};

for i = 2:length(NormStressCell)
    NormStress = outerjoin(NormStress,NormStressCell{i},'MergeKeys',true);
end
%change the pink stuff in quotes to whatever you want the filename to be
%writetable(NormStress,'XXXNormalizedStresses.csv')
title('Stress vs. SOC for various C-rates')
'C-rates are:'
string(unique(rates))
