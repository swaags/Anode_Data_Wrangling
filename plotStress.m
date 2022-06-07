

figure(2)
hold on
for i = 1:length(chargeCycles)
    thisChg = chargeCycles{i};
    hsi = find(thisChg(:,6)<.013);%find the row corresponding to the start of the bottom hold
    %CHARGE
    startStrs = thisChg(1,2);%build stress differences for each section of cycle
    chgHoldStartStrs = thisChg(hsi(1),2);
    chgHoldEndStrs = thisChg(end,2);
    if isnan(chgHoldStartStrs)%use last available data if final stress value is NaN
        lastIndex = sum(~isnan(thisChg(:,2)));
        if lastIndex > 0
            chgHoldStartStrs = thisChg(lastIndex,2);
        end
    end
    chgStrsChangeStart = abs(chgHoldStartStrs-startStrs);
    chgStrsChangeEnd = abs(chgHoldEndStrs-startStrs);
    lowHoldStrs = abs(chgHoldStartStrs - chgHoldEndStrs);
    crate = round(abs(thisChg(2,7)/C_rate),2);
    t = thisChg(1,1);
    times = (thisChg(:,1)-t).*crate;
    stress = thisChg(:,2)-startStrs;
    CR = round(crate);
    CR = CR+1;
    clr = {'r-','g-','b-','g-'};
    plot(times(1:hsi),stress(1:hsi),clr{CR})
    
end

    