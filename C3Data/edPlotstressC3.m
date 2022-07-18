%run <x>reallign first

figure(2)
hold on
rates = [];
cnt=1;
for i = 1:length(chargeCycles)
    thisChg = chargeCycles{i};
    if min(thisChg(:,6)) < 0
        %hsi = find(thisChg(:,6)< -0.097) %find the row corresponding to the start of the bottom hold (plating)
        hsi = length(thisChg);
    else
        hsi = find(thisChg(:,6)<0.013) %find the row corresponding to the start of the bottom hold
    end
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
    times = (thisChg(:,1)-t).*crate;
    stress = thisChg(:,2)-startStrs;
    rates(i) = round(crate,1);
    if i == 1
        valhold(cnt)=[rates(i)];
        cntvar=cnt;
    else
        fnd=find(valhold==rates(i));
        if length(fnd) == 0
            cnt=cnt+1;
            valhold(cnt)=[rates(i)];
            cntvar=cnt;
        else
            cntvar=fnd(1);
        end
    end
    clr = {'r-','g-','b-','c-'};
    plot(times,stress,clr{cntvar})
end

title('Stress vs. SOC for various C-rates')
'C-rates are:'
string(unique(rates));

