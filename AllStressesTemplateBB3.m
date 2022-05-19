
%plot things from each half cycle (can just copy and paste this section
%once this has been run once)
sz = [length(chargeCycles) 13];
varNames = {'cycle_number','chg_stress_start','chg_stress_holdstart','chg_stress_change','dischg_stress_start','dischg_stress_holdend',...
    'dischg_stress_change','chgcap_start_mAhg','chgcap_end_mAhg','dischgcap_end_mAhg','dischgcap_start_mAhg','charge_current','c_rate'};
chronoArray = zeros(sz);
for i = 1:length(chargeCycles)
    thisChg = chargeCycles{i};
    thisDis = dischargeCycles{i+1};
    hsi = find(thisChg(:,6)<.011);%find the row corresponding to the start of the bottom hold
    dhsi = find(thisDis(:,6)>1.49);% top hold
    chgHoldStart = thisChg(hsi(1),2);
    if isnan(chgHoldStart)
        lastIndex = sum(~isnan(thisChg(:,2)));
        if lastIndex > 0
            chgHoldStart = thisChg(lastIndex,2);
        end
    end
    startStrs = thisChg(1,2);
    chgStrsChange = abs(chgHoldStart-startStrs);
    dischgStart = thisDis(1,2);
    endStrs = thisDis(end,2);
    if isnan(endStrs)
        lastIndex = sum(~isnan(thisDis(:,2)));
        if lastIndex > 0
            endStrs = thisDis(lastIndex,2);
        end
    end
    dischgStrsChange = abs(endStrs - dischgStart);
    crate = round(abs(thisChg(2,7)/C_rate),2);
    chronoArray(i,:) = [i, startStrs,chgHoldStart, chgStrsChange, dischgStart, endStrs, dischgStrsChange, thisChg(hsi(1),9)/mass,...
        thisChg(end,9)/mass, thisDis(dhsi(1),8)/mass, thisDis(end,8)/mass,thisChg(2,7),crate];
end

chrono = array2table(chronoArray,'VariableNames',varNames)
%writetable(chrono,fullfile(projdir,'BaselineB3Chronology.csv'))


figure(6)
yyaxis left
plot(chrono.cycle_number,chrono.chg_stress_start,'-bo')
hold on
plot(chrono.cycle_number,chrono.chg_stress_holdstart,'-bs')
plot(chrono.cycle_number,chrono.chg_stress_change,'-b*')
plot(chrono.cycle_number,chrono.dischg_stress_start,'-ro')
plot(chrono.cycle_number,chrono.dischg_stress_holdend,'-rs')
plot(chrono.cycle_number,chrono.dischg_stress_change,'-r*')
ylabel('Stress MPa')
%{
yyaxis right
plot(chrono.cycle_number,chrono.chgcap_end_mAhg,'-o')
plot(chrono.cycle_number,chrono.dischgcap_end_mAhg,'-s')
hold off
%}
legend('Charge Stress Start','Charge Stress Hold Start','Charge Stress Change','Discharge Stress Start','Discharge Stress End','Discharge Stress Change')
xlabel('Cycle Number')
%ylabel('Capacity mAh/g')
title('Baseline B3 Cycles')
