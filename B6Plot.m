%MUST RUN B6reallign FIRST
%plot things from each half cycle (can just copy and paste this section
%once this has been run once)
sz = [length(chargeCycles) 12];
varNames = {'cycle_number','chg_stress_change_start','chg_stress_change_end','dischg_stress_change_start',...
    'dischg_stress_change_end','low_hold_stress_change','chgcap_start_mAhg','chgcap_end_mAhg','dischgcap_start_mAhg',...
    'dischgcap_end_mAhg','charge_current','c_rate','chg_stress_change_start_norm','chg_stress_change_end_norm',...
    'dischg_stress_change_start_norm','dischg_stress_change_end_norm'};
chronoArray = zeros(sz);
for i = 1:length(chargeCycles)
    thisChg = chargeCycles{i};
    thisDis = dischargeCycles{i+1};
    hsi = find(thisChg(:,6)<.011);%find the row corresponding to the start of the bottom hold
    dhsi = find(thisDis(:,6)>1.49);% top hold
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
    %DISCHARGE
    dischgStartStrs = thisDis(1,2);
    dischgHoldStartStrs = thisDis(dhsi(1),2);
    dischgHoldEndStrs = thisDis(end,2);
    if isnan(dischgHoldEndStrs)%use last available data if final stress value is NaN
        lastIndex = sum(~isnan(thisDis(:,2)));
        if lastIndex > 0
            dischgHoldEndStrs = thisDis(lastIndex,2);
        end
    end
    dischgStrsChangeStart = abs(dischgHoldStartStrs - dischgStartStrs);
    dischgStrsChangeEnd = abs(dischgHoldEndStrs - dischgStartStrs);
    crate = round(abs(thisChg(2,7)/C_rate),2);
    chronoArray(i,:) = [i, chgStrsChangeStart,chgStrsChangeEnd, dischgStrsChangeStart,dischgStrsChangeEnd,lowHoldStrs,...
        thisChg(hsi(1),9)/mass, thisChg(end,9)/mass, thisDis(dhsi(1),8)/mass, thisDis(end,8)/mass,thisChg(2,7),crate];
end
%make normalized stress collumns
chronoArray(:,13:16) = chronoArray(:,2:5)./chronoArray(:,7:10);
chrono = array2table(chronoArray,'VariableNames',varNames)
%writetable(chrono,fullfile(projdir,'B6Chronology.csv'))


figure(6)
yyaxis left
plot(chrono.cycle_number,chrono.chg_stress_change_start,'-bo')
hold on
plot(chrono.cycle_number,chrono.chg_stress_change_end,'-bs')
plot(chrono.cycle_number,chrono.dischg_stress_change_start,'-ro')
plot(chrono.cycle_number,chrono.dischg_stress_change_end,'-rs')
plot(chrono.cycle_number,chrono.low_hold_stress_change,'-go')
ylabel('Stress MPa')
yyaxis right
plot(chrono.cycle_number,chrono.chgcap_start_mAhg,'-b*')
plot(chrono.cycle_number,chrono.chgcap_end_mAhg,'-b+')
plot(chrono.cycle_number,chrono.dischgcap_start_mAhg,'-r*')
plot(chrono.cycle_number,chrono.dischgcap_end_mAhg,'-r+')
hold off

legend('Charge Stress Change St','Charge Stress Change End','Discharge Stress Change St','Discharge Stress Change End',...
    'Low Hold Stress Change','Charge Capacity Start','Charge Capacity End','Discharge Capacity Start','Discharge Capacity End')
xlabel('Cycle Number')
ylabel('Capacity mAh/g')
title('Baseline 4 Cycles')
%{
figure (7)
plot(chrono.cycle_number,chrono.chg_stress_change_start_norm,'-bo')
hold on
plot(chrono.cycle_number,chrono.chg_stress_change_end_norm,'-bs')
plot(chrono.cycle_number,chrono.dischg_stress_change_start_norm,'-ro')
plot(chrono.cycle_number,chrono.dischg_stress_change_end_norm,'-rs')
legend('Charge Stress Change St','Charge Stress Change End','Discharge Stress Change St','Discharge Stress Change End')
%}