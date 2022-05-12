
%plot things from each half cycle (can just copy and paste this section
%once this has been run once)
sz = [length(chargeCycles) 10];
varNames = {'cycle_number','max_stress','start_stress','stress_change','chgcap_start_mAhg','chgcap_end_mAhg','dischgcap_end_mAhg','dischgcap_start_mAhg','charge_current','c_rate'};
chronoArray = zeros(sz);
for i = 1:length(chargeCycles)
    hsi = find(chargeCycles{i}(:,6)<.011);%find the row corresponding to the start of the bottom hold
    dhsi = find(dischargeCycles{i+1}(:,6)>1.49);% top hold
    maxStrs = chargeCycles{i}(hsi(1),2);
    startStrs = chargeCycles{i}(1,2);
    strsChg = maxStrs-startStrs;
    crate = round(abs(chargeCycles{i}(2,7)/C_rate),2);
    chronoArray(i,:) = [i, maxStrs, startStrs, strsChg, chargeCycles{i}(hsi(1),9)/mass, chargeCycles{i}(end,9)/mass, dischargeCycles{i+1}(dhsi(1),8)/mass, dischargeCycles{i+1}(end,8)/mass,chargeCycles{i}(2,7),crate];
end

chrono = array2table(chronoArray,'VariableNames',varNames)
writetable(chrono,fullfile(projdir,'BaselineB3Chronology.csv'))


figure(6)
yyaxis left
plot(chrono.cycle_number,chrono.start_stress,'-bo')
hold on
plot(chrono.cycle_number,chrono.max_stress,'-bs')
plot(chrono.cycle_number,chrono.stress_change,'-b*')
ylabel('Stress MPa')
yyaxis right
plot(chrono.cycle_number,chrono.chgcap_end_mAhg,'-o')
plot(chrono.cycle_number,chrono.dischgcap_end_mAhg,'-s')
hold off
legend('Stress at start of charge','Stress at end of charge','Stress Change','Charge Capacity (end)','Discharge Capacity (end)')
xlabel('Cycle Number')
ylabel('Capacity mAh/g')
title('Baseline B3 Cycles')
