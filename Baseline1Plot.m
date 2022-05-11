%MUST RUN Baseline1reallign FIRST
%plot things from each half cycle (can just copy and paste this section
%once this has been run once)
sz = [length(chargeCycles) 8];
varNames = {'cycle_number','max_stress','start_stress','stress_change','capacity_start_mAhg','capacity_end_mAhg','charge_current','c_rate'};
chronoArray = zeros(sz);
for i = 1:length(chargeCycles)
    hsi = find(chargeCycles{i}(:,6)<.011);
    maxStrs = chargeCycles{i}(hsi(1),2);
    startStrs = chargeCycles{i}(1,2);
    strsChg = maxStrs-startStrs;
    crate = round(abs(chargeCycles{i}(2,7)/C_rate),2);
    chronoArray(i,:) = [i, maxStrs, startStrs, strsChg, chargeCycles{i}(hsi(1),9)/mass, chargeCycles{i}(end,9)/mass,chargeCycles{i}(2,7),crate];
end

chrono = array2table(chronoArray,'VariableNames',varNames)
writetable(chrono,fullfile(projdir,'Baseline1Chronology.csv'))


figure(6)
yyaxis left
plot(chrono.cycle_number,chrono.start_stress,'-bo')
hold on
plot(chrono.cycle_number,chrono.max_stress,'-bs')
plot(chrono.cycle_number,chrono.stress_change,'-b*')
ylabel('Stress MPa')
yyaxis right
plot(chrono.cycle_number,chrono.capacity_start_mAhg,'-o')
plot(chrono.cycle_number,chrono.capacity_end_mAhg,'-s')
hold off
legend('Stress at start of charge','Stress at end of charge','Stress Change','Capacity at hold start','Capacity at hold end')
xlabel('Cycle Number')
ylabel('Capacity mAh/g')
title('Baseline 1 Cycles')

