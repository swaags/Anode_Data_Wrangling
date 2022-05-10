
%plot things from each half cycle (can just copy and paste this section
%once this has been run once)
sz = [length(chargeCycles) 7];
varNames = {'cycle_number','max_stress','start_stress','stress_change','capacity_start_mAhg','capacity_end_mAhg','charge_current'};
chronoArray = zeros(sz);
for i = 1:length(chargeCycles)
    hsi = find(chargeCycles{i}(:,6)<.011);%find the row corresponding to the start of the hold
    maxStrs = chargeCycles{i}(hsi(1),2);
    startStrs = chargeCycles{i}(1,2);
    strsChg = maxStrs-startStrs;
    chronoArray(i,:) = [i, maxStrs, startStrs, strsChg, chargeCycles{i}(hsi(1),9)/mass, chargeCycles{i}(end,9)/mass,chargeCycles{i}(2,7)];
end

chrono = array2table(chronoArray,'VariableNames',varNames)
writetable(chrono,fullfile(projdir,'Pattern6Chronology.csv'))


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
title('Pattern 6 Cycles')

%{
figure(5)
hold on
for i = hCycIndex'
    if halfCycles{i}.ox_red(2) == 0
        plot(halfCycles{i}.QDischarge_mA_h,halfCycles{i}.stress_MPa)
        %EC labs definition of charge is the opposite of ours
        title('Pattern 6')
        xlabel('Charge mAh')
        ylabel('Stress MPa')
    end
end

figure(4)
hold on
for j = hCycIndex'
    if halfCycles{j}.ox_red(2) == 1
        plot(halfCycles{j}.QCharge_mA_h,halfCycles{j}.stress_MPa)
        %EC labs definition of charge is the opposite of ours
        title('Pattern 6')
        xlabel('Discharge mAh')
        ylabel('Stress MPa')
    end
end





figure(1)
plot(P6master.time_s/3600,P6master.stress_MPa);
set(gca,'FontSize',26)
xlabel('Time (Hrs)'); ylabel('Stress (MPa)');
hold on;
plot(P6master.time_s/3600,P6master.x_I__mA)
plot(P6master.time_s/3600,P6master.Ewe_V)
hold off;

figure(2)
plot(P6master.stress_MPa,P6master.Ewe_V)

figure(3)
plot(P6master.x_Q_Qo__mA_h,P6master.stress_MPa)
%}