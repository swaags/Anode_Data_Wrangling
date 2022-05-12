
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
%writetable(chrono,fullfile(projdir,'Baseline4Chronology.csv'))


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
title('Baseline 4 Cycles')

%{
figure(3)
hold on
for i = cycIndex'
    if halfCycles{i}.ox_red(2) == 0
        plot(halfCycles{i}.QDischarge_mA_h,halfCycles{i}.stress_MPa)
        %EC labs definition of charge is the opposite of ours
        title('Baseline 4')
        xlabel('Charge mAh')
        ylabel('Stress MPa')
    end
end

figure(4)
hold on
for j = cycIndex'
    if halfCycles{j}.ox_red(2) == 1
        plot(halfCycles{j}.QCharge_mA_h,halfCycles{j}.stress_MPa)
        %EC labs definition of charge is the opposite of ours
        title('Baseline 4')
        xlabel('Discharge mAh')
        ylabel('Stress MPa')
    end
end

figure(1)
plot(B4master.time_s/3600,B4master.stress_MPa);
set(gca,'FontSize',26)
xlabel('Time (Hrs)'); ylabel('Stress (MPa)');
hold on;
plot(B4master.time_s/3600,B4master.x_I__mA)
plot(B4master.time_s/3600,B4master.Ewe_V)
hold off;

figure(2)
plot(B4master.stress_MPa,B4master.Ewe_V)
%}