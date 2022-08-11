sample = "C7";
figure(6)
yyaxis left
plot(chrono.cycle_number(1:7),chrono.irrev_stress(1:7),'-bo')
ylabel('Stress MPa')
yyaxis right
plot(chrono.cycle_number(1:7),chrono.irrev_cap(1:7),'-ro')
grid on
legend('Irreversible Stress','Irreversible Capacity')
xlabel('Cycle Number')
ylabel('Capacity mAh/g')
title(strcat(sample,' Cycle Summary'))
saveas(gcf,strcat(sample,' Irrev'),'png');

figure(7)
c = chrono.c_rate(1:7);
scatter(chrono.irrev_cap(1:7),chrono.irrev_stress(1:7),[],c,'filled')
xlabel('Irreversible Capacity (mAh/g)')
ylabel('Irreversible Stress (MPa)')
title(strcat(sample,' Cycles colored by C-Rate'))
cb = colorbar;
saveas(gcf,strcat(sample,' Scatter'),'png');