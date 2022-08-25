%run <x>reallign first
clear hystArray
%input range of cycles you want to plot:
showCycs = [1:3];
names = cell(1,4*length(chargeCycles));
for i = 1:length(chargeCycles)
    i
    thisChg = chargeCycles{i};
    
    %rezero passed current and stress
    startSOC = thisChg(1,5); 
    startStress = thisChg(1,2);
    thisChg(:,5) = -(thisChg(:,5)-startSOC);
    thisChg(:,2) = thisChg(:,2)-startStress;
    colRange = [4*i-3,4*i-2];
    hystArray(1:numel(thisChg(:,5)),colRange) = [thisChg(:,5),thisChg(:,2)];
    try
        thisDis = dischargeCycles{i+1};
        thisDis(:,5) = -(thisDis(:,5)-startSOC);
        thisDis(:,2) = thisDis(:,2)-startStress;
        %
        %SOC = [thisChg(:,5);thisDis(:,5)];
        %stress = [thisChg(:,2);thisDis(:,2)];
        hystArray(1:numel(thisDis(:,2)),colRange+2) = [thisDis(:,5),thisDis(:,2)];
        names([colRange,colRange+2]) = {sprintf('Charge %d SOC',i),sprintf('Charge %d Stress',i),sprintf('Dischg %d SOC',i),sprintf('Dischg %d Stress',i)};
    catch
        names(colRange) = {sprintf('Charge %d SOC',i),sprintf('Charge %d Stress',i)};
    end
end

for i = showCycs
    last = find(hystArray(:,4*i-3),1,'last');
    figure(i)
    plot(hystArray(1:last,4*i-3),hystArray(1:last,4*i-2))
    hold on
    plot(hystArray(:,4*i-1),hystArray(:,4*i))
    hold off
    xlabel('Current Passed (mAh)');
    ylabel('Stress Change (MPa)');
    legend('Charge','Discharge');
    title(sprintf('Cycle %d',i));
end


%%

hystTable = array2table(hystArray,'VariableNames',names);


%This outputs a CSV of stress and SOC for each cycle as adjacent columns
%change the XXX to avoid overwriting other files
%writetable(hystTable,'XXXHysteresisStress.csv')


