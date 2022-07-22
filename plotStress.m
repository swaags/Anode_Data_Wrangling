%run <x>reallign first

%choose what SOC you want to see stress values at
SOCInq = .2;
timeInq = SOCInq*3600;

figure(3)
hold on
rates = []; noStress = []; slice = zeros(length(chargeCycles),5); 
for i = 1:length(chargeCycles)
    i
    [~, ia, ~] = unique(chargeCycles{i}(:,1),'stable');
    thisChg = chargeCycles{i}(ia,:);
    %find the row corresponding to the start of the bottom hold
    if min(thisChg(:,6)) < 0
        hsi = length(thisChg);%include hold in plating case
    else
        hsi = find(thisChg(:,6)<0.013);
        hsi = hsi(1);
    end
    crate = chrono.c_rate(i);
    adjCRrate = chrono.actual_c_rate(i);
    t = thisChg(1,1);
    %scale times by c rate (equivalent to SOC at 1C)
    adjTimes = (thisChg(:,1)-t).*adjCRrate;
    adjTimes = adjTimes(1:hsi);
    startStrs = thisChg(1,2);%build stress differences for each section of cycle
    stress = thisChg(:,2)-startStrs;
    stress = stress(1:hsi);
    try
        stressInq = interp1(adjTimes,smoothdata(stress,'Gaussian',40),timeInq);
    catch
        noStress = [noStress, i];
        stressInq = 0;
    end
    slice(i,:) = [i,SOCInq,stressInq,chrono.adjusted_cap(i),chrono.actual_c_rate(i)];
    if hsi(1) == 1 || sum(isnan(stress)) == length(stress)
        continue
    end
    %{
    %downsample the most scrunched runs if desired
    timeStep = adjTimes(2) - adjTimes(1);
    if timeStep < 5
        adjTimes = round(downsample(adjTimes,round(10/timeStep)),-1);
        stress = downsample(stress,round(10/timeStep));
    end
    %}
    CR = round(crate*2);
    CR = CR+1;
    rates(i) = round(crate,2);
    clr = {'r-','g-','b-','c-','m-','y-','k-'};
    
    plot(adjTimes/3600,stress,clr{CR})
    names = ["times",strcat("stress ",string(i)," C: ",string(crate))];
    SOCStressCell{i} = table(adjTimes,stress,'VariableNames',names);
end
SOCStressCell=SOCStressCell(~cellfun('isempty',SOCStressCell));
NormStress = SOCStressCell{1};

sliceTable = array2table(slice,'VariableNames',{'cycle number','SOC','Stress','norm. capacity','norm. C rate'});

for i = 2:length(SOCStressCell)
    NormStress = outerjoin(NormStress,SOCStressCell{i},'MergeKeys',true);
end
%This outputs the entirety of the normalized stress data curves
%change the pink stuff in quotes to whatever you want the filename to be
%writetable(NormStress,'XXXNormalizedStresses.csv')

%This outputs just the stress values at the specified SOC above:
%writetable(sliceTable,'XXXSlice.csv')

line([SOCInq, SOCInq], [2 -6], 'LineStyle','--')
title('Stress vs. SOC for various C-rates')
ylabel('Stress (MPa)')
xlabel('SOC')
'Adjusted C-rates are:'
string(unique(rates))
'No Stress data for the following cycles at that SOC:'
string(noStress)
