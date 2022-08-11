close all;clear
format long;

%initialize sample-specific constants
hs = 668e-6; %microns current substrate 
C_angle = 0.78173; %[ cos(a)/ 2L ]
Bs = 115; %GPa for copper substrate (from McMaster)
hf = 22e-6; %film thickness
projdir = 'B5Data';
mass = .0158*.92;
OneCRate = 4.36;

MOSSFileList = dir(fullfile(projdir,'*.txt'));
M = size(MOSSFileList,1);
ECFileList = dir(fullfile(projdir,'*.mpt'));
N = size(ECFileList,1);
MOSSlineups = csvread(fullfile(projdir,'MOSStimes.csv'));
EClineups = csvread(fullfile(projdir,'ECtimes.csv'));
MList = cell(1,M);
EList = cell(1,N);
lastlines = zeros(N+1,10);

for j = 1:N
    opts = detectImportOptions(fullfile(projdir,ECFileList(j).name),'FileType','text');
    opts.SelectedVariableNames = {'ox_red','NsChanges','time_s','x_Q_Qo__mA_h','Ewe_V','x_I__mA','QCharge_mA_h','QDischarge_mA_h','halfCycle','cycleNumber'};%select whatever data you want
    EList{j} = readmatrix(fullfile(projdir,ECFileList(j).name),opts);
    EList{j}(:,3) = round(EList{j}(:,3)+EClineups(j,2),-1); %add time offsets for each file
    EList{j}(:,[4 9 10]) = EList{j}(:,[4 9 10])+lastlines(j,[4 9 10]); %increment cumulative columns by previous file value
    lastlines(j+1,:) = EList{j}(end,:);%record final line to adjust subsequent files
    lastlines(j+1,9) = lastlines(j+1,9) +1;%update half cycle numbering
    EList{j} = array2table(EList{j},'VariableNames',opts.SelectedVariableNames);
end
Etable = vertcat(EList{1:N});

for k = 1:M
    MList{k} = readmatrix(fullfile(projdir,MOSSFileList(k).name)); %load data
    MList{k}(:,2) = (Bs*hs^2*MList{k}(:,2)*10*C_angle)/(6*hf);%Stoney
    MList{k}(:,1) = round(MList{k}(:,1)+MOSSlineups(k,2),-1); %add time offsets for each file
    [~,uid] = unique(MList{k}(:,1),'stable');%remove rows with same time entry
    MList{k} = array2table(MList{k}(uid,:),'VariableNames',{'time_s','stress_MPa'});
end
Mtable = vertcat(MList{1:M});

B5master = outerjoin(Mtable,Etable,'MergeKeys',true);
% UNCOMMENT THIS TO MAKE A CSV FILE
%writetable(B5master,fullfile(projdir,'B5Stress.csv'));

totHC = max(B5master.halfCycle);
halfCycles = cell(1,totHC);
chargeCycles = cell(1,ceil(totHC/2));
%make list of half cycle indices
[l hCycIndex] = findgroups(B5master.halfCycle);
if hCycIndex(1)==0
    hCycIndex = hCycIndex(2:end);
end
%build subtables from each half cycle and just charge cycles
for i = hCycIndex'
    thisCycle = find(B5master.halfCycle==i);
    thisArr = table2array(B5master(thisCycle,:));
    startVals = thisArr(1,3);
    halfCycles{i} = thisArr;
    if startVals == 0
        chargeCycs{i} = thisArr;
    elseif startVals == 1
        dischargeCycs{i} = thisArr;
    end
end
empInd = cellfun(@isempty, chargeCycs) == 0;
chargeCycles = chargeCycs(empInd);
DempInd = cellfun(@isempty, dischargeCycs) == 0;
dischargeCycles = dischargeCycs(DempInd);

%make list of FULL cycle indices, defined by EC labs as a climb and fall of
%voltage, so opposite of what we are considering. 
[l fCycIndex] = findgroups(B5master.cycleNumber);
if fCycIndex(1)==0
    fCycIndex = fCycIndex(2:end);
end
%build subtables from each FULL cycle and rezero time
for i = fCycIndex'
    thisCycle = find(B5master.cycleNumber==i);
    startVals = B5master(thisCycle(1),:);
    Cycles{i} = B5master(thisCycle,:);
    %Cycles{i}.time_s = Cycles{i}.time_s - startVals.time_s;
end

figure(1)
plot(B5master.time_s/3600,B5master.stress_MPa);
xlabel('Time (Hrs)');
hold on;
plot(B5master.time_s/3600,B5master.x_I__mA)
plot(B5master.time_s/3600,B5master.Ewe_V)
title(projdir)
legend('MOS Stress (MPa)','Current (mA)','Voltage (V)')
hold off;


%plot things from each half cycle (can just copy and paste this section
%once this has been run once)
sz = [length(chargeCycles) 14];
varNames = {'cycle_number','chg_stress_change_start','chg_stress_change_end','dischg_stress_change_start',...
    'dischg_stress_change_end','low_hold_stress_change','chgcap_start_mAhg','chgcap_end_mAhg','dischgcap_start_mAhg',...
    'dischgcap_end_mAhg','charge_current','c_rate','irrev_stress','irrev_cap',...
    'chg_stress_change_start_norm_MPa_mAhg','chg_stress_change_end_norm_MPa_mAhg',...
    'dischg_stress_change_start_norm_MPa_mAhg','dischg_stress_change_end_norm_MPa_mAhg'};
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
    dischgHoldEndStrs = thisDis(end-5,2);
    if isnan(dischgHoldEndStrs)%use last available data if final stress value is NaN
        lastIndex = sum(~isnan(thisDis(:,2)));
        if lastIndex > 0
            dischgHoldEndStrs = thisDis(lastIndex,2);
        end
    end 
    dischgStrsChangeStart = abs(dischgHoldStartStrs - dischgStartStrs);
    dischgStrsChangeEnd = abs(dischgHoldEndStrs - dischgStartStrs);
    %IRREVERSIBLE STRESS
    irrevStrs = dischgHoldEndStrs - startStrs;
    crate = round(abs(thisChg(2,7)/OneCRate),2);
    chronoArray(i,:) = [i, chgStrsChangeStart,chgStrsChangeEnd, dischgStrsChangeStart,dischgStrsChangeEnd,lowHoldStrs,...
        thisChg(hsi(1),9)/mass, thisChg(end,9)/mass, thisDis(dhsi(1),8)/mass, thisDis(end,8)/mass,thisChg(2,7),crate,...
        irrevStrs, (thisChg(end,9) - thisDis(end,8))/mass];
end
%make normalized stress collumns
chronoArray(:,15:18) = chronoArray(:,2:5)./chronoArray(:,7:10);
chrono = array2table(chronoArray,'VariableNames',varNames);


%normalize capacity by fitting c/10 discharge capacity fade curve
figure(5)
cTenths = chrono(find(chrono.c_rate<.11),{'cycle_number','dischgcap_end_mAhg'});
plot(cTenths.cycle_number,cTenths.dischgcap_end_mAhg,'ko')
xlabel('cycle number');
ylabel('capacity (mAh/g)');
modelFun = @(b,x) b(1) + b(2) * exp(-b(3)*x(:, 1));
beta0 = [220,100,.2]; 
mdl = fitnlm(cTenths,modelFun,beta0);
yfit = modelFun(mdl.Coefficients{:,1}',chrono.cycle_number);
hold on
%plot results to check fit
plot(chrono.cycle_number,yfit)
title(strcat(projdir,' Curve Fit of C/10 Cycles'));
legend('data','fit');
hold off

%add normalized capacity and C rate column to chrono table
chrono.adjusted_cap  = yfit;
chrono.actual_c_rate = abs(chrono.charge_current./chrono.adjusted_cap/mass);

% UNCOMMENT THIS TO MAKE A CSV FILE
%writetable(chrono,'Pattern6Chronology.csv')

figure(6)
yyaxis left
plot(chrono.cycle_number,chrono.chg_stress_change_end,'-bo')
hold on
plot(chrono.cycle_number,chrono.dischg_stress_change_end,'-rs')
plot(chrono.cycle_number,chrono.irrev_stress,'-go')
ylabel('Stress MPa')
yyaxis right
plot(chrono.cycle_number,chrono.chgcap_end_mAhg,'-b*')
plot(chrono.cycle_number,chrono.dischgcap_end_mAhg,'-r*')
plot(chrono.cycle_number,chrono.irrev_cap,'-g+')
hold off
grid on
legend('Charge Stress Change End','Discharge Stress Change End','Irreversible Stress',...
    'Charge Capacity End','Discharge Capacity End','Irreversible Capacity')
xlabel('Cycle Number')
ylabel('Capacity mAh/g')
title(strcat(projdir,' Cycle Summary'))