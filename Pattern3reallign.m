close all;clear
format long;

%taken from code Catherine provided me
hs = 500e-6; %microns current substrate 
C_angle = 0.78173; %[ cos(a)/ 2L ]
Bs = 86.38; %GPa for copper substrate (from McMaster)
hf = 43e-6; %film thickness
projdir = 'Pattern3';
mass = .0123;
C_rate = 3.69;

MOSSFileList = dir(fullfile(projdir,'MOSS','*.txt'));
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
    MList{k} = readmatrix(fullfile(projdir,'MOSS',MOSSFileList(k).name)); %load data
    MList{k}(:,2) = (Bs*hs^2*MList{k}(:,2)*10*C_angle)/(6*hf);%Stoney
    MList{k}(:,1) = round(MList{k}(:,1)+MOSSlineups(k,2),-1); %add time offsets for each file
    [~,uid] = unique(MList{k}(:,1),'stable');%remove rows with same time entry
    MList{k} = array2table(MList{k}(uid,:),'VariableNames',{'time_s','stress_MPa'});
end
Mtable = vertcat(MList{1:M});

P3master = outerjoin(Mtable,Etable,'MergeKeys',true);
% UNCOMMENT THIS TO MAKE A CSV FILE
writetable(P3master,fullfile(projdir,'Pattern3Stress.csv'));

%make list of half cycle indices
[l hCycIndex] = findgroups(P3master.halfCycle);
if hCycIndex(1)==0
    hCycIndex = hCycIndex(2:end);
end
%build subtables from each half cycle and just charge cycles
for i = hCycIndex'
    thisCycle = find(P3master.halfCycle==i);
    startVals = P3master(thisCycle(1),:);
    halfCycles{i} = P3master(thisCycle,:);
    %halfCycles{i}.time_s = halfCycles{i}.time_s - startVals.time_s;
    if startVals.ox_red == 0
        chargeCycs{i} = table2array(P3master(thisCycle,:));
    elseif startVals.ox_red == 1
        dischargeCycs{i} = table2array(P3master(thisCycle,:));
    end
end
empInd = cellfun(@isempty, chargeCycs) == 0;
chargeCycles = chargeCycs(empInd);
DempInd = cellfun(@isempty, dischargeCycs) == 0;
dischargeCycles = dischargeCycs(DempInd);

%DELETE THIRD ELEMENT BECAUSE REASONS
chargeCycles(3)=[];

totC = max(P3master.cycleNumber)+1;
Cycles = cell(1,totC);
%make list of FULL cycle indices, defined by EC labs as a climb and fall of
%voltage, so opposite of what we are considering. 
[l fCycIndex] = findgroups(P3master.cycleNumber);
if fCycIndex(1)==0
    fCycIndex = fCycIndex(2:end);
end
%build subtables from each FULL cycle and rezero time
for i = fCycIndex'
    thisCycle = find(P3master.cycleNumber==i);
    startVals = P3master(thisCycle(1),:);
    Cycles{i} = P3master(thisCycle,:);
    %Cycles{i}.time_s = Cycles{i}.time_s - startVals.time_s;
end

figure(1)
plot(P3master.time_s/3600,P3master.stress_MPa);
set(gca,'FontSize',26)
xlabel('Time (Hrs)'); ylabel('Stress (MPa)');
hold on;
plot(P3master.time_s/3600,P3master.x_I__mA)
plot(P3master.time_s/3600,P3master.Ewe_V)
hold off;

