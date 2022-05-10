close all;clear
format long;

%initialize sample-specific vars
hs = 513e-6;%microns current substrate 
%THIS NEEDS TO BE CHANGED FOR JUNY'S MOSS MEASUREMENTS!!!
C_angle = 0.78173;%[ cos(a)/ 2L ]
Bs = 115;%GPa for copper substrate (from McMaster)
hf = 24e-6;%film thickness
projdir = 'Pattern6Data';
mass = .0155*.92;
C_rate = 4.2;

MOSSFileList = dir(fullfile(projdir,'*.txt'));
M = size(MOSSFileList,1);
ECFileList = dir(fullfile(projdir,'*.mpt'));
N = size(ECFileList,1);
MOSSlineups = csvread(fullfile(projdir,'P6MOSStimes.csv'));%load time offset files
EClineups = csvread(fullfile(projdir,'P6ECtimes.csv'));
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
    MList{k}(:,2) = (Bs*hs^2*MList{k}(:,2)*10*C_angle)/(6*hf);%Stoney /100 for percent *1000 to get MPa
    MList{k}(:,1) = round(MList{k}(:,1)+MOSSlineups(k,2),-1); %add time offsets for each file, round to nearest 10 s
    [~,uid] = unique(MList{k}(:,1),'stable');%remove rows with same time entry
    MList{k} = array2table(MList{k}(uid,:),'VariableNames',{'time_s','stress_MPa'});
end
Mtable = vertcat(MList{1:M});

P6master = outerjoin(Mtable,Etable,'MergeKeys',true); %join on time collumn, hence rounding to nearest 10 s
% UNCOMMENT THIS TO MAKE A CSV FILE
%writetable(P6master,fullfile(projdir,'Pattern6Stress.csv'))

totHC = max(P6master.halfCycle);
halfCycles = cell(1,totHC);
chargeCycles = cell(1,ceil(totHC/2));
%make list of half cycle indices
[l hCycIndex] = findgroups(P6master.halfCycle);
if hCycIndex(1)==0
    hCycIndex = hCycIndex(2:end);
end
%build subtables from each half cycle and just charge cycles
for i = hCycIndex'
    thisCycle = find(P6master.halfCycle==i);
    startVals = P6master(thisCycle(1),:);
    halfCycles{i} = P6master(thisCycle,:);
    %halfCycles{i}.time_s = halfCycles{i}.time_s - startVals.time_s;
    if startVals.ox_red == 0
        chargeCycs{i} = table2array(P6master(thisCycle,:));
    end
end
empInd = cellfun(@isempty, chargeCycs) == 0;
chargeCycles = chargeCycs(empInd);

totC = max(P6master.cycleNumber)+1;
Cycles = cell(1,totC);
%make list of FULL cycle indices, defined by EC labs as a climb and fall of
%voltage, so opposite of what we are considering. 
[l fCycIndex] = findgroups(P6master.cycleNumber);
if fCycIndex(1)==0
    fCycIndex = fCycIndex(2:end);
end
%build subtables from each FULL cycle and rezero time
for i = fCycIndex'
    thisCycle = find(P6master.cycleNumber==i);
    startVals = P6master(thisCycle(1),:);
    Cycles{i} = P6master(thisCycle,:);
    %Cycles{i}.time_s = Cycles{i}.time_s - startVals.time_s;
end


figure(1)
plot(P6master.time_s/3600,P6master.stress_MPa);
set(gca,'FontSize',26)
xlabel('Time (Hrs)'); ylabel('Stress (MPa)');
hold on;
plot(P6master.time_s/3600,P6master.x_I__mA)
plot(P6master.time_s/3600,P6master.Ewe_V)
hold off;

