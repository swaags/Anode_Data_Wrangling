slice = C7master(15645:20461,[1,2,6]);
%%
figure(1)
slice.strsGauss = smoothdata(slice.stress_MPa,'gaussian',5000,'SamplePoints',slice.time_s);
slice.voltGauss = smoothdata(slice.Ewe_V,'gaussian',150,'SamplePoints',slice.time_s);
plot(slice.time_s,slice.stress_MPa,'.',slice.time_s,slice.strsGauss)
hold on
plot(slice.time_s,slice.Ewe_V,'.',slice.time_s,slice.voltGauss)
evenTs = [min(slice.time_s):10:max(slice.time_s)]';
evenStress = interp1(slice.time_s,slice.strsGauss,evenTs);
evenVolt= interp1(slice.time_s,slice.Ewe_V,evenTs);
dFirst = diff([evenStress evenVolt]);  

dFirstSmooth = smoothdata(dFirst,'gaussian',50);
figure(3)
plot(evenTs(1:length(dFirstSmooth)),dFirstSmooth(:,1),evenTs(1:length(dFirstSmooth)),dFirstSmooth(:,2))

dSecond = diff(dFirstSmooth);
figure(2)
plot(evenTs(1:length(dSecond)),dSecond(:,1),evenTs(1:length(dSecond)),dSecond(:,2))
