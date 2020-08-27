clc;
clear;
close all;

fund_consts;
DataLoad= 'Off';
if strcmp(DataLoad,'Off')
    wl = (1:0.1:10)*1e-6;
    Theta = (0:1:3)*pi/180;

    d = [0;500e-6;0];
    n = ones(length(d),length(wl));
    n(2,:) = linspace(2.2,3.3,length(wl));
else
    load('Testv1.mat');
end
for i=1:length(wl)
    for j =1:length(Theta)
        [rp(i,j),rs(i,j)] =TMatrix(n(:,i),d,Theta(j),wl(i));
        
    end
end

AngleWeightingFunc = @(x) sin(x).*cos(x);
   
AngleWeight = AngleWeightingFunc(Theta);
Rp = abs(rp).^2;
Rs = abs(rs).^2;   
    
    
AngleAveragedRef = 0.5*(trapz(Theta,Rp.*AngleWeight,2)+trapz(Theta,Rs.*AngleWeight,2))./(trapz(Theta,AngleWeight,2));



wn = 1./(wl*1e2);
Data= zeros(length(wl),2);
Data(:,1)=wn;
Data(:,2)=  AngleAveragedRef;

save('TestGenData.mat','Data');