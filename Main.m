clc;
clear;
close all;


% defining the structure 

d = [0 550e-6 0];

Data=load('substrate_R.txt');

wl = 1./Data(:,1)*1e-2;
R = Data(:,2);
Rexp = R;
ThetaMin = 9.8;
ThetaMax = 23.6;
DelTheta = 10;
Theta = (ThetaMin:DelTheta:ThetaMax)*pi/180;
d =[0;500e-6;0];

n = ones(length(d),length(wl));
DesignLayer =2; % layer# for which to calculate n

save('FixedData.mat','Theta','wl','R','d','n','DesignLayer');

% initial value of the index
nInit = zeros(1,length(wl)*2);
nInit(1:length(wl)) = 3.5;%linspace(2.2,3.3,length(wl));% real
nInit((length(wl)+1):end) = 0;%wl/(4*pi)*(-log(0.1)/d(DesignLayer));%imaginary

% lower bound on the index
nlb = zeros(1,length(wl)*2);
nlb(1:length(wl)) = 1.0;%real
nlb((length(wl)+1):end) = 0;%imaginary

% upper bound on the index
nub = zeros(1,length(wl)*2);
nub(1:length(wl)) = 5.0;%real 
nub((length(wl)+1):end) = 10;%imaginary

SolverType = 'manual'; % could be 'MATLAB' or 'manual'
if strcmp(SolverType,'MATLAB')
    options = optimoptions(@fmincon,'Display','iter','MaxFunctionEvaluations',10000,'ScaleProblem',true,'UseParallel',true,'Algorithm','sqp');
    parpool('local')
    f =fmincon(@FOM,nInit,[],[],[],[],nlb,nub,[],options);
    delete(gcp)
else
    IterationNos = 100;
    alpha = 1e-4;
    deln = 0.01+1j*1e-4;
    nStart = nInit;
    for i =1:IterationNos
        i
        [f(i),rp,rs] =MeritFunc(nStart);
        dfdn = GradientFOM(nStart,deln,rp,rs);
        tempDat =zeros(1,length(wl)*2);
        tempDat(1,1:length(wl)) = real(dfdn);
        tempDat(1,(length(wl)+1):end) = imag(dfdn);
        nStart = nStart-alpha*tempDat;
    end
    
    plot(wl/1e-6,nStart(1:length(wl)),'b')
    figure
    plot(wl/1e-6,nStart(1,(length(wl)+1):end),'r')
end
%save('Testv1.mat','f','Theta','wl','R','d','n','DesignLayer');





