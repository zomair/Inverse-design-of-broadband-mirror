clc;
clear;
close all;

% plot parameters
set(0, 'DefaultAxesFontSize',16)
set(0, 'defaultfigurecolor',[1 1 1]);
MyLineWidth = 1.5;
MyGridWidth = 1.2;

% defining the structure 

d = [0 550e-6 0];

Data=load('substrate_R.txt');


wl = 1./Data(:,1)*1e-2;

R = Data(:,2);
% R((wl<9.5e-6) | (wl>10e-6))=[];
% wl((wl<9.5e-6) | (wl>10e-6))=[];
Rexp = R;
T = R;
T = zeros(size(R));
ThetaMin = 9.6;
ThetaMax = 28.5;
DelTheta =10;
Theta = (ThetaMin:DelTheta:ThetaMax)*pi/180;
d =[0;500e-6;0];

n = ones(length(d),length(wl));
DesignLayer =2; % layer# for which to calculate n

Kmin = wl/(4*pi)*(-log(0.1)/d(DesignLayer));
save('FixedData.mat','Theta','wl','R','T','d','n','Kmin','DesignLayer');

% initial value of the index


nInit = zeros(1,length(wl)*2);
nInit((length(wl)+1):end) = Kmin;%imaginary
nInit(1:length(wl)) = (1+sqrt(R))./(1-sqrt(R));% real



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
    
    % defining lagrange parameters for constraint violation
    alpha =0; % penalty for real(index)>4
    beta = 0; % penalty for real(index)<1
    gamma = 0;%1j*1e-3; % penalty for imaginary(index)>10
    delta = 0;%1j*1e2; % penalty for imaginary(index)<0
    save('LargrangCoeff','alpha','beta','gamma','delta');
    
    IterationNos = 500;
    StepSize = 1e-1;
    deln = 1e-8;
    nStart = nInit;
    for i =1:IterationNos
       
        [loss,rp,rs,tp,ts,~,~] =MeritFunc(nStart);
        f(i) = mean(loss);
        sprintf('Step %d,mean loss: %f',i,f(i))
        if ((i>1) && (f(i)>f(i-1)))
            StepSize = StepSize/2;
%         elseif ((i>1) && (f(i)<=f(i-1)))
%             StepSize = StepSize*2;
%        
        end
       
        dfdn = GradientFOM(nStart,deln,rp,rs,tp,ts);
        tempDat =zeros(1,length(wl)*2);
        tempDat(1,1:length(wl)) = real(dfdn);
        tempDat(1,(length(wl)+1):end) = imag(dfdn);
        nStart = (nStart-StepSize*tempDat);
        
        
    end
    [~,~,~,~,~,Rfinal,Tfinal] =MeritFunc(nStart);
    figure
    plot(wl/1e-6,nStart(1:length(wl)),'b','linewidth',MyLineWidth)
    hold on;
    plot(wl/1e-6,nInit(1:length(wl)),'r','linewidth',MyLineWidth)
    xlabel('Wavelength (\mu m)')
    ylabel('Real refractive index');
    datacursurmode off
    figure
    plot(wl/1e-6,nStart(1,(length(wl)+1):end),'b','linewidth',MyLineWidth)
    hold on;
    plot(wl/1e-6,nInit(1,(length(wl)+1):end),'r','linewidth',MyLineWidth)
    xlabel('Wavelength (\mu m)')
    ylabel('Imaginary refractive index');
    datacursurmode off
    
    figure
    plot(wl/1e-6,Rfinal*100,'b','linewidth',MyLineWidth)
    hold on;
    plot(wl/1e-6,R*100,'r','linewidth',MyLineWidth)
    legend('Theory','Measured')
    xlabel('Wavelength (\mu m)')
    ylabel('Reflectivity (%)');
    datacursurmode off
    
    
    figure
    plot(wl/1e-6,Tfinal*100,'b','linewidth',MyLineWidth)
    hold on;
    plot(wl/1e-6,T*100,'r','linewidth',MyLineWidth)
    legend('Theory','Measured')
    xlabel('Wavelength (\mu m)')
    ylabel('Transmissivity (%)');
    datacursurmode off
    figure
    plot(1:1:IterationNos,f,'b','linewidth',MyLineWidth)
    xlabel('Iteration Nos')
    ylabel('Mean loss');
    datacursurmode off
    n(DesignLayer,:) = nStart(1:length(wl))+1j*nStart(1,(length(wl)+1):end);
end
save('Testv1.mat','f','Theta','wl','R','d','n','DesignLayer');





