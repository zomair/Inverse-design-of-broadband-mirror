% clc;
% clear;
% close all;
% 
% 
% % defining the structure 
% n =[1 3.5 1.6 1];
% d = [0 500e-9 300e-9 0];
% Wavelength = 4e-6;
% Angle = 0;

% calculating the first t-matrix coefficients

function[rp,rs] = TMatrix(n,d,Angle,Wavelength)

sM = [1 0
    0 1];
pM = [1 0
    0 1];
for i=2:length(n)-1
    sM = sM*sTransmissionMatrixBuilder(n(i-1),n(i),Wavelength,Angle,n(1));
    sM = sM*propagationMatrixBuilder(n(i),Wavelength,d(i),Angle,n(1));

    pM = pM*pTransmissionMatrixBuilder(n(i-1),n(i),Wavelength,Angle,n(1));
    pM = pM*propagationMatrixBuilder(n(i),Wavelength,d(i),Angle,n(1));
    
    
end

sM = sM*sTransmissionMatrixBuilder(n(end-1),n(end),Wavelength,Angle,n(1)); 


pM = pM*pTransmissionMatrixBuilder(n(end-1),n(end),Wavelength,Angle,n(1)); 


rs = sM(2,1)/sM(1,1);
rp = pM(2,1)/pM(1,1);

end