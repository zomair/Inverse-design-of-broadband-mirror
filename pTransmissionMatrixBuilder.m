function[M] = pTransmissionMatrixBuilder(IndexFirst,IndexSecond,Wavelength,IncidentAngle,IndexFirstLayer)

    fund_consts;
    AngleFirst = asin(IndexFirstLayer*sin(IncidentAngle)/IndexFirst);
    AngleSecond = asin(IndexFirstLayer*sin(IncidentAngle)/IndexSecond);    
    
    KxFirst = 2*pi*IndexFirst./Wavelength.*cos(AngleFirst);
    KxSecond = 2*pi*IndexSecond./Wavelength.*cos(AngleSecond);
    
    L = KxSecond/KxFirst;
    
    U = [L*IndexFirst/IndexSecond+IndexSecond/IndexFirst;L*IndexFirst/IndexSecond-IndexSecond/IndexFirst];
    M = 0.5*[U flipud(U)];
end


