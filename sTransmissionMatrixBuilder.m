function[M] = sTransmissionMatrixBuilder(IndexFirst,IndexSecond,Wavelength,IncidentAngle,IndexFirstLayer)

    fund_consts;
    AngleFirst = asin(IndexFirstLayer*sin(IncidentAngle)/IndexFirst);
    AngleSecond = asin(IndexFirstLayer*sin(IncidentAngle)/IndexSecond);    
    
    KxFirst = 2*pi*IndexFirst./Wavelength.*cos(AngleFirst);
    KxSecond = 2*pi*IndexSecond./Wavelength.*cos(AngleSecond);
    
    L = KxSecond/KxFirst;
    
    U = [1+L;1-L];
    M = 0.5*[U flipud(U)];
end



