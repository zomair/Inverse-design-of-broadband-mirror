function[P] = propagationMatrixBuilder(Index,Wavelength,Thickness,IncidentAngle,IndexFirstLayer)

    fund_consts;
    Angle = asin(IndexFirstLayer*sin(IncidentAngle)/Index);
        
    
    Kx = 2*pi*Index./Wavelength.*cos(Angle);
    P = [exp(-1j*Kx*Thickness) 0
        0 exp(1j*Kx*Thickness)];
end


