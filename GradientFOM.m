function [delfdeln] =GradientFOM(TrialIndexColumn,deln,rpOrig,rsOrig)

    fund_consts;
    
    
    load('FixedData.mat');
    
    %% loading the data in a column form
    
    TrialIndex = TrialIndexColumn(1:length(wl))+1j*TrialIndexColumn((length(wl)+1):end);
    
    n(DesignLayer,:) = TrialIndex+deln;
    
    rs = size(length(wl),length(Theta));
    rp = size(length(wl),length(Theta));
    
    
    for i =1:length(wl)
        for j=1:length(Theta)
            [rp(i,j),rs(i,j)] = TMatrix(n(:,i),d,Theta(j),wl(i));
        end
    end
    
    delrpdeln = (rp-rpOrig)/deln;
    delrsdeln = (rs-rsOrig)/deln;
    
    AngleWeightingFunc = @(x) sin(x).*cos(x);
    AngleWeight = AngleWeightingFunc(Theta);
    
    delrpdelnAngleAveraged = trapz(Theta,delrpdeln.*AngleWeight,2)./(trapz(Theta,AngleWeight,2));
    delrsdelnAngleAveraged = trapz(Theta,delrsdeln.*AngleWeight,2)./(trapz(Theta,AngleWeight,2));
    
    rsOrigAngleAverage = trapz(Theta,rsOrig.*AngleWeight,2)./(trapz(Theta,AngleWeight,2)); 
    rpOrigAngleAverage = trapz(Theta,rpOrig.*AngleWeight,2)./(trapz(Theta,AngleWeight,2)); 
    
    delRdeln = (conj(rsOrigAngleAverage).*delrsdelnAngleAveraged+rsOrigAngleAverage.*conj(delrsdelnAngleAveraged)+...
        conj(rpOrigAngleAverage).*delrpdelnAngleAveraged+rpOrigAngleAverage.*conj(delrpdelnAngleAveraged));
    
    
    Rp = abs(rpOrig).^2;
    Rs = abs(rsOrig).^2;
    
    
   
    
    
    AngleAveragedRef = 0.5*(trapz(Theta,Rp.*AngleWeight,2)+trapz(Theta,Rs.*AngleWeight,2))./(trapz(Theta,AngleWeight,2));
    
    delfdeln = -2*(R-AngleAveragedRef).*delRdeln;
 
end

