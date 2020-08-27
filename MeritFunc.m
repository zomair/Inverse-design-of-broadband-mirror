function [Loss,rp,rs,AngleAveragedRef] =MeritFunc(TrialIndexColumn)

    fund_consts;
    
    
    load('FixedData.mat');
    
    %% loading the data in a column form
    
    TrialIndex = TrialIndexColumn(1:length(wl))+1j*TrialIndexColumn((length(wl)+1):end);
    
    n(DesignLayer,:) = TrialIndex;
    
    rs = size(length(wl),length(Theta));
    rp = size(length(wl),length(Theta));
    
    
    for i =1:length(wl)
        for j=1:length(Theta)
            [rp(i,j),rs(i,j)] = TMatrix(n(:,i),d,Theta(j),wl(i));
        end
    end
    Rp = abs(rp).^2;
    Rs = abs(rs).^2;
    
    AngleWeightingFunc = @(x) sin(x).*cos(x);
   
    AngleWeight = AngleWeightingFunc(Theta);
   
    
    
    AngleAveragedRef = 0.5*trapz(Theta,(Rs+Rp).*AngleWeight,2)./trapz(Theta,AngleWeight,2);
    
    load('LargrangCoeff.mat','alpha','beta','gamma','delta');
    
    Loss = ((R-AngleAveragedRef).').^2+alpha*(real(n(DesignLayer,:))>4)-beta*(real(n(DesignLayer,:))<=1)...
        +gamma*(imag(n(DesignLayer,:))>10)-delta*(imag(n(DesignLayer,:))<0);
 
end

