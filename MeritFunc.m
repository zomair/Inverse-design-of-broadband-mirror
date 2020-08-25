function [f,rp,rs] =MeritFunc(TrialIndexColumn)

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
   
    
    
    AngleAveragedRef = 0.5*(trapz(Theta,Rp.*AngleWeight,2)+trapz(Theta,Rs.*AngleWeight,2))./(trapz(Theta,AngleWeight,2));
    
    f = mean((R-AngleAveragedRef).^2);
 
end

