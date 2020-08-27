function [delfdeln] =GradientFOM(TrialIndexColumn,deln,rpOrig,rsOrig,tpOrig,tsOrig)

    fund_consts;
    
    
    load('FixedData.mat');
    
    %% loading the data in a column form
    
    TrialIndex = TrialIndexColumn(1:length(wl))+1j*TrialIndexColumn((length(wl)+1):end);
    
    n(DesignLayer,:) = TrialIndex+deln;
    
    rs = size(length(wl),length(Theta));
    rp = size(length(wl),length(Theta));
    ts = size(length(wl),length(Theta));
    tp = size(length(wl),length(Theta));
    
    
    for i =1:length(wl)
        for j=1:length(Theta)
            [rp(i,j),rs(i,j),tp(i,j),ts(i,j)] = TMatrix(n(:,i),d,Theta(j),wl(i));
        end
    end
    
    delrpdeln = (rp-rpOrig)/deln;
    delrsdeln = (rs-rsOrig)/deln;

    deltpdeln = (tp-tpOrig)/deln;
    deltsdeln = (ts-tsOrig)/deln;
    
    
    delrpdelnConj = conj(rp-rpOrig)/deln;
    delrsdelnConj = conj(rs-rsOrig)/deln;

    
    deltpdelnConj = conj(tp-tpOrig)/deln;
    deltsdelnConj = conj(ts-tsOrig)/deln;
    
    
    AngleWeightingFunc = @(x) sin(x).*cos(x);
    AngleWeight = AngleWeightingFunc(Theta);
    
    Delrp =rpOrig.*delrpdelnConj+delrpdeln.*conj(rpOrig);
    Delrs =rsOrig.*delrsdelnConj+delrsdeln.*conj(rsOrig);

    Deltp =tpOrig.*deltpdelnConj+deltpdeln.*conj(tpOrig);
    Delts =tsOrig.*deltsdelnConj+deltsdeln.*conj(tsOrig);
    
    DelrpPlusDelrs = Delrp+Delrs;
    DeltpPlusDelts = Deltp+Delts;    
    
    delRdeln = 0.5*trapz(Theta,DelrpPlusDelrs.*AngleWeight,2)./(trapz(Theta,AngleWeight,2));
    delTdeln = 0.5*trapz(Theta,DeltpPlusDelts.*AngleWeight,2)./(trapz(Theta,AngleWeight,2));    
    
    Rp = abs(rpOrig).^2;
    Rs = abs(rsOrig).^2;
    Tp = abs(tpOrig).^2;
    Ts = abs(tsOrig).^2;

    
    AngleAveragedRef =  0.5*trapz(Theta,(Rs+Rp).*AngleWeight,2)./trapz(Theta,AngleWeight,2);
    AngleAveragedTrans =  0.5*trapz(Theta,(Ts+Tp).*AngleWeight,2)./trapz(Theta,AngleWeight,2);    
    %delfdeln = -2*(R-AngleAveragedRef).*delRdeln;
    load('LargrangCoeff.mat','alpha','beta','gamma','delta');
   
    
    
    delfdeln = -2*(delRdeln.').*((R-AngleAveragedRef).')+-2*(delTdeln.').*((T-AngleAveragedTrans).')+...
        alpha*(real(n(DesignLayer,:))>4)-beta*(real(n(DesignLayer,:))<=1)...
        +gamma*(imag(n(DesignLayer,:))>10)-delta*(imag(n(DesignLayer,:))<0);
    
end

