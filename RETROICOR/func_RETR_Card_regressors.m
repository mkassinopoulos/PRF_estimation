function [Regr]=func_RETR_Card_regressors(time,PPGlocs,M)
% RETROICOR (Cardiac regressors)

NV = length(time);

Phi = zeros(NV,1);
for i = 1:NV    
    t = time(i);   
    [~,minI] = min(abs(PPGlocs-t));
    
    minOnLeft = t-PPGlocs(minI)>0;
    if (minI == 1 && ~minOnLeft)
        t2 = PPGlocs(minI);
        t1 = t2-1;        
    elseif (minI == length(PPGlocs) && minOnLeft)
        t1 = PPGlocs(minI);
        t2 = t1+1;
    elseif minOnLeft       
        t1 = PPGlocs(minI);
        t2 = PPGlocs(minI+1);
    else
        t1 = PPGlocs(minI-1);
        t2 = PPGlocs(minI);
    end  
    
    Phi(i) = 2*pi*(t-t1)/(t2-t1);    
end
    
Regr = zeros(NV,M*2);
for i = 1:M
    Regr(:,(i-1)*2+1) = cos(i*Phi);
    Regr(:,i*2) = sin(i*Phi);
end


%%
