function [Regr]=func_RETR_Resp_regressors(resp_f,M,Fs)
% RETROICOR (Respiratory regressors)

NT=length(resp_f);
Phi=zeros(NT,1);

resp_f = smooth(resp_f, 1*Fs);
resp_der=diff(resp_f); resp_der=[resp_der;resp_der(end)];

resp_der = filloutliers(resp_der,'linear','movmedian',0.5*Fs);

NB=500;
[Val,edges]=histcounts(resp_f,NB);
for i=1:NT    
    v=resp_f(i);    
    [~,edge]=min(abs(v-edges));
    if (edge-v)>0, edge=edge-1; end
    area=sum(Val(1:edge));
    sign_resp=sign(resp_der(i)); if  sign_resp==0, sign_resp=1; end
    Phi(i)=pi*area*sign_resp/NT;    
end    

Regr=zeros(NT,M*2);
for i=1:M
    Regr(:,(i-1)*2+1)=cos(i*Phi);
    Regr(:,i*2)=sin(i*Phi);
end



%%



