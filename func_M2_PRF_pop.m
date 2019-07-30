function [obj_function,CRF_pop,RRF_pop,r_PRF_pop  ]=func_M2_PRF_pop(P,Ts_10,HR_all,RF_all,ind_BOLD_10,GS_all)


t1c=P(1) ; d1c=P(2);
t2c=P(3);  d2c=P(4);
t1r=P(5); d1r = P(6);
t2r=P(7); d2r=P(8);
R_CRF = P(9); R_RRF = P(10);

[NV,nScans] = size(GS_all);
t_win= 0 :Ts_10:60;

a1= sqrt(t1c)/d1c; a2= sqrt(t1c)*d1c ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardF=IR/max(IR);  
a1= sqrt(t2c)/d2c; a2= sqrt(t2c)*d2c ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardS=IR/max(IR);     

a1= sqrt(t1r)/d1r; a2= sqrt(t1r)*d1r ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respF=IR/max(IR);   
a1= sqrt(t2r)/d2r; a2= sqrt(t2r)*d2r ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respS=IR/max(IR);  

CRF_pop = IR_cardF + R_CRF*IR_cardS;     
RRF_pop = IR_respF + R_RRF*IR_respS;     

r_PRF_pop = zeros(nScans,1);

parfor sc = 1:nScans
    GS = zscore(GS_all(:,sc)); HR = zscore(HR_all(:,sc)); RF = zscore(RF_all(:,sc));    
    
    HR_conv = conv(HR,CRF_pop); HR_conv_MR = HR_conv(ind_BOLD_10);
    RF_conv = conv(RF,RRF_pop); RF_conv_MR = RF_conv(ind_BOLD_10); 
    regr = [HR_conv_MR,RF_conv_MR];     regr = detrend(regr,'linear');
    regr = [ones(NV,1), regr];
    
    B = regr\GS;     yPred = regr*B;
    r_PRF_pop(sc) = corr(GS, yPred);     
    
end
obj_function =1 - mean(r_PRF_pop);



%%   ----------------------------------------------------