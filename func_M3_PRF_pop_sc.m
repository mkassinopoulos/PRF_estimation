function [obj_function,CRF_pop_1, CRF_pop_2, RRF_pop_1, RRF_pop_2 ,r_PRF_pop ]=func_M3_PRF_pop_sc(P,Ts_10,HR_all,RF_all,ind_BOLD_10,GS_all)


t1c=P(1) ; d1c=P(2);
t2c=P(3);  d2c=P(4);
t1r=P(5); d1r = P(6);
t2r=P(7); d2r=P(8);

[NV,nScans] = size(GS_all);
t_win= 0 :Ts_10:60;

a1= sqrt(t1c)/d1c; a2= sqrt(t1c)*d1c ;
IR = t_win.^a1.*exp(-t_win/a2); CRF_pop_1=IR/max(IR);  
a1= sqrt(t2c)/d2c; a2= sqrt(t2c)*d2c ;
IR = t_win.^a1.*exp(-t_win/a2); CRF_pop_2=IR/max(IR);     

a1= sqrt(t1r)/d1r; a2= sqrt(t1r)*d1r ;
IR = t_win.^a1.*exp(-t_win/a2); RRF_pop_1=IR/max(IR);   
a1= sqrt(t2r)/d2r; a2= sqrt(t2r)*d2r ;
IR = t_win.^a1.*exp(-t_win/a2); RRF_pop_2=IR/max(IR);  

r_PRF_pop = zeros(nScans,1);
parfor sc = 1:nScans
    GS = zscore(GS_all(:,sc)); HR = zscore(HR_all(:,sc)); RF = zscore(RF_all(:,sc));    
    
    HR_conv_1 = conv(HR,CRF_pop_1); HR_conv_1_MR = HR_conv_1(ind_BOLD_10);
    HR_conv_2 = conv(HR,CRF_pop_2); HR_conv_2_MR = HR_conv_2(ind_BOLD_10);
    RF_conv_1 = conv(RF,RRF_pop_1); RF_conv_1_MR = RF_conv_1(ind_BOLD_10); 
    RF_conv_2 = conv(RF,RRF_pop_2); RF_conv_2_MR = RF_conv_2(ind_BOLD_10); 
    
    regr = [HR_conv_1_MR, HR_conv_2_MR, RF_conv_1_MR, RF_conv_2_MR];     regr = detrend(regr,'linear');
    regr = [ones(NV,1), regr];
    
    B = regr\GS;     yPred = regr*B;
    r_PRF_pop(sc) = corr(GS, yPred);    
    
end
obj_function =1 - mean(r_PRF_pop);




%%   ----------------------------------------------------