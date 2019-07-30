function [obj_function,CRF_sc,RRF_sc,HR_conv,RF_conv,r_PRF_sc,yPred, HR_Fconv_MR, HR_Sconv_MR, RF_Fconv_MR, RF_Sconv_MR]=func_M5_PRF_sc_vxl(P,Ts_10,HR,RF,ind_BOLD_10,GS,only_OF)

t1c=P(1) ; d1c=P(2);
t2c=P(3);  d2c=P(4);
t1r=P(5); d1r = P(6);
t2r=P(7); d2r=P(8);

NV = length(GS);
t_win= 0 :Ts_10:60;

a1= sqrt(t1c)/d1c; a2= sqrt(t1c)*d1c ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardF=IR/max(IR);  
a1= sqrt(t2c)/d2c; a2= sqrt(t2c)*d2c ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardS=IR/max(IR);     

a1= sqrt(t1r)/d1r; a2= sqrt(t1r)*d1r ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respF=IR/max(IR);   
a1= sqrt(t2r)/d2r; a2= sqrt(t2r)*d2r ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respS=IR/max(IR);  

r_PRF_sc=zeros(3,1);

HR = zscore(HR); RF = zscore(RF);
HR_Fconv=conv(HR,IR_cardF); HR_Fconv_MR=HR_Fconv(ind_BOLD_10);
HR_Sconv=conv(HR,IR_cardS); HR_Sconv_MR=HR_Sconv(ind_BOLD_10);
RF_Fconv=conv(RF,IR_respF); RF_Fconv_MR=RF_Fconv(ind_BOLD_10);
RF_Sconv=conv(RF,IR_respS); RF_Sconv_MR=RF_Sconv(ind_BOLD_10);
regr = [HR_Fconv_MR,HR_Sconv_MR,RF_Fconv_MR,RF_Sconv_MR];  regr = detrend(regr,'linear');
regr = [ones(NV,1),regr];

B = regr\GS;     yPred = regr*B;
obj_function = 1 - corr(yPred,GS) ;


if only_OF == 0
    CRF_sc = B(2) * IR_cardF + B(3) * IR_cardS;     CRF_sc = CRF_sc/max(abs(CRF_sc));
    RRF_sc = B(4)*IR_respF + B(5)*IR_respS; RRF_sc = RRF_sc/max(abs(RRF_sc));
    
    HR_Fconv=conv(HR,IR_cardF);
    HR_Sconv=conv(HR,IR_cardS);
    
    HR_conv = B(2)*HR_Fconv + B(3)*HR_Sconv;    HR_conv = HR_conv(1:length(HR));
    RF_conv = B(4)*RF_Fconv + B(5)*RF_Sconv;     RF_conv = RF_conv(1:length(RF));
    
    r_PRF_sc(1) = corr(yPred,GS);
    yPred_card = regr(:,2:3)*B(2:3);  r_PRF_sc(2) = corr(yPred_card,GS);
    yPred_resp = regr(:,4:5)*B(4:5);  r_PRF_sc(3) = corr(yPred_resp,GS);
    
    CRF_sc = B(2) * IR_cardF + B(3) * IR_cardS;     CRF_sc = CRF_sc/max(abs(CRF_sc));
    RRF_sc = B(4)*IR_respF + B(5)*IR_respS; RRF_sc = RRF_sc/max(abs(RRF_sc));
               
end



%%   ----------------------------------------------------



