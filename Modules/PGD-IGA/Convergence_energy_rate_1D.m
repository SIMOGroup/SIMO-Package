clear all;
format long;
L=1;

h=[L/20 L/40 L/60 L/80 L/100 L/160]; % h
   
dispErrorNormp1=[0.0023 5.4057e-4 2.3624e-4 1.3177e-4 8.3912e-5 3.2532e-5];
dispErrorNormp2=[9.5086e-5 7.9024e-6 1.8663e-6 6.7347e-7 3.0659e-7 6.0977e-8];
dispErrorNormp3=[2.3545e-4 2.0521e-5 4.8952e-6 1.7724e-6 8.0664e-7 1.5492e-7];

%----------
 for i=1:length(h)
   logh(i)=log10(h(i));
   logDisp_p1(i)=log10(dispErrorNormp1(i));
   logDisp_p2(i)=log10(dispErrorNormp2(i));
   logDisp_p3(i)=log10(dispErrorNormp3(i));   
 end
fitDisp_p1=polyfit(logh,logDisp_p1,1);
nDp1=fitDisp_p1(1)
fitDisp_p2=polyfit(logh,logDisp_p2,1);
nDp2=fitDisp_p2(1)
fitDisp_p3=polyfit(logh,logDisp_p3,1);
nDp3=fitDisp_p3(1)

figure
plot(logh,logDisp_p1,'bo-','LineWidth',1.5);
hold on;
plot(logh,logDisp_p2,'bd-','LineWidth',1.5);
plot(logh,logDisp_p3,'bs-','LineWidth',1.5);

legend('p=1','p=2','p=3');
