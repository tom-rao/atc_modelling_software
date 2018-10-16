% ԭʼ����
% XactΪ�����tetRŨ�ȣ�XtΪtetRŨ�ȣ�SxΪaTCŨ��
% (Xact) = (Xt)/(1+(Sx/Kx)^n);
% (A_mRNA) = (Ptet)*(beta)/(1+Xact/Kd); mRNAת¼����
% (dmRNA)/dt = (A_mRNA)-(Kdeg_mRNA)*(mRNA); mRNA�ۻ�����
% (dtetR)/dt = (Ktrans)*mRNA-(Kdeg_tetR)*(tetR); tetR�ۻ�����
% aTCΪaTCŨ��, betaΪPtet��������������, KxΪaTC��tetR�Ľ��볣��, KdΪtetR��DNA�Ľ��볣��
% KtransΪmRNAת¼����, Kdeg_mRNAΪmRNA��������, Kdeg_tetRΪtetR��������, PtetΪ�����ӿ�����
% ln(aTC(t))=ln(aTC(0))-kdeg_atc*t; aTC���⣬kdeg_atcΪatc�������ʳ���

%�����趨
f1=figure();
title('tetR(nM)-t(s)')
xlabel('t(s)')
ylabel('tetR(nM)')
hold on;
f2=figure();
title('GSDMD(nM)-t(s)')
xlabel('t(s)')
ylabel('GSDMD(nM)')
hold on;
f3=figure();
title('GSDMDMax(nM)-ATc(nM)')
xlabel('ATc(nM)')
ylabel('Max concentration of GSDMD(nM)')
hold on;
beta=0.0023;
Kx=0.36;
Kd=0.1;
Ktrans=240;
Kdeg_mRNA=0.009;
Kdeg_tetR=0.631;
Ptet=4;
atc_down=1;
atc_up=10;
Kdeg_atc=7.6e-4;
Ktrans_GSDMD=200;
Kdeg_GSDMD=0.8;

i=1;
insight=zeros(1,int16((atc_up-atc_down)/0.5+1));
aTC_Con=zeros(1,int16((atc_up-atc_down)/0.5+1));


for aTC=atc_down:0.5:atc_up
    aTC_Con(i)=aTC;
    insight(i)=parameter(aTC, beta, Kx, Kd, Ktrans, Kdeg_mRNA, Kdeg_tetR, Ptet, Ktrans_GSDMD, Kdeg_GSDMD, Kdeg_atc, f1, f2);
    i=i+1;
end
y_curve=log(1./insight);
x_curve=log(1./aTC_Con);
[p,s]=polyfit(x_curve, y_curve,1);


figure(f3)
plot(aTC_Con, insight, 'red', 'linewidth', 2);
plot(x_curve,y_curve,'color','blue','linewidth',1)
hold on;

%���볣����΢�ַ��̣������ƶ���ѧͼ��
function equali=parameter(aTC, beta, Kx, Kd, Ktrans, Kdeg_mRNA, Kdeg_tetR, Ptet,  Ktrans_GSDMD, Kdeg_GSDMD, Kdeg_atc, f1, f2)
[t,y]=ode45(@odefun,[0,300],[0,0,0],[],aTC, beta, Kx, Kd, Ktrans, Kdeg_mRNA, Kdeg_tetR, Ptet, Ktrans_GSDMD, Kdeg_GSDMD, Kdeg_atc);
figure(f1)
plot(t,y(:,2),'red')
figure(f2)
plot(t,y(:,3),'green')
equali=max(y(:,3));
end

%΢�ַ�����
function dydt=odefun(t,y,aTC, beta, Kx, Kd, Ktrans, Kdeg_mRNA, Kdeg_tetR, Ptet, Ktrans_GSDMD, Kdeg_GSDMD, Kdeg_atc)
dydt=zeros(3,1);
%y1��mRNAת¼����y2��tetRת¼��
dydt(1)=Ptet*beta/(1+(y(2)/(1+(aTC*exp(-Kdeg_atc*t).^2/Kx))).^2/Kd)-Kdeg_mRNA*y(1);
dydt(2)=Ktrans*y(1)-Kdeg_tetR*y(2);
dydt(3)=Ktrans_GSDMD*y(1)-Kdeg_GSDMD*y(3);
end
