clear;clc;

c0 = 299792458;

lamb0 = 1555.02e-9;    %pump, calibrated w/ OSA
w10 = c0./lamb0.*2*pi;
Pin_min_unconv = 108.5e-6; %W
Pump_conv = 8.098;
Pin_min = Pin_min_unconv.*Pump_conv;
kappa1_in = 5.3233e6.*2*pi; %rad
Trans = 0.2;
kappa1_ex = kappa1_in.*(1-sqrt(Trans))./(1+sqrt(Trans));

a1sq0 = kappa1_ex.*Pin_min./((kappa1_ex+kappa1_in)./2).^2;
d_lamb = (1555.03-1554.91).*1e-9;
d_w1 = 2*pi.*c0./lamb0.^2.*d_lamb;
B1 = d_w1./a1sq0;

% delta a1sq = r.*delta voltage
%1V pz =2.1175e-3nm = 15.6264uW pump (converted)
% r = -4.603e-6.*kappa1_ex./((kappa1_ex+kappa1_in)./2).^2;
r = -1./B1.*(2*pi*278.19e6);

%delta w2 = B2r.*delta volt
% 1V pz = 5.749e-4nm lambda2 change = 278.19MHz f2 change
B2r = 2*pi.*518e6;  %assumed

%delta w1 = r0.*delta volt
%1V pz =2.1175e-3nm = 262.53MHz f1 change
r0 = 2*pi.*262.53e6;

%% fitting SH power-detuning
load('SHpower_Volt');
% SHpower_n = SHpower_n2;
[max_P, max_n] = max(SHpower_n);
xdata = (1:length(SHpower_n)).'.*incre_V;
xdata = xdata-xdata(max_n);

options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-12);
fun = @(x, xdata)(1+r./a1sq0.*xdata).^2./((2.*r0.*xdata-B2r.*xdata+x(2)).^2./(x(1).^2)+1).*1e9;
x0 = [2*pi.*5e6,0];
ydata = SHpower_n.*1e9;
x = lsqcurvefit(fun,x0,xdata,ydata);
figure; plot(xdata, ydata,'o');
hold on;
xfit = linspace(xdata(1),xdata(length(xdata)),1000);
yfit = (1+r./a1sq0.*xfit).^2./((2.*r0.*xfit-B2r.*xfit+x(2)).^2./(x(1).^2)+1).*1e9;
plot(xfit, yfit);

%% P2max
P2max = 384.31867e-12;  %W
ke2gsq_k2sq = P2max./(kappa1_ex.^2.*Pin_min.^2./((kappa1_ex+kappa1_in)./2).^4);
k2_2_assume = x(1);     %(kp2_ex+kp2_in)/2
ke2gsq_assume = ke2gsq_k2sq.*k2_2_assume.^2;
ke2_assume = 2*pi.*0.7934e6;
gsq_assume = ke2gsq_assume./ke2_assume;

%% Pin, critical
% chasing = (2w10-w20)/(-B2+2B1)
chasing = Pin_min.*kappa1_ex./((kappa1_in+kappa1_ex)./2).^2;
B2_assume = -B2r./r;
D_w = chasing.*(-B2_assume+B1.*2);

%% power line 
% g*a2a1* << B1|a1|^2*a1, ignored in actual calculation
% options = optimoptions('fsolve','MaxFunEval', 400,'MaxIter', 1000, 'algorithm', 'levenberg-marquardt');
% fun1 = @(win, Pin, xa1sq)xa1sq.*((win-w10+B1.*xa1sq).^2+((kappa1_in+kappa1_ex)./2).^2)-kappa1_ex.*Pin;
w20 = 2.*w10-D_w;
Pini = Pump_conv.*150e-6;
lambin_pool = (lamb0+linspace(-0.1, 0.25, 2000).*1e-9).';
win_pool = w20+[-2:0.0001:2].'.*1e10;%c0./lambin_pool.*2*pi;
% x0 = kappa1_ex.*Pini./((kappa1_ex+kappa1_in)./2).^2;
% a1sq_solu = zeros(length(win_pool),1);
A = B1.^2;
B = 2.*B1.*(win_pool-w10);
C = ((win_pool-w10).^2+((kappa1_ex+kappa1_in)./2).^2);
D = -kappa1_ex.*Pini;
Delt0 = B.^2-3.*A.*C;
Delt1 = 2.*B.^3-9.*A.*B.*C+27.*A.^2.*D;
CC1 = ((Delt1+sqrt(Delt1.^2-4.*Delt0.^3))./2).^(1/3);
CC2 = CC1.*(-1+sqrt(3).*1i)./2;
CC3 = CC1.*(-1-sqrt(3).*1i)./2;
x1 = -1./(3.*A).*(B+CC1+Delt0./CC1);
x2 = -1./(3.*A).*(B+CC2+Delt0./CC2);
x3 = -1./(3.*A).*(B+CC3+Delt0./CC3);
a1sq_solu = max([real(x1),real(x2),real(x3)].').';
a1sq_solu_min = min([real(x1),real(x2),real(x3)].').';
a1sq_solu(real(x2)==real(x3)) = a1sq_solu_min(real(x2)==real(x3));


A2 = B2_assume.^2;
B2 = 2.*B2_assume.*(2.*win_pool-w20);
C2 = ((2.*win_pool-w20).^2+(k2_2_assume).^2);
D2 = -gsq_assume.*a1sq_solu.^2;
Delt20 = B2.^2-3.*A2.*C2;
Delt21 = 2.*B2.^3-9.*A2.*B2.*C2+27.*A2.^2.*D2;
CC21 = ((Delt21+sqrt(Delt21.^2-4.*Delt20.^3))./2).^(1/3);
CC22 = CC21.*(-1+sqrt(3).*1i)./2;
CC23 = CC21.*(-1-sqrt(3).*1i)./2;
x21 = -1./(3.*A).*(B+CC21+Delt0./CC21);
x22 = -1./(3.*A).*(B+CC22+Delt0./CC22);
x23 = -1./(3.*A).*(B+CC23+Delt0./CC23);


