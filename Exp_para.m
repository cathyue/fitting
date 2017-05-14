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
ke2gsq_assume = ke2gsq_k2sq.*x(1).^2;
ke2_assume = 2*pi.*0.7934e6;
gsq_assume = ke2gsq_assume./ke2_assume;

%% Pin, critical
% chasing = (2w10-w20)/(B2-2B1)
chasing = Pin_min.*kappa1_ex./((kappa1_in+kappa1_ex)./2).^2;
B2_assume = -B2r./r;
D_w = chasing.*(B2_assume-B1.*2);

%% power line 
% g*a2a1* << B1|a1|^2*a1, ignored in actual calculation
fun1 = @(win, Pin, xa1sq)xa1sq.*((win-w10-B1.*xa1sq).^2+((kappa1_in+kappa1_ex)./2).^2)-kappa1_ex.*Pin;
Pini = Pump_conv.*150e-6;
lambin_pool = (lamb0+linspace(0, 0.25, 200).*1e-9).';
win_pool = c0./lambin_pool.*2*pi;
x0 = kappa1_ex.*Pini./((kappa1_ex+kappa1_in)./2).^2;
a1sq_solu = zeros(length(win_pool),1);
for k = 1:length(win_pool)
    a1sq_solu = fsolve(fun1(win_pool(k), Pini))
end


