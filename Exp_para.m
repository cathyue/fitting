clear;clc;

Pin_min_unconv = 103.2e-6; %W
Pump_conv = 8.098;
Pin_min = Pin_min_unconv.*Pump_conv;
kappa1_in = 5.3233e6.*2*pi; %rad
Trans = 0.2;
kappa1_ex = kappa1_in.*(1-sqrt(Trans))./(1+sqrt(Trans));

a1sq0 = kappa1_ex.*Pin_min./((kappa1_ex+kappa1_in)./2).^2;

% delta a1sq = r.*delta voltage
%1V pz =2.1175e-3nm = 4.603uW pump (converted)
r = -4.603e-6.*kappa1_ex./((kappa1_ex+kappa1_in)./2).^2;

%delta w2 = B2r.*delta volt
% 1V pz = 5.749e-4nm lambda2 change = 278.19MHz f2 change
B2r = 2*pi.*521e6;

%delta w1 = r0.*delta volt
%1V pz =2.1175e-3nm = 2.1175MHz f1 change
r0 = 2*pi.*262.53e6;

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

