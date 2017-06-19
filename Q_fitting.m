clear;clc;
load exp_data;  %choose between data and data_f2;

%%
ydata = Trans_wo;    %need to change between w, wo, f2
xdata = (1:length(ydata)).';

%% general fitting
% fun = @(x, xdata)level.*(1-(x(1)-x(2))./((xdata-x(3)).^2+x(1)));
% x0 = [20.^2, 0, 650];
% x = lsqcurvefit(fun,x0,xdata,ydata);
% x1 = x;
% ydata_fit = level.*(1-(x(1)-x(2))./((xdata-x(3)).^2+x(1)));
% figure; plot(xdata, ydata); 
% hold on;
% plot(xdata, ydata_fit);

%% with splitting
fun = @(x,xdata)level.*((2.*xdata-x(2)-x(4)).^2.*x(3).^2+((xdata-x(2)).*(xdata-x(4))-x(1).^2+2.*(x(1)-x(3)).*x(1)).^2)...
    ./((xdata-x(2)).^2+x(1).^2)./((xdata-x(4)).^2+x(1).^2);
x0 = [97, 500, 50, 900];%x = x0;

% ydata_fit = level.*((2.*xdata-x(2)-x(4)).^2.*x(3).^2+((xdata-x(2)).*(xdata-x(4))-x(1).^2+2.*(x(1)-x(3)).*x(1)).^2)...
%     ./((xdata-x(2)).^2+x(1).^2)./((xdata-x(4)).^2+x(1).^2);
figure; plot(xdata, ydata); 
hold on;
% plot(xdata, ydata_fit);
x = lsqcurvefit(fun,x0,xdata,ydata);
ydata_fit = level.*((2.*xdata-x(2)-x(4)).^2.*x(3).^2+((xdata-x(2)).*(xdata-x(4))-x(1).^2+2.*(x(1)-x(3)).*x(1)).^2)...
    ./((xdata-x(2)).^2+x(1).^2)./((xdata-x(4)).^2+x(1).^2);
plot(xdata, ydata_fit);

%% convert to normal units
%for 1555
% Full_range = 1./freq_sweep./2;
% unit_wl = increm_wo./Full_range.*sweep_Vpp.*Vpp_nm;  %need to change between w, wo
% c0 = 299792458;
% unit_omega = c0./(lambda.*1e-9).^2.*(unit_wl.*1e-9);
% kap_ex = (x(1)-x(3)).*unit_omega./1e6; %MHz
% kap_in = (x(1)+x(3)).*unit_omega./1e6;

%for 778
Full_range = 1./freq_sweep./2;
unit_omega = increm_f2./Full_range.*Vpp_GHz.*1e9;
kap_ex = (sqrt(x(1))-sqrt(x(2))).*unit_omega./1e6;  %MHz
kap_in = (sqrt(x(1))+sqrt(x(2))).*unit_omega./1e6;
