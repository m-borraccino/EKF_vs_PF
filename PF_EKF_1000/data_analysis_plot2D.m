% Analisi e confronto di EKF e PF
% - Marco Borraccino
% Università di Pisa, Identificazione Sistemi incerti 2020


%% plot 2D reale-EKF
figure(7)
hold on
title('EKF Trajectory 2D')

%disegno la traiettoria della posizione reale NED dal Vehicle Model
real_tr = plot(out.Lon_ts.data,out.Lat_ts.data,'black-');
%disegno la traiettoria della posizione stimata NED dal Navigation System
EKF_tr = plot(out.EKF_Lon_es.data,out.EKF_Lat_es.data,'r.');

%creazione del bacino
plot(nav_A(2),nav_A(1),'r*')
text(nav_A(2),nav_A(1),' A')
plot(nav_B(2),nav_B(1),'r*')
text(nav_B(2),nav_B(1),' B')
plot(nav_C(2),nav_C(1),'r*')
text(nav_C(2),nav_C(1),' C')
plot(nav_D(2),nav_D(1),'r*')
text(nav_D(2),nav_D(1),' D')
plot([nav_A(2) nav_B(2)],[nav_A(1) nav_B(1)],'g')
plot([nav_B(2) nav_C(2)],[nav_B(1) nav_C(1)],'g')
plot([nav_C(2) nav_D(2)],[nav_C(1) nav_D(1)],'g')
plot([nav_D(2) nav_A(2)],[nav_D(1) nav_A(1)],'g')
xlabel('yEast [m]')
ylabel('xNorth [m]')

for i=1:1:length(out.EKF_Lat_es.data)-1
    xell = 3*out.EKF_devx.data(i);
    yell = 3*out.EKF_devy.data(i);
    [x,y] = getEllipse(yell,xell,[out.EKF_Lon_es.data(i) out.EKF_Lat_es.data(i)]);
    plot(x,y);
end

legend([real_tr EKF_tr],'Real trajectory','Estimated trajectory')

%% plot 2D reale-PF

figure(8)
hold on
title('PF Trajectory 2D')

%disegno la traiettoria della posizione reale NED dal Vehicle Model
real_tr = plot(out.Lon_ts.data,out.Lat_ts.data,'black-');
%disegno la traiettoria della posizione stimata NED dal Navigation System
PF_tr = plot(out.PF_Lon_es.data,out.PF_Lat_es.data,'b.');

%creazione del bacino
plot(nav_A(2),nav_A(1),'r*')
text(nav_A(2),nav_A(1),' A')
plot(nav_B(2),nav_B(1),'r*')
text(nav_B(2),nav_B(1),' B')
plot(nav_C(2),nav_C(1),'r*')
text(nav_C(2),nav_C(1),' C')
plot(nav_D(2),nav_D(1),'r*')
text(nav_D(2),nav_D(1),' D')
plot([nav_A(2) nav_B(2)],[nav_A(1) nav_B(1)],'g')
plot([nav_B(2) nav_C(2)],[nav_B(1) nav_C(1)],'g')
plot([nav_C(2) nav_D(2)],[nav_C(1) nav_D(1)],'g')
plot([nav_D(2) nav_A(2)],[nav_D(1) nav_A(1)],'g')
xlabel('yEast [m]')
ylabel('xNorth [m]')

for i=1:1:length(out.PF_Lat_es.data)-1
    xell = 3*out.PF_devx.data(i);
    yell = 3*out.PF_devy.data(i);
    [x,y] = getEllipse(yell,xell,[out.PF_Lon_es.data(i) out.PF_Lat_es.data(i)]);
    plot(x,y);
end

legend([real_tr PF_tr],'Real trajectory','Estimated trajectory')

%% plot 2D reale-EKF-PF

figure(9)
hold on
title('PF-EKF Trajectory 2D')

%disegno la traiettoria della posizione reale NED dal Vehicle Model
real_tr = plot(out.Lon_ts.data,out.Lat_ts.data,'black-');
%disegno la traiettoria della posizione stimata NED dal Navigation System
EKF_tr = plot(out.EKF_Lon_es.data,out.EKF_Lat_es.data,'r.');
%disegno la traiettoria della posizione stimata NED dal Navigation System
PF_tr = plot(out.PF_Lon_es.data,out.PF_Lat_es.data,'b.');

%creazione del bacino
plot(nav_A(2),nav_A(1),'r*')
text(nav_A(2),nav_A(1),' A')
plot(nav_B(2),nav_B(1),'r*')
text(nav_B(2),nav_B(1),' B')
plot(nav_C(2),nav_C(1),'r*')
text(nav_C(2),nav_C(1),' C')
plot(nav_D(2),nav_D(1),'r*')
text(nav_D(2),nav_D(1),' D')
plot([nav_A(2) nav_B(2)],[nav_A(1) nav_B(1)],'g')
plot([nav_B(2) nav_C(2)],[nav_B(1) nav_C(1)],'g')
plot([nav_C(2) nav_D(2)],[nav_C(1) nav_D(1)],'g')
plot([nav_D(2) nav_A(2)],[nav_D(1) nav_A(1)],'g')
xlabel('yEast [m]')
ylabel('xNorth [m]')

legend([real_tr EKF_tr PF_tr],'Real trajectory','EKF Estimated trajectory','PF Estimated trajectory')



%% ------------------------------------------------------------------------
function [x,y] = getEllipse(r1,r2,C)
    beta = linspace(0,2*pi,100);
    x = r1*cos(beta) - r2*sin(beta);
    y = r1*cos(beta) + r2*sin(beta);
    x = x + C(1,1);
    y = y + C(1,2);
end