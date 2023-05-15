% Analisi e confronto di EKF e PF
% - Marco Borraccino
% Università di Pisa, Identificazione Sistemi incerti 2020

%% grafico Deviazione standard xy

dev_std_xy_EKF = sqrt(out.EKF_devx.Data.^2+out.EKF_devy.Data.^2);
dev_std_xy_PF = sqrt(out.PF_devx.Data.^2+out.PF_devy.Data.^2);

figure
hold on
title('Standard Deviation xy')
axis on
grid minor
xlabel('Time [s]')
xlim([-10 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('dev_{std} [m]')
ylim([0 2])

%disegno la deviazione standard xy relativa alla stima del filtro EKF
plot(dev_std_xy_EKF,'r-','MarkerSize',2);
%disegno la deviazione standard xy relativa alla stima del filtro EKF
plot(dev_std_xy_PF,'b-','MarkerSize',2);

legend('Standard Deviation XY EKF','Standard Deviation XY PF')
%% grafico Deviazione standard z

figure
hold on
title('Standard Deviation z')
axis on
grid minor
xlabel('Time [s]')
xlim([-10 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('dev_{std} [m]')
ylim([0 0.3])

%disegno la deviazione standard z relativa alla stima del filtro EKF
plot(out.EKF_devz.Data,'r-','MarkerSize',2);
%disegno la deviazione standard z relativa alla stima del filtro PF
plot(out.PF_devz.Data,'b-','MarkerSize',2);

legend('Standard Deviation z EKF','Standard Deviation z PF')