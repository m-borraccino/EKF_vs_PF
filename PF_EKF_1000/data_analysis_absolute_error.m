% Analisi e confronto di EKF e PF
% - Marco Borraccino
% Università di Pisa, Identificazione Sistemi incerti 2020

%% grafico xNorth
figure
hold on
title('Position component [1]')
axis on
grid on
xlabel('Time [s]')
xlim([0 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('xNorth [m]')
ylim([min(out.Lat_ts.data)-5 max(out.Lat_ts.data)+5])

%disegno la traiettoria della posizione reale NED dal Vehicle Model
real_tr_xNorth = plot(out.Lat_ts.data,'black-','MarkerSize',2);
%disegno la traiettoria della posizione stimata NED dal Navigation System
EKF_tr_xNorth = plot(out.EKF_Lat_es.data,'r-','MarkerSize',2);
%disegno la traiettoria della posizione stimata NED dal Navigation System
PF_tr_xNorth = plot(out.PF_Lat_es.data,'b-','MarkerSize',2);

legend('xNorth Real','xNorth EKF Estimated','xNorth PF Estimated')

figure
hold on
title('Error Position component [1]')
axis on
grid on
xlabel('Time [s]')
xlim([0 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('Absolute Error [m]')
error_Lat_EKF = out.Lat_ts.data-out.EKF_Lat_es.data;
error_Lat_PF = out.Lat_ts.data-out.PF_Lat_es.data;
ylim([min(error_Lat_EKF)-0.5 max(error_Lat_EKF)+0.5])
plot(error_Lat_EKF,'r-','MarkerSize',2);
plot(error_Lat_PF,'b-','MarkerSize',2);
legend('xNorth EKF Error','xNorth PF Error')

%% grafico yEast
figure
hold on
title('Position component [2]')
axis on
grid on
xlabel('Time [s]')
xlim([0 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('yEast [m]')
ylim([min(out.Lon_ts.data)-5 max(out.Lon_ts.data)+5])

%disegno la traiettoria della posizione reale NED dal Vehicle Model
real_tr_yEast = plot(out.Lon_ts.data,'black-','MarkerSize',2);
%disegno la traiettoria della posizione stimata NED dal Navigation System
EKF_tr_yEast = plot(out.EKF_Lon_es.data,'r-','MarkerSize',2);
%disegno la traiettoria della posizione stimata NED dal Navigation System
PF_tr_yEast = plot(out.PF_Lon_es.data,'b-','MarkerSize',2);

legend('yEast Real','yEast EKF Estimated','yEast PF Estimated')

figure
hold on
title('Error Position component [2]')
axis on
grid on
xlabel('Time [s]')
xlim([0 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('Absolute Error [m]')
error_Lon_EKF = out.Lon_ts.data-out.EKF_Lon_es.data;
error_Lon_PF = out.Lon_ts.data-out.PF_Lon_es.data;
ylim([min(error_Lon_EKF)-0.5 max(error_Lon_EKF)+0.5])
plot(error_Lon_EKF,'r-','MarkerSize',2);
plot(error_Lon_PF,'b-','MarkerSize',2);
legend('yEast EKF Error','yEast PF Error')

%% grafico Depth
figure 
hold on
title('Position component [3]')
axis on
grid on
xlabel('Time [s]')
xlim([0 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('Depth [m]')
ylim([min(out.Depth_ts.data)-1 max(out.Depth_ts.data)+1])

%disegno la traiettoria della posizione reale NED dal Vehicle Model
real_tr_Depth = plot(out.Depth_ts.data,'black-','MarkerSize',2);
%disegno la traiettoria della posizione stimata NED dal Navigation System
EKF_tr_Depth = plot(out.EKF_Depth_es.data,'r-','MarkerSize',2);
%disegno la traiettoria della posizione stimata NED dal Navigation System
PF_tr_Depth = plot(out.PF_Depth_es.data,'b-','MarkerSize',2);

legend('Depth Real','Depth EKF Estimated','Depth PF Estimated')

figure
hold on
title('Error Position component [3]')
axis on
grid on
xlabel('Time [s]')
xlim([0 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('Absolute Error [m]')
error_Depth_EKF = out.Depth_ts.data-out.EKF_Depth_es.data;
error_Depth_PF = out.Depth_ts.data-out.PF_Depth_es.data;
ylim([min(error_Depth_EKF)-0.5 max(error_Depth_EKF)+0.5])
plot(error_Depth_EKF,'r-','MarkerSize',2);
plot(error_Depth_PF,'b-','MarkerSize',2);
legend('Depth EKF Error','Depth PF Error')