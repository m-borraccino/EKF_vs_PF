% Analisi e confronto di EKF e PF
% - Marco Borraccino
% Università di Pisa, Identificazione Sistemi incerti 2020

clear all
close all
clc

init_all;
%scelta del filtro con il quale eseguire la simulazione
nav_filter = -1;             %  1 = EKF
                             % -1 = PF
                            
out = sim('sim_all',2000);
close_system

%% Sezione per il calcolo della norma dell'errore e confronto con indice RMSE

error_EKF = [out.Lat_ts.data'-out.EKF_Lat_es.data';...
             out.Lon_ts.data'-out.EKF_Lon_es.data';...
             out.Depth_ts.data'-out.EKF_Depth_es.data'];

error_PF = [out.Lat_ts.data'-out.PF_Lat_es.data';...
            out.Lon_ts.data'-out.PF_Lon_es.data';...
            out.Depth_ts.data'-out.PF_Depth_es.data'];

error_norm_EKF = zeros(1,length(error_EKF));
error_norm_PF = zeros(1,length(error_PF));
for k=1:length(error_EKF)
    error_norm_EKF(k) = norm(error_EKF(:,k));
end
for k=1:length(error_PF)
    error_norm_PF(k) = norm(error_PF(:,k));
end

figure
hold on
title(['Error Position norm'])
axis on
grid on
xlabel('Time [s]')
xlim([0 out.SimulationMetadata.ModelInfo.StopTime])
ylabel('Absolute Error norm [m]')
ylim([min([error_norm_PF,error_norm_EKF])-0.5 max([error_norm_PF,error_norm_EKF])+0.5])
plot(error_norm_EKF,'r-','MarkerSize',2);
plot(error_norm_PF ,'b-','MarkerSize',2);  
legend('EKF Error','PF Error')


RMSE_EKF_x = sqrt(mean((error_EKF(1,:)).^2));
RMSE_EKF_y = sqrt(mean((error_EKF(2,:)).^2));
RMSE_EKF_z = sqrt(mean((error_EKF(3,:)).^2));
RMSE_PF_x = sqrt(mean((error_PF(1,:)).^2));
RMSE_PF_y = sqrt(mean((error_PF(2,:)).^2));
RMSE_PF_z = sqrt(mean((error_PF(3,:)).^2));
RMSE_EKF = norm([RMSE_EKF_x RMSE_EKF_y RMSE_EKF_z]);
RMSE_PF  = norm([RMSE_PF_x RMSE_PF_y RMSE_PF_z]);
VarNames = {'RMSEx', 'RMSEy', 'RMSEz','RMSE total'};
RowNames = {'EKF', 'PF'};
T = table([RMSE_EKF_x; RMSE_PF_x],[RMSE_EKF_y ;RMSE_PF_y],[RMSE_EKF_z; RMSE_PF_z],[RMSE_EKF;RMSE_PF], ...
    'VariableNames',VarNames,'RowNames',RowNames)

%% eseguire uno di questi script per altre analisi e confronti

% data_analysis_plot2D

% data_analysis_particles

% data_analysis_dev_standard


% data_analysis_absolute_error

