% Analisi e confronto di EKF e PF
% - Marco Borraccino
% Università di Pisa, Identificazione Sistemi incerti 2020

%%
n=nav_N_part;
% definire l'intervallo temporale per disegnare le particelle
t_1_init=1;       %[secondi]
t_2_final=80;      %[secondi]
time_step=2;        %[secondi]
step_particle=10;       %particelle disegnate= n/step_particle

%% grafico xNorth con particelle
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
% %disegno la traiettoria della posizione stimata NED dal Navigation System
EKF_tr_xNorth = plot(out.EKF_Lat_es.data,'r-','MarkerSize',2);
% %disegno la traiettoria della posizione stimata NED dal Navigation System
PF_tr_xNorth = plot(out.PF_Lat_es.data,'b-','MarkerSize',2);

%disegno le particelle
for i=t_1_init:time_step:t_2_final
    plot(i,out.PF_s.data(1,1:step_particle:n,i),'b.','MarkerSize',5);
end
    
legend('xNorth Real','xNorth EKF Estimated','xNorth PF Estimated')

%% grafico yEast con particelle
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
%disegno le particelle
for i=t_1_init:time_step:t_2_final
    plot(i,out.PF_s.data(2,1:step_particle:n,i),'b.','MarkerSize',5);
end

legend('yEast Real','yEast EKF Estimated','yEast PF Estimated')

%% grafico Depth con particelle
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

%disegno le particelle
for i=t_1_init:time_step:t_2_final
    plot(i,out.PF_s.data(3,1:step_particle:n,i),'b.','MarkerSize',5);
end

legend('Depth Real','Depth EKF Estimated','Depth PF Estimated')

%% Plot 2d della posizione stimata con particelle

step_particle=10;       %particelle disegnate= n/step_particle
time_step=1;

figure
hold on
title('PF Trajectory 2D with particle')

%disegno la traiettoria della posizione reale NED dal Vehicle Model
real_tr = plot(out.Lon_ts.data,out.Lat_ts.data,'black-');

%disegno le particelle
for i=t_1_init:time_step:t_2_final
   Particle_PF = plot(out.PF_s.data(2,1:step_particle:n,i),out.PF_s.data(1,1:step_particle:n,i),'b.','MarkerSize',5);
end

%disegno la traiettoria della posizione stimata NED dal Navigation System
PF_tr = plot(out.PF_Lon_es.data,out.PF_Lat_es.data,'r*','MarkerSize',5);

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

legend([real_tr PF_tr Particle_PF],'Real trajectory','Estimated trajectory','PF Particles')
