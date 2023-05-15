%% PARAMETRI CONTROLLORE
ctrl_clockwise = false;     
if strcmp(inspectionDirection, 'clockwise') == 1
    ctrl_clockwise = true;
else
    ctrl_clockwise = false;
end

%Dati generali
ctrl_Ts = 0.1;          % tempo di campionamento
ctrl_r = 0.2851;        % raggio veicolo
ctrl_dc = 0.1;          % distanza CM e CB
ctrl_T_max = 45;        % spinta massima di un thruster
ctrl_M_max = ctrl_T_max*ctrl_r;     % coppia massima di un thruster, saranno usati a coppie
ctrl_TAM = [  1     1           0           0               0       0;...
              0     0           1           1               0       0;...
              0     0           0           0               1       1;...
              0     0  -(ctrl_r+ctrl_dc)  ctrl_r+ctrl_dc    0       0;...
              0     0           0           0            -ctrl_r    ctrl_r;...
          -ctrl_r ctrl_r        0           0               0       0];


%Controllore orientazione
% equivale a PI con Kp = 10 e Ki = 0.1
ctrl_Cva_td3 = zpk([0.999], [1], 10, ctrl_Ts);      % tuning con sisotool

%Parametri simulink PID
ctrl_k_pos = 1/3;       % guadagno proporzionale pos->vel
ctrl_k_or = 1;          % guadagno proporzionale or->va
ctrl_Ki_pos = 0.05;     % guadagno integrale pos->vel
ctrl_Ki_or = 0.01;      % guadagno integrale or->va
ctrl_int_bound_pos = 1;     % [m] soglia attivazione integratore posizione, 0 disattiva integrale
ctrl_int_bound_or = 0.1;    % [rad] soglia attivazione integratore orientazione, 0 disattiva integrale
ctrl_Kp_v_pos = 200;    % azione proporzionale del PI di velocità nel controllore di posizione
ctrl_Ki_v_pos = 1;      % guadagno azione integrale del PI di velocità nel controllore di posizione
ctrl_Kp_vel = 200;      % azione proporzionale del PI del controllore di velocità
ctrl_Ki_vel = 100;      % azione integrale del PI del controllore di velocità
ctrl_bound_sat_or = 1;          % [rad] soglia della massima integrazione sull'errore angolare (max vel integratore 0.05 rad/s)
ctrl_bound_sat_vel = 1;         % [m/s] soglia della massima integrazione sulla velocità (max azione integrale 100N)
ctrl_bound_sat_v_pos = 10;      % [m/s] soglia della massima integrazione sulla velocità (max azione integrale 10N)
ctrl_bound_sat_pos = 10;        % [m] soglia della massima integrazione sulla velocità (max vel integratore 0.5 m/s)
ctrl_vel_uncontrolled = 0.5;    % [m/s] velocità impostata al di fuori dello spazio di frenata
ctrl_va_uncontrolled = 0.15;    % [rad] velocità angolare impostata al di fuori dello spazio di frenata