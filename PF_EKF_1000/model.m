%% Vehicle Model Mission Data

% Origine terna NED (cornerC), espressa in coordinate ECEF:
mod_lat0 = 43.780189;                          % [decimal degrees]
mod_lon0 = 11.282698;                          % [decimal degrees]
mod_h0 = 0;                                    % [m]

% Punto inizio missione, espresso in coordinate NED:
mod_wgs84 = wgs84Ellipsoid('meter');
[mod_xNorth,mod_yEast,mod_zDown] = geodetic2ned(initPoint(1), initPoint(2), initPoint(3), mod_lat0, mod_lon0, mod_h0, mod_wgs84);

% Posa iniziale del veicolo, espressa in terna NED:
mod_eta_init = [mod_xNorth; mod_yEast; -mod_zDown; 0; 0; 0];

% Velocità iniziale del veicolo, espressa in terna body:
mod_ni_init = [0; 0; 0; 0; 0; 0];

%% Vehicle Model Parameters

% Il veicolo è stato supposto sferico.

% Parametri fondamentali:
mod_m = 100;                                   % Massa del veicolo [kg]
mod_rho_AUV = rho;                             % Densità del veicolo [kg/m^3]
mod_Vol = mod_m/mod_rho_AUV;                   % Volume del veicolo [m^3]
mod_R = ((3*mod_Vol)/(4*pi))^(1/3);            % Raggio del veicolo [m]
mod_g = 9.81;                                  % Accelerazione di gravità [m/s^2]
mod_W = mod_m*mod_g;                           % Forza peso [N] 

% Parametri per la matrice di massa aggiunta e per la matrice di Coriolis aggiunta:
mod_Xu_dot = (2/3)*mod_rho_AUV*pi*(mod_R^3);   % Termine di massa aggiunta lungo la direzione di surge [kg]
mod_Yv_dot = mod_Xu_dot;                       % Termine di massa aggiunta lungo la direzione di sway [kg]
mod_Zw_dot = mod_Xu_dot;                       % Termine di massa aggiunta lungo la direzione di heave [kg]
mod_Kp_dot = 0;                                % Termine di inerzia aggiunta per una rotazione in rollio [kg*m^2]
mod_Mq_dot = mod_Kp_dot;                       % Termine di inerzia aggiunta per una rotazione in beccheggio [kg*m^2]
mod_Nr_dot = mod_Kp_dot;                       % Termine di inerzia aggiunta per una rotazione in imbardata [kg*m^2]

% Parametri per la matrice di damping:
mod_Rm = (4*mod_R)/(3*pi);                     % Raggio medio semisfera [m]
mod_A = pi*mod_R*mod_R;                        % Area della sezione trasversale della sfera [m^2]
mod_Cd_F = 0.45;                               % Coefficiente di drag per la forza []
mod_Cd_M = 0.4;                                % Coefficiente di drag per il momento []

%% Inertia Tensor Computation

% Si definiscono i momenti di inerzia della sfera calcolati assumendo il 
% centro di massa coincidente con il centro di simmetria (cioè con il 
% centro di galleggiamento) e una massa distribuita in modo uniforme:
mod_Mxx = (2/5)*mod_m*(mod_R^2);               % [kg*m^2]
mod_Myy = (2/5)*mod_m*(mod_R^2);               % [kg*m^2]
mod_Mzz = (2/5)*mod_m*(mod_R^2);               % [kg*m^2]

% Si definisce il tensore di inerzia sotto le ipotesi di cui sopra:
mod_icm = [mod_Mxx mod_Myy mod_Mzz];
mod_Icm = diag(mod_icm);           

% Si definisce il vettore distanza tra il centro di massa abbassato di 10 
% cm lungo l'asse z del veicolo e il centro di simmetria:
mod_r_cm = [0 0 0.10];

mod_x_cm = mod_r_cm(1);                        % [m]
mod_y_cm = mod_r_cm(2);                        % [m]
mod_z_cm = mod_r_cm(3);                        % [m]

% Si definisce la matrice Skew-simmetrica e la si applica al vettore
% distanza:
mod_S = [    0     -mod_z_cm   mod_y_cm;
          mod_z_cm     0      -mod_x_cm;
         -mod_y_cm  mod_x_cm      0     ];
 
% Si applica il teorema di Huygens-Steiner per calcolare il nuovo tensore 
% di inerzia
mod_Io = mod_Icm - mod_m*(mod_S*mod_S);

% Si ricavano le nuove componenti del tensore di inerzia:
mod_Ixx = mod_Io(1,1);                         % [kg*m^2]
mod_Iyy = mod_Io(2,2);                         % [kg*m^2]
mod_Izz = mod_Io(3,3);                         % [kg*m^2]

%% Thrusters Data

% Thruster position:
mod_p1 = [0; mod_R; 0];
mod_p2 = [0; -mod_R; 0];
mod_p3 = [0; 0; (mod_R + mod_z_cm)];
mod_p4 = [0; 0; -(mod_R + mod_z_cm)];
mod_p5 = [mod_R; 0; 0];
mod_p6 = [-mod_R; 0; 0];

% Thruster direction:
mod_t1 = [1; 0; 0];
mod_t2 = [1; 0; 0];
mod_t3 = [0; 1; 0];
mod_t4 = [0; 1; 0];
mod_t5 = [0; 0; 1];
mod_t6 = [0; 0; 1];
