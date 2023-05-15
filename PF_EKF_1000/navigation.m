%% Questo script deve essere eseguito solo una volta durante la simulazione
% va eseguito solo una volta che sono stati caricati i dati del bacino nel
% file missionC.m

nav_filter = -1;             %  1 = EKF
                            % -1 = PF

%% TEMPO DI CAMPIONAMENTO DEL FILTRO DI KALMAN

nav_DT = 0.1;

res_enable_n_eff=1; %ricampionamento per N eff

nav_resampling=10000;  %periodo [s] per ricampionamento sistematico
                       %posto un valore alto verrà sfruttata la tecnica
                       %delle particelle effettive, posto a 0.1 verrà
                       %eseguito il ricampionamento ad ogni iterazione

nav_N_part=1000;

nav_contatore = [0 0 0];  %[DVL AHRS Gyro]

%% 1. Calcolo coefficenti dei piani delle pareti del bacino

% Definizione terna NED con origine (lat0,long0,h0) coincidente al
% cornerC

nav_lat0 = cornerC(1);
nav_lon0 = cornerC(2);
nav_h0 = 0;

% Conversione ECEF to NED dei corner

[nav_A(1),nav_A(2),nav_A(3)]= geodetic2ned(cornerA(1),cornerA(2),0,nav_lat0,nav_lon0,nav_h0,wgs84Ellipsoid);
[nav_B(1),nav_B(2),nav_B(3)]= geodetic2ned(cornerB(1),cornerB(2),0,nav_lat0,nav_lon0,nav_h0,wgs84Ellipsoid);
[nav_C(1),nav_C(2),nav_C(3)]= geodetic2ned(cornerC(1),cornerC(2),0,nav_lat0,nav_lon0,nav_h0,wgs84Ellipsoid);
[nav_D(1),nav_D(2),nav_D(3)]= geodetic2ned(cornerD(1),cornerD(2),0,nav_lat0,nav_lon0,nav_h0,wgs84Ellipsoid);

%Calcolo coefficienti dei piani delle pareti
nav_plane_AB = coefficient_plane(nav_A,nav_B);
nav_plane_BC = coefficient_plane(nav_B,nav_C);
nav_plane_CD = coefficient_plane(nav_C,nav_D);
nav_plane_DA = coefficient_plane(nav_D,nav_A);

%inizializzazione variabile per la direzione di pattugliamento
if strcmp(inspectionDirection,'clockwise')
    nav_giro = 1;           %orario
else
    nav_giro = 2;           %antiorario
end
%%



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo coefficienti di un piano verticale passanti per due punti X1 e X2
% in ingresso, in uscita si ha un vettore riga (1x4) con i coefficienti dell'equazione del piano
% cioè y = [a b c d] con equazione ax + by + cz + d = 0

function  y = coefficient_plane(x1,x2)     
	y(1) = 1;                               %coefficiente a = 1
    y(2) = (x1(1)-x2(1))/(x2(2)-x1(2));     %coefficiente b
    y(3) = 0;                               %coefficiente c = 0 (perché il piano è verticale)
    y(4) = -x1(1)-y(2)*x1(2);               %coefficiente d
end
     
