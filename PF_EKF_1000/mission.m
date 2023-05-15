%% Mission Supervisor & Reference Generator parameters

% Inizializzazione variabili

    fine_missione = 0;   % Viene settata ad 1 nello stato finale di fine missione

    stato = 100;         % Variabile utilizzata nella gestione dell'attivazione degli stati 
                         % e nell'aggiornamento del current_state

    contatore = 0;       % Variabile utilizzata nel ciclo di pattugliamento e 
                         % nell'aggiornamento della parete

    Ts_supervisor = 0.1; % Tempo di campionamento

    SOS=0;               % Utilizzata per segnalare eventuali guasti o situazioni di errore
         
    timer = 40*2;        % [s] Attesa massima per la ricezione dei dati del GPS

    

% Soglie per l'emergenza
 
    soglia_SOS_pos = [0.00005, 0.00005, 1]; % [deg, deg, m] Errore tollerato nel mantenimento della profondità 

    soglia_SOS_or = [0.25, 0.25, 0.25]; % [rad] Errore tollerato nel raggiungimento dell'orientazione iniziale

    soglia_SOS_orp = 0.25; % [rad] Errore tollerato nel mantenimento dell'orientazione desiderata

    soglia_SOS_parete = 1; % [m] Errore tollerato durante il pattugliamento 

    soglia_SOS_pos_ris = 2; % [m] Errore tollerato durante la risalita

    soglia_SOS_collisione = 2; % [m] Errore tollerato durante il rientro

    soglia_SOS_vel = 0.2; % [m/s] Errore tollerato nel mantenimento di velocità



% Soglie per le condizioni di switch

    soglia_pos = [0.00001, 0.00001, 0.5]; %[deg,deg,m] Intervallo ammesso nel raggiungimento di posizione

    soglia_or = [0.1, 0.1, 0.1]; % [rad] Intervallo ammesso nel raggiungimento dell'orientazione iniziale

    soglia_orp = 0.1; % [rad] Intervallo ammesso nel raggiungimento dell'orientazione desiderata

    soglia_parete = 0.5; %  [m] Intervallo ammesso per la distanza dalla parete

    soglia_parete_frenata = 0.5; %  [m] Intervallo ammesso nel raggiungimento del punto iniziale di frenata
    
    