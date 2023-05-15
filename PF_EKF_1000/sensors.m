%% Sensors Parameters 
sens_Variance_GPS = 3;      %[m^2]
sens_Sample_time_GPS = 1;   %[s]
sens_Resolution_GPS = 0.1;  %[m]
sens_threshold = 0.485;     %[m]         % [m] distanza dal centro di massa alla
                                    % superficie (0.385) + 0.1 (antenna
                                    % GPS)

sens_Variance_depth1 = (0.2)^2;    % [m^2]
sens_Variance_depth2 = (0.4)^2;    % [m^2] 
sens_Sample_time_depth = 0.1;  % [s]
sens_Resolution_depth = 0.012; % [m]

sens_Variance_DVL = (speed*0.03)^2;    % [(m/s)^2]
sens_Sample_time_dvl = 0.2;            % [s]
sens_Resolution_dvl = 0.001;           % [m/s]

sens_Variance_AHRS_pitch = (0.0087)^2;   % [rad^2]
sens_Variance_AHRS_roll = (0.0087)^2;    % [rad^2]
sens_Variance_AHRS_yaw = (0.0017)^2;     % [rad^2]
sens_Sample_time_AHRS = 0.1;         % [s]
sens_Resolution_AHRS = 0.0008726;    % [rad]

sens_Variance_GYRO = 2.5133e-06;   % [(rad/s)^2]
sens_Sample_time_GYRO = 0.1;       % [s]
sens_Resolution_GYRO = 0.0008726;  % [rad/s]
sens_Bias_GYRO= 3.8785e-05;        

sens_Variance_sonar = 0.05;   % [m^2]
sens_Resolution_sonar = 0.02; % [m]
sens_Sample_time_sonar = 0.5; % [s]
sens_distCM=0.267;            % [m]
sens_BeamWidth=10;            %[gradi]
