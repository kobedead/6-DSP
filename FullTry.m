t = readtable("10Ax1.tsv", "FileType","text",'Delimiter', '\t');
% Extract data
data = t{:,:};
[nRows, nCols] = size(data);



%First we need to filter the data to reduce noice.
%For this we can use a low-pass Butterworth filter with zero-phase
%distortion (filtfilt function) with good cutoff freq -> need to search
%explaination

% Example: Applying a 4th order Butterworth filter with 10Hz cutoff
fs = 300; % Original sampling frequency
fc = 10;  % Cutoff frequency (example, adjust as needed)
Wn = fc / (fs/2); % Normalized cutoff frequency
[b, a] = butter(4, Wn, 'low'); % 4th order low-pass Butterworth filter

% Assuming 'data' is your nRows x nCols matrix of marker coordinates
filteredData = zeros(size(data));
for i = 1:nCols
    % Apply filter column by column
    filteredData(:, i) = filtfilt(b, a, data(:, i));
end

% Convert filtered data back to a table 
T_filt = array2table(filteredData, 'VariableNames', t.Properties.VariableNames);



% PREALLOCATION of attitude matrices 
U_upperArm = zeros(3, 3, nRows);
F_forearm = zeros(3, 3, nRows);
T_thorax = zeros(3, 3, nRows);
P_pelvis = zeros(3, 3, nRows);
TL_thighL = zeros(3, 3, nRows); 
SL_shankL = zeros(3, 3, nRows); 

% PREALLOCATION of rotaion matrices 
R_rel_UT = zeros(3, 3, nRows); % Shoulder
R_rel_FU = zeros(3, 3, nRows); % Elbow
R_rel_TP = zeros(3, 3, nRows); % Core
R_rel_STL = zeros(3, 3, nRows);% Left Knee

shoulder_angles_deg = zeros(nRows, 3);
elbow_angles_deg = zeros(nRows, 3);
core_angles_deg = zeros(nRows, 3);
pelvis_angles_deg = zeros(nRows, 3);
thorax_angles_deg = zeros(nRows, 3);
kneeL_angles_deg = zeros(nRows, 3);
% LOOP THROUGH FRAMES 
for i = 1:nRows
    % Extract all necessary markers for this frame
  
    AR = [T_filt.ARX(i), T_filt.ARY(i), T_filt.ARZ(i)];
    ELR = [T_filt.ELRX(i), T_filt.ELRY(i), T_filt.ELRZ(i)];
    EMR = [T_filt.EMRX(i), T_filt.EMRY(i), T_filt.EMRZ(i)]; 
    PLR = [T_filt.PLRX(i), T_filt.PLRY(i), T_filt.PLRZ(i)];
    PMR = [T_filt.PMRX(i), T_filt.PMRY(i), T_filt.PMRZ(i)];

    MS = [T_filt.MSX(i), T_filt.MSY(i), T_filt.MSZ(i)];
    PX = [T_filt.PXX(i), T_filt.PXY(i), T_filt.PXZ(i)];
    C7 = [T_filt.C7X(i), T_filt.C7Y(i), T_filt.C7Z(i)];
    T7 = [T_filt.T7X(i), T_filt.T7Y(i), T_filt.T7Z(i)];
    SIPSL = [T_filt.SIPSLX(i), T_filt.SIPSLY(i), T_filt.SIPSLZ(i)];
    SIPSR = [T_filt.SIPSRX(i), T_filt.SIPSRY(i), T_filt.SIPSRZ(i)];
    SIASL = [T_filt.SIASLX(i), T_filt.SIASLY(i), T_filt.SIASLZ(i)];
    SIASR = [T_filt.SIASRX(i), T_filt.SIASRY(i), T_filt.SIASRZ(i)];

    % Extract markers for Left Leg (assuming HL, CLL, CML, MLL, MML are available)
    % HL = ... CLL = ... CML = ... MLL = ... MML = ...

    CLL = [T_filt.CLLX(i), T_filt.CLLY(i), T_filt.CLLZ(i)];
    CML = [T_filt.CMLX(i), T_filt.CMLY(i), T_filt.CMLZ(i)];
    CLR = [T_filt.CLRX(i), T_filt.CLRY(i), T_filt.CLRZ(i)];
    CMR = [T_filt.CMRX(i), T_filt.CMRY(i), T_filt.CMRZ(i)];
    MLL = [T_filt.MLLX(i), T_filt.MLLY(i), T_filt.MLLZ(i)];
    MML = [T_filt.MMLX(i), T_filt.MMLY(i), T_filt.MMLZ(i)];
    MLR = [T_filt.MLRX(i), T_filt.MLRY(i), T_filt.MLRZ(i)];
    MMR = [T_filt.MMRX(i), T_filt.MMRY(i), T_filt.MMRZ(i)];



    %                       --- Calculate Upper Arm (U) ---
    Mid_ELR_EMR = (ELR + EMR) / 2;
    % Define Y_h2: from the midpoint (elbow) to the shoulder (AR)
    Y_h = AR - Mid_ELR_EMR;
    Y_h = Y_h / norm(Y_h);
    % Define Y_f: from PMR to the midpoint (this vector helps define the plane)
    Y_f_temp = Mid_ELR_EMR - PMR;       
    Y_f_temp = Y_f_temp / norm(Y_f_temp); 
    % Calculate Z_h2: perpendicular to Y_h2 and Y_f
    Z_h = cross(Y_h, Y_f_temp);
    Z_h = Z_h / norm(Z_h);
    % Calculate X_h2: perpendicular to Z_h2 and Y_h2 (to complete right-handed system)
    X_h = cross(Y_h, Z_h); 
    X_h = X_h / norm(X_h);
    
    % Store in attitude matrix 
    U_upperArm(:,:,i) = [X_h(:), Y_h(:), Z_h(:)]; % Store as columns


    %                       --- Calculate Forearm (F) ---
    Mid_EP = (ELR + EMR) / 2;  % Midpoint of the epicondyles
    v_1 = PLR - PMR;
    v_2 = Mid_EP - PMR;
    X_f = cross(v_1, v_2);
    X_f = X_f / norm(X_f);
    Y_f = Mid_EP - PMR;
    Y_f = Y_f / norm(Y_f);
    Z_f = cross(X_f, Y_f);
    Z_f = Z_f / norm(Z_f);
    F_forearm(:,:,i) = [X_f(:), Y_f(:), Z_f(:)];

    %                    --- Calculate Thorax (T) ---
    Mid_PX_T7 = (PX + T7) / 2;
    Mid_MS_C7 = (MS + C7) / 2;

    Y_t = Mid_MS_C7 - Mid_PX_T7;
    Y_t = Y_t / norm(Y_t); 
    
    % Vector roughly pointing left->right or right->left
    temp_X = MS - C7; 
    
    %Zt perpendicular to MS->C7 and Mid_PX_and_T7 (pointing to right -> check)
    Z_t = cross(Y_t, temp_X);
    Z_t = Z_t / norm(Z_t); 
   
    %Xt perpendicular to both, point forwards -> check
    X_t = cross(Y_t, Z_t);  
    X_t = X_t / norm(X_t);
       
    %store attitude
    T_thorax(:,:,i) = [X_t(:), Y_t(:), Z_t(:)];

    %                           --- Calculate Pelvis (P) ---

    % Z-axis -> line between SIAS to the right(negative) 
    Z_p = SIASR-SIASL;
    Z_p = Z_p / norm(Z_p);
    
    
    %X-axis (negative)
    Midpoint_SIAS = (SIASR+SIASL)/2;
    Midpoint_SIPS = (SIPSR+SIPSL)/2;
    
    X_p = Midpoint_SIAS - Midpoint_SIPS;
    X_p = X_p * norm(X_p);
    
    %Y-axis (check orientation)
    Y_p = cross(Z_p , X_p);
    Y_p = Y_p / norm(Y_p);
    
    
    P_pelvis(:,:,i) = [X_p(:) , Y_p(:) , Z_p(:)];


    % --- Calculate Left Thigh (TL)
    knee_midpoint = (CLL + CLR) / 2;
    Y_tl = SIASL - knee_midpoint;
    Y_tl = Y_tl / norm(Y_tl);
    v1 = CLL - SIASL; %
    v2 = CLR - SIASL; %
    plane_normal = cross(v1, v2);
    Z_tl = cross(Y_f, plane_normal); 
    Z_tl = Z_tl / norm(Z_tl); % normalized z_tl
    X_tl = cross(Y_tl, Z_tl); % x_tl is perpendicular to y_tl and z_tl
    X_tl = X_tl / norm(X_tl); % Normalized x_tl
    TL_thighL(:,:,i) = [X_tl(:), Y_tl(:), Z_tl(:)]; % Attitude matrix


    % --- Calculate Left Shank (SL) --- (Example - Adapt based on ISB/Wu)
    % Assuming MLL = Malleolus Lateralis Left, MML = Malleolus Medialis Left
    Mid_Ankle_L = (MLL + MML) / 2;
    Mid_Knee_L = (CLL + CML) / 2; % Needs CML/CLL markers from above
    Y_sl = Mid_Knee_L - Mid_Ankle_L; Y_sl = Y_sl / norm(Y_sl); % Points proximally
    temp_Z = MLL - MML; % Vector points medially
    X_sl = cross(Y_sl, temp_Z); X_sl = X_sl / norm(X_sl); % Points anteriorly
    Z_sl = cross(X_sl, Y_sl); Z_sl = Z_sl / norm(Z_sl); % Points Left (Lateral)
    SL_shankL(:,:,i) = [X_sl(:), Y_sl(:), Z_sl(:)];






        % --- Calculate Relative Rotation Matrices ---
    % Shoulder: Upper Arm relative to Thorax (Distal rel Proximal)
    R_rel_UT(:,:,i) =   T_thorax(:,:,i)' * U_upperArm(:,:,i);
    % Elbow: Forearm relative to Upper Arm
    R_rel_FU(:,:,i) = U_upperArm(:,:,i)' * F_forearm(:,:,i);
    % Core: Thorax relative to Pelvis
    R_rel_TP(:,:,i) = P_pelvis(:,:,i)' * T_thorax(:,:,i);
    % Left Knee: Shank relative to Thigh
    R_rel_STL(:,:,i) = TL_thighL(:,:,i)' * SL_shankL(:,:,i);

    

    % --- Calculate Euler/Cardan Angles ---
    % SEQUENCE CHOICE IS CRITICAL - Based on ISB JCS definitions / common practice
    % Ensure results match clinical definitions (e.g., Flexion positive) - may need sign changes

    % Pelvis angles relative to Global (ZXY: Sagittal rot, Frontal rot, Transverse rot)
    pelvis_angles = rotm2eul(P_pelvis(:,:,i), 'ZXY');
    pelvis_angles_deg(i,:) = rad2deg(pelvis_angles);


    % Thorax angles relative to Global (ZXY: Sagittal rot, Frontal rot, Transverse rot)
    thorax_angles = rotm2eul(T_thorax(:,:,i), 'ZXY');
    thorax_angles_deg(i,:) = rad2deg(thorax_angles);

    % Shoulder (Thorax -> Humerus) - ISB: Y(Thorax)-Y(Humerus)-Floating Z -> 'YXY' common
    % Rotation order: Plane of Elev (Y), Elevation (X'), Axial Rotation (Y'')
    shoulder_angles = rotm2eul(R_rel_UT(:,:,i), 'YXY');
    shoulder_angles_deg(i,:) = rad2deg(shoulder_angles);


    % Elbow (Humerus -> Forearm) - ISB: Z(Hum)-Y(Forearm)-Floating X -> 'ZYX' common
    % Rotation order: Flex/Ext (Z), Varus/Valgus (Y'), Pro/Sup (X'') - Check ISB elbow details
    % Wu Part II Sec 3.5 seems different: e1=Zhum(Flex/Ext), e3=Yforearm(Pro/Sup) -> ZYX common interpretation
    elbow_angles = rotm2eul(R_rel_FU(:,:,i), 'ZYX');
    elbow_angles_deg(i,:) = rad2deg(elbow_angles);

    % Core (Pelvis -> Thorax) - ISB Spine: Z(Prox)-Y(Dist)-Floating X -> 'ZXY' interpretation
    % Rotation order: Flex/Ext (Z), Lat Bend (X'), Axial Rot (Y'')
    core_angles = rotm2eul(R_rel_TP(:,:,i), 'ZXY');
    core_angles_deg(i,:) = rad2deg(core_angles);

    % Left Knee (Thigh -> Shank) - Grood & Suntay Knee: Z(Femur)-Y(Tibia)-Floating X -> 'ZXY' common
    % Rotation order: Flex/Ext (Z), Add/Abd (X'), Int/Ext Rot (Y'')
    kneeL_angles = rotm2eul(R_rel_STL(:,:,i), 'ZXY');
    core_angles_deg(i,:) = rad2deg(kneeL_angles);



end




% --- 6. CALCULATE ROTATIONAL VELOCITIES AND ACCELERATIONS ---
% Using numerical differentiation (gradient function is suitable for unevenly spaced data, but time is even here)
% Or use diff and divide by dt. For centered difference, adjust indexing.
dt = 1/fs; % Time step

% Angular Velocities (degrees/s)
pelvis_velocities_deg_s = diff(pelvis_angles_deg, 1, 1) / dt;
thorax_velocities_deg_s = diff(thorax_angles_deg, 1, 1) / dt;
shoulder_velocities_deg_s = diff(shoulder_angles_deg, 1, 1) / dt;
elbow_velocities_deg_s = diff(elbow_angles_deg, 1, 1) / dt;
core_velocities_deg_s = diff(core_angles_deg, 1, 1) / dt;
% kneeL_velocities_deg_s = diff(kneeL_angles_deg, 1, 1) / dt; % Uncomment

% Append a NaN row or repeat last value to match original time vector length for plotting
pelvis_velocities_deg_s = [pelvis_velocities_deg_s; nan(1,3)];
thorax_velocities_deg_s = [thorax_velocities_deg_s; nan(1,3)];
shoulder_velocities_deg_s = [shoulder_velocities_deg_s; nan(1,3)];
elbow_velocities_deg_s = [elbow_velocities_deg_s; nan(1,3)];
core_velocities_deg_s = [core_velocities_deg_s; nan(1,3)];
% kneeL_velocities_deg_s = [kneeL_velocities_deg_s; nan(1,3)]; % Uncomment

% Angular Accelerations (degrees/s^2)
% Differentiate velocities. Apply a filter to velocities before differentiation if noisy.
pelvis_accelerations_deg_s2 = diff(pelvis_velocities_deg_s(1:end-1,:), 1, 1) / dt; % Exclude the NaN row for diff
thorax_accelerations_deg_s2 = diff(thorax_velocities_deg_s(1:end-1,:), 1, 1) / dt;
shoulder_accelerations_deg_s2 = diff(shoulder_velocities_deg_s(1:end-1,:), 1, 1) / dt;
elbow_accelerations_deg_s2 = diff(elbow_velocities_deg_s(1:end-1,:), 1, 1) / dt;
core_accelerations_deg_s2 = diff(core_velocities_deg_s(1:end-1,:), 1, 1) / dt;
% kneeL_accelerations_deg_s2 = diff(kneeL_velocities_deg_s(1:end-1,:), 1, 1) / dt; % Uncomment

% Append two NaN rows or repeat last values to match original time vector length for plotting
pelvis_accelerations_deg_s2 = [pelvis_accelerations_deg_s2; nan(2,3)];
thorax_accelerations_deg_s2 = [thorax_accelerations_deg_s2; nan(2,3)];
shoulder_accelerations_deg_s2 = [shoulder_accelerations_deg_s2; nan(2,3)];
elbow_accelerations_deg_s2 = [elbow_accelerations_deg_s2; nan(2,3)];
core_accelerations_deg_s2 = [core_accelerations_deg_s2; nan(2,3)];
% kneeL_accelerations_deg_s2 = [kneeL_accelerations_deg_s2; nan(2,3)]; % Uncomment

% --- 7. CALCULATE CONTINUOUS RELATIVE PHASE (CRP) ---
% CRP is calculated between two oscillating signals.
% Example: CRP between shoulder flexion/extension and elbow flexion/extension
% Requires selection of specific angle components.
% Shoulder Angle 1 (Plane of Elev) vs Elbow Angle 1 (Flex/Ext)

% Ensure signals are zero-mean or appropriately detrended if necessary for Hilbert transform
signal1_angle = shoulder_angles_deg(:,1); % Example: Shoulder Plane of Elevation
signal2_angle = elbow_angles_deg(:,1);   % Example: Elbow Flexion/Extension

% Optional: Filter angles again if CRP is noisy (use a bandpass if a specific frequency range is expected)
% [b_crp, a_crp] = butter(2, [0.5 5]/(fs/2), 'bandpass'); % Example: 0.5-5 Hz activity
% signal1_angle_filt = filtfilt(b_crp, a_crp, signal1_angle);
% signal2_angle_filt = filtfilt(b_crp, a_crp, signal2_angle);

% Calculate phase using Hilbert transform
phase_signal1 = angle(hilbert(detrend(signal1_angle))); % Phase in radians
phase_signal2 = angle(hilbert(detrend(signal2_angle))); % Phase in radians

% Calculate CRP
crp_shoulder_elbow_rad = phase_signal1 - phase_signal2;

% Unwrap CRP to avoid jumps (optional, depends on interpretation)
% crp_shoulder_elbow_rad_unwrapped = unwrap(crp_shoulder_elbow_rad);

% Convert CRP to degrees if preferred
crp_shoulder_elbow_deg = rad2deg(crp_shoulder_elbow_rad);

% Normalize CRP to be between -180 and 180 or 0 and 360
crp_shoulder_elbow_deg_norm = mod(crp_shoulder_elbow_deg + 180, 360) - 180;


% --- 8. VISUALIZATION ---
% Plot Euler/Cardan angles, velocities, and accelerations for each segment/joint
% Example for Shoulder joint


% Shoulder Angles
figure('Name', 'Shoulder Kinematics', 'NumberTitle', 'off', 'WindowState', 'maximized');
sgtitle('Shoulder Joint Kinematics'); % Super title for the figure

timeVector = 1 : 1  : nRows;


subplot(3,1,1);
plot(timeVector, shoulder_angles_deg(:,1));
hold on;
plot(timeVector, shoulder_angles_deg(:,2));
plot(timeVector, shoulder_angles_deg(:,3));
hold off;
title('Shoulder Angles (YXY: Plane of Elev, Elev, Axial Rot)');
xlabel('Time (s)');
ylabel('Angle (degrees)');
legend('Plane of Elev (Y)', 'Elevation (X'')', 'Axial Rot (Y'''')', 'Location', 'best');
grid on;

% Shoulder Velocities
subplot(3,1,2);
plot(timeVector, shoulder_velocities_deg_s(:,1), 'r', 'LineWidth', 1.5);
hold on;
plot(timeVector, shoulder_velocities_deg_s(:,2), 'g', 'LineWidth', 1.5);
plot(timeVector, shoulder_velocities_deg_s(:,3), 'b', 'LineWidth', 1.5);
hold off;
title('Shoulder Angular Velocities');
xlabel('Time (s)');
ylabel('Velocity (degrees/s)');
legend('Vel Plane of Elev', 'Vel Elevation', 'Vel Axial Rot', 'Location', 'best');
grid on;

% Shoulder Accelerations
subplot(3,1,3);
plot(timeVector, shoulder_accelerations_deg_s2(:,1), 'r', 'LineWidth', 1.5);
hold on;
plot(timeVector, shoulder_accelerations_deg_s2(:,2), 'g', 'LineWidth', 1.5);
plot(timeVector, shoulder_accelerations_deg_s2(:,3), 'b', 'LineWidth', 1.5);
hold off;
title('Shoulder Angular Accelerations');
xlabel('Time (s)');
ylabel('Acceleration (degrees/s^2)');
legend('Acc Plane of Elev', 'Acc Elevation', 'Acc Axial Rot', 'Location', 'best');
grid on;


% Example for Elbow joint (similar structure)
figure('Name', 'Elbow Kinematics', 'NumberTitle', 'off', 'WindowState', 'maximized');
sgtitle('Elbow Joint Kinematics');

subplot(3,1,1);
    plot(timeVector, elbow_angles_deg(:,1), 'r', 'LineWidth', 1.5); hold on;
    plot(timeVector, elbow_angles_deg(:,2), 'g', 'LineWidth', 1.5);
    plot(timeVector, elbow_angles_deg(:,3), 'b', 'LineWidth', 1.5); hold off;
    legend('Flex/Ext (Z)', 'Pro/Sup (X'')', 'Varus/Valgus (Y'''')', 'Location','best');

title('Elbow Angles (ZXY: Flex/Ext, Pro/Sup, Varus/Valgus)');
xlabel('Time (s)'); ylabel('Angle (degrees)');
grid on;

% Subplot 2: Elbow Velocities
subplot(3,1,2);
    plot(timeVector, elbow_velocities_deg_s(:,1), 'r', 'LineWidth', 1.5); hold on;
    plot(timeVector, elbow_velocities_deg_s(:,2), 'g', 'LineWidth', 1.5);
    plot(timeVector, elbow_velocities_deg_s(:,3), 'b', 'LineWidth', 1.5); hold off;
    legend('Vel Flex/Ext', 'Vel Pro/Sup', 'Vel Varus/Valgus', 'Location','best');
title('Elbow Angular Velocities');
xlabel('Time (s)'); ylabel('Velocity (degrees/s)');
grid on;

% Subplot 3: Elbow Accelerations
subplot(3,1,3);
    plot(timeVector, elbow_accelerations_deg_s2(:,1), 'r', 'LineWidth', 1.5); hold on;
    plot(timeVector, elbow_accelerations_deg_s2(:,2), 'g', 'LineWidth', 1.5);
    plot(timeVector, elbow_accelerations_deg_s2(:,3), 'b', 'LineWidth', 1.5); hold off;
    legend('Acc Flex/Ext', 'Acc Pro/Sup', 'Acc Varus/Valgus', 'Location','best');

title('Elbow Angular Accelerations');
xlabel('Time (s)'); ylabel('Acceleration (degrees/s^2)');
grid on;

% Plot CRP
figure('Name', 'Continuous Relative Phase (CRP)', 'NumberTitle', 'off', 'WindowState', 'maximized');
    plot(timeVector, crp_shoulder_elbow_deg_norm, 'm', 'LineWidth', 1.5);
    ylim([-180 180]); 
    yticks(-180:60:180);
title('CRP: Shoulder (Plane of Elev) vs Elbow (Flex/Ext)');
xlabel('Time (s)');
ylabel('CRP (degrees)');
grid on;

















