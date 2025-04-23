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

shoulder_angles = zeros(nRows, 3);
elbow_angles = zeros(nRows, 3);
core_angles = zeros(nRows, 3);
pelvis_angles = zeros(nRows, 3);
thorax_angles = zeros(nRows, 3);
kneeL_angles = zeros(nRows, 3);


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
    Mid_EP = (ELR + EMR) / 2;
    Mid_Styloids = (PLR + PMR) / 2; % Wrist joint center approx
    
    % Y_f axis is the unit vector from US (PMR) to Mid_EP
    Y_f = Mid_EP - PMR; 
    Y_f = Y_f / norm(Y_f);
    % Vector defining wrist plane
    v_elbow_plane = PLR - PMR; 
    
    X_f = cross(Y_f, v_elbow_plane);
    X_f = Z_f / norm(Z_f); 
    
    Z_f = cross(X_f, Y_f);
    Z_f = X_f / norm(X_f); % Points anteriorly (check ISB)
    
    % Store attidute matrix
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
    
    
    P_pelvis(:,:,i) = [X_p(:) ; Y_p(:) ; Z_p(:)];


    % --- Calculate Left Thigh (TL) --- (Example - Adapt based on ISB/Wu)
    % Assuming HL = Hip Left, CLL = Condylus Lateralis Left, CML = Condylus Medialis Left
    % Mid_Knee_L = (CLL + CML) / 2;
    % Y_tl = HL - Mid_Knee_L; Y_tl = Y_tl / norm(Y_tl); % Points proximally
    % temp_Z = CLL - CML; % Vector points medially
    % X_tl = cross(Y_tl, temp_Z); X_tl = X_tl / norm(X_tl); % Points anteriorly
    % Z_tl = cross(X_tl, Y_tl); Z_tl = Z_tl / norm(Z_tl); % Points Left (Lateral)
    % TL_thighL(:,:,i) = [X_tl(:), Y_tl(:), Z_tl(:)];

    % --- Calculate Left Shank (SL) --- (Example - Adapt based on ISB/Wu)
    % Assuming MLL = Malleolus Lateralis Left, MML = Malleolus Medialis Left
    % Mid_Ankle_L = (MLL + MML) / 2;
    % Mid_Knee_L = (CLL + CML) / 2; % Needs CML/CLL markers from above
    % Y_sl = Mid_Knee_L - Mid_Ankle_L; Y_sl = Y_sl / norm(Y_sl); % Points proximally
    % temp_Z = MLL - MML; % Vector points medially
    % X_sl = cross(Y_sl, temp_Z); X_sl = X_sl / norm(X_sl); % Points anteriorly
    % Z_sl = cross(X_sl, Y_sl); Z_sl = Z_sl / norm(Z_sl); % Points Left (Lateral)
    % SL_shankL(:,:,i) = [X_sl(:), Y_sl(:), Z_sl(:)];



        % --- Calculate Relative Rotation Matrices ---
    % Shoulder: Upper Arm relative to Thorax (Distal rel Proximal)
    R_rel_UT(:,:,i) = T_thorax(:,:,i)' * U_upperArm(:,:,i);
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
    pelvis_angles(i,:) = rotm2eul(P_pelvis(:,:,i), 'ZXY');

    % Thorax angles relative to Global (ZXY: Sagittal rot, Frontal rot, Transverse rot)
    thorax_angles(i,:) = rotm2eul(T_thorax(:,:,i), 'ZXY');

    % Shoulder (Thorax -> Humerus) - ISB: Y(Thorax)-Y(Humerus)-Floating Z -> 'YXY' common
    % Rotation order: Plane of Elev (Y), Elevation (X'), Axial Rotation (Y'')
    shoulder_angles(i,:) = rotm2eul(R_rel_UT(:,:,i), 'YXY');

    % Elbow (Humerus -> Forearm) - ISB: Z(Hum)-Y(Forearm)-Floating X -> 'ZYX' common
    % Rotation order: Flex/Ext (Z), Varus/Valgus (Y'), Pro/Sup (X'') - Check ISB elbow details
    % Wu Part II Sec 3.5 seems different: e1=Zhum(Flex/Ext), e3=Yforearm(Pro/Sup) -> ZYX common interpretation
    elbow_angles(i,:) = rotm2eul(R_rel_FU(:,:,i), 'ZYX');

    % Core (Pelvis -> Thorax) - ISB Spine: Z(Prox)-Y(Dist)-Floating X -> 'ZXY' interpretation
    % Rotation order: Flex/Ext (Z), Lat Bend (X'), Axial Rot (Y'')
    core_angles(i,:) = rotm2eul(R_rel_TP(:,:,i), 'ZXY');

    % Left Knee (Thigh -> Shank) - Grood & Suntay Knee: Z(Femur)-Y(Tibia)-Floating X -> 'ZXY' common
    % Rotation order: Flex/Ext (Z), Add/Abd (X'), Int/Ext Rot (Y'')
    kneeL_angles(i,:) = rotm2eul(R_rel_STL(:,:,i), 'ZXY');






end















