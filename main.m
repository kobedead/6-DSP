t = readtable("10Ax1.tsv", "FileType", "text", 'Delimiter', '\t');
data = table2array(t);
upsampling_factor = 2;  % Upsampling factor
nOriginal = size(data, 1);
xOriginal = 1:nOriginal;
xInterp = linspace(1, nOriginal, nOriginal * upsampling_factor);
dataInterp = interp1(xOriginal, data, xInterp, 'spline');
t = array2table(dataInterp, 'VariableNames', t.Properties.VariableNames);
    
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

% Store missing data for all frames
all_missing_data = struct(); 

% Calculate the average jump size
avg_jump = calculateAverageJumpSize(T_filt);
disp(['Average jump size in the data: ' num2str(avg_jump)]);

% Now you can use avg_jump to set your allowed_jump_threshold
allowed_jump_threshold = avg_jump * 5 * upsampling_factor;  % Tweak this parameter
disp(['Jump threshold set to: ' num2str(allowed_jump_threshold)]);


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
    CLL = [T_filt.CLLX(i), T_filt.CLLY(i), T_filt.CLLZ(i)];
    CLR = [T_filt.CLRX(i), T_filt.CLRY(i), T_filt.CLRZ(i)];
    CML = [T_filt.CMLX(i), T_filt.CMLY(i), T_filt.CMLZ(i)];
    MLL = [T_filt.MLLX(i), T_filt.MLLY(i), T_filt.MLLZ(i)];
    MML = [T_filt.MMLX(i), T_filt.MMLY(i), T_filt.MMLZ(i)];

    % Test invalid data

    % Check whether this frame has missing marker data (0 values)
    missing_data_info = findMissingMarkerData(T_filt, i, allowed_jump_threshold);
    if ~isempty(fieldnames(missing_data_info))
        if ~isfield(all_missing_data, ['frame_' num2str(i)])
            all_missing_data.(['frame_' num2str(i)]) = struct();
        end
        % Merge missing_data_info into all_missing_data for the current frame
        current_frame_data = all_missing_data.(['frame_' num2str(i)]);
        fields_to_add = fieldnames(missing_data_info);
        for f = 1:length(fields_to_add)
            current_frame_data.(fields_to_add{f}) = missing_data_info.(fields_to_add{f});
        end
        all_missing_data.(['frame_' num2str(i)]) = current_frame_data;
    end


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

    %midpoint PX and T7 (below)
    Mid_PX_and_T7 = (PX+T7)/2;
    %midpoint  MS and C7 (above)
    Mid_MS_and_C7 = (MS + C7)/2 ; 
    
    Y_t = Mid_MS_and_C7 - Mid_PX_and_T7 ; 
    Y_t = Y_t / norm(Y_t);
    
    %Zt perpendicular to MS->C7 and Mid_PX_and_T7 (pointing to right -> check)
    v1 = C7 - MS;
    v2 = Mid_PX_and_T7 - MS;
    
    %needs to be perpendicular to the 2 above vectors
    Z_t = cross(v1,v2);
    Z_t = Z_t / norm(Z_t);
    
    
    %Xt perpendicular to both, point forwards -> check
    X_t = cross(Y_t , Z_t);
    X_t = X_t / norm(X_t);

    %store attitude
    T_thorax(:,:,i) = [X_t(:), Y_t(:), Z_t(:)];

    %                           --- Calculate Pelvis (P) ---

        
    % === Z-axis: from left to right ASIS ===
    Z_p = SIASR - SIASL;
    Z_p = Z_p / norm(Z_p);
    
    % === Midpoints ===
    Mid_SIAS = (SIASR + SIASL) / 2;
    Mid_PSIS = (SIPSR + SIPSL) / 2;
    
    % === Vector from ASISL to PSIS midpoint ===
    w = Mid_PSIS - SIASL;
    
    % === Projection of w onto Z_p ===
    proj_length = dot(w, Z_p);
    proj_vector = proj_length * Z_p;
    
    % === X-axis: vector orthogonal to Z_p, pointing anterior ===
    X_p = w - proj_vector;
    X_p = -X_p / norm(X_p);
    
    
    % === Y-axis: cross product to ensure right-handed system ===
    Y_p = cross(Z_p, X_p);
    Y_p = Y_p / norm(Y_p);
    
    % === Re-orthogonalize just in case ===
    X_p = cross(Y_p, Z_p);
    X_p = X_p / norm(X_p);
    
    P_pelvis(:,:,i) = [X_p(:) , Y_p(:) , Z_p(:)];


    % --- Calculate Left Thigh (TL)
    knee_midpoint = (CLL + CML) / 2;
    Y_tl = SIASL - knee_midpoint;
    Y_tl = Y_tl / norm(Y_tl);
    v1 = CLL - SIASL; %
    v2 = CML - SIASL; %
    plane_normal = cross(v1, v2);
    Z_tl = cross(Y_tl, plane_normal); 
    Z_tl = Z_tl / norm(Z_tl); % normalized z_tl
    X_tl = cross(Y_tl, Z_tl); % x_tl is perpendicular to y_tl and z_tl
    X_tl = X_tl / norm(X_tl); % Normalized x_tl
    TL_thighL(:,:,i) = [X_tl(:), Y_tl(:), Z_tl(:)]; % Attitude matrix


    % --- Calculate Left Shank (SL) --- (Example - Adapt based on ISB/Wu)
    Mid_Ankle_L = (MLL + MML) / 2;
    Mid_Knee_L  = (CLL + CML) / 2;
    % Define coordinate system for left shank (SL)
    Y_sl = Mid_Knee_L - Mid_Ankle_L; Y_sl = Y_sl / norm(Y_sl);  % Proximal
    temp_Z = MLL - MML;  % Vector pointing medial
    X_sl = cross(Y_sl, temp_Z); X_sl = X_sl / norm(X_sl);       % Anterior
    Z_sl = cross(X_sl, Y_sl); Z_sl = Z_sl / norm(Z_sl);         % Lateral
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
    elbow_angles = rotm2eul(R_rel_FU(:,:,i), 'ZXY');
    elbow_angles_deg(i,:) = rad2deg(elbow_angles);

    % Core (Pelvis -> Thorax) - ISB Spine: Z(Prox)-Y(Dist)-Floating X -> 'ZXY' interpretation
    % Rotation order: Flex/Ext (Z), Lat Bend (X'), Axial Rot (Y'')
    core_angles = rotm2eul(R_rel_TP(:,:,i), 'ZXY');
    core_angles_deg(i,:) = rad2deg(core_angles);

    % Left Knee (Thigh -> Shank) - Grood & Suntay Knee: Z(Femur)-Y(Tibia)-Floating X -> 'ZXY' common
    % Rotation order: Flex/Ext (Z), Add/Abd (X'), Int/Ext Rot (Y'')
    kneeL_angles = rotm2eul(R_rel_STL(:,:,i), 'ZXY');
    kneeL_angles_deg(i,:) = rad2deg(kneeL_angles);
end

% --- CHECK FOR MISSING MARKER DATA ---
% Prints out for which body parts at what frame missing marker data 
% was found. Uses a custom function called frames_with_missing_data

disp('Missing Data Summary (including big jumps):');
if exist('all_missing_data', 'var') == 1
    frames_with_missing_data = fieldnames(all_missing_data);
    if isempty(frames_with_missing_data)
        disp('No missing data (0-valued or big jumps) found in any frame.');
    else
        for k = 1:length(frames_with_missing_data)
            frame_name = frames_with_missing_data{k};
            frame_number = str2double(frame_name(7:end)); % Extract frame number
            missing_info = all_missing_data.(frame_name);
            missing_markers = fieldnames(missing_info);
            if ~isempty(missing_markers)
                fprintf('Frame %d: ', frame_number);
                for m = 1:length(missing_markers)
                    marker_problem = missing_markers{m};
                    fprintf('%s ', marker_problem);
                end
                fprintf('\n');
            end
        end
    end
else
    disp('Error: all_missing_data was not created. Missing data detection may have failed.');
end


% --- VERIFY CALCULATION OF ROTATION MATRIXES ---
% Check orthogonality and determinant of the rotation matrixes

dt = 1/fs; % Time step
tolerance = 1e-10;
invalidFrames = [];

for i = 1:nRows
    % List of matrices to check
    names = {'R_rel_UT', 'R_rel_FU', 'R_rel_TP', 'R_rel_STL'};
    matrices = {R_rel_UT(:,:,i), R_rel_FU(:,:,i), R_rel_TP(:,:,i), R_rel_STL(:,:,i)};
    
    for j = 1:length(matrices)
        R = matrices{j};
        orthogonalityError = norm(R' * R - eye(3));
        determinantError = abs(det(R) - 1);
        
        if orthogonalityError > tolerance || determinantError > tolerance
            fprintf('Frame %d: %s failed validation | OrthoErr = %.2e | DetErr = %.2e\n', ...
                i, names{j}, orthogonalityError, determinantError);
            invalidFrames = [invalidFrames; i];
        end
    end
end

if isempty(invalidFrames)
    disp('All rotation matrices are valid.');
else
    disp('Some rotation matrices failed validation. See messages above.');
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
kneeL_velocities_deg_s = diff(kneeL_angles_deg, 1, 1) / dt; % Uncomment

% Append a NaN row or repeat last value to match original time vector length for plotting
pelvis_velocities_deg_s = [pelvis_velocities_deg_s; nan(1,3)];
thorax_velocities_deg_s = [thorax_velocities_deg_s; nan(1,3)];
shoulder_velocities_deg_s = [shoulder_velocities_deg_s; nan(1,3)];
elbow_velocities_deg_s = [elbow_velocities_deg_s; nan(1,3)];
core_velocities_deg_s = [core_velocities_deg_s; nan(1,3)];
kneeL_velocities_deg_s = [kneeL_velocities_deg_s; nan(1,3)]; % Uncomment

% Angular Accelerations (degrees/s^2)
% Differentiate velocities. Apply a filter to velocities before differentiation if noisy.
pelvis_accelerations_deg_s2 = diff(pelvis_velocities_deg_s(1:end-1,:), 1, 1) / dt; % Exclude the NaN row for diff
thorax_accelerations_deg_s2 = diff(thorax_velocities_deg_s(1:end-1,:), 1, 1) / dt;
shoulder_accelerations_deg_s2 = diff(shoulder_velocities_deg_s(1:end-1,:), 1, 1) / dt;
elbow_accelerations_deg_s2 = diff(elbow_velocities_deg_s(1:end-1,:), 1, 1) / dt;
core_accelerations_deg_s2 = diff(core_velocities_deg_s(1:end-1,:), 1, 1) / dt;
kneeL_accelerations_deg_s2 = diff(kneeL_velocities_deg_s(1:end-1,:), 1, 1) / dt; % Uncomment

% Append two NaN rows or repeat last values to match original time vector length for plotting
pelvis_accelerations_deg_s2 = [pelvis_accelerations_deg_s2; nan(2,3)];
thorax_accelerations_deg_s2 = [thorax_accelerations_deg_s2; nan(2,3)];
shoulder_accelerations_deg_s2 = [shoulder_accelerations_deg_s2; nan(2,3)];
elbow_accelerations_deg_s2 = [elbow_accelerations_deg_s2; nan(2,3)];
core_accelerations_deg_s2 = [core_accelerations_deg_s2; nan(2,3)];
kneeL_accelerations_deg_s2 = [kneeL_accelerations_deg_s2; nan(2,3)]; % Uncomment

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

%%%%%%%%%%%%

% --- Angle-Velocity Method for Shoulder–Elbow CRP ---
% Get angular velocities
% vel_shoulder = gradient(signal1_angle, 1/fs);
% vel_elbow = gradient(signal2_angle, 1/fs);
% 
% % Normalize
% normalize = @(x) (std(x) == 0).*zeros(size(x)) + (std(x) ~= 0).* ((x - mean(x)) / std(x));
% angle1_norm = normalize(signal1_angle);
% angle2_norm = normalize(signal2_angle);
% vel1_norm = normalize(vel_shoulder);
% vel2_norm = normalize(vel_elbow);
% 
% % Compute phase angles
% phi1 = atan2(vel1_norm, angle1_norm);
% phi2 = atan2(vel2_norm, angle2_norm);
% crp_shoulder_elbow_av = rad2deg(wrapToPi(phi1 - phi2));


% --- CRP: Shoulder vs Elbow using Angle-Velocity Method ---

% 1. Extract signals
angle_shoulder = shoulder_angles_deg(:,1); % Plane of Elevation
angle_elbow = elbow_angles_deg(:,1);       % Flexion/Extension

% 2. Compute velocities
vel_shoulder = gradient(angle_shoulder, 1/fs);
vel_elbow = gradient(angle_elbow, 1/fs);

% 3. Normalize both angles and velocities (z-score with zero std protection)
normalize = @(x) (std(x) == 0).*zeros(size(x)) + (std(x) ~= 0).* ((x - mean(x)) / std(x));
shoulder_angle_norm = normalize(angle_shoulder);
elbow_angle_norm = normalize(angle_elbow);
shoulder_vel_norm = normalize(vel_shoulder);
elbow_vel_norm = normalize(vel_elbow);

% 4. Compute phase angles
phi_shoulder_av = atan2(shoulder_vel_norm, shoulder_angle_norm);
phi_elbow_av = atan2(elbow_vel_norm, elbow_angle_norm);

% 5. Compute CRP (Angle-Velocity Method)
crp_shoulder_elbow_av = rad2deg(wrapToPi(phi_shoulder_av - phi_elbow_av));

%%%%%%


%  BETWEEN THORAX AND PELVIS (ZXY -> Flexion/Extension component) ---
% Extract angle components for flexion/extension (Z-axis from 'ZXY')
angle_thorax = thorax_angles_deg(:,1);  % Thorax flexion/extension
angle_pelvis = pelvis_angles_deg(:,1);  % Pelvis flexion/extension

% Time vector
fs = 300;
t = (0:length(angle_thorax)-1) / fs;

% Smooth and compute velocities
vel_thorax = gradient(angle_thorax, 1/fs);
vel_pelvis = gradient(angle_pelvis, 1/fs);

% Normalize (z-score method)
normalize = @(x) (std(x) == 0) .* zeros(size(x)) + (std(x) ~= 0) .* ((x - mean(x)) / std(x));
angle_thorax_norm = normalize(angle_thorax);
angle_pelvis_norm = normalize(angle_pelvis);
vel_thorax_norm = normalize(vel_thorax);
vel_pelvis_norm = normalize(vel_pelvis);

% Phase angles (A/V method)
phi_thorax = atan2(vel_thorax_norm, angle_thorax_norm);
phi_pelvis = atan2(vel_pelvis_norm, angle_pelvis_norm);
crp_core_av = rad2deg(wrapToPi(phi_thorax - phi_pelvis));

% Hilbert transform
phi_thorax_h = angle(hilbert(detrend(angle_thorax)));
phi_pelvis_h = angle(hilbert(detrend(angle_pelvis)));
crp_core_hilbert = rad2deg(wrapToPi(phi_thorax_h - phi_pelvis_h));


% CRP BETWEEN PELVIS AND LEFT KNEE (new) ---
% Use flexion/extension component from each
pelvis_flex = pelvis_angles_deg(:,1);
kneeL_flex = kneeL_angles_deg(:,1);

% Time vector
t_crp_pelvis_knee = (0:length(pelvis_flex)-1) / fs;

% Compute and smooth velocities
vel_pelvis = gradient(pelvis_flex, 1/fs);
vel_knee = gradient(kneeL_flex, 1/fs);

% Normalize
pelvis_angle_norm = normalize(pelvis_flex);
knee_angle_norm = normalize(kneeL_flex);
pelvis_vel_norm = normalize(vel_pelvis);
knee_vel_norm = normalize(vel_knee);

% Angle-Velocity phase
phi_pelvis = atan2(pelvis_vel_norm, pelvis_angle_norm);
phi_knee = atan2(knee_vel_norm, knee_angle_norm);
crp_pk_av = rad2deg(wrapToPi(phi_knee - phi_pelvis));

% Hilbert
phi_pelvis_h = angle(hilbert(detrend(pelvis_flex)));
phi_knee_h = angle(hilbert(detrend(kneeL_flex)));
crp_pk_hilbert = rad2deg(wrapToPi(phi_knee_h - phi_pelvis_h));



% --- 8. VISUALIZATION ---
% Plot Euler/Cardan angles, velocities, and accelerations for each segment/joint
% Shoulder
figure('Name', 'Shoulder Kinematics', 'NumberTitle', 'off');
sgtitle('Shoulder Joint Kinematics');

subplot(3,1,1);
plot(t, shoulder_angles_deg(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, shoulder_angles_deg(:,2), 'g', 'LineWidth', 1.5);
plot(t, shoulder_angles_deg(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Shoulder Angles (YXY: Plane of Elev, Elev, Axial Rot)');
xlabel('Time (s)'); ylabel('Angle (degrees)');
legend('Plane of Elev', 'Elevation', 'Axial Rot'); grid on;

subplot(3,1,2);
plot(t, shoulder_velocities_deg_s(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, shoulder_velocities_deg_s(:,2), 'g', 'LineWidth', 1.5);
plot(t, shoulder_velocities_deg_s(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Shoulder Angular Velocities'); xlabel('Time (s)'); ylabel('Velocity (deg/s)'); grid on;

subplot(3,1,3);
plot(t, shoulder_accelerations_deg_s2(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, shoulder_accelerations_deg_s2(:,2), 'g', 'LineWidth', 1.5);
plot(t, shoulder_accelerations_deg_s2(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Shoulder Angular Accelerations'); xlabel('Time (s)'); ylabel('Accel (deg/s^2)'); grid on;

% Elbow
figure('Name', 'Elbow Kinematics', 'NumberTitle', 'off');
sgtitle('Elbow Joint Kinematics');

subplot(3,1,1);
plot(t, elbow_angles_deg(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, elbow_angles_deg(:,2), 'g', 'LineWidth', 1.5);
plot(t, elbow_angles_deg(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Elbow Angles (ZXY: Flex/Ext, Pro/Sup, Varus/Valgus)');
xlabel('Time (s)'); ylabel('Angle (degrees)'); grid on;

subplot(3,1,2);
plot(t, elbow_velocities_deg_s(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, elbow_velocities_deg_s(:,2), 'g', 'LineWidth', 1.5);
plot(t, elbow_velocities_deg_s(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Elbow Angular Velocities'); xlabel('Time (s)'); ylabel('Velocity (deg/s)'); grid on;

subplot(3,1,3);
plot(t, elbow_accelerations_deg_s2(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, elbow_accelerations_deg_s2(:,2), 'g', 'LineWidth', 1.5);
plot(t, elbow_accelerations_deg_s2(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Elbow Angular Accelerations'); xlabel('Time (s)'); ylabel('Accel (deg/s^2)'); grid on;

% Pelvis
figure('Name', 'Pelvis Kinematics', 'NumberTitle', 'off');
sgtitle('Pelvis Kinematics');

subplot(3,1,1);
plot(t, pelvis_angles_deg(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, pelvis_angles_deg(:,2), 'g', 'LineWidth', 1.5);
plot(t, pelvis_angles_deg(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Pelvis Angles (ZXY)'); xlabel('Time (s)'); ylabel('Angle (deg)'); grid on;

subplot(3,1,2);
plot(t, pelvis_velocities_deg_s(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, pelvis_velocities_deg_s(:,2), 'g', 'LineWidth', 1.5);
plot(t, pelvis_velocities_deg_s(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Pelvis Angular Velocities'); xlabel('Time (s)'); ylabel('Velocity (deg/s)'); grid on;

subplot(3,1,3);
plot(t, pelvis_accelerations_deg_s2(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, pelvis_accelerations_deg_s2(:,2), 'g', 'LineWidth', 1.5);
plot(t, pelvis_accelerations_deg_s2(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Pelvis Angular Accelerations'); xlabel('Time (s)'); ylabel('Accel (deg/s^2)'); grid on;

% Thorax
figure('Name', 'Thorax Kinematics', 'NumberTitle', 'off');
sgtitle('Thorax Kinematics');

subplot(3,1,1);
plot(t, thorax_angles_deg(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, thorax_angles_deg(:,2), 'b', 'LineWidth', 1.5);
plot(t, thorax_angles_deg(:,3), 'g', 'LineWidth', 1.5); hold off;
title('Thorax Angles (ZXY)'); xlabel('Time (s)'); ylabel('Angle (deg)'); grid on;

subplot(3,1,2);
plot(t, thorax_velocities_deg_s(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, thorax_velocities_deg_s(:,2), 'g', 'LineWidth', 1.5);
plot(t, thorax_velocities_deg_s(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Thorax Angular Velocities'); xlabel('Time (s)'); ylabel('Velocity (deg/s)'); grid on;

subplot(3,1,3);
plot(t, thorax_accelerations_deg_s2(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, thorax_accelerations_deg_s2(:,2), 'g', 'LineWidth', 1.5);
plot(t, thorax_accelerations_deg_s2(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Thorax Angular Accelerations'); xlabel('Time (s)'); ylabel('Accel (deg/s^2)'); grid on;

% Core
figure('Name', 'Core Kinematics', 'NumberTitle', 'off');
sgtitle('Core (Thorax rel. Pelvis) Kinematics');

subplot(3,1,1);
plot(t, core_angles_deg(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, core_angles_deg(:,2), 'g', 'LineWidth', 1.5);
plot(t, core_angles_deg(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Core Angles (ZXY)'); xlabel('Time (s)'); ylabel('Angle (deg)'); grid on;

subplot(3,1,2);
plot(t, core_velocities_deg_s(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, core_velocities_deg_s(:,2), 'g', 'LineWidth', 1.5);
plot(t, core_velocities_deg_s(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Core Angular Velocities'); xlabel('Time (s)'); ylabel('Velocity (deg/s)'); grid on;

subplot(3,1,3);
plot(t, core_accelerations_deg_s2(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, core_accelerations_deg_s2(:,2), 'g', 'LineWidth', 1.5);
plot(t, core_accelerations_deg_s2(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Core Angular Accelerations'); xlabel('Time (s)'); ylabel('Accel (deg/s^2)'); grid on;

% Left Knee
figure('Name', 'Left Knee Kinematics', 'NumberTitle', 'off');
sgtitle('Left Knee Kinematics');

subplot(3,1,1);
plot(t, kneeL_angles_deg(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, kneeL_angles_deg(:,2), 'g', 'LineWidth', 1.5);
plot(t, kneeL_angles_deg(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Left Knee Angles (ZXY)'); xlabel('Time (s)'); ylabel('Angle (deg)'); grid on;

subplot(3,1,2);
plot(t, kneeL_velocities_deg_s(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, kneeL_velocities_deg_s(:,2), 'g', 'LineWidth', 1.5);
plot(t, kneeL_velocities_deg_s(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Left Knee Angular Velocities'); xlabel('Time (s)'); ylabel('Velocity (deg/s)'); grid on;

subplot(3,1,3);
plot(t, kneeL_accelerations_deg_s2(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, kneeL_accelerations_deg_s2(:,2), 'g', 'LineWidth', 1.5);
plot(t, kneeL_accelerations_deg_s2(:,3), 'b', 'LineWidth', 1.5); hold off;
title('Left Knee Angular Accelerations'); xlabel('Time (s)'); ylabel('Accel (deg/s^2)'); grid on;


% --- Plot CRP Comparison: Shoulder vs Elbow ---
figure('Name', 'CRP: Shoulder vs Elbow', 'NumberTitle', 'off');
plot(t, crp_shoulder_elbow_av, 'b', 'LineWidth', 1.5); hold on;
plot(t, crp_shoulder_elbow_deg_norm, 'r--', 'LineWidth', 1.5);
legend('Angle-Velocity CRP', 'Hilbert CRP');
title('CRP: Shoulder (Plane of Elevation) vs Elbow (Flexion)');
xlabel('Time (s)');
ylabel('CRP (degrees)');
ylim([-180 180]); yticks(-180:60:180);
grid on;

sample_frame = 150;
fprintf('[Shoulder–Elbow CRP] t = %.2fs | A/V = %.1f° | Hilbert = %.1f°\n', ...
        sample_frame/fs, crp_shoulder_elbow_av(sample_frame), crp_shoulder_elbow_deg_norm(sample_frame));


% Plot comparison: THORAX AND PELVIS
figure('Name','CRP Thorax vs Pelvis','NumberTitle','off');
plot(t, crp_core_av, 'b', 'LineWidth', 1.5); hold on;
plot(t, crp_core_hilbert, 'r--', 'LineWidth', 1.5);
legend('Angle-Velocity CRP','Hilbert CRP');
title('CRP: Thorax vs Pelvis (Flexion/Extension)');
xlabel('Time (s)'); ylabel('CRP (degrees)');
ylim([-180 180]); yticks(-180:60:180);
grid on;

% --- Output sample CRP value for Thorax–Pelvis at sample_frame
fprintf('[Thorax-Pelvis CRP] t = %.2fs | A/V = %.1f deg | Hilbert = %.1f deg\n', ...
    t(sample_frame), crp_core_av(sample_frame), crp_core_hilbert(sample_frame));


% Plot Comparision: CPR Pelvis and Knee
figure('Name','CRP Pelvis vs Knee','NumberTitle','off');
plot(t_crp_pelvis_knee, crp_pk_av, 'b', 'LineWidth', 1.5); hold on;
plot(t_crp_pelvis_knee, crp_pk_hilbert, 'r--', 'LineWidth', 1.5);
legend('Angle-Velocity CRP','Hilbert CRP');
title('CRP: Pelvis vs Left Knee (Flexion)');
xlabel('Time (s)'); ylabel('CRP (degrees)');
ylim([-180 180]); yticks([-180:60:180]);
grid on;

% Sample print
fprintf('[Pelvis-Knee CRP] t = %.2fs | A/V = %.1f deg | Hilbert = %.1f deg\n', ...
    t_crp_pelvis_knee(sample_frame), crp_pk_av(sample_frame), crp_pk_hilbert(sample_frame));



% --- Interactive Phase Plane Visualization: Shoulder vs Elbow (Plane of Elev vs Flexion) ---

% Extract angle and velocity
angle_s = shoulder_angles_deg(:,1);       % Plane of Elevation
vel_s = shoulder_velocities_deg_s(:,1);   % Angular velocity

% Time vector
timeVec = (0:length(angle_s)-1) / fs;

% Initialize figure
fig = figure('Name','Interactive Phase Plane: Shoulder vs Elbow','NumberTitle','off');
ax = axes('Parent', fig);
plot(ax, angle_s, vel_s, 'b.'); % Trajectory
hold on;
startPt = plot(ax, angle_s(1), vel_s(1), 'go', 'MarkerSize', 10, 'DisplayName', 'Start');
endPt = plot(ax, angle_s(end), vel_s(end), 'ro', 'MarkerSize', 10, 'DisplayName', 'End');
markerPt = plot(ax, angle_s(1), vel_s(1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k');
title(ax, 'Shoulder Angle-Velocity Phase Plane');
xlabel(ax, 'Angle (deg)');
ylabel(ax, 'Velocity (deg/s)');
legend(ax, 'Trajectory', 'Start', 'End', 'Selected Point');
grid on;

% Slider
uicontrol('Style', 'text', 'String', 'Frame:', 'Position', [20, 20, 50, 20]);
slider = uicontrol('Style', 'slider',...
    'Min',1,'Max',length(angle_s),'Value',1,...
    'Position', [70, 20, 300, 20],...
    'SliderStep', [1/(length(angle_s)-1) , 10/(length(angle_s)-1)]);

% Label for info
infoText = uicontrol('Style','text', 'Position', [380, 20, 250, 20], ...
    'String', sprintf('t = %.2fs | Angle = %.1f° | Vel = %.1f°/s', timeVec(1), angle_s(1), vel_s(1)));

% Define callback function inline with access to outer variables
slider.Callback = @(src,~) updateMarker(round(src.Value), angle_s, vel_s, timeVec, markerPt, infoText);
updateMarker(1, angle_s, vel_s, timeVec, markerPt, infoText);

% Update function (moved outside to work with explicit args)
function updateMarker(frameIdx, angle_s, vel_s, timeVec, markerPt, infoText)
    markerPt.XData = angle_s(frameIdx);
    markerPt.YData = vel_s(frameIdx);
    infoText.String = sprintf('t = %.2fs | Angle = %.1f° | Vel = %.1f°/s', ...
                              timeVec(frameIdx), angle_s(frameIdx), vel_s(frameIdx));
end


% --- Interactive Phase Plane: Thorax vs Pelvis (Flexion/Extension) ---

% Extract thorax and pelvis flexion/extension angles (ZXY → 1st column = Flex/Ext)
angle_tp = thorax_angles_deg(:,1);       % Thorax flex/ext
vel_tp = thorax_velocities_deg_s(:,1);   % Thorax angular velocity

% Time vector
time_tp = (0:length(angle_tp)-1) / fs;

% Initialize figure
fig_tp = figure('Name','Interactive Phase Plane: Thorax vs Pelvis','NumberTitle','off');
ax_tp = axes('Parent', fig_tp);
plot(ax_tp, angle_tp, vel_tp, 'b.'); % Trajectory
hold on;
start_tp = plot(ax_tp, angle_tp(1), vel_tp(1), 'go', 'MarkerSize', 10, 'DisplayName', 'Start');
end_tp = plot(ax_tp, angle_tp(end), vel_tp(end), 'ro', 'MarkerSize', 10, 'DisplayName', 'End');
marker_tp = plot(ax_tp, angle_tp(1), vel_tp(1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k');
title(ax_tp, 'Thorax Angle-Velocity Phase Plane');
xlabel(ax_tp, 'Angle (deg)');
ylabel(ax_tp, 'Velocity (deg/s)');
legend(ax_tp, 'Trajectory', 'Start', 'End', 'Selected Point');
grid on;

% Slider
uicontrol('Style', 'text', 'String', 'Frame:', 'Position', [20, 20, 50, 20]);
slider_tp = uicontrol('Style', 'slider',...
    'Min',1,'Max',length(angle_tp),'Value',1,...
    'Position', [70, 20, 300, 20],...
    'SliderStep', [1/(length(angle_tp)-1) , 10/(length(angle_tp)-1)]);

% Info Text
infoText_tp = uicontrol('Style','text', 'Position', [380, 20, 300, 20], ...
    'String', sprintf('t = %.2fs | Angle = %.1f° | Vel = %.1f°/s', time_tp(1), angle_tp(1), vel_tp(1)));

% Callback function
slider_tp.Callback = @(src,~) updateMarkerTP(round(src.Value), angle_tp, vel_tp, time_tp, marker_tp, infoText_tp);
updateMarkerTP(1, angle_tp, vel_tp, time_tp, marker_tp, infoText_tp);

% Callback definition
function updateMarkerTP(frameIdx, angle_tp, vel_tp, time_tp, marker_tp, infoText_tp)
    marker_tp.XData = angle_tp(frameIdx);
    marker_tp.YData = vel_tp(frameIdx);
    infoText_tp.String = sprintf('t = %.2fs | Angle = %.1f° | Vel = %.1f°/s', ...
                                 time_tp(frameIdx), angle_tp(frameIdx), vel_tp(frameIdx));
end




% --- Interactive Phase Plane: Pelvis vs Left Knee (Flexion/Extension) ---

% Extract pelvis and left knee flexion/extension angles
angle_pk = kneeL_angles_deg(:,1);           % Left knee flex/ext
vel_pk = kneeL_velocities_deg_s(:,1);       % Left knee angular velocity

% Time vector
time_pk = (0:length(angle_pk)-1) / fs;

% Initialize figure
fig_pk = figure('Name','Interactive Phase Plane: Pelvis vs Left Knee','NumberTitle','off');
ax_pk = axes('Parent', fig_pk);
plot(ax_pk, angle_pk, vel_pk, 'b.'); % Trajectory
hold on;
start_pk = plot(ax_pk, angle_pk(1), vel_pk(1), 'go', 'MarkerSize', 10, 'DisplayName', 'Start');
end_pk = plot(ax_pk, angle_pk(end), vel_pk(end), 'ro', 'MarkerSize', 10, 'DisplayName', 'End');
marker_pk = plot(ax_pk, angle_pk(1), vel_pk(1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k');
title(ax_pk, 'Pelvis vs Left Knee Angle-Velocity Phase Plane');
xlabel(ax_pk, 'Angle (deg)');
ylabel(ax_pk, 'Velocity (deg/s)');
legend(ax_pk, 'Trajectory', 'Start', 'End', 'Selected Point');
grid on;

% Slider
uicontrol('Style', 'text', 'String', 'Frame:', 'Position', [20, 20, 50, 20]);
slider_pk = uicontrol('Style', 'slider',...
    'Min',1,'Max',length(angle_pk),'Value',1,...
    'Position', [70, 20, 300, 20],...
    'SliderStep', [1/(length(angle_pk)-1) , 10/(length(angle_pk)-1)]);

% Info Text
infoText_pk = uicontrol('Style','text', 'Position', [380, 20, 300, 20], ...
    'String', sprintf('t = %.2fs | Angle = %.1f° | Vel = %.1f°/s', time_pk(1), angle_pk(1), vel_pk(1)));

% Callback function
slider_pk.Callback = @(src,~) updateMarkerPK(round(src.Value), angle_pk, vel_pk, time_pk, marker_pk, infoText_pk);
updateMarkerPK(1, angle_pk, vel_pk, time_pk, marker_pk, infoText_pk);

% Callback definition
function updateMarkerPK(frameIdx, angle_pk, vel_pk, time_pk, marker_pk, infoText_pk)
    marker_pk.XData = angle_pk(frameIdx);
    marker_pk.YData = vel_pk(frameIdx);
    infoText_pk.String = sprintf('t = %.2fs | Angle = %.1f° | Vel = %.1f°/s', ...
                                 time_pk(frameIdx), angle_pk(frameIdx), vel_pk(frameIdx));
end






% --- 9. TIME POINTS ---
% Calculate time points
% 1. Foot contact left leg (FC)
% 2. Ball release (BR)
% 3. Maximal external rotation of the right shoulder (MER)

% 1. Foot contact left leg (FC)
avgMMLX = mean(T_filt.MMLX(200:end));

xMMLX = T_filt.MMLX(200:end); % x pos of MMLX after frame 200
v_xMMLX = gradient(xMMLX, dt); % Velocity
ax_xMMLX = gradient(v_xMMLX, dt); % Acceleration

vx_smooth = smoothdata(v_xMMLX, 'movmean', 5); % Smoothing velocity
ax_smooth = smoothdata(ax_xMMLX, 'movmean', 5); % Smoothing accel

[~, idx_min_acc] = min(ax_smooth); % min acceleration

% Map back to original frame index
frame_contact = 199 + idx_min_acc;  % because xMMLX starts at frame 200

% Also identify peak velocity for comparison (optional)
[~, idx_peak_v] = max(vx_smooth);
frame_peak_velocity = 199 + idx_peak_v;

% 2. Ball release (BR)
frame_FC = frame_contact;

% Extract PLR components from FC onward
PLR = [T_filt.PLRX(frame_FC:end), T_filt.PLRY(frame_FC:end), T_filt.PLRZ(frame_FC:end)];

% Compute velocity in each direction
v_PLR = gradient(PLR, dt);  % rows: frames, cols: X/Y/Z

% Compute magnitude of 3D velocity
v_PLR_mag = sqrt(sum(v_PLR.^2, 2));

% Smooth velocity magnitude (optional)
v_PLR_mag_smooth = smoothdata(v_PLR_mag, 'movmean', 5);

% Find peak velocity (ball release)
[~, idx_BR] = max(v_PLR_mag_smooth);
frame_BR = frame_FC - 1 + idx_BR;  % map to global frame index

% Full PLR data for plotting
PLR_full = [T_filt.PLRX, T_filt.PLRY, T_filt.PLRZ];
v_full = gradient(PLR_full, dt);
a_full = gradient(v_full, dt);

v_mag_full = sqrt(sum(v_full.^2, 2));
a_mag_full = sqrt(sum(a_full.^2, 2));


% 3. Maximal external rotation of the right shoulder (MER)
axial_rotation = shoulder_angles_deg(frame_FC:end, 3);  % external/internal rotation

% Min rotation shows maximum external rotation
% Max rotation show maxmimum internal rotation
[max_ext_rot, idx_max_rot] = min(axial_rotation);

frame_MER = frame_FC - 1 + idx_max_rot;

% -------------------------------
% Timeline Visualization
% -------------------------------

% Unified Timeline Visualization + Kinematics Verification
t_rel = ((1:height(T_filt)) - frame_MER) / fs;
x_limits = [-nRows/fs, nRows/fs];

figure('Name','Kinematic Overview - Aligned to MER', 'NumberTitle','off');

% --- Subplot 1: Left Malleoli  (X) ---
subplot(4,2,1);
plot(t_rel(200:end), xMMLX, 'k');
hold on;
xline((frame_FC - frame_MER)/fs, '--r');
xline(0, '--m');
xline((frame_BR - frame_MER)/fs, '--g');
text((frame_FC - frame_MER)/fs, max(xMMLX)*1.05, 'FC', 'Color','r', 'HorizontalAlignment','center');
text(0, max(xMMLX)*1.05, 'MER', 'Color','m', 'HorizontalAlignment','center');
text((frame_BR - frame_MER)/fs, max(xMMLX)*1.05, 'BR', 'Color','g', 'HorizontalAlignment','center');
xlabel('Time (s)')
ylabel('Position (mm)');
title('X Position of Left Malleoli');
xlim(x_limits)
grid on;

% --- Subplot 2: Left Malleoli Velocity (X) ---
subplot(4,2,3);
plot(t_rel(200:end), vx_smooth, 'k');
hold on;
xline((frame_FC - frame_MER)/fs, '--r');
xline(0, '--m');
xline((frame_BR - frame_MER)/fs, '--g');
xlabel('Time (s)')
ylabel('Vel. (mm/s)');
title('Left Malleoli Velocity');
xlim(x_limits)
grid on;

% --- Subplot 3: Left Malleoli Acceleration (X) ---
subplot(4,2,5);
plot(t_rel(200:end), ax_smooth, 'k');
hold on;
xline((frame_FC - frame_MER)/fs, '--r');
xline(0, '--m');
xline((frame_BR - frame_MER)/fs, '--g');
xlabel('Time (s)')
ylabel('Acc. (mm/s²)');
title('Left Malleoli Acceleration');
xlim(x_limits)
grid on;

% --- Subplot 4: PLR X Position ---
subplot(4,2,2);
plot(t_rel, T_filt.PLRX);
hold on;
xline((frame_FC - frame_MER)/fs, '--r');
xline(0, '--m');
xline((frame_BR - frame_MER)/fs, '--g');
text((frame_FC - frame_MER)/fs, max(T_filt.PLRX)*1.05, 'FC', 'Color','r', 'HorizontalAlignment','center');
text(0, max(T_filt.PLRX)*1.05, 'MER', 'Color','m', 'HorizontalAlignment','center');
text((frame_BR - frame_MER)/fs, max(T_filt.PLRX)*1.05, 'BR', 'Color','g', 'HorizontalAlignment','center');
title('PLR x Position ');
xlabel('Time (s)')
ylabel('Position (mm)');
xlim(x_limits)
grid on;

% --- Subplot 2: PLR Velocity Magnitude ---
subplot(4,2,4);
plot(t_rel, v_mag_full, 'b');
hold on;
xline((frame_FC - frame_MER)/fs, '--r');
xline(0, '--m');
xline((frame_BR - frame_MER)/fs, '--g');
xlabel('Time (s)')
ylabel('Velocity (mm/s)');
title('PLR Velocity Magnitude');
xlim(x_limits)
grid on;

% --- Subplot 3: Shoulder Axial Rotation ---
subplot(4,2,6);
plot(t_rel, shoulder_angles_deg(:,3), 'b');
hold on;
xline((frame_FC - frame_MER)/fs, '--r');
xline(0, '--m');
xline((frame_BR - frame_MER)/fs, '--g');
xlabel('Time (s)')
ylabel('Angle (°)');
title('Shoulder Axial Rotation');
xlim(x_limits)
grid on;

% --- 3D graph of PLR trajectory ---
figure('Name','PLR 3D Trajectory', 'NumberTitle','off');
plot3(PLR_full(:,1), PLR_full(:,2), PLR_full(:,3), 'b');
hold on;
plot3(PLR_full(frame_FC,1), PLR_full(frame_FC,2), PLR_full(frame_FC,3), 'ro', 'DisplayName','Foot Contact');
plot3(PLR_full(frame_BR,1), PLR_full(frame_BR,2), PLR_full(frame_BR,3), 'go', 'DisplayName','Ball Release');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('PLR Marker 3D Trajectory');
legend;
grid on;
axis equal;
view(3);
rotate3d on;


% Summary of time points
disp(['--- Time points ---'])

time_FC = frame_FC / upsampling_factor / fs;
time_MER = frame_MER / upsampling_factor / fs;
time_BR = frame_BR / upsampling_factor / fs;

disp(['The left foot contacted the ground at frame: ' num2str(frame_FC / upsampling_factor) ...
      ' (Time: ' num2str(time_FC, '%.3f') ' s)']);

disp(['Maximum external rotation occurred at frame: ' num2str(frame_MER / upsampling_factor) ...
      ' (Time: ' num2str(time_MER, '%.3f') ' s)']);

disp(['The ball was released at frame: ' num2str(frame_BR / upsampling_factor) ...
      ' (Time: ' num2str(time_BR, '%.3f') ' s)']);








