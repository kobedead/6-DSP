% Read and filter data
t = readtable("10Ax1.tsv", "FileType","text",'Delimiter', '\t');
data = t{:,:};
[nRows, nCols] = size(data);

% Butterworth filter (4th order, low-pass at 10Hz)
fs = 300; % Sampling frequency
fc = 10;  % Cutoff
Wn = fc / (fs/2); 
[b, a] = butter(4, Wn, 'low');

% Filter data column-wise
filteredData = zeros(size(data));
for i = 1:nCols
    filteredData(:, i) = filtfilt(b, a, data(:, i));
end

% Convert filtered matrix back to table
T_filt = array2table(filteredData, 'VariableNames', t.Properties.VariableNames);

% Extract filtered marker positions for frame i
i = 1;

% Left leg markers
CLL = [T_filt.CLLX(i), T_filt.CLLY(i), T_filt.CLLZ(i)];
CML = [T_filt.CMLX(i), T_filt.CMLY(i), T_filt.CMLZ(i)];
MLL = [T_filt.MLLX(i), T_filt.MLLY(i), T_filt.MLLZ(i)];
MML = [T_filt.MMLX(i), T_filt.MMLY(i), T_filt.MMLZ(i)];

% Midpoints
Mid_Ankle_L = (MLL + MML) / 2;
Mid_Knee_L  = (CLL + CML) / 2;

% Define coordinate system for left shank (SL)
Y_sl = Mid_Knee_L - Mid_Ankle_L; Y_sl = Y_sl / norm(Y_sl);  % Proximal
temp_Z = MLL - MML;  % Vector pointing medial
X_sl = cross(Y_sl, temp_Z); X_sl = X_sl / norm(X_sl);       % Anterior
Z_sl = cross(X_sl, Y_sl); Z_sl = Z_sl / norm(Z_sl);         % Lateral

F_shankL = [X_sl(:), Y_sl(:), Z_sl(:)];

% Visualization
figure;
hold on; grid on; axis equal;
scale_factor = 1000;

% Plot key anatomical points
scatter3(MLL(1), MLL(2), MLL(3), 100, 'b', 'filled'); % Malleolus Lateralis
scatter3(MML(1), MML(2), MML(3), 100, 'r', 'filled'); % Malleolus Medialis
scatter3(CLL(1), CLL(2), CLL(3), 100, 'g', 'filled'); % Lateral Knee
scatter3(CML(1), CML(2), CML(3), 100, 'c', 'filled'); % Medial Knee

% Draw shank coordinate system at ankle midpoint
quiver3(Mid_Ankle_L(1), Mid_Ankle_L(2), Mid_Ankle_L(3), ...
    X_sl(1)*scale_factor, X_sl(2)*scale_factor, X_sl(3)*scale_factor, 0, 'r', 'LineWidth', 2); % X_sl (anterior)
quiver3(Mid_Ankle_L(1), Mid_Ankle_L(2), Mid_Ankle_L(3), ...
    Y_sl(1)*scale_factor, Y_sl(2)*scale_factor, Y_sl(3)*scale_factor, 0, 'g', 'LineWidth', 2); % Y_sl (proximal)
quiver3(Mid_Ankle_L(1), Mid_Ankle_L(2), Mid_Ankle_L(3), ...
    Z_sl(1)*scale_factor, Z_sl(2)*scale_factor, Z_sl(3)*scale_factor, 0, 'b', 'LineWidth', 2); % Z_sl (lateral)

legend('MLL','MML','CLL','CML','X_{shankL}','Y_{shankL}','Z_{shankL}');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Left Shank Coordinate System');
view(3);
hold off;
