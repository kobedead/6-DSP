function average_jump_size = calculateAverageJumpSize(T_filt)
    % CALCULATEAVERAGEJUMPSIZE Calculates the average jump size for marker data in a table.
    %
    %   average_jump_size = calculateAverageJumpSize(T_filt)
    %
    %   Input:
    %       T_filt: Table containing marker data.
    %
    %   Output:
    %       average_jump_size: Average jump size across all markers and frames.
    
    marker_cols = {'ARX', 'ARY', 'ARZ', 'ELRX', 'ELRY', 'ELRZ', 'EMRX', 'EMRY', 'EMRZ', ...
                   'PLRX', 'PLRY', 'PLRZ', 'PMRX', 'PMRY', 'PMRZ', 'MSX', 'MSY', 'MSZ', ...
                   'PXX', 'PXY', 'PXZ', 'C7X', 'C7Y', 'C7Z', 'T7X', 'T7Y', 'T7Z', ...
                   'SIPSLX', 'SIPSLY', 'SIPSLZ', 'SIPSRX', 'SIPSRY', 'SIPSRZ', ...
                   'SIASLX', 'SIASLY', 'SIASLZ', 'SIASRX', 'SIASRY', 'SIASRZ', ...
                   'CLLX', 'CLLY', 'CLLZ', 'CLRX', 'CLRY', 'CLRZ', 'CMLX', 'CMLY', 'CMLZ', ...
                   'MLLX', 'MLLY', 'MLLZ', 'MMLX', 'MMLY', 'MMLZ'};
    
    total_jump_size = 0;
    num_jumps = 0;
    
    for j = 1:length(marker_cols) / 3 % Iterate through each marker (X, Y, Z components)
        for i = 2:height(T_filt) % Start from the second frame to calculate differences
            prev_marker_data = [T_filt.(marker_cols{3*(j-1)+1})(i-1), ...
                                T_filt.(marker_cols{3*(j-1)+2})(i-1), ...
                                T_filt.(marker_cols{3*j})(i-1)];
            marker_data = [T_filt.(marker_cols{3*(j-1)+1})(i), ...
                           T_filt.(marker_cols{3*(j-1)+2})(i), ...
                           T_filt.(marker_cols{3*j})(i)];
            
            diff_x = abs(marker_data(1) - prev_marker_data(1));
            diff_y = abs(marker_data(2) - prev_marker_data(2));
            diff_z = abs(marker_data(3) - prev_marker_data(3));
            
            jump_size = sqrt(diff_x^2 + diff_y^2 + diff_z^2);
            
            total_jump_size = total_jump_size + jump_size;
            num_jumps = num_jumps + 1;
        end
    end
    
    if num_jumps > 0
        average_jump_size = total_jump_size / num_jumps;
    else
        average_jump_size = 0; % Or NaN, depending on how you want to handle no jumps
    end
end