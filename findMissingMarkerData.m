function missing_data_info = findMissingMarkerData(T_filt, i, allowedJump)
    % FINDMISSINGMARKERDATA Finds missing marker data (0-valued) and big jumps for a given frame.
    %
    %   missing_data_info = findMissingMarkerData(T_filt, i, allowedJump)
    %
    %   Input:
    %       T_filt:      Table containing marker data.
    %       i:           Frame number to check.
    %       allowedJump: Maximum allowed absolute difference between consecutive
    %                    frames for a data point to be considered valid.
    %
    %   Output:
    %       missing_data_info: Structure containing names of markers with
    %                          missing data (0-valued or big jumps) and their
    %                          corresponding frame numbers.

    missing_data_info = struct();
    marker_names = {'AR', 'ELR', 'EMR', 'PLR', 'PMR', 'MS', 'PX', 'C7', 'T7', ...
                    'SIPSL', 'SIPSR', 'SIASL', 'SIASR', 'CLL', 'CLR', 'CML', ...
                    'MLL', 'MML'};
    marker_cols = {'ARX', 'ARY', 'ARZ', 'ELRX', 'ELRY', 'ELRZ', 'EMRX', 'EMRY', 'EMRZ', ...
                   'PLRX', 'PLRY', 'PLRZ', 'PMRX', 'PMRY', 'PMRZ', 'MSX', 'MSY', 'MSZ', ...
                   'PXX', 'PXY', 'PXZ', 'C7X', 'C7Y', 'C7Z', 'T7X', 'T7Y', 'T7Z', ...
                   'SIPSLX', 'SIPSLY', 'SIPSLZ', 'SIPSRX', 'SIPSRY', 'SIPSRZ', ...
                   'SIASLX', 'SIASLY', 'SIASLZ', 'SIASRX', 'SIASRY', 'SIASRZ', ...
                   'CLLX', 'CLLY', 'CLLZ', 'CLRX', 'CLRY', 'CLRZ', 'CMLX', 'CMLY', 'CMLZ', ...
                   'MLLX', 'MLLY', 'MLLZ', 'MMLX', 'MMLY', 'MMLZ'};

    for j = 1:length(marker_names)
        marker_data = [T_filt.(marker_cols{3*j-2})(i), ...
                       T_filt.(marker_cols{3*j-1})(i), ...
                       T_filt.(marker_cols{3*j})(i)];

        % Check for zero values (missing data)
        if any(marker_data == 0)
            if ~isfield(missing_data_info, marker_names{j})
                missing_data_info.(marker_names{j}) = [];
            end
            missing_data_info.(marker_names{j}) = ...
                [missing_data_info.(marker_names{j}); i];
        end

        % Check for big jumps (if not the first frame)
        if i > 1
            prev_marker_data = [T_filt.(marker_cols{3*j-2})(i-1), ...
                                T_filt.(marker_cols{3*j-1})(i-1), ...
                                T_filt.(marker_cols{3*j})(i-1)];

            diff_x = abs(marker_data(1) - prev_marker_data(1));
            diff_y = abs(marker_data(2) - prev_marker_data(2));
            diff_z = abs(marker_data(3) - prev_marker_data(3));

            if diff_x > allowedJump || diff_y > allowedJump || diff_z > allowedJump
                if ~isfield(missing_data_info, [marker_names{j} '_jump'])
                    missing_data_info.([marker_names{j} '_jump']) = [];
                end
                missing_data_info.([marker_names{j} '_jump']) = ...
                    [missing_data_info.([marker_names{j} '_jump']); i];
            end
        end
    end
end