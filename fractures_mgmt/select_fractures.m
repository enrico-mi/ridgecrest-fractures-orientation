function [selected_fractures, flag_fractures] = select_fractures(fractures, selection, flag_compute_angles)
% select fractures read from shp file based on some characteristics
% fractures: structures of cells containing fractures
% selection: structure containing information about selection to perform
%
% selection.type (with type=verified,coseismic,resolution or angles): substructure containing parameters needed for each selection type
% selection.type.boolean: if True, perform selection 'type'
% selection types and parameters:
% selection.verified: select fractures from Ponti's dataset based on whether they've been verified by the authors or partially verified. Note that read_fractures_shp does this too and might be faster.
% selection.verified.condition: 'Yes', 'Partial' or 'All' (not that 'All' simply keeps the same selection as provided).
% selection.coseismic: select fractures that Thompson Jobe et al (2020) did not identify as existing before 2019 sequence; it does so by comparing the geometrical center of each fracture with the geometrical center of each fracture in Thompson Jobe dataset, if it is closer than a threshold distance, the fractures is discarded.
% selection.coseismic.distance: threshold distance (defaults to 50m).
  
% selection.method: select fractures from Ponti's dataset based on the observation method
    
% selection.method.list: vector of booleans to select or not the 4 possible observation methods; the order of the entries is ["Field", "Imagery", "Inferred", "Remote Sensing"]. All four values need to be specified, and at least one has to be 1 (true). Example: [0 1 1 0] select fractures observed with imagery or remote sensing, discarding field and inferred observations.
    
    selected_fractures = fractures;

    % flag needed to avoid continuing computation when selected subset is empty
    flag_fractures = check_fractures(selected_fractures);
    
    if selection.method.boolean && flag_fractures
        selected_fractures = select_by_method(selected_fractures,selection.method.list);
        flag_fractures = check_fractures(selected_fractures);
    end
    
    if selection.verified.boolean && flag_fractures
        selected_fractures = select_by_verified(selected_fractures,selection.verified.condition);
        flag_fractures = check_fractures(selected_fractures);
    end

    if selection.coseismic.boolean  && flag_fractures
        selected_fractures = select_by_coseismic(selected_fractures,selection.coseismic);
        flag_fractures = check_fractures(selected_fractures);
    end
    
    if selection.angles.boolean && flag_fractures
        selected_fractures = select_by_angles(selected_fractures,selection.angles);
        flag_fractures = check_fractures(selected_fractures);
    end
    
    if flag_compute_angles && flag_fractures
        selected_fractures = compute_angles_shp(selected_fractures);
    end
    
    if selection.resolution.boolean && flag_fractures
        if ~(isfield(selection.resolution,'ref_length'))
            error('No reference length provided for resolution selection.');
        elseif selection.resolution.segments
            selected_fractures  = select_segments_by_resolution(selected_fractures, selection.resolution.ref_length, selection.resolution.negative,false);
        elseif ~selection.resolution.negative
            [selected_fractures,~] = lengthscale_fractures(selected_fractures,selection.resolution.ref_length);
        elseif selection.resolution.negative
            [~, selected_fractures] = lengthscale_fractures(selected_fractures,selection.resolution.ref_length);
        end
        flag_fractures = check_fractures(selected_fractures);
    end

    if ~flag_fractures
        warning('Selected subset is empty.');
    end
end

function flag_fractures = check_fractures(selected_fractures)
    if size(selected_fractures,2) == 0
        flag_fractures = false;
    else
        flag_fractures = true;        
    end
end