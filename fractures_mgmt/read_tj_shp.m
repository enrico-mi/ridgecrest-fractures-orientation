function [tj] = read_tj_shp(lon_lims, lat_lims)
    
    filepath = fullfile('fractures_mgmt', 'shp_files', 'Ridgecrest_PreEQ_Features.shp');

    if (exist('lat_lims','var')) && (~isempty(lat_lims))
        mask = true;
    else
        mask = false;
    end
    
    switch mask
        
      case true
        
        %% Importing fracture lines from Thompson Jobe et al. (2020)
        % Fractures are filtered by class (to avoid artificial and subtle type)
        
        boundaries = [lon_lims(1), lat_lims(1); lon_lims(2), lat_lims(2)];
        
        tj = shaperead(filepath,...
                       'UseGeoCoords',true,...
                       'BoundingBox',boundaries,...
                       'Selector', ...
                       {@(class)...
                        ((strcmp(class,'geologic')) || (isempty(class))),...
                        'CLASS'});

      case false

        
        %% Importing fracture lines from Ponti et al. (2020)
        % Fractures are filtered by origin (tectonic, to avoid uncertain or shaking
        % origins) and event source (foreshock or mainshock)
        
        tj = shaperead(filepath,...
                       'UseGeoCoords',true,...
                       'Selector', ...
                       {@(class)...
                        ((strcmp(class,'geologic')) || (isempty(class))),...
                        'CLASS'});

    end

    % there's a fracture that has Nan values
    for f=length(tj):-1:1
        frac = tj(f);
        if isnan(frac.BoundingBox(1,1))
            tj(f) = [];
        end
    end
    tj = calc_UTM(tj);
end