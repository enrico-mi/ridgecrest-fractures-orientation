function selected_fractures = select_by_faultdistance(fractures,limit_distance)
% select fractures whose center is more than 'distance' away from 2019
% mainshock fault (at northern termination)

% limit_distance is provided in m

    flag_explode_only = true;
    fractures = calc_center(fractures);
    selected_fractures = select_segments_by_resolution(fractures, 0, false, flag_explode_only);
    tmp_fractures = selected_fractures;

    % Add UTM x and y coordinates
    zone_utm = utmzone(35.8,-117.7);
    disp("Using UTM zone" + string(zone_utm));
    [ellipsoid,~] = utmgeoid(zone_utm);
    utmstruct = defaultm('utm');
    utmstruct.zone = zone_utm;
    utmstruct.geoid = ellipsoid;
    utmstruct = defaultm(utmstruct);
    
    % fault trace, y=mx+q;
    lat1 = 35.8214;
    lon1 = -117.6304;
    lat2 = 35.9380;
    lon2 = -117.7701;
    [x1,y1] = projfwd(utmstruct,lat1,lon1);
    [x2,y2] = projfwd(utmstruct,lat2,lon2);
    q = (x2*y1-x1*y2)/(x2-x1);
    m = (y1-q)/x1;

    for f=length(tmp_fractures):-1:1
        
        x0 = tmp_fractures(f).XCenters;
        y0 = tmp_fractures(f).YCenters;
        distance = distance_from_fault(x0,y0,m,q);        
                
        if distance < limit_distance 
            selected_fractures(f) = [];
        end
        
    end
    
end
    
function fractures = calc_center(fractures)

    for k=1:length(fractures)
        x = fractures(k).UTMx;
        y = fractures(k).UTMy;
        
        xc = (x(2:end)+x(1:end-1))/2;
        yc = (y(2:end)+y(1:end-1))/2;
        
        fractures(k).XCenters = xc;
        fractures(k).YCenters = yc;
    end

end

function distance = distance_from_fault(x0,y0,m,q);
    
    distance = abs(y0-m*x0-q)/sqrt(m*m+1);
    
end