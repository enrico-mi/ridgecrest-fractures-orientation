function selected_fractures = select_by_coseismic(fractures,coseismic)
% select fractures whose center is more than 'distance' away from any fracture
% that Thompson Jobe et al (2020) re-defined antecedent 2019 earthquake

    if isfield(coseismic,'distance')
        distance = coseismic.distance;
    else
        distance = 30;
    end
    selected_fractures = fractures;
    tmp_fractures = calc_center(fractures);
    
    thompsonjobe = read_tj_shp([-117.74 -117.63],[35.82 35.94]); 
    %thompsonjobe = shaperead('shp_files/Ridgecrest_PreEQ_Features.shp',...
    %                         'UseGeoCoords',true);

    thompsonjobe = calc_UTM(thompsonjobe);
    thompsonjobe = calc_center(thompsonjobe);
    
    maxt = length(thompsonjobe);

    for f=length(tmp_fractures):-1:1
        
        t = 1;
        checked = false;
        while (checked == false) && (t <= maxt)
            
            for ss=1:length(thompsonjobe(t).XCenters)
                
                delta_x = thompsonjobe(t).XCenters(ss) - tmp_fractures(f).XCenters;
                delta_y = thompsonjobe(t).YCenters(ss) - tmp_fractures(f).YCenters;
                
                if any( sqrt(delta_x.^2 + delta_y.^2) < distance )
                    selected_fractures(f) = [];
                    checked = true;
                    break;
                end

            end
            t = t + 1;
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

