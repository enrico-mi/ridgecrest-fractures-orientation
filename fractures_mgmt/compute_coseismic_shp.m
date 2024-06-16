function fractures = compute_coseismic_shp(fractures)

    %% Add coseismic y/n flag to shp file
    %  It uses a reference distance of 30m, hardcoded;
    %  Use select_by_coseismic to use different length values.
    
    distance = 30;
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
                % check if any fracture in Thompson Jobe dataset is closer than
                % 30m to the current fracture from Ponti et al dataset
                delta_x = thompsonjobe(t).XCenters(ss) - tmp_fractures(f).XCenters;
                delta_y = thompsonjobe(t).YCenters(ss) - tmp_fractures(f).YCenters;
                
                if any( sqrt(delta_x.^2 + delta_y.^2) < distance )
                    fractures(f).coseismic = 'No';
                    checked = true;
                    break;
                else
                    fractures(f).coseismic = 'Yes';
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
