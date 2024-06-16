function fractures = compute_angles_shp(fractures)

    %% Computing angles
    
    if ~isfield(fractures(1),'UTMx')
        warning('UTM coords not present, computing them now.')
        fractures = calc_UTM(fractures);
    end

    for k=1:length(fractures)

        frac_angle_utm = atan( diff(fractures(k).UTMy) ./ diff(fractures(k).UTMx) );
        frac_angle_ave_utm = atan( (fractures(k).UTMy(end) - fractures(k).UTMy(1))/(fractures(k).UTMx(end) - fractures(k).UTMx(1)) );
        sin_arg = sum(diff(fractures(k).UTMy));
        cos_arg = sum(diff(fractures(k).UTMx));
        angle_ave_utm = atan(sin_arg/cos_arg);
        fractures(k).AngleUTM = frac_angle_utm;
        fractures(k).AngleAveUTM = angle_ave_utm;
        
        if (angle_ave_utm - frac_angle_ave_utm) > 1e-15
            error("Different average computations give different results!");
        end
        
    end

end

