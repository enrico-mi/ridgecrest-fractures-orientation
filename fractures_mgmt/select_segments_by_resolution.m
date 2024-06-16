function selected_segments = select_segments_by_resolution(fractures, ref_length, negative_resolution, flag_explode_only)
% get input fractures and re-organize them in segments instead, to which the
% length threshold is applied

    selected_fractures = fractures;
    
    selected_segments = [];
    idx_segments = 1;

    if flag_explode_only
        seg_length = 1e6;
        ref_length = 0;
        negative_resolution = false;
    end               
    
    if ref_length < 1e-10
        warning("ref_length = 0, no selection or rediscretization applied, I am just converting the fractures structure to a segment one.");
    end
    
    for f=length(fractures):-1:1
        frac = selected_fractures(f);
        if isfield(frac,'Lat')
            frac_segments = length(frac.Lat)-2; %nb of segments, last entry of Lat is Nan
        elseif isfield(frac,'theta_I') && flag_explode_only
            frac_segments = length(frac.theta_I); %nb of segments
        else
            error('Fractures are not in a supported format');
        end
        
        for s = 1:frac_segments
            if isfield(frac,'UTMx')
                dx2 = (frac.UTMx(s+1)-frac.UTMx(s))^2;
                dy2 = (frac.UTMy(s+1)-frac.UTMy(s))^2;
                seg_length = sqrt(dx2+dy2);
            end
            
            if ( (seg_length > ref_length) && (~negative_resolution) )
                selected_segments(idx_segments).fid = frac.fid;
                if isfield(frac,'AngleUTM')
                    selected_segments(idx_segments).AngleUTM = frac.AngleUTM(s);
                end
                if isfield(frac,'UTMx')
                    selected_segments(idx_segments).UTMx = frac.UTMx(s:s+1);
                end
                if isfield(frac,'UTMy')
                    selected_segments(idx_segments).UTMy = frac.UTMy(s:s+1);
                end
                if isfield(frac,'XCenters')
                    selected_segments(idx_segments).XCenters = frac.XCenters(s);
                end
                if isfield(frac,'YCenters')
                    selected_segments(idx_segments).YCenters = frac.YCenters(s);
                end
                if isfield(frac,'mechanism')
                    selected_segments(idx_segments).mechanism = frac.mechanism;
                end
                if isfield(frac,'theta_I')
                    selected_segments(idx_segments).theta_I = frac.theta_I(s);
                end
                if isfield(frac,'sigma_I')
                    selected_segments(idx_segments).sigma_I = frac.sigma_I(s);
                end                
                if isfield(frac,'sigma_II')
                    selected_segments(idx_segments).sigma_II = frac.sigma_II(s);
                end                
                selected_segments(idx_segments).sid = idx_segments;
                idx_segments = idx_segments+1;
                
            elseif ( (seg_length <= ref_length) && negative_resolution )
                selected_segments(idx_segments).fid = frac.fid;
                selected_segments(idx_segments).AngleUTM = frac.AngleUTM(s);
                selected_segments(idx_segments).UTMx = frac.UTMx(s:s+1);
                selected_segments(idx_segments).UTMy = frac.UTMy(s:s+1);
                if isfield(frac,'mechanism')
                    selected_segments(idx_segments).mechanism = frac.mechanism;
                end
                if isfield(frac,'theta_I')
                    selected_segments(idx_segments).theta_I = frac.theta_I(s);
                end
                selected_segments(idx_segments).sid = idx_segments;
                idx_segments = idx_segments+1;

            end
        end
    end

end