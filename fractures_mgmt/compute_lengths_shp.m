function fractures = compute_lengths_shp(fractures)

    %% Add segment lengths to shp file

    if ~isfield(fractures(1),'UTMx')
        warning('UTM coords not present, computing them now.')
        fractures = calc_UTM(fractures);
    end
    
    for f=1:length(fractures)
        frac_segments = length(fractures(f).UTMx)-1;
        seg_lengths = zeros(frac_segments,1);

        for k=1:frac_segments
            dx2 = (fractures(f).UTMx(k+1)-fractures(f).UTMx(k))^2;
            dy2 = (fractures(f).UTMy(k+1)-fractures(f).UTMy(k))^2;
            seg_lengths(k) = sqrt(dx2+dy2);
        end
        
        fractures(f).Lengths = seg_lengths;
    end
end
