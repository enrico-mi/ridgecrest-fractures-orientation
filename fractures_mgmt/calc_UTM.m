function fractures = calc_UTM(fractures)

    % Add UTM x and y coordinates
    zone_utm = utmzone(35.8,-117.7);
    disp("Using UTM zone" + string(zone_utm));
    [ellipsoid,estr] = utmgeoid(zone_utm);
    utmstruct = defaultm('utm');
    utmstruct.zone = zone_utm;
    utmstruct.geoid = ellipsoid;
    utmstruct = defaultm(utmstruct);

    for f=1:length(fractures)
        [UTMx, UTMy] = projfwd(utmstruct,fractures(f).Lat(1:end-1),fractures(f).Lon(1:end-1));
        fractures(f).UTMx = UTMx;
        fractures(f).UTMy = UTMy;
    end

end