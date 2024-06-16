function SRw = getSRw(string_stress_angles, mu_idx, magn_factor_idx, slip_model, realization, randang, flag_mechanisms, res_ext, shear_ext, dist_ext, smallerk_ext)
    
    flag_mu = true;
    if strcmp(res_ext,"")
        string_resolution = "0000";
    else
        string_resolution = res_ext;
    end
    SSR_type = "errW";
    SSR_exp = 1;

    % only getting 100m depth with magn_factors
    [csv_filename, ~, ~] = generateFilename(flag_mu,flag_mechanisms,SSR_type,SSR_exp,string_stress_angles,string_resolution,slip_model,realization,randang,shear_ext,dist_ext,smallerk_ext);
    
    SRw_matrix = readmatrix(csv_filename);

    % SRw_matrix has dimensions mu (0, 0.1, 0.2, ..., 1, inf, rand) and k
    % (0:0.1:1) for all orientations;
    shapeSRw = size(SRw_matrix);
    debugFlag = 0;
    if ~(shapeSRw(1)==18 && shapeSRw(2) == 11)
        debugFlag = 1;
    end
    SRw = SRw_matrix(mu_idx,magn_factor_idx);

end