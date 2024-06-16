function [csv_filename, mat_filename, flag_skip] = generateFilename(flag_mu,flag_mechanisms,SSR_type,SSR_exp,string_stress_angles,string_resolution,slip_model,realization,randang,shear_ext,dist_ext,smallerk_ext)

    if flag_mu
        filename = "./output/angles_analysis/" + slip_model +...
            "_mu_N" + string_stress_angles +...
            "E_res" + string_resolution +...
            "_" + SSR_type + string(SSR_exp);
        
    else
        filename = "./output/angles_analysis/" + slip_model +...
            "_beta_N" + string_stress_angles +...
            "E_res" + string_resolution +...
            "_" + SSR_type + string(SSR_exp);
    end
    
    if flag_mechanisms
        filename = filename + "_mechs";
    end
    
    filename = filename + "_randang" + num2str(randang,'%02.f');
    
    %filename = filename + "_fracture";

    filename = filename + "_realization" + num2str(realization,'%03.f');
    
    filename = filename + shear_ext;
    
    filename = filename + dist_ext;
    
    filename = filename + smallerk_ext;

    csv_filename = filename + "_degrees.csv";
    mat_filename = filename + "_degrees.mat";

    if isfile(csv_filename)
        flag_skip = true;
    else
        flag_skip = false;
    end
    
end