function mu_matrices_sigma_I_NW(mechanisms_data, realization, randangle, shearFlag, dist_ext, smallerk)
% perform all test cases
    flag_mu = true;
    
    if randangle > 15
        sigma_H_orientations = -15:40;
    else
        sigma_H_orientations = -10:25;
    end
    mu_range = [0.:0.1:1.5 1e16]; % 1e16 == beta=0 == sigma_I
    
    if smallerk
        magn_factors = 0.:0.01:0.1;
    else
        magn_factors = 0.:0.1:1;
    end
    
    resolutions = [0,10];
    
    slip_models = ["JB2";"RB2";"XB2"];
    if shearFlag
        slip_models = ["JB2";"RB2";"XB2";"JB2";"JB2"];
    end
    if ~strcmp(dist_ext,"")
        slip_models = ["JB2"];
    end
    
    for mech = 1:2
        
        if smallerk && mech>1
            continue
        end
        
        if mech==1
            flag_mechanisms = true;
            disp('Analyzing known slip mechanisms');
        else
            flag_mechanisms = false;
            disp('Analyzing all slip mechanisms indifferently');
        end
        
        for r=1:2
            if r==1
                resNeg_selection = false;
                disp('Analyzing fractures longer than provided resolution');
            else
                resNeg_selection = true;
                
                disp('Analyzing negative resolution, i.e. fractures shorter than provided value');
            end
            
            if smallerk && r>1
                continue
            end
            
            angles_analysis_NW(flag_mu, mu_range, sigma_H_orientations, magn_factors,...
                               resolutions, slip_models, flag_mechanisms,...
                               resNeg_selection, mechanisms_data, realization, randangle, dist_ext);
            
        end
    end
end