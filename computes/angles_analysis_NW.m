function angles_analysis_NW(flag_mu, mu_beta, stress_angles, magn_factors, resolutions, slip_models, flag_mechanisms, resNeg_selection, mechanisms_data, realization, randangle, dist_ext)

% files_list: list of output files from run_calc_stresses; first element 
% flag_mu: determines if second parameter is mu (true) or beta (false)
% mu_beta: rock internal friction or critical angle beta from sigma_H
% [beta: failure angle with respect to sigma_1 direction: 
%       beta = pi/4 - atan(mu)/2 ]
% stress_angles: regional stress orientation in degrees with respect from North,
%             positive if clockwise (e.g. stress_angles=5 corresponds to N5E)
% depths: depths at which to compute angles, in m (e.g. 100 for 100m)
% selection: structure for select_fractures function

    selection_params.resolution.negative = resNeg_selection;

    SSR_type = "errW";
    SSR_exp = 1;

    nb_cases = length(stress_angles)*length(magn_factors)*length(resolutions)*length(slip_models);
    if resNeg_selection && ~all(resolutions)        
        nb_cases = length(stress_angles)*length(magn_factors)*(length(resolutions)-1)*length(slip_models);
    end
    cases_counter = 1;

    % Loop #0 over slip models
    for sm = 1:length(slip_models)
        slip_model_name = slip_models(sm);
        
        if sm == 4
            shear_ext = "_10GPa";
        elseif sm == 5
            shear_ext = "_20GPa";
        else
            shear_ext = '';
        end
        
        if sm > 1
            if res > 0 || resNeg_selection
                continue
            end
        end       
        
        % Loop #1 over resolutions - SSR_matrix created inside here
        for r = 1:length(resolutions)
            
            res = resolutions(r);
            
            if res > 0 && ( mean(magn_factors) > 0 && mean(magn_factors) < 0.1 )
                continue
            end
            
            
            %skip negative for zero resolution (it does not apply)
            if res < 1e-10 && selection_params.resolution.negative
                disp("negative resolution with zero reference length does not exist, skipping case");
                continue
            end
            
            string_res = num2str(res,'%04.f');
            if resNeg_selection
                string_res = "Neg" + string_res;
            end
            
            selection_params.resolution.ref_length = res;
            % Loop #2 over stress orientation
            for sa = 1:length(stress_angles)
                % SSR_matrix created with dimensions mu_values, magnifying factors
                SSR_matrix = zeros( length(mu_beta)+1, length(magn_factors) );
                
                string_stress_angles = num2str(stress_angles(sa),'%02.f');
                %disp("Analyzing angle " + string_stress_angles);
                
                % generating file name (one per stress orientation and
                % resolution + selection params, but not separate for
                % magn_factors or mu which are the next two loops)
    
                lower_k = "";
                if mean(magn_factors) > 0 && mean(magn_factors) < 0.1
                    lower_k = "_lowerk";
                end
                [csv_filename, mat_filename, flag_skip] = generateFilename(flag_mu,flag_mechanisms,SSR_type,SSR_exp,string_stress_angles,string_res,slip_model_name,realization,randangle,shear_ext,dist_ext,lower_k);
                
                if flag_skip
                    disp("File named " + csv_filename + " already exists, skipping.");
                    continue;
                end
                
                % Loop #3 over range of magnifying factors
                for mf = 1:length(magn_factors)
                    magn_factor = magn_factors(mf); % this is called k in the manuscript
                    if magn_factor > 0 && magn_factor < 0.1
                        string_mf = num2str(magn_factor*100,'%03.f');
                    else
                        string_mf = num2str(magn_factor*10,'%02.f');
                    end

                    sim_name = slip_model_name + "_N" + string_stress_angles + "E_d0100_res" + string_res + "_mf" + string_mf + "_NW";
                    
                    if sm > 3
                        sim_name = sim_name + shear_ext;
                    end
                    
                    postEQ_analysis = load("./output/"+sim_name+"/"+sim_name+"_fractures_data.mat");
                    fractures_data = postEQ_analysis.fractures_data;
                    clearvars postEQ_analysis;
                    
                    % reminding structure of fractures_data
                    % {1} = analysis name (string)
                    % {2} = slip model name (string)
                    % {3} = angle of regional stress (string)
                    % {4} = coseisimic intensity / magnifying factor (float)
                    % {5} = resolution (string)
                    % {6} = selection parameters for fracture datasets (struct, for select_fractures.m)
                    % {7} = fractures with their id and theta_I (struct)
                    
                    if ~strcmp(sim_name, fractures_data{1})
                        error("name of analysis asked for and of analysis retrieved do not match");
                    end
                    fractures = fractures_data{7};
                    % "explode" fracture into segments
                    %fractures = select_segments_by_resolution(fractures, 0, false, true);
                    
                    % mechanisms_data has been filtered to contain only the
                    % desired fractures (e.g. the coseismic subset), so now we
                    % use that subset to remove the extra fractures. In this
                    % first loop we are interested in having UTM coords and
                    % angle values, so we don't care about the mechanisms per
                    % se and we can use this also when not investigating the
                    % faults with known mechanims (mechanisms_data still
                    % contains the fractures with unknown mechanisms, they're
                    % labeled 'U')
                    for f=1:length(mechanisms_data)
                        targetIdx = find( [fractures.fid] == mechanisms_data(f).fid );
                        if ~isempty(targetIdx)
                            fractures(targetIdx).UTMx = mechanisms_data(f).UTMx;
                            fractures(targetIdx).UTMy = mechanisms_data(f).UTMy;
                            fractures(targetIdx).AngleUTM = mechanisms_data(f).AngleUTM;
                        end
                    end

                    % add mechanism (left-lateral, right-lateral, dip-slip) if asked for
                    if flag_mechanisms
                        % mechanisms_data has been filtered to contain only the
                        % desired fractures (e.g. the coseismic subset), so now
                        % we use that subset to assign the slip mechanism to
                        % each fractures in the dataset with theta_I and remove
                        % the extra fractures
                        for f=1:length(mechanisms_data)
                            targetIdx = find( [fractures.fid] == mechanisms_data(f).fid );
                            if ~isempty(targetIdx)
                                fractures(targetIdx).mechanism = mechanisms_data(f).mechanism;
                            end
                        end
                        % removing fractures with unknown/no slip mechanism
                        for f=length(fractures):-1:1
                            if isempty(fractures(f).mechanism) || strcmp(fractures(f).mechanism,"U")
                                fractures(f) = [];
                            end
                        end
                    end
                    % removing extra fractures
                    for f=length(fractures):-1:1
                        if isempty(fractures(f).AngleUTM)
                            fractures(f) = [];
                        end
                    end

                    % Loop #4 over mu values
                    %disp("Analyzing case #" + num2str(cases_counter) + " of " + num2str(nb_cases) );
                    cases_counter = cases_counter + 1;
                    
                    for m_b=1:length(mu_beta)
                        
                        if flag_mu
                            mu = mu_beta(m_b);
                            beta = pi/4 - atan(mu_beta(m_b))/2;
                        else
                            beta = mu_beta(m_b);
                            mu = tan(pi/2-2*beta);
                        end

                        % compute error
                        [err_abs, err_rel, err_w] = angles_error_shear_degrees_NW(fractures,beta,randangle,flag_mechanisms);
                        
                        if strcmp(SSR_type,'errAbs')
                            err_SSR = sum(abs(err_abs).^SSR_exp);
                        elseif strcmp(SSR_type,'errRel')
                            err_SSR = sum(abs(err_rel).^SSR_exp);
                        elseif strcmp(SSR_type,'errW')
                            err_SSR = sum(abs(err_w).^SSR_exp);
                        else
                            error('Wrong error type for SSR compute.');
                        end

                        SSR_matrix(m_b,mf) = err_SSR;
                        
                    end

                    % compute error assuming for each fracture a random
                    % value of mu between 0. and 1.5
                    nb_segments = 0;
                    for s=1:length(fractures)
                        nb_segments = nb_segments + length(fractures(s).AngleUTM);
                    end
                    if flag_mu
                        mu_rand = rand(nb_segments,1)*1.5; %interval [0,1.5]
                        beta_rand = pi/4 - atan(mu_rand)/2;
                    else
                        beta_rand = deg2rad( 22.5 + rand(nb_segments,1)*22.5 ); %beta randomly between 45 and 22.5deg (mu=0 and mu=1);
                        mu_rand = tan(pi/2-2*beta);
                    end
                    
                    [rand_err_abs, rand_err_rel, rand_err_w] = angles_error_shear_degrees_NW(fractures,beta_rand,randangle,flag_mechanisms);
                    
                    if strcmp(SSR_type,'errAbs')
                        rand_err_SSR = sum(abs(rand_err_abs).^SSR_exp);
                    elseif strcmp(SSR_type,'errRel')
                        rand_err_SSR = sum(abs(rand_err_rel).^SSR_exp);
                    elseif strcmp(SSR_type,'errW')
                        rand_err_SSR = sum(abs(rand_err_w).^SSR_exp);
                    else
                        error('Wrong error type for SSR compute.');
                    end
                    
                    SSR_matrix(m_b+1,mf) = rand_err_SSR;
                    
                end

                % We have done all loop #4 (over mu) and we are the end of the
                % current iteration over loop #3 (over depths). We now save to
                % disk a whole heatmap (SSR for given depth and mu) 
                disp("Saving to file " + csv_filename);
                writematrix(SSR_matrix,csv_filename);
                
                %save(mat_filename,'SSR_matrix','stress_angles','mu','magn_factors','selection_params', 'SSR_type','beta');
            
            end
        end
    end
end
