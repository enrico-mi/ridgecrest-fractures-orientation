function run_calc_stresses_NW(stress_angle,shearFlag)

% stress_angle: regional stress orientation in degrees with respect from North,
%               positive if clockwise (e.g. stress_angle=5 corresponds to N5E,
%               -10 to N10W)
% shearFlag: if true, parametrical study includes different shear moduli

    if (~exist('stress_angle','var'))
        stress_angle = 0.;
        string_stress_angle = num2str(stress_angle,'%02.f');
    else
        string_stress_angle = num2str(stress_angle,'%02.f');
        stress_angle = deg2rad(-stress_angle);
    end

    close all;
    clearvars -except stress_angle string_stress_angle shearFlag;
    addpath('../computes');
    addpath('../plots');
    addpath('./tools');

    %% Background stress
    % three quantities are passed:
    % sigma_bg.sigma_xx: stress along x
    % sigma_bg.sigma_yy: stress along y
    % sigma_bg.sigma_xy: deviatoric stress

    % Lythostatic stress as in Fialko (2021), varying depth

    rho_water = 1.0e3; %kg/m3
    rho_rock  = 2.7e3;
    gravity   = 9.8; %m/s2
                     % note: kg/m3*m/s2*m = Pa, consistent with Okada code output
    depth = 100; %m

    sigma_H = (rho_rock - rho_water)*gravity*depth;
    sigma_h = sigma_H/3;

    sigma_xx = sigma_h;
    sigma_yy = sigma_H;
    sigma_xy =       0.;
    sigma_zz = sigma_H;

    sigma_bg = rotate_stress_2D(sigma_xx, sigma_yy, sigma_xy, stress_angle);
    sigma_bg.sigma_zz = sigma_zz;

    %% Parametrical study

    magn_factors = 0.:0.1:0.1;
    resolutions = [0,10];

    % Strike = 261.75 ;  Dip = 82.52 ;  Rake = -179.64 ;
    theta_xy = deg2rad(-(136.89-90)); %136.89 is measured from N
                                      %warning('During visualization phase (e.g. plots for sigma_ij), the frame of reference is rotated by Ridgecrest global strike angle (136.89 deg from N, i.e. ~45 deg clockwise with respect to x-axis).');

    [~, ponti_mainshock_NW, ~] = read_fractures_shp([-117.74 -117.63],[35.82 35.94],'all');

    % selection needed only for resolutions, but the other "constraints" should
    % be applied in the angles analysis phase
    selection.angles.boolean = false;
    selection.verified.boolean = false;
    selection.coseismic.boolean = false;
    selection.resolution.boolean = true;
    selection.resolution.negative = false;
    selection.method.boolean = false;
    selection.resolution.segments = false;

    string_depth = num2str(depth,'%04.f');

    slip_models = ["JB2";"RB2";"XB2"];
    
    files_saved = [];
    counter = 1;
    
    if shearFlag
        sm_max = 5;
        slip_models = ["JB2";"RB2";"XB2";"JB2";"JB2"];
    else
        sm_max = 3;
    end
    
    smallerk = false;
    if ( mean(magn_factors) > 0 && mean(magn_factors) < 0.1 )
        smallerk = true;
    end
    
    for sm = 1:sm_max
        
        % sm will be the slip model later;
        slip_model_name = slip_models(sm);

        if sm > 1 && smallerk
            continue
        end
        
        %% Mechanical properties
        % Some of them are specified/derived in compute_calc_stresses_shp.m

        if sm == 1
            string_foreshock = "01JINx";
            string_mainshock = "02JINx";
            shearModulus = 33e9; %Pa
        elseif sm == 2
            string_foreshock = "01ROSS";
            string_mainshock = "02ROSS";
            shearModulus = 36e9; %Pa
        elseif sm == 3
            string_foreshock = "03XUxx";
            string_mainshock = "04XUxx";
            shearModulus = 10e9; %Pa
        elseif sm == 4
            string_foreshock = "01JINx";
            string_mainshock = "02JINx";
            shearModulus = 10e9; %Pa
            extFile = "_10GPa";
        elseif sm == 5
            string_foreshock = "01JINx";
            string_mainshock = "02JINx";
            shearModulus = 20e9; %Pa
            extFile = "_20GPa";
        end            
        
        
        for sel = 1:2
            
            if sel==1
                selection.resolution.negative = false;
            else
                selection.resolution.negative = true;
            end
            
            if sel > 1 && smallerk
                continue
            end
            
            for res = 1:length(resolutions)
                
                if resolutions(res) > 0 && (sm > 1 || smallerk)
                    continue
                end
                
                selection.resolution.ref_length = resolutions(res);
                selected_fractures = select_fractures(ponti_mainshock_NW,selection,true,false);
                string_res = num2str(resolutions(res),'%04.f');
                
                for mf = 1:length(magn_factors)
                    magn_factor = magn_factors(mf); % this is called k in the manuscript
                    if magn_factor > 0 && magn_factor < 0.1
                        string_mf = num2str(magn_factor*100,'%03.f');
                    else
                        string_mf = num2str(magn_factor*10,'%02.f');
                    end
                    
                    fig_folder =  "./output/" + slip_model_name + "_N" + string_stress_angle + "E_d" + string_depth + "_res" + string_res + "_mf" + string_mf + "_NW";
                    if selection.resolution.negative
                        fig_folder =  "./output/" + slip_model_name + "_N" + string_stress_angle + "E_d" + string_depth + "_resNeg" + string_res + "_mf" + string_mf + "_NW";
                    end
                    
                    if sm > 3
                        fig_folder = fig_folder + extFile;
                    end
                    
                    if ~(exist(fig_folder, 'dir') == 7)
                        mkdir(fig_folder);

                        filename_output_mains  = char("./coseismic-input/output2019RIDGEC" + string_mainshock + "_shp_NW_"+string_depth+"m_res"+string_res+"m");
                        filename_output_fores  = char("./coseismic-input/output2019RIDGEC" + string_foreshock + "_shp_NW_"+string_depth+"m_res"+string_res+"m");
                        if selection.resolution.negative
                            filename_output_mains  = char("./coseismic-input/output2019RIDGEC" + string_mainshock + "_shp_NW_"+string_depth+"m_resNeg"+string_res+"m");
                            filename_output_fores  = char("./coseismic-input/output2019RIDGEC" + string_foreshock + "_shp_NW_"+string_depth+"m_resNeg"+string_res+"m");
                        end
                        
                        if sm > 3
                            filename_output_mains = filename_output_mains + extFile;
                            filename_output_fores = filename_output_fores + extFile;
                        end
                        
                        disp('data and figures folder: ' + string(fig_folder));
                        disp('mainshock stress: ' + string(filename_output_mains));
                        disp('foreshock stress: ' + string(filename_output_fores));

                        [~, prefix, ~] = fileparts(fig_folder);
                        fractures_data_name = string(fig_folder) + '/' + string(prefix) + '_fractures_data.mat';

                        fractures_data{1} = prefix; % coincides with id of analysis
                        fractures_data{2} = slip_model_name;
                        fractures_data{3} = string(string_stress_angle);
                        fractures_data{4} = magn_factor;
                        fractures_data{5} = string(string_res);
                        fractures_data{6} = selection; % preserve selection parameters to filter fractures
                        if magn_factor > 0
                            magn_factor_fore = 1.;
                        else
                            magn_factor_fore = 0.;
                        end
                        fractures_data{7} = compute_calc_stresses_shp(filename_output_mains, fig_folder, sigma_bg, shearModulus, magn_factor, theta_xy, filename_output_fores, magn_factor_fore, selected_fractures);
                        

                        save(fractures_data_name,"fractures_data");
                        clearvars fractures_data;
                        
                        files_saved{counter} = fractures_data_name;
                        counter = counter + 1;
                        
                        disp('Done!');
                    else
                        disp('Folder already exists, case not computed.');
                    end
                end
            end
        end
    end
    
    disp("Files produced:");
    for f = 1:length(files_saved)
        disp(files_saved{f});
    end
end