function selected_fractures = compute_calc_stresses_shp(filename_output, fig_folder, sigma_bg, shearModulus, magn_factor_main, theta_xy, filename_output2, magn_factor_fore, selected_fractures)
    
    %% NB: if stresses compute in filename_output are not on a regular grid, not all code capacities can be used (mostly plotting)
    
% filename_output: output file from calc_stress 
% extension is '.pscmp_', then Plot_faulttrace script is used
% fig_folder: figures destination folder (the path must already exist)
% capped (optional): max absolute value to cap stresses to (for better plot
% visualization)
% filename_output2 (optional): second event to include (it only adds stress,
% not geometry of the fault); grid has to be the same as the first file

    % plot figures in the background
    set(groot,'defaultFigureVisible','off')

    %% Mechanical properties as in Xu et al (2020, Science)

    G = shearModulus; %in Pa (e.g 10*1e9 for 10 GPa)
    
    % Further assumed parameter by me to convert stress to strain
    nu = 0.25;

    % The following follow from previous assumptions
    Gprime = G/(1-nu); %for plane strain, I assume sigma_zz -> 0 at shallow depths
    E = 2*G*(1+nu);

    % these are not used in this function, but can come handy for debugging
    % when saving the workspace
    [~, prefix, ~] = fileparts(fig_folder);
    bname = string(fig_folder) + '/' + string(prefix) + '_';
    workspace_name = string(fig_folder) + '/' + string(prefix) + '_workspace.mat';

    %% I/O
    delimiterIn = ' ';
    headerlinesIn = 1;
    stress_data = importdata(filename_output,delimiterIn,headerlinesIn);
    for k=4:9
        stress_data.data(:,k) = stress_data.data(:,k)*magn_factor_main;
    end
    
    if (exist('filename_output2','var'))
        stress_data2 = importdata(filename_output2,delimiterIn,headerlinesIn);
        check1 = isequal(stress_data.data(:,1),stress_data2.data(:,1));
        check2 = isequal(stress_data.data(:,2),stress_data2.data(:,2));
        check3 = isequal(stress_data.data(:,3),stress_data2.data(:,3));
        if ~check1 && ~check2 && ~check3
            error('Input stress grids do not have same x,y,z coordinates grid');
        end
        if (~exist('magn_factor_fore','var'))
            magn_factor_fore = 1.;
        end
        for k=4:9
            stress_data.data(:,k) = stress_data.data(:,k) + stress_data2.data(:,k)*magn_factor_fore;
        end
    end
    
    %% Sigma plots - % test plots with scatter from raw, calc_stress data
% $$$     scatter_stress(stress_data,dims(1),dims(2),s_cols(1),bname+'sigma'+sfx(1),'sigma'+sfx(1));
% $$$     scatter_stress(stress_data,dims(1),dims(2),s_cols(2),bname+'sigma'+sfx(2),'sigma'+sfx(2));
% $$$     scatter_stress(stress_data,dims(1),dims(2),s_cols(3),bname+'sigma'+sfx(3),'sigma'+sfx(3));
% $$$     scatter_stress(stress_data,dims(1),dims(2),s_cols(4),bname+'sigma'+sfx(4),'sigma'+sfx(4));
% $$$     scatter_stress(stress_data,dims(1),dims(2),s_cols(5),bname+'sigma'+sfx(5),'sigma'+sfx(5));
% $$$     scatter_stress(stress_data,dims(1),dims(2),s_cols(6),bname+'sigma'+sfx(6),'sigma'+sfx(6));

    % NB: Okada uses x for NS and y for EW directions; he also uses
    % tension-positive for normal stresses; not sure yet about shear stresses
    % at the moment if I don't change signs for sigma_xy, I get a negative stress drop for a left-dim2eral fault, which should be consistent with a positive left-dim2eral imposed shear stress outside the crack (the crack goes agains the positive imposed shear)

    sigma_xx = -stress_data.data(:,5); % swapping x-y (EW-NS)
    sigma_yy = -stress_data.data(:,4); % swapping x-y (EW-NS)
    sigma_zz = -stress_data.data(:,6);
    sigma_xy = -stress_data.data(:,7);
    sigma_yz =  stress_data.data(:,9); % swapping x-y (EW-NS)
    sigma_xz =  stress_data.data(:,8); % swapping x-y (EW-NS)

    sigma_xy = sigma_xy + sigma_bg.sigma_xy;
    sigma_xx = sigma_xx + sigma_bg.sigma_xx;
    sigma_yy = sigma_yy + sigma_bg.sigma_yy;
    sigma_zz = sigma_zz + sigma_bg.sigma_zz;

    % sigma_zz_2D = 0.25 * (sigma_xx + sigma_yy);

    % Headache inducing note on sigma_ij & co structure.
    
    % When using contourf(x,y,Z), Z is a matrix where the n-th row corresponds
    % to the n-th y-coordinate, and the m-th column to the m-th x-coordinate.
    
    % When reshaping sigma_ij, the first size (rows) has to be the nb of
    % components along the dimension that in the output file from calc-stress
    % increases first (e.g. it's Lat in outputTest1) and the second size (cols)
    % the nb of components along the other dimension (e.g. Lon in outputTest1).
    
    % We want to call contourf(x,y,Z) with x=dim1 and y=dim2, so that our dim1
    % (dim2) is the horizontal (vertical) axis in the plots. If dim1 does not
    % coincide with the second dimension to increase in the output file from
    % calc-stress, it means that the reshaped sigma_ij is such that its m-th
    % row corresponds to the m-th x-coordinate, and the n-th column to the n-th
    % y-coordinate, i.e. the opposite of what we need for Z in contourf, so
    % need to transpose sigma_ij.
    
    % Finally, when one of the two dimensions is depth (z), we want to flip all
    % the values along that direction, to plot the results in a frame of
    % reference where z is positive upward (NB: this is done in plot_stress by
    % calling plane_vars_plot).

    % working in compression-positive convention, then sigma_1 is the greatest
    % compressive stress, and theta_1pxy is the angle between local x-axis and
    % sigma_1
    flag_3D = true;
    [sigma_1, sigma_3, theta_p1, theta_p3, sigma_2, sigma_1pxy, sigma_2pxy, theta_1pxy] = principal_sigmas_eigen(sigma_xx, sigma_yy, sigma_xy, flag_3D, sigma_zz, sigma_xz, sigma_yz);

    flag_3D = false;
    [sigma_I, sigma_II, theta_I, theta_II] = principal_sigmas_eigen(sigma_xx, sigma_yy, sigma_xy, flag_3D);

    %% Time for plots
    
    % Rotating horizontal sigma_ij for visualization purposes
    if theta_xy > 0
        sense = ' counterclockwise';
    else
        sense = ' clockwise';
    end
    disp('sigma_ij are rotated by ' + string(rad2deg(abs(theta_xy))) + sense + ' to align with mainshock strike for visualization purposes');
     
    %Angles are computed with respect to fault strike angle, so I am rotating them by ' +string(rad2deg(-theta_xy)) + ' degrees to plot them as angles with respect to global x-axis.');
    rotated_sigmas = rotate_frame_2D(sigma_xx, sigma_yy, sigma_xy, theta_xy);

    lat = stress_data.data(:,1);
    lon = stress_data.data(:,2);

    theta_counter = 0;
    for f=1:length(selected_fractures)
        selected_fractures(f).theta_I = zeros(length(selected_fractures(f).AngleUTM),1);
        for segment = 1:length(selected_fractures(f).Lat)-2 %last entry is NaN by construction, no segment
            lonc = (selected_fractures(f).Lon(segment+1) + selected_fractures(f).Lon(segment))/2;
            latc = (selected_fractures(f).Lat(segment+1) + selected_fractures(f).Lat(segment))/2;
            theta_counter = theta_counter + 1;
            if abs(lonc-lon(theta_counter)) < 5e-4 && abs(latc-lat(theta_counter)) < 1e-5
                selected_fractures(f).theta_I(segment) = theta_I(theta_counter);
            else
                error('it does not look like the right coordinates for this theta_I point');
            end
        end
    end
    if ~(theta_counter==length(sigma_I))
        warning('we did not store all theta_I points!')
    end
    
    sigma_xy_rot = rotated_sigmas.sigma_xy;
    sigma_xx_rot = rotated_sigmas.sigma_xx;
    sigma_yy_rot = rotated_sigmas.sigma_yy;

    % reducing file size by removing extra fields
    fields = fieldnames(selected_fractures);
    for f = 1:length(fields)
        if ~strcmp(fields(f),"fid") && ~strcmp(fields(f),"theta_I")
            selected_fractures = rmfield(selected_fractures, fields(f));
        end
    end

end