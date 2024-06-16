%% Batch compute of error matrices (selecting fractures further than 1 or 2km from main fault)

% runs 100 random realizations to apply random error to original angles
addpath('./computes');
addpath('./fractures_mgmt');

rinit = 1;

randangles = [0,10];

distances_from_fault = [1000,2000];

for l=1:length(distances_from_fault)
    
    d_from_fault = distances_from_fault(l);

    ponti_mechanisms_data = load("ponti_mainshock_NW_mechs.mat");
    mechanisms_data = ponti_mechanisms_data.ponti_mainshock_NW;
    mechanisms_data = compute_angles_shp(mechanisms_data);
    % restrain to coseismic fractures (the other selection parameters are needed for compatibility)
    selection.coseismic.boolean = true;
    selection.angles.boolean = false;
    selection.verified.boolean = false;
    selection.resolution.boolean = false;
    selection.resolution.negative = false;
    selection.method.boolean = false;
    selection.resolution.segments = false;
    mechanisms_data = select_fractures(mechanisms_data,selection,false,false);
    mechanisms_data = select_by_faultdistance(mechanisms_data,d_from_fault);

    clearvars ponti_mechanisms_data;

    dist_ext = "_dist" + num2str(d_from_fault,'%04.f');

    for a=1:length(randangles)
        randangle = randangles(a);
        
        if randangle == 0
            rmax = 1;
        else
            rmax = 100;
        end

        parfor r=rinit:rmax
            mu_matrices_sigma_I_NW(mechanisms_data,r,randangle,true,dist_ext,false);
            disp("Realization " + num2str(r));
        end
    end
end