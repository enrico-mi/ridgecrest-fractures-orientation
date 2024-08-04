%% Batch compute of error matrices

% runs 100 random realizations to apply random error to original angles
addpath( fullfile('.','computes') );
addpath( fullfile('.','fractures_mgmt') );

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
mechanisms_data = select_fractures(mechanisms_data,selection,false);

clearvars ponti_mechanisms_data;

rinit = 1;
randangles = [0,5,10,20];
smallerk = [true, false];

for a=1:length(randangles)
    randangle = randangles(a);
    if randangle > 0
        rmax = 100;
    else
        rmax = 1;
    end
    for sk = 1:length(smallerk)
        if smaller(sk)
            if ~(randangle==0 || randangle==10)
                continue
            end
        end
        parfor r=rinit:rmax
            mu_matrices_sigma_I_NW(mechanisms_data,r,randangle,false,'',smallerk(sk));
            disp("Realization " + num2str(r));
        end
    end
end
