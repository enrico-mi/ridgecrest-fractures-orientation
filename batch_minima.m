%% Batch compute of best fits

addpath('./computes');

rinit = 1;

randangles = [0,5,10,20];

slip_models = ["JB2";"RB2";"XB2";"JB2";"JB2"];
sm_max = length(slip_models);
dist_ext = "";
res_ext = "";
mech_flags=[true,false];
smallerk = [true, false];

for a=1:length(randangles)
    randangle = randangles(a);
    if randangle == 0
        rmax = 1;
    else
        rmax = 100;
    end
    for f=1:2
        flag_mechanisms = mech_flags(f);
        if (randangle ~= 10 && ~flag_mechanisms)
            continue
        end

        string_ang = num2str(randangle,'%02.f');
        for sm=1:sm_max
            if (sm > 1 && ~flag_mechanisms)
                continue
            end

            extFile = '';
            if sm == 4
                extFile = "_10GPa";
            elseif sm == 5
                extFile = "_20GPa";
            end            
            minima = zeros(rmax,4);
            
            for sk = 1:length(smallerk)
                flag_smallerk = smallerk(sk);
                if flag_smallerk
                    if (sm > 1) || (~flag_mechanisms) || (randangle > 10)
                        continue
                    end
                    smallerk_ext = "_lowerk";
                else
                    smallerk_ext = "";
                end
                
                parfor r=rinit:rmax
                    disp("Minima " + slip_models(sm) + " - Realization " + num2str(r));
                    minima(r,:) = minSRw(r,slip_models(sm),randangle,flag_mechanisms,res_ext,extFile,dist_ext,smallerk_ext);
                end
                if flag_mechanisms
                    save("./minima_" + slip_models(sm) + "_" + string_ang + extFile + smallerk_ext, "minima");
                else
                    save("./minima_all_mechs" + slip_models(sm) + "_" + string_ang + extFile + smallerk_ext, "minima");
                end
            end
        end
    end
end
clearvars dist_ext;

%%%% for faults by length

randangles = [0,10];
resolution = [10];

slip_model = "JB2";
shear_ext = "";
res_exts = ["0010","Neg0010"];
dist_ext = "";
mech_flags=[true,false];

for a=1:length(randangles)

    randangle = randangles(a);    
    if randangle == 0
        rmax = 1;
    else
        rmax = 100;
    end

    for f=1:1
        flag_mechanisms = mech_flags(f);
        if (randangle ~= 10 && ~flag_mechanisms)
            continue
        end

        for res = 2:2
            res_ext = res_exts(res);


            string_ang = num2str(randangle,'%02.f');
            extFile = '';
            minima = zeros(rmax,4);
            for r=rinit:rmax
                disp("Minima " + slip_model + " - Distance from fault: " + dist_ext + " - Realization " + num2str(r));
                minima(r,:) = minSRw(r,slip_model,randangle,flag_mechanisms,res_ext,shear_ext,dist_ext);
            end
            if flag_mechanisms
                save("./minima_" + slip_model + "_" + string_ang + "_" + res_ext + "m","minima");
            else
                save("./minima_all_mechs_red_" + slip_model + "_" + string_ang + "_" + res_ext + "m","minima");
            end
        end
    end
end


%%%% for distance from faults

randangles = [0,10];
distances_from_fault = [1000,2000,3000];

slip_model = "JB2";
shear_ext = "";
res_ext = "";

for l=1:length(distances_from_fault)

    d_from_fault = distances_from_fault(l);
    dist_ext = "_dist" + num2str(d_from_fault,'%04.f');

    for a=1:length(randangles)

        randangle = randangles(a);    
        if randangle == 0
            rmax = 1;
        else
            rmax = 100;
        end

        string_ang = num2str(randangle,'%02.f');
        extFile = '';
        minima = zeros(rmax,4);
        parfor r=rinit:rmax
            disp("Minima " + slip_model + " - Distance from fault: " + dist_ext + " - Realization " + num2str(r));
            minima(r,:) = minSRw(r,slip_model,randangle,res_ext,shear_ext,dist_ext);
        end
        save("./minima_" + slip_model + "_" + string_ang + dist_ext,"minima");
    end
end