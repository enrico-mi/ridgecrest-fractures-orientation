function minimum = minSRw(realization, slip_model_name, randangle, flag_mechanisms, res_ext, shear_ext, dist_ext, smallerk_ext)
    
    if randangle > 10
        sigmaH_angles = -15:1:40;
    elseif strcmp(res_ext,'Neg0010')
        sigmaH_angles = 12:1:16;
    else
        sigmaH_angles = -10:1:25;
    end

    if randangle > 10 && ~flag_mechanisms
        error("no all mechanisms computation for uncertainty angles larger than 10 degrees")
    end
    
    if ~flag_mechanisms
        sigmaH_angles = 0:1:15;
    end
    
    if strcmp(smallerk_ext,"_lowerk")
        magn_factors = 0.:0.01:0.1;
    else 
        magn_factors = 0.:0.1:1.;
    end

    mus = {};
    for j=1:16
        mus{j} = (j-1)/10;%*0.3;
    end
    mus{17} = 1e6;
    mus{18} = "rand";
    mu_idx = 1:18;
    
    SRw_min = 1e10*ones(length(magn_factors),1);
    thetaH_min = -30*ones(length(magn_factors),1); %random value outside range to
                                                 %check that a value is always
                                                 %assigned
    mu_min = {};
    
    for mf = 1:length(magn_factors)
        for sa = 1:length(sigmaH_angles)
            string_stress_angles = num2str(sigmaH_angles(sa),'%02.f');
            SRw_tmp = getSRw(string_stress_angles,mu_idx,mf,slip_model_name,realization,randangle,flag_mechanisms,res_ext,shear_ext,dist_ext,smallerk_ext);
            [SRw_min_mu, idx] = min(SRw_tmp);
            if SRw_min_mu < SRw_min(mf)
                SRw_min(mf) = SRw_min_mu;
                mu_min{mf} = mus{idx};
                thetaH_min(mf) = sigmaH_angles(sa);
            end
        end
    end
    
    [SRw_min_abs1, idx1] = min(SRw_min);
% $$$     SRw_min2 = SRw_min;
% $$$     SRw_min2(idx1) = [];
% $$$     [SRw_min_abs2, idx2] = min(SRw_min2);
% $$$     if idx2 >= idx1
% $$$         idx2 = idx2 + 1;
% $$$     end
    
    disp("Minimum1 SRw=" + num2str(SRw_min_abs1));
    disp("for        k=" + num2str(magn_factors(idx1)));
    disp("at        mu=" + num2str(mu_min{idx1}));
    disp("and   thetaH=N" + num2str(thetaH_min(idx1)) + "E");
            
% $$$     disp("Minimum2 SRw=" + num2str(SRw_min_abs2));
% $$$     disp("for        k=" + num2str(magn_factors(idx2)));
% $$$     disp("at        mu=" + num2str(mu_min{idx2}));
% $$$     disp("and   thetaH=N" + num2str(thetaH_min(idx2)) + "E");
    
    minimum = [SRw_min_abs1, magn_factors(idx1), mu_min{idx1}, thetaH_min(idx1)];
% $$$     minimum2 = [SRw_min_abs2, magn_factors(idx2), mu_min{idx2}, thetaH_min(idx2)];
            
end
