function [err_abs_array, err_rel_array, err_scaled_array] = angles_error_shear_degrees_NW(fractures, beta, randangle, flag_mechanisms)

% NB on angle conversions.  Geophysics measures angles from N, positive if
% clockwise; in my previous computations, all angles are in classical
% trigonometric space, with angles measured from y=0 and positive if
% counter-clockwise. The latter is valid for both principal stresses angles and
% for fractures' angles (which are computed with atan), and they are all in the
% interval [-pi/2;pi/2]. In the following code I do two types of
% conversion. The first one is to place trigonometrical angles in the
% geographic once; the angles are then in the interval [0, pi/2], with 0 being
% N0E, pi/2 N90E, and pi S0E. The second conversion is to represent the angles
% symmetrically with respect to the angle of local largest stresses. (As a S10E
% angle is equivalent to N10W in these cases, because I take into account slip
% direction separetly when I have it.) The angles are then expressed in
% radiants in the domain [thetaH-pi/2;thetaH+pi/2]. Only at the end I convert
% radians to degrees.    
    
    nb_segments = 0;
    longest_segment = 0.;
    cumulative_length = 0.;
    for f=1:length(fractures)
        frac_segments = length(fractures(f).AngleUTM);
        seg_lengths = zeros(frac_segments,1);
        nb_segments = nb_segments + frac_segments;

        for k=1:frac_segments
            dx2 = (fractures(f).UTMx(k+1)-fractures(f).UTMx(k))^2;
            dy2 = (fractures(f).UTMy(k+1)-fractures(f).UTMy(k))^2;
            seg_lengths(k) = sqrt(dx2+dy2);
            if seg_lengths(k) > longest_segment
                longest_segment = seg_lengths(k);
            end
        end
        
        fractures(f).Lengths = seg_lengths;
        cumulative_length = cumulative_length + sum(seg_lengths);
    end
    
    if ~flag_mechanisms
        [err_abs_array, err_rel_array, err_scaled_array] = error_no_mechanisms(fractures, beta, cumulative_length, nb_segments, randangle);
    else
        [err_abs_array, err_rel_array, err_scaled_array] = error_mechanisms(fractures, beta, cumulative_length, nb_segments, randangle);
    end
end

function [err_abs_array, err_rel_array, err_scaled_array] = error_no_mechanisms(fractures, beta, cumulative_length, nb_segments, randangle)
    
    delta_min_array = zeros(nb_segments,1);
    err_abs_array = zeros(nb_segments,1);
    err_rel_array = zeros(nb_segments,1);
    err_scaled_array = zeros(nb_segments,1);

    segment_idx = 0;

    for f=1:length(fractures)
        
        for k=1:length(fractures(f).AngleUTM)
            
            segment_idx = segment_idx + 1;

            sim_sigma_1 = fractures(f).theta_I(k);

            if length(beta)==1
                sim_LL = sim_sigma_1 - beta;
                sim_RL = sim_sigma_1 + beta;
            else
                sim_LL = sim_sigma_1 - beta(segment_idx);
                sim_RL = sim_sigma_1 + beta(segment_idx);
            end
            
            frac_angle = fractures(f).AngleUTM(k);
            % perturb angle by +/- 5 degrees
            randnum = rand(1);
            randang = randnum*deg2rad(randangle); %re-scaling angle to 0-5 degrees
            randsign = sign(rand(1)-0.5); % rand(1) is between 0 and 1; not using randnum or sign and magnitude are correlated
            frac_angle = fractures(f).AngleUTM(k);
            frac_angle = frac_angle+randang*randsign;
            if frac_angle < -pi/2
                frac_angle = frac_angle + pi;
            elseif frac_angle > pi/2
                frac_angle = frac_angle - pi;
            end
            
            % structure: [angle_LL, angle_RL];
            sim_angles_shear = [sim_LL, sim_RL];
            
            % two measures of distances because angles are in -pi/2:pi/2
            % domain, but, e.g., if a fracture is at -pi/2, and theta_1 = pi/4,
            % theta_2 = pi/2, theta_2 is the segment with an orientation closer
            % to the fracture (they're actually aligned)
            distance_1 = abs(sim_angles_shear - frac_angle);
            distance_2 = abs( (sim_angles_shear - sign(sim_angles_shear)*pi) - frac_angle);
            [delta_min_1, ~] = min(distance_1);
            [delta_min_2, ~] = min(distance_2);
            if delta_min_1 < delta_min_2
                delta_min = delta_min_1;
            else
                delta_min = delta_min_2;
            end
            delta_min_array(segment_idx) = delta_min;
            
            % remember positions: 1=LL, 2=RL
            if length(beta)==1
                max_error = pi/2 - beta; % only one direction (LL)
            else
                max_error = pi/2 - beta(segment_idx);
            end
            delta_min = rad2deg(delta_min);
            max_error = rad2deg(max_error);
            err_rel = delta_min/max_error;
            err_weight = err_rel * (fractures(f).Lengths(k)/cumulative_length);
            err_abs_array(segment_idx) = delta_min;
            err_rel_array(segment_idx) = err_rel;
            err_scaled_array(segment_idx) = err_weight;
            
        end
        
    end
    
end

function [err_abs_array, err_rel_array, err_scaled_array] = error_mechanisms(fractures, beta, cumulative_length, nb_segments, randangle)

    delta_min_array = zeros(nb_segments,1);
    err_abs_array = zeros(nb_segments,1);
    err_rel_array = zeros(nb_segments,1);
    err_scaled_array = zeros(nb_segments,1);

    segment_idx = 0;
    for f=1:length(fractures)
        
        colour_plot_LL = zeros(length(fractures(f).AngleUTM),4);
        colour_plot_RL = zeros(length(fractures(f).AngleUTM),4);
        colour_plot_shearall = zeros(length(fractures(f).AngleUTM),4);
        colour_plot_shearall_grey = zeros(length(fractures(f).AngleUTM),3);
        
        if isempty(fractures(f).mechanism) || strcmp(fractures(f).mechanism,"U")
            continue
        elseif strcmp(fractures(f).mechanism,"RL")
            direction = +1;
        elseif strcmp(fractures(f).mechanism,"LL")
            direction = -1;
        elseif strcmp(fractures(f).mechanism,"V")
            direction =  0;
        end
        
        % perturb angle by +/- randangle
        % randnum = rand(1);
        % randang = randnum*deg2rad(randangle); %re-scaling angle to 0-5 degrees
        % randsign = sign(rand(1)-0.5); % rand(1) is between 0 and 1; not using randnum or sign and ma    
        
        for k=1:length(fractures(f).AngleUTM)
            
            segment_idx = segment_idx + 1;

            sim_sigma_1 = fractures(f).theta_I(k);

            if length(beta)==1
                sim_angle = sim_sigma_1 + direction*beta;
            else
                sim_angle = sim_sigma_1 + direction*beta(segment_idx);
            end
            
            if sim_angle > pi/2
                sim_angle = sim_angle - pi;
            elseif sim_angle < -pi/2
                sim_angle = sim_angle + pi;
            end
            
            frac_angle = fractures(f).AngleUTM(k);
            % perturb angle by +/- 5 degrees
            randnum = rand(1);
            randang = randnum*deg2rad(randangle); %re-scaling angle to 0-5 degrees
            randsign = sign(rand(1)-0.5); % rand(1) is between 0 and 1; not using randnum or sign and magnitude are correlated
            frac_angle = fractures(f).AngleUTM(k);
            frac_angle = frac_angle+randang*randsign;
            if frac_angle < -pi/2
                frac_angle = frac_angle + pi;
            elseif frac_angle > pi/2
                frac_angle = frac_angle - pi;
            end
            
            % two errors because angles are in [-pi/2,pi/2], but pi/4 is as
            % close to -pi/2 as to pi/2 and 0;

            if frac_angle > sim_angle
                error_one = frac_angle - sim_angle;
                error_two = pi - (frac_angle - sim_angle);
            else
                error_one = sim_angle - frac_angle;
                error_two = pi - (sim_angle - frac_angle);
            end
            
            if error_one < 0
                error('error_one < 0');
            elseif error_two < 0
                error('error_two < 0');
            end
            
            delta_min = min([error_one, error_two]);
            delta_min = rad2deg(delta_min);
            
            delta_min_array(segment_idx) = delta_min;

            max_error = rad2deg(pi/2);

            err_rel = delta_min/max_error;
            
            err_weight = delta_min * (fractures(f).Lengths(k)/cumulative_length);
            err_colour = 1 - err_weight;
            err_abs_array(segment_idx) = delta_min;
            err_rel_array(segment_idx) = err_rel;
            err_scaled_array(segment_idx) = err_weight;
        end
    end
end
