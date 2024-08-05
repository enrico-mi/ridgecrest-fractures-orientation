function [selected_fractures, negative_selection] = lengthscale_fractures(fractures, ref_length)
% select and re-discretize fractures at resolution ref_length and above

% concept: selecting fractures discretized at a scale larger than ref_length;
% if the overall fracture is smaller than ref_length, it is discarded; if
% re-discretization of the fracture would lead to shorter segments
% (e.g. refined), it is kept as is (we want to avoid adding information that is
% not there); if fracture is discretized at scale smaller than provided
% ref_length, it is re-descritized at coarser scale (to look at information
% only at that scale and above).

% negative_selection is the complementary selection, i.e. for a length scale L, it contains all fractures smaller than L. No re-discretization is performed.

    selected_fractures = fractures;
    negative_selection = fractures;
    if ref_length > 1e-10
        for f=length(fractures):-1:1
            frac = selected_fractures(f);
            UTMx = frac.UTMx;
            UTMy = frac.UTMy;
            frac_length = sqrt( (UTMx(end) - UTMx(1))^2 + (UTMy(end) - UTMy(1))^2);
            nb_seg = floor(frac_length/ref_length);
            if nb_seg == 0
                % discarding fractures shorter than provided resolution
                selected_fractures(f) = [];
            elseif nb_seg == 1
                negative_selection(f) = [];
                % for selected_fractures, skip and keep fracture as is
            elseif nb_seg > length(UTMx)-1
                % preserving fractures coarser than provided resolution (L fracture segments > L resolution)
                selected_fractures(f) = frac;
                negative_selection(f) = [];
            else
                negative_selection(f) = [];
                % re-discretizing fractures with segment smaller than resolution (L fracture segments < L resolution)
                UTMxi = linspace(UTMx(1), UTMx(end), nb_seg+1);
                UTMyi = linspace(UTMy(1), UTMy(end), nb_seg+1);
                
                % distance looked for for intermediate nodes where we're discretizing
                UTMx_l = zeros(length(UTMxi)-2,1);
                UTMy_l = zeros(length(UTMxi)-2,1);
                UTMxi = UTMxi(2:end-1);
                UTMyi = UTMyi(2:end-1);

                for j=1:length(UTMxi)
                    % we rediscretize by finding the points along the original
                    % fracture's segments whose projection on the line linking the
                    % first and last point of the fracture coincides with the new
                    % discretization points of such linking line
                    
                    % linking-line slope and intercept
                    m1 = (UTMy(end) - UTMy(1))/(UTMx(end) - UTMx(1));
                    q1 = UTMy(1) - m1*UTMx(1);
                    
                    % slope and intercept of line orthogonal to linking line
                    % through UTMxi(j)
                    m1_p = -1/m1;
                    q1_p = UTMyi(j) - m1_p*UTMxi(j);
                    
                    % projection of each original point to linking line
                    q_orig = UTMy - m1_p*UTMx;
                    UTMx_perp = (q1 - q_orig)/(m1_p - m1);
                    
                    % for each new point UTM_xi(j), find segment of reference to find
                    % point whose projection lies on linking line
                    dist_pos = UTMx_perp;
                    dist_pos(UTMx_perp > UTMxi(j)) = Inf;
                    [~,ind_pos] = min(abs(dist_pos - UTMxi(j)));
                    dist_neg = UTMx_perp;
                    dist_neg(UTMx_perp < UTMxi(j)) = Inf;
                    [~,ind_neg] = min((dist_neg - UTMxi(j)));
                    ind_segm = sort([ind_pos, ind_neg]);
                    
                    % find slope and intercept of line where segment lies
                    m2 = (UTMy(ind_segm(2))-UTMy(ind_segm(1)))/(UTMx(ind_segm(2))-UTMx(ind_segm(1)));
                    q2 = UTMy(ind_segm(1)) - m2*UTMx(ind_segm(1));
                    
                    % find x and y of new interpolated point
                    x_target = (q1_p - q2)/(m2 - m1_p);
                    y_target = m2*x_target + q2;
                    
                    delta_lines = y_target - (m1_p*x_target + q1_p);
                    if delta_lines > 1e-3 %m
                        warning('Orthogonal to linking-line and segment line predict two different points ' + string(delta_lines) + 'm apart.');
                    end
                    
                    UTMx_l(j) = x_target;
                    UTMy_l(j) = y_target;

                end
                frac.UTMx = [UTMx(1), UTMx_l', UTMx(end)];
                frac.UTMy = [UTMy(1), UTMy_l', UTMy(end)];

                %updating angles
                frac_angle_utm = atan( diff(frac.UTMy) ./ diff(frac.UTMx) );
                frac_angle_ave_utm = atan( (frac.UTMy(end) - frac.UTMy(1))/(frac.UTMx(end) - frac.UTMx(1)) );
                sin_arg = sum(diff(frac.UTMy));
                cos_arg = sum(diff(frac.UTMx));
                angle_ave_utm = atan(sin_arg/cos_arg);
                frac.AngleUTM = frac_angle_utm;
                frac.AngleAveUTM = angle_ave_utm;
                
                if (angle_ave_utm - frac_angle_ave_utm) > 1e-15
                    error("Different average computations give different results!");
                end
                
                zone_utm = utmzone(35.8,-117.7);
                disp("Using UTM zone" + string(zone_utm));
                [ellipsoid,estr] = utmgeoid(zone_utm);
                utmstruct = defaultm('utm');
                utmstruct.zone = zone_utm;
                utmstruct.geoid = ellipsoid;
                utmstruct = defaultm(utmstruct);
                [lat, lon] = projinv(utmstruct,frac.UTMx,frac.UTMy);
                frac.Lat = lat;
                frac.Lon = lon;
                
            end
        end
    else
        warning("ref_length = 0, no selection or rediscretization applied.");
    end

end