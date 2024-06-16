function sigma = rotate_frame_2D(sigma_xx, sigma_yy, sigma_xy, theta_xy)
    
% Rotate frame of reference by theta_xy, stresses do not rotate.
% Normally used to visualize stresses in a frame of reference whose x-axis is
% aligned with the fault strike
    
    Q = [  cos(theta_xy)  sin(theta_xy);...
          -sin(theta_xy)  cos(theta_xy)];
    
    for p = 1:size(sigma_xx,1)
        for q = 1:size(sigma_xx,2)
            sigma_tensor = [sigma_xx(p,q) sigma_xy(p,q);...
                            sigma_xy(p,q) sigma_yy(p,q)];

            sigma_rotated = Q*sigma_tensor*Q';
            
            sigma_xx(p,q) = sigma_rotated(1,1);
            sigma_yy(p,q) = sigma_rotated(2,2);
            sigma_xy(p,q) = sigma_rotated(1,2);
        end
    end
    
    sigma.sigma_xx = sigma_xx;
    sigma.sigma_yy = sigma_yy;
    sigma.sigma_xy = sigma_xy;
    
end