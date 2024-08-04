function sigma = rotate_stress_2D(sigma_xx, sigma_yy, sigma_xy, theta_xy)
    
% Rotate stresses by theta_xy, frame of reference is fixed. 
% Normally used to apply regional stresses at different orientation in frame of
% reference whose x-axis is aligned with EW direction
    
    R = [  cos(theta_xy) -sin(theta_xy);...
           sin(theta_xy)  cos(theta_xy)];
    
    for p = 1:size(sigma_xx,1)
        for q = 1:size(sigma_xx,2)
            sigma_tensor = [sigma_xx(p,q) sigma_xy(p,q);...
                            sigma_xy(p,q) sigma_yy(p,q)];

            sigma_rotated = R*sigma_tensor*R';
            
            sigma_xx(p,q) = sigma_rotated(1,1);
            sigma_yy(p,q) = sigma_rotated(2,2);
            sigma_xy(p,q) = sigma_rotated(1,2);
        end
    end
    
    sigma.sigma_xx = sigma_xx;
    sigma.sigma_yy = sigma_yy;
    sigma.sigma_xy = sigma_xy;
    
end