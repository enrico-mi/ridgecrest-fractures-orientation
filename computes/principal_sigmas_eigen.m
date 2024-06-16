function [sigma_1, sigma_3, theta_1, theta_3, varargout] = principal_sigmas_eigen(sigma_xx, sigma_yy, sigma_xy, flag_3D, sigma_zz, sigma_xz, sigma_yz)
    
    % NB on the sign convention.

    % If you are working in *compression-positive* convention, then this script
    % returns:
    % sigma_1, defined as the greatest *compressive* principal stress
    % sigma_3, defined as the greatest *tensile* principal stress
        
    % If you are working in *tension-positive* convention, then this script
    % returns:
    % sigma_1, defined as the greatest *tensile* principal stress
    % sigma_3, defined as the greatest *compressive* principal stress

    % The angles are always defined as follows.
    
    % theta_1 is the angle that measures by how much you need to rotate your
    % infinitesimal element so that its x-axis is aligned with
    % sigma_1. theta_1 > 0 if counterclockwise.
    
    % theta_3 is the angle that measures by how much you need to rotate your
    % infinitesimal element so that its x-axis is aligned with
    % sigma_3. theta_3 > 0 if counterclockwise.
        
    % I adopted sigma_1, sigma_3 in 2D too to easily switch between 2D and 3D.
        
    if (~exist('flag_3D','var'))
        flag_3D = false;
    end
        
    % NB on 3D:
    
    % varargout contains sigma_2, theta_2, where sigma_3 < sigma_2 < sigma_1
    % and sigma_1 (sigma_3) are always the maximum (minimum) principal
    % stresses. This is a bit uncomfortable and definitely not elegant, but
    % it's staying like this for the moment, as the 3D is a patch (the original
    % function was written for plane strain).
    
    % Angles kind of lose meaning: they are all with respect to the <xy> plane
    % as that coincides with the free surface for my Ridgecrest study and it is
    % the plane the fault is (roughly) perpendicular to. It should be a good
    % enough approximation close to the surface, but be aware that the
    % direction is actually 3D (in particular at deeper depths).
        
    switch flag_3D
        
      case false
        %% 2D case (default for backward compatibility)

        if (exist('sigma_zz','var')) || (exist('sigma_xz','var')) || (exist('sigma_yz','var'))
        warning(['2D case is called for, but 3D stress tensor is passed: note that z components will be discarded.)']);
        end
        
        sigma_1 = zeros(size(sigma_xx));
        sigma_3 = zeros(size(sigma_xx));
        theta_1 = zeros(size(sigma_xx));
        theta_3 = zeros(size(sigma_xx));
        
        xsize = size(sigma_xx,1);
        ysize = size(sigma_xx,2);
        
        for j=1:xsize
            for k=1:ysize
                
                sigma_tensor = [sigma_xx(j,k) sigma_xy(j,k); sigma_xy(j,k) sigma_yy(j,k)];
                
                [V,D] = eig(sigma_tensor);
                [d,ind] = sort(diag(D),'descend');
                Ds = D(ind,ind);
                Vs = V(:,ind);
                
                sigma_1(j,k) = Ds(1,1);
                sigma_3(j,k) = Ds(2,2);
                n_1 = Vs(:,1);
                n_3 = Vs(:,2);
                
                theta_1(j,k) = atan(n_1(2)/n_1(1));
                theta_3(j,k) = atan(n_3(2)/n_3(1));
                
            end
        end
        
      case true
        %% 3D case
        sigma_1 = zeros(size(sigma_xx));
        sigma_2 = zeros(size(sigma_xx));
        sigma_3 = zeros(size(sigma_xx));
        theta_1 = zeros(size(sigma_xx));
        theta_2 = zeros(size(sigma_xx));
        theta_3 = zeros(size(sigma_xx));
        
        sigma_1pxy = zeros(size(sigma_xx));
        sigma_2pxy = zeros(size(sigma_xx));
        theta_1pxy = zeros(size(sigma_xx));
        theta_2pxy = zeros(size(sigma_xx));
        
        xsize = size(sigma_xx,1);
        ysize = size(sigma_xx,2);
        
        warningc = 0;
        for j=1:xsize
            for k=1:ysize
                
                sigma_tensor = [sigma_xx(j,k) sigma_xy(j,k) sigma_xz(j,k);...
                                sigma_xy(j,k) sigma_yy(j,k) sigma_yz(j,k);...
                                sigma_xz(j,k) sigma_yz(j,k) sigma_zz(j,k)];
                
                [V,D] = eig(sigma_tensor);
                [d,ind] = sort(diag(D),'descend');
                Ds = D(ind,ind);
                Vs = V(:,ind);
                
                sigma_1(j,k) = Ds(1,1);
                sigma_2(j,k) = Ds(2,2);
                sigma_3(j,k) = Ds(3,3);
                n_1 = Vs(:,1);
                n_2 = Vs(:,2);
                n_3 = Vs(:,3);
                
%                 if n_1(3,1) < 0.3
% $$$                 speremoben1(j,k) = n_1(3,1);
% $$$                 speremoben2(j,k) = n_2(3,1);
% $$$                 speremoben3(j,k) = n_3(3,1);
%                 end
                
                theta_1(j,k) = atan(n_1(2)/n_1(1));
                theta_2(j,k) = atan(n_2(2)/n_2(1));
                theta_3(j,k) = atan(n_3(2)/n_3(1));
                
                if abs(n_1(3)) < 0.1
                    sigma_1pxy(j,k) = sigma_1(j,k);
                    theta_1pxy(j,k) = theta_1(j,k);
                    if abs(n_2(3)) < 0.1
                        sigma_2pxy(j,k) = sigma_2(j,k);
                        theta_2pxy(j,k) = theta_2(j,k);
                    elseif abs(n_3(3)) < 0.1
                        sigma_2pxy(j,k) = sigma_3(j,k);
                        theta_2pxy(j,k) = theta_3(j,k);
                    else
%                         warning("something went wrong with the principal angles normal vectors");
                        warningc = warningc + 1;
                    end
                elseif abs(n_2(3)) < 0.1
                    sigma_1pxy(j,k) = sigma_2(j,k);
                    sigma_2pxy(j,k) = sigma_3(j,k);
                    theta_1pxy(j,k) = theta_2(j,k);
                    theta_2pxy(j,k) = theta_3(j,k);
                else
%                     warning("something went wrong with the principal angles normal vectors");
                    warningc = warningc + 1;
                end
                
            end
        end
        
        warning("Total warnings = " + string(warningc/(xsize*ysize)*100) + "%.");
        %varargout{1} = sigma_2;
        %varargout{2} = theta_2;

        varargout{1} = sigma_2;
        varargout{2} = sigma_1pxy;
        varargout{3} = sigma_2pxy;
        varargout{4} = theta_1pxy;
    
    end

end