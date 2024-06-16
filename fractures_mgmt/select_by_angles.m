function selected_fractures = select_by_angles(fractures, angles)
% select fractures by their average orientation; interval given in radians

    selected_fractures = fractures;
    for f=length(fractures):-1:1
        angle = rad2deg(fractures(f).AngleAveUTM);
        if angle < angles.range(1) || angle > angles.range(2)
            selected_fractures(f) = [];
        end
    end

end
