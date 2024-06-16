function selected_fractures = select_by_method(fractures,list)
% select Ponti's database fractures based on their verification status
% list = [b,b,b,b]: 4 booleans 1 or 0, where 1 includes methodology, 0
% excludes it

    methodologies = ["Field","Imagery","Inferred","Remote Sensing"];
    methodologies(list==1) = [];

    selected_fractures = fractures;
    for f=length(fractures):-1:1
        frac = fractures(f);
        for m = 1:length(methodologies)
            if strcmp(frac.Class,methodologies(m))
                selected_fractures(f) = [];
                break
            end
        end
    end

end

