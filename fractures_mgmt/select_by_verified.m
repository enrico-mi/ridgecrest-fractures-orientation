function selected_fractures = select_by_verified(fractures,condition)
% select Ponti's database fractures based on their verification status

    selected_fractures = fractures;
    if strcmp(condition,'Yes') || strcmp(condition,'Partial')
        for f=length(fractures):-1:1
            frac = fractures(f);
            if ~(strcmp(frac.Verified,condition))
                selected_fractures(f) = [];
            end
        end
    elseif strcmp(condition,'All')
        % keep current selection
    else
        error('selection.verified.condition not valid');
    end
end

