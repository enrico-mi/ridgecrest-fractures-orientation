% Batch compute of stresses at fractures' centers

% angles_list: range of angles to investigate in degrees from North (CW positive); for uncertainty analyses with large uncertainty angles, it is recommended to expand the range to -20:35

addpath('./computes');

angles_list = -10:1:25;

parfor sa=1:length(angles_list)
    t = getCurrentTask();
    logname = "logTask" + num2str(t.ID) + "angle" + num2str(angles_list(sa));
    diary logname
    diary on
    disp("Angle " + num2str(angles_list(sa)) + " is performed by task " + num2str(t.ID)); 
    run_calc_stresses_NW(angles_list(sa),true);
    diary off
end