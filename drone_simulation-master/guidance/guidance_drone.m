function [ desired_position,desired_psi ] = guidance_drone()
%GUIDANCE_DRONE Summary of this function goes here
%   Detailed explanation goes here
global time guidance pointer
step_distance = 1;
switch guidance
    
    case 'STEP'
        if time(pointer)<10
            desired_position = step_distance*[3 0 -1.5]; %*time(pointer)]';
%         elseif time(pointer) < 10
%             desired_position = step_distance*[0 1 -1*time(pointer)]';
%         elseif time(pointer) < 15
%             desired_position = step_distance*[-1 0 -1*time(pointer)]';
%         elseif time(pointer)<20
%             desired_position = step_distance*[0 -1 -1*time(pointer)]';
%         else
%             desired_position = step_distance*[1 0 -1*time(pointer)]';
        end
        desired_psi = 0;
end

end

