function [] = main()
close all;
clear;
clc;
%MAIN ????????????
%   ????????
global simulation_time step pointer;
global drone_states actuator_states;
global i time desired_omega desired_angular_velocity
global desired_angle desired_velocity_body desired_position m g
global controller guidance

controller = 'PID';  % INDI or PID
guidance = 'STEP'; % 

global_parameters();
%desired_angular_velocity = zeros(3,simulation_time/step);
for i = 0:step:simulation_time
    %desired_position(:,pointer)=[3 0 -1.5]';
    [ desired_position(:,pointer),desired_angle(3,pointer) ] = guidance_drone();
    %desired_position(:,pointer)=[3 0 -1.5]';
    
 %% sensor model
     sensor_model();
 %% position controller
    [ desired_velocity_body(:,pointer) ] = controller_position(desired_position(:,pointer));

   %% velocity controller
   %desired_velocity_body(:,pointer) = [1 1 -1]';
    [ desired_angle(1:2,pointer),F ] = controller_velocity_body( desired_velocity_body(:,pointer) );

 %% angle controller
  %desired_angle(:,pointer) = [5/180*pi 0 0]';
    desired_angular_velocity(:,pointer) = controller_angle( desired_angle(:,pointer));
    %temp = desired_angular_velocity(:,pointer)
 %% angular velocity controller 
 %desired_angular_velocity(:,pointer) = [0 0 0.5];

    M = controller_angular_velocity( desired_angular_velocity(:,pointer));

    desired_omega(:,pointer) = convert_force_2_omega(F,M);
    %% actuator dynamics
    actuator_states(:,pointer+1) = runge_kutta_actuator(actuator_states(:,pointer),desired_omega(:,pointer));
    
    %% drone dynamics
    drone_states(:,pointer+1) = runge_kutta_drone(drone_states(:,pointer),actuator_states(:,pointer+1));
    %%
    pointer = pointer+1;
    time(pointer) = i;
end

plot_actuator();
end

