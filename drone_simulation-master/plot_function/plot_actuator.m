function [ ] = plot_actuator()
%PLOT_ACTUATOR ????????????
%   ????????
global actuator_states time drone_states desired_omega pointer desired_angular_velocity
global desired_angle desired_velocity_body desired_position
figure(1);
desired_omega(:,pointer) = desired_omega(:,end);
subplot(2,2,1)
plot(time,actuator_states(1,:),time,desired_omega(1,:));
xlabel('t/s');
ylabel('rotor1 speed[rpm]');
legend('real','ref');
subplot(2,2,2)
plot(time,actuator_states(2,:),time,desired_omega(2,:));
xlabel('t/s');
ylabel('rotor2 speed[rpm]');
legend('real','ref');
subplot(2,2,3)
plot(time,actuator_states(3,:),time,desired_omega(3,:));
xlabel('t/s');
ylabel('rotor3 speed[rpm]');
legend('real','ref');
subplot(2,2,4)
plot(time,actuator_states(4,:),time,desired_omega(4,:));
xlabel('t/s');
ylabel('rotor4 speed[rpm]');
legend('real','ref');


figure(2)
desired_angular_velocity(:,pointer) = desired_angular_velocity(:,end);
desired_angle(:,pointer) = desired_angle(:,end);
desired_velocity_body(:,pointer) = desired_angle(:,end);
desired_position(:,pointer) =desired_position(:,end);
subplot(4,3,1)
plot(time,drone_states(1,:),time,desired_position(1,:));
xlabel('t/s');
ylabel('x earth[m]');
legend('real','ref');
subplot(4,3,2)
plot(time,drone_states(2,:),time,desired_position(2,:));
xlabel('t/s');
ylabel('y earth[m]');
legend('real','ref');
subplot(4,3,3)
plot(time,drone_states(3,:),time,desired_position(3,:));
xlabel('t/s');
ylabel('z earth[m]');
legend('real','ref');
subplot(4,3,4)
plot(time,drone_states(4,:),time,desired_velocity_body(1,:));
xlabel('t/s');
ylabel('x body[m/s]');
legend('real','ref');
subplot(4,3,5)
plot(time,drone_states(5,:),time,desired_velocity_body(2,:));
xlabel('t/s');
ylabel('y body[m/s]');
legend('real','ref');
subplot(4,3,6)
plot(time,drone_states(6,:),time,desired_velocity_body(3,:));
xlabel('t/s');
ylabel('z body[m/s]');
legend('real','ref');
subplot(4,3,7)
plot(time,drone_states(7,:)/pi*180,time,desired_angle(1,:)/pi*180);
xlabel('t/s');
ylabel('Phi[degree]');
legend('real','ref');
subplot(4,3,8)
plot(time,drone_states(8,:)/pi*180,time,desired_angle(2,:)/pi*180);
xlabel('t/s');
ylabel('Theta[degree]');
legend('real','ref');
subplot(4,3,9)
plot(time,drone_states(9,:)/pi*180,time,desired_angle(3,:)/pi*180);
xlabel('t/s');
ylabel('Psi[degree]');
legend('real','ref');
subplot(4,3,10)
plot(time,drone_states(10,:)/pi*180,time,desired_angular_velocity(1,:)/pi*180);
xlabel('t/s');
ylabel('p[degree/s]');
legend('real','ref');
subplot(4,3,11)
plot(time,drone_states(11,:)/pi*180,time,desired_angular_velocity(2,:)/pi*180);
xlabel('t/s');
ylabel('q[degree/s]');
legend('real','ref');
subplot(4,3,12)
plot(time,drone_states(12,:)/pi*180,time,desired_angular_velocity(3,:)/pi*180);
xlabel('t/s');
ylabel('r[degree/s]');
legend('real','ref');

figure(3)
plot3(drone_states(1,:),drone_states(2,:),-drone_states(3,:),'b','LineWidth',2);
hold on;
plot3(desired_position(:,1),desired_position(:,2),-desired_position(:,3),...
    'r --','LineWidth',2);
xlabel('x[m]');
ylabel('y[m]');
zlabel('z[m]');
legend('real','ref');
grid on;
end

