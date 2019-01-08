close all

figure;



subplot(2,3,1)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,2), 'r')
title('x');

subplot(2,3,2)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,3), 'r')
title('y');

subplot(2,3,3)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,4), 'r')
title('z');

subplot(2,3,4)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,5), 'r')
title('v_x');

subplot(2,3,5)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,6), 'r')
title('v_y');

subplot(2,3,6)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,7), 'r')
title('v_z');


figure;

subplot(2,3,1)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.CONTROLS(:,2), 'r')
title('\delta_{\phi}');

subplot(2,3,2)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.CONTROLS(:,3), 'r')
title('\delta_{\theta}');

subplot(2,3,3)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.CONTROLS(:,4), 'r')
title('\delta_{\psi}');

subplot(2,3,4)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,8), 'r')
title('\phi');

subplot(2,3,5)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,9), 'r')
title('\theta');

subplot(2,3,6)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,10), 'r')
title('\psi');


figure
subplot(2,1,1)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.CONTROLS(:,5), 'r')
title('\delta_{Thrust}');
subplot(2,1,2)
plot(linspace(0, out.PARAMETERS(1,2), length(out.PARAMETERS(:,2))), out.STATES(:,11), 'r')
title('Thrust');
% 
% subplot(1,5,2)
% plot(out.CONTROLS(:,1), out.CONTROLS(:,3), 'r')
% title('omega_2');
% 
% subplot(1,5,3)
% plot(out.CONTROLS(:,1), out.CONTROLS(:,4), 'r')
% title('omega_3');
% 
% subplot(1,5,4)
% plot(out.CONTROLS(:,1), out.CONTROLS(:,5), 'r')
% title('omega_4');
% 
% subplot(1,5,5)
% plot(out.CONTROLS(:,1), out.CONTROLS(:,2)-out.CONTROLS(:,4), 'r')
% title('q-index');
