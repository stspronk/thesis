close all

figure;

subplot(2,3,1)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,2), 'r')
title('x');

subplot(2,3,2)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,3), 'r')
title('y');

subplot(2,3,3)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,4), 'r')
title('z');

subplot(2,3,4)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,5), 'r')
title('v_x');

subplot(2,3,5)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,6), 'r')
title('v_y');

subplot(2,3,6)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,7), 'r')
title('v_z');


figure;

subplot(2,3,1)
plot(out.CONTROLS(:,1), out.CONTROLS(:,2), 'r')
<<<<<<< HEAD
title('\delta_{\phi}');

subplot(2,3,2)
plot(out.CONTROLS(:,1), out.CONTROLS(:,3), 'r')
title('\delta_{\theta}');

subplot(2,3,3)
plot(out.CONTROLS(:,1), out.CONTROLS(:,4), 'r')
title('\delta_{\psi}');

subplot(2,3,4)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,8), 'r')
title('\phi');

subplot(2,3,5)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,9), 'r')
title('\theta');

subplot(2,3,6)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,10), 'r')
title('\psi');
=======
title('delta_phi');

subplot(2,3,2)
plot(out.CONTROLS(:,1), out.CONTROLS(:,3), 'r')
title('delta_theta');

subplot(2,3,3)
plot(out.CONTROLS(:,1), out.CONTROLS(:,4), 'r')
title('delta_psi');

subplot(2,3,4)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,8), 'r')
title('phi');

subplot(2,3,5)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,9), 'r')
title('theta');

subplot(2,3,6)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,10), 'r')
title('psi');
>>>>>>> 2c74ee34512ed2c9377a395dee51bf258fe49c80


figure
subplot(2,1,1)
plot(out.CONTROLS(:,1), out.CONTROLS(:,5), 'r')
title('\delta_{Thrust}');
subplot(2,1,2)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,11), 'r')
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
