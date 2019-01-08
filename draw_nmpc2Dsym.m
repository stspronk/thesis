close all

figure;

subplot(2,3,1)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,2), 'r')
title('q');

subplot(2,3,2)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,3), 'r')
title('theta');

subplot(2,3,3)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,4), 'r')
title('x');

subplot(2,3,4)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,5), 'r')
title('z');

subplot(2,3,5)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,6), 'r')
title('v_x');

subplot(2,3,6)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,7), 'r')
title('v_z');


figure;

subplot(1,5,1)
plot(out.CONTROLS(:,1), out.CONTROLS(:,2), 'r')
title('omega_1');

subplot(1,5,2)
plot(out.CONTROLS(:,1), out.CONTROLS(:,3), 'r')
title('omega_2');

subplot(1,5,3)
plot(out.CONTROLS(:,1), out.CONTROLS(:,4), 'r')
title('omega_3');

subplot(1,5,4)
plot(out.CONTROLS(:,1), out.CONTROLS(:,5), 'r')
title('omega_4');

subplot(1,5,5)
plot(out.CONTROLS(:,1), out.CONTROLS(:,2)-out.CONTROLS(:,4), 'r')
title('q-index');