figure;

subplot(2,3,1)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,1), 'r')
title('q');

subplot(2,3,2)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,2), 'r')
title('theta');

subplot(2,3,3)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,3), 'r')
title('x');

subplot(2,3,4)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,4), 'r')
title('z');

subplot(2,3,5)
plot(out.STATES_SAMPLED(:,1), out.STATES_SAMPLED(:,5), 'r')
title('v_x');

subplot(2,3,6)
plot(out.STATES_SAMPLED(:,4), out.STATES_SAMPLED(:,6), 'r')
title('v_z');


figure;

subplot(1,3,1)
plot(out.CONTROLS(:,1), out.CONTROLS(:,2), 'r')
title('CONTROL 1 DDR0');

subplot(1,3,2)
plot(out.CONTROLS(:,1), out.CONTROLS(:,3), 'r')
title('CONTROL 2 DPSI');

subplot(1,3,3)
plot(out.CONTROLS(:,1), out.CONTROLS(:,4), 'r')
title('CONTROL 3 DCL');
