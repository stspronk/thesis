clear all, close all, clc;
% sysdeps/unix/sysv/linux/x86_64/pthread_cond_wait
% mex -v GCC='/usr/bin/gcc-4.7' timestwo.c
BEGIN_ACADO;                                % Always start with "BEGIN_ACADO". 
    
    acadoSet('problemname', 'nmpc2Dsym');     % Set your problemname. If you 
                                           % skip this, all files will
                                          % be named "myAcadoProblem"
    
    
    %DifferentialState       p;              % roll rate
    DifferentialState        q theta x z v_x v_z;              % pitch rate     
    %DifferentialState       r;              % yaw rate
    %DifferentialState       phi;
    %DifferentialState       theta;
    %DifferentialState       psi;
    %DifferentialState       x;              % x-position(inertial frame)
    %DifferentialState       y;              % y-postion(inertial frame)
    %DifferentialState       z;              % z-position(inertial frame)
    %DifferentialState       v_x;            % x-velocity(body frame)
    %DifferentialState       v_y;            % y-velocity(body frame)
    %DifferentialState       v_z;            % z-velocity(body frame)

%     DifferentialState       Omega;
%     DifferentialState       position;
%     DifferentialState       velocity_body;
%     DifferentialState       angle;
%     Parameter T;

    Control                 omega1 omega2 omega3 omega4;              % 
     
    
%     Disturbance             d;              % 
%     static Parameters;
    SAMPLING_TIME=0.05;
    OCP_HORIZON=3;
    SIMULATION_TIME=7.0;
    k_F = 6.11*10^-8;
    k_M = 1.5*10^-9;
    L = 0.175;
    I_xx = 2.32*10^-3;
    I_yy = 2.32*10^-3;
    I_zz = 4.00*10^-3;
    m = 0.5;
    g = 9.8;
    I = [I_xx 0 0; 0 I_yy 0; 0 0 I_zz];
    %omega=[omega1 omega2 omega3 omega4];
    %F = (omega(1)^2+omega(2)^2+omega(3)^2+omega(4)^2)*(-k_F);
    %M = [(-omega(2)^2+omega(4)^2)*k_F*L (omega(1)^2-omega(3)^2)*k_F*L ...
    %(-k_M*omega(1)^2+k_M*omega(2)^2-k_M*omega(3)^2+k_M*omega(4)^2)]';


    %R_d_angle = [1 tan(theta)*sin(phi) tan(theta)*cos(phi);...
                    %0 cos(phi) -sin(phi);...
                    %0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
                
                    
    %R_B_E = [cos(psi)*cos(theta) 
    %cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi) ...
      %cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);...
      
      %sin(psi)*cos(theta) 
      %sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi) ...
      %sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);...
      
      % -sin(theta)
      %cos(theta)*sin(phi)
      %cos(theta)*cos(phi)];
      
      
   % R_E_B = [cos(theta)*cos(psi)
   %cos(theta)*sin(psi)
   %-sin(theta);...
   
    % sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)...
   % sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi)
   %sin(phi)*cos(theta);...
   
    % cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi)...
    % cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)
    %cos(phi)*cos(theta)];
 
    
    
   % position = [x y z]';
   % velocity_body = [v_x v_y v_z]';
   % angle = [phi theta psi]';
   % Omega = [p q r]';

%     d_position = R_B_E * velocity_body;
%     d_velocity_body = (R_E_B*[0 0 g]'*m + [0 0 F]')/m -cross(Omega,velocity_body);
%     d_angle = R_d_angle*Omega;
%     d_Omega = inv(I)*(M-cross(Omega,I*Omega));
    
    
%     current_angular_velocity_accel = d_Omega;

%     d_states = [d_position;d_velocity_body;d_angle;d_Omega];
    %% Model Differential Equations
    f = acado.DifferentialEquation();       % Set the differential equation object
    
    f.add(dot(q) == 1/I_yy * ((omega1^2-omega3^2)*k_F*L));
    f.add(dot(theta) == q);

    % f.add(dot(position) == R_B_E * velocity_body);   
    f.add(dot(x) == cos(theta)*v_x+sin(theta)*v_z);
    %f.add(dot(y) == sin(psi)*cos(theta)*v_x+sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi)*v_y...
       % +sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi)*v_z);
    f.add(dot(z) == -sin(theta)*v_x+cos(theta)*v_z);
    
    %f.add(dot(velocity_body) == (R_E_B*[0 0 g]'*m + [0 0 ((omega1^2+omega2^2+omega3^2+omega4^2)...
        %*(-k_F))]')/m -cross(Omega,velocity_body));                                      
    f.add(dot(v_x) == (-sin(theta)*g*m)/m-(q*v_z));
    %f.add(dot(v_y) == (sin(phi)*cos(theta)*g*m)/m-(p*v_z-r*v_x));
    f.add(dot(v_z) == (cos(theta)*g*m+(omega1^2+omega2^2+omega3^2+omega4^2)*-k_F)/m -(-q*v_x));
    
    %f.add(dot(angle) == R_d_angle*Omega);                                                                    
    %f.add(dot(phi) == p+q*tan(theta)*sin(phi)+r*tan(theta)*cos(phi));
    %f.add(dot(psi) == q*(sin(phi)/cos(theta))+r*(cos(phi)/cos(theta)));
   
    %f.add(dot(Omega) == inv(I)*(([(-omega2^2+omega4^2)*k_F*L (omega1^2-omega3^2)*k_F*L...
        %(-k_M*omega1^2+k_M*omega2^2-k_M*omega3^2+k_M*omega4^2)]')-cross(Omega,I*Omega)));  
    %f.add(dot(p) == 1/I_xx * ((-omega2^2+omega4^2)*k_F*L -(q*r*I_zz-r*q*I_yy)));
    
    %f.add(dot(r) == 1/I_zz * ((-k_M*omega1^2+k_M*omega2^2-k_M*omega3^2+k_M*omega4^2) - (p*q*I_yy - q*p*I_xx)));
    
    %% reference trajectory
    ref = [ 0.00      0.00       0.00        0.00      -1.50        0.00        0.00;      % Set up a given reference trajectory
            2.0       0.00       0.00        0.00      -2.00        0.00        0.00;        
            4.0       0.00       0.00        1.00      -2.00        0.00        0.00;
            6.0       0.00       0.00        2.00      -1.50        0.00        0.00];
            %7.0       0.00       0.00        0.00      -1.50        0.00        0.00];
            %30.0       0.00       0.00        6.00      -1.50        0.00        0.00]; 
    %      TIME      q_REF      theta_REF    x_REF     z_REf      v_x_REF     v_z_REF 
    
    %% Optimal Control Problem
    ocp = acado.OCP(0.0, OCP_HORIZON ,30.0);    % Set up the Optimal Control Problem (OCP)
                                            % Start at 0s, control in 20
                                            % intervals upto 4s
                                            
    %ocp.minimizeMayerTerm(T);  
    h={q, theta, x,z, v_x, v_z};   % the LSQ-Function

    %Q = eye(6);                             % The weighting matrix UPDATE!!!
    r = zeros(1, 6);                         % The reference
        
    ocp.minimizeLSQ( h, r );             % Minimize this Least Squares Term
    %ocp.minimizeLSQ( h, r );               % (Other possibility)
    %ocp.minimizeLSQ( h );                  % (Other possibility)
    
    ocp.subjectTo( f );                     % Your OCP is always subject to your 
                                            % differential equation
  
    %ocp.subjectTo(z < 0)
    %control initial conditions
    %ocp.subjectTo( 'AT_START', omega1 == 0.0 ); 
    %ocp.subjectTo( 'AT_START', omega2 == 0.0 );     
    %ocp.subjectTo( 'AT_START', omega3 == 0.0 );     
    %ocp.subjectTo( 'AT_START', omega4 == 0.0 );     
    
    %states initial conditions
    %ocp.subjectTo( 'AT_START', x == 0 );
    %ocp.subjectTo( 'AT_START', z == 0 );
    %ocp.subjectTo( 'AT_START', q == 0 );
    %ocp.subjectTo( 'AT_START', v_x == 0 );
    %ocp.subjectTo( 'AT_START', v_z == 0 );
    %ocp.subjectTo( 'AT_START', theta == 0 );

    %states end conditions
    %ocp.subjectTo( 'AT_END', x == 0 );
    %ocp.subjectTo( 'AT_END', z == -1 );
    %ocp.subjectTo( 'AT_END', q == 0 );
    %ocp.subjectTo( 'AT_END', v_x == 0 );
    %ocp.subjectTo( 'AT_END', v_z == 0 );
    %ocp.subjectTo( 'AT_END', theta == 0 );
    
    %constraints
    %ocp.subjectTo(-0.4*pi <= theta <= 0.4*pi)
    ocp.subjectTo(z <= 0.0);
    ocp.subjectTo( -10 <= v_z <= 10);
    ocp.subjectTo( -10 <= v_x <= 10);
    ocp.subjectTo( 0.0 <= omega1 <= 8960.0 );
    ocp.subjectTo( 0.0 <= omega2 <= 8960.0 );
    ocp.subjectTo( 0.0 <= omega3 <= 8960.0 );
    ocp.subjectTo( 0.0 <= omega4 <= 8960.0 );
    %ocp.subjectTo(    0 <= T <= 15  ); 

    %% Process
    outputFunction = acado.OutputFcn();
    dynamicSystem = acado.DynamicSystem(f, outputFunction);    % f is ODE
    process = acado.Process(dynamicSystem, 'INT_RK45');
    
    %d=[0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    %   1.0 0.0 0.0 0.0 0.0 0.1 0.0];
    
    %process.setProcessDisturbance(d)           % matrix analogous to measurement matric of opc
                                               % rows are different
                                               % measurements, collumns
                                               % contain different measured
                                               % values of the disturbances
     %process.initializeAlgebraicStates(v)       % 

    
    %% Controller
    algo = acado.RealTimeAlgorithm(ocp, SAMPLING_TIME);
    %algo.set('HESSIAN_APPROXIMATION', 'GAUSS_NEWTON' );
    algo.set('INTEGRATOR_TYPE', 'INT_RK45');
    algo.set( 'INTEGRATOR_TOLERANCE',   1e-5);    
    algo.set( 'ABSOLUTE_TOLERANCE',     1e-4 );
    algo.set('MAX_NUM_ITERATIONS', 5.0 );

    %algo.initializeDifferentialStates(powerkite_states);
    %algo.initializeControls          (powerkite_controls);
    
    reference =acado.StaticReferenceTrajectory(ref);
    controller = acado.Controller(algo , reference );
    
    %reference = acado.StaticReferenceTrajectory(r);            % static reference
%     reference = acado.PeriodicReferenceTrajectory([r])          % periodic reference
    
%     zeroReference = acado.StaticReferenceTrajectory();            %zero reference
        
    
    %% Simulation Environment
    sim = acado.SimulationEnvironment(0.0, SIMULATION_TIME, process, controller);
    
    x0=zeros(1,6);
    sim.init( x0 )              %starting values of all states
    
    

    
    
END_ACADO;           % Always end with "END_ACADO".
                     % This will generate a file problemname_ACADO.m. 
                     % Run this file to get your results. You can
                     % run the file problemname_ACADO.m as many
                     % times as you want without having to compile again.



out = nmpc2Dsym_RUN();                % Run the test. The name of the RUN file
                                            % is problemname_RUN, so in
                                            % this case getting_started_RUN
draw_nmpc2Dsym                                            
