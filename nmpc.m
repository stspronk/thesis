clear all, clc;

BEGIN_ACADO;                                % Always start with "BEGIN_ACADO". 
    
    acadoSet('problemname', 'nmpc');        % Set your problemname. If you 
    acadoSet('results_to_file', true);     % skip this, all files will
                                            % be named "myAcadoProblem"
    
    DifferentialState       x;              % x-position(inertial frame)
    DifferentialState       y;              % y-postion(inertial frame)
    DifferentialState       z;              % z-position(inertial frame)
    DifferentialState       v_x;            % x-velocity(body frame)
    DifferentialState       v_y;            % y-velocity(body frame)
    DifferentialState       v_z;            % z-velocity(body frame)
    DifferentialState       phi;
    DifferentialState       theta;
    DifferentialState       psi;
    DifferentialState       p;              % roll rate q r phi theta psi x y z v_x v_y v_z
    DifferentialState       q;              % pitch rat     
    DifferentialState       r;              % yaw rate
    
    Control                 omega1 omega2 omega3 omega4;
     
    
%     Disturbance             d;              % 
%   Parameters;
    SAMPLING_TIME=0.05;
    OCP_HORIZON=2;
    SIMULATION_TIME=6.0;
    k_F = 6.11*10^-8;
    k_M = 1.5*10^-9;
    L = 0.175;
    I_xx = 2.32*10^-3;
    I_yy = 2.32*10^-3;
    I_zz = 4.00*10^-3;
    m = 0.5;
    g = 9.81;
    %omega=[omega1 omega2 omega3 omega4];
    %F = (omega(1)^2+omega(2)^2+omega(3)^2+omega(4)^2)*(-k_F);
    %M = [(-omega(2)^2+omega(4)^2)*k_F*L (omega(1)^2-omega(3)^2)*k_F*L ...
    %(-k_M*omega(1)^2+k_M*omega(2)^2-k_M*omega(3)^2+k_M*omega(4)^2)]';


    %R_d_angle = [1 tan(theta)*sin(phi) tan(theta)*cos(phi);...
                    %0 cos(phi) -sin(phi);...
                    %0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
                
    %R_B_E = [cos(psi)*cos(theta) cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi) ...
      %cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);...
      %sin(psi)*cos(theta) sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi) ...
      %sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);...
      % -sin(theta) cos(theta)*sin(phi) cos(theta)*cos(phi)];
   % R_E_B = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);...
    % sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)...
   % sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);...
    % cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi)...
    % cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];
 
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

    %% reference trajectory
 ref = [0.0       0.00       0.00     -1.50      0.00        0.00        0.00     0.00       0.00      0.00        0.00        0.00        0.00;      % Set up a given reference trajectory
        2.0       0.00       0.00	  -1.50      0.00        0.00        0.00     0.00       0.00     0.00        0.00        0.00        0.00;
        3.0       1.00       -2.00	  -1.50      0.00        0.00        0.00     0.00       0.00     0.00        0.00        0.00        0.00];

            
    %   TIME      x_REF      y_REF    z_REF     v_x_REf     v_y_REF     v_z_REF  phi_REF    theta_REF psi_REF     p_REF       q_REF       r_REF 
    %% Model Differential Equations
    f = acado.DifferentialEquation();       % Set the differential equation object
    
    % f.add(dot(position) == R_B_E * velocity_body);   
    f.add(dot(x) == (cos(psi)*cos(theta))*v_x+(cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi))*v_y+(cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi))*v_z);
    f.add(dot(y) == (sin(psi)*cos(theta))*v_x+(sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi))*v_y+(sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi))*v_z);
    f.add(dot(z) == -sin(theta)*v_x+(cos(theta)*sin(phi))*v_y+(cos(theta)*cos(phi))*v_z);
    
    %f.add(dot(velocity_body) == (R_E_B*[0 0 g]'*m + [0 0 ((omega1^2+omega2^2+omega3^2+omega4^2)...
        %*(-k_F))]')/m -cross(Omega,velocity_body));                                      
    f.add(dot(v_x) == (-sin(theta)*g*m)/m-(q*v_z-r*v_y));
    f.add(dot(v_y) == (sin(phi)*cos(theta)*g*m)/m-(p*v_z-r*v_x));
    f.add(dot(v_z) == (cos(phi)*cos(theta)*g*m+(omega1^2+omega2^2+omega3^2+omega4^2)*-k_F)/m -(p*v_y-q*v_x));
    
    %f.add(dot(angle) == R_d_angle*Omega);                                                                    
    f.add(dot(phi) == p+q*tan(theta)*sin(phi)+r*tan(theta)*cos(phi));
    f.add(dot(theta) == q*cos(phi)+r*-sin(phi));
    f.add(dot(psi) == q*(sin(phi)/cos(theta))+r*(cos(phi)/cos(theta)));
   
    %f.add(dot(Omega) == inv(I)*(([(-omega2^2+omega4^2)*k_F*L (omega1^2-omega3^2)*k_F*L...
        %(-k_M*omega1^2+k_M*omega2^2-k_M*omega3^2+k_M*omega4^2)]')-cross(Omega,I*Omega)));  
    f.add(dot(p) == 1/I_xx * ((-omega2^2+omega4^2)*k_F*L -(q*r*I_zz-r*q*I_yy)));
    f.add(dot(q) == 1/I_yy * ((omega1^2-omega3^2)*k_F*L - (p*r*I_zz - r*p*I_xx)));
    f.add(dot(r) == 1/I_zz * ((-k_M*omega1^2+k_M*omega2^2-k_M*omega3^2+k_M*omega4^2) - (p*q*I_yy - q*p*I_xx)));
    
    
    
    %% Optimal Control Problem
    ocp = acado.OCP(0.0, OCP_HORIZON ,20.0);    % Set up the Optimal Control Problem (OCP)
                                            % Start at 0s, control in 100
                                            % intervals upto 5s
                                            
    
    h={x, y, z, v_x, v_y, v_z, phi, theta, psi, p, q, r};   % the LSQ-Function

    %Q = eye(3);                             % The weighting matrix 
    ref_h = zeros(1,12);                         % The reference
        
    ocp.minimizeLSQ( h, ref_h );             % Minimize this Least Squares Term
    %ocp.minimizeLSQ( h, r );               % (Other possibility)
    %ocp.minimizeLSQ( h );                  % (Other possibility)
    
    ocp.subjectTo( f );                     % Your OCP is always subject to your 
    ocp.subjectTo(z <= 0.0);                 % differential equation
    %ocp.subjectTo( -10 <= v_x <= 10);
    %ocp.subjectTo( -10 <= v_y <= 10);
    %ocp.subjectTo( -10 <= v_z <= 10);
    ocp.subjectTo( 0.0 <= omega1 <= 8960.0 );
    ocp.subjectTo( 0.0 <= omega2 <= 8960.0 );
    ocp.subjectTo( 0.0 <= omega3 <= 8960.0 );
    ocp.subjectTo( 0.0 <= omega4 <= 8960.0 );
    

    
    
    %% Process
    outputFunction = acado.OutputFcn();
    dynamicSystem = acado.DynamicSystem(f, outputFunction);    % f is ODE
    process = acado.Process(dynamicSystem, 'INT_RK45');
    
    
%     process.setProcessDisturbance(d)           % matrix analogous to measurement matric of opc
                                               % rows are different
                                               % measurements, collumns
                                               % contain different measured
                                               % values of the disturbances
%     process.initializeAlgebraicStates(v)       % vector of initial values of states

    
    %% Controller
    algo = acado.RealTimeAlgorithm(ocp, SAMPLING_TIME);
    algo.set('INTEGRATOR_TYPE', 'INT_RK78');
    algo.set( 'INTEGRATOR_TOLERANCE',   1e-8);    
    algo.set( 'ABSOLUTE_TOLERANCE',     1e-8);
    algo.set('MAX_NUM_ITERATIONS', 5.0 );
   
    reference =acado.StaticReferenceTrajectory(ref);
    controller = acado.Controller(algo , reference );
%     reference = acado.StaticReferenceTrajectory([r]);            % static reference
%     reference = acado.PeriodicReferenceTrajectory([r])          % periodic reference
    
%     zeroReference = acado.StaticReferenceTrajectory();            %zero reference
        
    %% Simulation Environment
    sim = acado.SimulationEnvironment(0.0, SIMULATION_TIME, process, controller);

    x0=zeros(1,12);
    sim.init( x0 )              %starting values of all states
    
    

    
    
END_ACADO;           % Always end with "END_ACADO".
                     % This will generate a file problemname_ACADO.m. 
                     % Run this file to get your results. You can
                     % run the file problemname_ACADO.m as many
                     % times as you want without having to compile again.



out = nmpc_RUN();                % Run the test. The name of the RUN file
                                            % is problemname_RUN, so in
                                            % this case getting_started_RUN
                                            
draw_nmpc;