clear all, clc;

BEGIN_ACADO;                                % Always start with "BEGIN_ACADO". 
    
    acadoSet('problemname', 'nmpc2');        % Set your problemname. If you 
    acadoSet('results_to_file', false);     % skip this, all files will
                                            % be named "myAcadoProblem"
    
    DifferentialState       x y z v_x v_y v_z phi theta psi THRUST

%     Control                 phi theta psi THRUST
    
    Control                 delta_phi delta_theta delta_psi delta_THRUST

    OCP_HORIZON=2;
    HORIZON_INTERVALS=20;
    SIMULATION_TIME=6.0;
    SAMPLING_TIME = (OCP_HORIZON/HORIZON_INTERVALS)/1.5;
    k_F = 6.11*10^-8;
    k_M = 1.5*10^-9;
    L = 0.175;
    I_xx = 2.32*10^-3;
    I_yy = 2.32*10^-3;
    I_zz = 4.00*10^-3;
    m = 0.5;
    g = 9.81;
    C_d = 0.5;
    

    %% reference trajectory
 ref = [2.0       10.00       0.00     -1.50      0.00        0.00        0.00     0.00        0.00        0.00      -4.916];      % Set up a given reference trajectory
%        2.0       3.00       0.00	  -1.50      0.00        0.00        0.00     0.00        0.00        0.00      -4.916];
%         3         3.00       6.00	  -1.50      0.00        0.00        0.00     0.00        0.00        (0.5*pi);
%         4.0       1.00       8.00	  -1.50      0.00        0.00        0.00     0.00        0.00        pi
%         6         -1.00      8.00	  -1.50      0.00        0.00        0.00     0.00        0.00        pi];
%         1.0       1.00       0.00	  -1.50      0.00        0.00        0.00;
%         1.0       1.00       0.00	  -1.50      0.00        0.00        0.00;
%         1.0       1.00       0.00	  -1.50      0.00        0.00        0.00;
%          3.0       1.00       0.00	  -1.50      0.00        0.00        0.00];
        
            
    %   TIME      x_REF      y_REF    z_REF     v_x_REf     v_y_REF     v_z_REF  phi_REF    theta_REF   psi_REF     p_REF       q_REF       r_REF 
    %% Model Differential Equations
    f = acado.DifferentialEquation();       % Set the differential equation object
    
%     velocity in earth frame
    f.add(dot(x) == cos(psi)*v_x - sin(psi)*v_y);
    f.add(dot(y) == cos(psi)*v_y + sin(psi)*v_x);
    f.add(dot(z) == v_z);
        
%     acceleration in body frame
    f.add(dot(v_x) == sin(theta) * (THRUST/m) - C_d * v_x);
    f.add(dot(v_y) == -sin(phi)*(THRUST/m) - C_d * v_y);
    f.add(dot(v_z) == ((cos(phi) * cos(theta))*THRUST + g*m)/m - C_d * v_z ); 
    
    
    f.add(dot(phi) == delta_phi);
    f.add(dot(theta) == delta_theta);
    f.add(dot(psi) == delta_psi);
    f.add(dot(THRUST) == delta_THRUST);
    
    %% Optimal Control Problem
    ocp = acado.OCP(0.0, OCP_HORIZON , HORIZON_INTERVALS);    % Set up the Optimal Control Problem (OCP)
                                            % Start at 0s, control in 100
                                            % intervals upto 5s
                                            
    h={x, y, z, v_x, v_y, v_z, phi, theta,psi, delta_psi};   % the LSQ-Function

    Q = eye(10);                             % The weighting matrix 
    Q(1,1) = 10;
    Q(2,2) = 10;
    Q(3,3) = 10;

    Q(4,4) = 0;                             %constrain velocities to eliminate overshoot
    Q(5,5) = 0;
    Q(6,6) = 0;
    
    Q(7,7) = 10;
    Q(8,8) = 1;
    Q(9,9) = 2;
    
    Q(10,10) = 0;

    ref_h = zeros(1,10);                        % The reference
    ref_h(3)=-1.50;

    ocp.minimizeLSQ( Q, h, ref_h );             % Minimize this Least Squares Term
      
    ocp.subjectTo( f );                         % Your OCP is always subject to your 
    ocp.subjectTo( -1.6 <= z <= -1.4 );                  % differential equation
    ocp.subjectTo( -25*(pi/180) <= phi <= 25*(pi/180) );
    ocp.subjectTo( -25*(pi/180) <= theta <= 25*(pi/180) );
    ocp.subjectTo( -25*(pi/180) <= psi <= 25*(pi/180) );
    ocp.subjectTo( -23*SAMPLING_TIME <= delta_phi <= 23*SAMPLING_TIME );
    ocp.subjectTo( -23*SAMPLING_TIME <= delta_theta <= 23*SAMPLING_TIME );
    ocp.subjectTo( -23*SAMPLING_TIME <= delta_psi <= 23*SAMPLING_TIME );
    ocp.subjectTo( -40*SAMPLING_TIME <= delta_THRUST <= 40*SAMPLING_TIME );

    ocp.subjectTo( -19.621 <= THRUST <= 0 );

    %% Process
    outputFunction = acado.OutputFcn();
    dynamicSystem = acado.DynamicSystem(f, outputFunction);    % f is ODE
    process = acado.Process(dynamicSystem, 'INT_RK78');

    
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
        
    %% Simulation Environment
    sim = acado.SimulationEnvironment(0.0, SIMULATION_TIME, process, controller);

    x0=zeros(1,10);
    x0(3)=-1.50;
    x0(10)=-4.916;
    sim.init( x0 )              %starting values of all states
    
    

    
    
END_ACADO;           % Always end with "END_ACADO".
                     % This will generate a file problemname_ACADO.m. 
                     % Run this file to get your results. You can
                     % run the file problemname_ACADO.m as many
                     % times as you want without having to compile again.



out = nmpc2_RUN();                % Run the test. The name of the RUN file
                                            % is problemname_RUN, so in
                                            % this case getting_started_RUN
                                            
draw_nmpc_improved