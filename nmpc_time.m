clear all, clc;

BEGIN_ACADO;                                % Always start with "BEGIN_ACADO". 
    
    acadoSet('problemname', 'nmpc_time');        % Set your problemname. If you 
    acadoSet('results_to_file', false);     % skip this, all files will
                                            % be named "myAcadoProblem"
    
    DifferentialState       x y z v_x v_y v_z phi theta psi THRUST
    
    Control                 delta_phi delta_theta delta_psi delta_THRUST

    Parameter               Time
    
    SAMPLING_TIME = 1;
    k_F = 6.11*10^-8;
    k_M = 1.5*10^-9;
    L = 0.175;
    I_xx = 2.32*10^-3;
    I_yy = 2.32*10^-3;
    I_zz = 4.00*10^-3;
    m = 0.5;
    g = 9.81;
    C_d = 0.5;
    %% Model Differential Equations
    f = acado.DifferentialEquation(0.0, Time);       % Set the differential equation object
    
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
    ocp = acado.OCP(0.0, Time);    % Set up the Optimal Control Problem (OCP)
                                   % Start at 0s, optimize for time
    ocp.minimizeMayerTerm(Time);
    
    ocp.subjectTo( f );
    ocp.subjectTo( 'AT_START', x ==  0.0 ); 
    ocp.subjectTo( 'AT_START', y ==  0.0 ); 
    ocp.subjectTo( 'AT_START', z ==  -1.5 ); 
    ocp.subjectTo( 'AT_START', v_x ==  0.0 ); 
    ocp.subjectTo( 'AT_START', v_y ==  0.0 ); 
    ocp.subjectTo( 'AT_START', v_z ==  0.0 ); 
    ocp.subjectTo( 'AT_START', phi ==  0.0 ); 
    ocp.subjectTo( 'AT_START', theta ==  0.0 ); 
    ocp.subjectTo( 'AT_START', psi ==  0.0 ); 
    ocp.subjectTo( 'AT_START', THRUST ==  -4.916 ); 
    
    ocp.subjectTo( 'AT_END', x ==  10.0 ); 
    ocp.subjectTo( 'AT_END', y ==  0.0 ); 
    ocp.subjectTo( 'AT_END', z ==  -1.5 ); 
    ocp.subjectTo( 'AT_END', v_x ==  0.0 ); 
    ocp.subjectTo( 'AT_END', v_y ==  0.0 ); 
    ocp.subjectTo( 'AT_END', v_z ==  0.0 ); 
    ocp.subjectTo( 'AT_END', phi ==  0.0 ); 
    ocp.subjectTo( 'AT_END', theta ==  0.0 ); 
    ocp.subjectTo( 'AT_END', psi ==  0.0 ); 
    ocp.subjectTo( 'AT_END', THRUST ==  -4.916  ); 
                              
    ocp.subjectTo( -1.6 <= z <= -1.4 );    
    ocp.subjectTo( -45*(pi/180) <= phi <= 45*(pi/180) );
    ocp.subjectTo( -45*(pi/180) <= theta <= 45*(pi/180) );
    ocp.subjectTo( -45*(pi/180) <= psi <= 45*(pi/180) );
    ocp.subjectTo( -23*SAMPLING_TIME <= delta_phi <= 23*SAMPLING_TIME );
    ocp.subjectTo( -23*SAMPLING_TIME <= delta_theta <= 23*SAMPLING_TIME );
    ocp.subjectTo( -23*SAMPLING_TIME <= delta_psi <= 23*SAMPLING_TIME );
    ocp.subjectTo( -40*SAMPLING_TIME <= delta_THRUST <= 40*SAMPLING_TIME );

    ocp.subjectTo( -19.621 <= THRUST <= 0 );

    ocp.subjectTo(0 <= Time <= 5 );
    
    
    %% Controller
    algo = acado.OptimizationAlgorithm(ocp);
    algo.set( 'KKT_TOLERANCE', 1e-10 );     % Set a custom KKT tolerance
   
        
END_ACADO;           % Always end with "END_ACADO".
                     % This will generate a file problemname_ACADO.m. 
                     % Run this file to get your results. You can
                     % run the file problemname_ACADO.m as many
                     % times as you want without having to compile again.



out = nmpc_time_RUN();                % Run the test. The name of the RUN file
                                            % is problemname_RUN, so in
                                            % this case getting_started_RUN
                                            
draw_nmpc_time