/*
*    This file is part of ACADO Toolkit.
*
*    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
*    Copyright (C) 2008-2009 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
*    Developed within the Optimization in Engineering Center (OPTEC) under
*    supervision of Moritz Diehl. All rights reserved.
*
*    ACADO Toolkit is free software; you can redistribute it and/or
*    modify it under the terms of the GNU Lesser General Public
*    License as published by the Free Software Foundation; either
*    version 3 of the License, or (at your option) any later version.
*
*    ACADO Toolkit is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public
*    License along with ACADO Toolkit; if not, write to the Free Software
*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*
*/


/**
*    Author David Ariens, Rien Quirynen
*    Date 2009-2013
*    http://www.acadotoolkit.org/matlab 
*/

#include <acado_optimal_control.hpp>
#include <acado_toolkit.hpp>
#include <acado/utils/matlab_acado_utils.hpp>

USING_NAMESPACE_ACADO

#include <mex.h>


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
 { 
 
    MatlabConsoleStreamBuf mybuf;
    RedirectStream redirect(std::cout, mybuf);
    clearAllStaticCounters( ); 
 
    mexPrintf("\nACADO Toolkit for Matlab - Developed by David Ariens and Rien Quirynen, 2009-2013 \n"); 
    mexPrintf("Support available at http://www.acadotoolkit.org/matlab \n \n"); 

    if (nrhs != 0){ 
      mexErrMsgTxt("This problem expects 0 right hand side argument(s) since you have defined 0 MexInput(s)");
    } 
 
    TIME autotime;
    DifferentialState x;
    DifferentialState y;
    DifferentialState z;
    DifferentialState v_x;
    DifferentialState v_y;
    DifferentialState v_z;
    DifferentialState omega1;
    DifferentialState omega2;
    DifferentialState omega3;
    DifferentialState omega4;
    Control phi;
    Control theta;
    Control psi;
    Control p;
    Control q;
    Control r;
    Function acadodata_f2;
    acadodata_f2 << x;
    acadodata_f2 << y;
    acadodata_f2 << z;
    acadodata_f2 << v_x;
    acadodata_f2 << v_y;
    acadodata_f2 << v_z;
    DVector acadodata_v1(6);
    acadodata_v1(0) = 0;
    acadodata_v1(1) = 0;
    acadodata_v1(2) = 0;
    acadodata_v1(3) = 0;
    acadodata_v1(4) = 0;
    acadodata_v1(5) = 0;
    DVector acadodata_v2(6);
    acadodata_v2(0) = 0;
    acadodata_v2(1) = 0;
    acadodata_v2(2) = 0;
    acadodata_v2(3) = 0;
    acadodata_v2(4) = 0;
    acadodata_v2(5) = 0;
    DMatrix acadodata_M1;
    acadodata_M1.read( "nmpc_data_acadodata_M1.txt" );
    DVector acadodata_v3(6);
    acadodata_v3(0) = 0;
    acadodata_v3(1) = 0;
    acadodata_v3(2) = 0;
    acadodata_v3(3) = 0;
    acadodata_v3(4) = 0;
    acadodata_v3(5) = 0;
    DifferentialEquation acadodata_f1;
    acadodata_f1 << dot(x) == ((-cos(phi)*sin(psi)+cos(psi)*sin(phi)*sin(theta))*v_y+(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))*v_z+cos(psi)*cos(theta)*v_x);
    acadodata_f1 << dot(y) == ((cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta))*v_y+(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi))*v_z+cos(theta)*sin(psi)*v_x);
    acadodata_f1 << dot(z) == ((-sin(theta))*v_x+cos(phi)*cos(theta)*v_z+cos(theta)*sin(phi)*v_y);
    acadodata_f1 << dot(v_x) == ((-sin(theta))*5.00000000000000000000e-01/5.00000000000000000000e-01*9.81000000000000049738e+00-q*v_z+r*v_y);
    acadodata_f1 << dot(v_y) == (5.00000000000000000000e-01/5.00000000000000000000e-01*9.81000000000000049738e+00*cos(theta)*sin(phi)-p*v_z+r*v_x);
    acadodata_f1 << dot(v_z) == (((-6.10999999999999978704e-08)*(pow(omega1,2.00000000000000000000e+00)+pow(omega2,2.00000000000000000000e+00)+pow(omega3,2.00000000000000000000e+00)+pow(omega4,2.00000000000000000000e+00))+5.00000000000000000000e-01*9.81000000000000049738e+00*cos(phi)*cos(theta))/5.00000000000000000000e-01-p*v_y+q*v_x);

    OCP ocp1(0, 2, 20);
    ocp1.minimizeLSQ(acadodata_f2, acadodata_v2);
    ocp1.subjectTo(acadodata_f1);
    ocp1.subjectTo(z <= 0.00000000000000000000e+00);
    ocp1.subjectTo(0.00000000000000000000e+00 <= omega1 <= 8.96000000000000000000e+03);
    ocp1.subjectTo(0.00000000000000000000e+00 <= omega2 <= 8.96000000000000000000e+03);
    ocp1.subjectTo(0.00000000000000000000e+00 <= omega3 <= 8.96000000000000000000e+03);
    ocp1.subjectTo(0.00000000000000000000e+00 <= omega4 <= 8.96000000000000000000e+03);


    OutputFcn acadodata_f3;

    DynamicSystem dynamicsystem1( acadodata_f1,acadodata_f3 );
    Process process2( dynamicsystem1,INT_RK45 );

    RealTimeAlgorithm algo1(ocp1, 0.05);
    algo1.set( INTEGRATOR_TYPE, INT_RK78 );
    algo1.set( INTEGRATOR_TOLERANCE, 1.000000E-08 );
    algo1.set( ABSOLUTE_TOLERANCE, 1.000000E-08 );
    algo1.set( MAX_NUM_ITERATIONS, 5 );

    StaticReferenceTrajectory referencetrajectory(acadodata_M1);
    Controller controller3( algo1,referencetrajectory );

    SimulationEnvironment algo2(0, 6, process2, controller3);
     algo2.init(acadodata_v3);
    returnValue returnvalue = algo2.run();


    VariablesGrid out_processout; 
    VariablesGrid out_feedbackcontrol; 
    VariablesGrid out_feedbackparameter; 
    VariablesGrid out_states; 
    VariablesGrid out_algstates; 
    algo2.getSampledProcessOutput(out_processout);
    algo2.getProcessDifferentialStates(out_states);
    algo2.getFeedbackControl(out_feedbackcontrol);
    out_processout.print( "nmpc_OUT_states_sampled.m","STATES_SAMPLED",PS_MATLAB ); 
    out_feedbackcontrol.print( "nmpc_OUT_controls.m","CONTROLS",PS_MATLAB ); 
    out_feedbackparameter.print( "nmpc_OUT_parameters.m","PARAMETERS",PS_MATLAB ); 
    out_states.print( "nmpc_OUT_states.m","STATES",PS_MATLAB ); 
    out_algstates.print( "nmpc_OUT_algebraicstates.m","ALGEBRAICSTATES",PS_MATLAB ); 
    const char* outputFieldNames[] = {"STATES_SAMPLED", "CONTROLS", "PARAMETERS", "STATES", "ALGEBRAICSTATES", "CONVERGENCE_ACHIEVED"}; 
    plhs[0] = mxCreateStructMatrix( 1,1,6,outputFieldNames ); 
    mxArray *OutSS = NULL;
    double  *outSS = NULL;
    OutSS = mxCreateDoubleMatrix( out_processout.getNumPoints(),1+out_processout.getNumValues(),mxREAL ); 
    outSS = mxGetPr( OutSS );
    for( int i=0; i<out_processout.getNumPoints(); ++i ){ 
      outSS[0*out_processout.getNumPoints() + i] = out_processout.getTime(i); 
      for( int j=0; j<out_processout.getNumValues(); ++j ){ 
        outSS[(1+j)*out_processout.getNumPoints() + i] = out_processout(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"STATES_SAMPLED",OutSS );
    mxArray *OutS = NULL;
    double  *outS = NULL;
    OutS = mxCreateDoubleMatrix( out_states.getNumPoints(),1+out_states.getNumValues(),mxREAL ); 
    outS = mxGetPr( OutS );
    for( int i=0; i<out_states.getNumPoints(); ++i ){ 
      outS[0*out_states.getNumPoints() + i] = out_states.getTime(i); 
      for( int j=0; j<out_states.getNumValues(); ++j ){ 
        outS[(1+j)*out_states.getNumPoints() + i] = out_states(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"STATES",OutS );
    mxArray *OutC = NULL;
    double  *outC = NULL;
    OutC = mxCreateDoubleMatrix( out_feedbackcontrol.getNumPoints(),1+out_feedbackcontrol.getNumValues(),mxREAL ); 
    outC = mxGetPr( OutC );
    for( int i=0; i<out_feedbackcontrol.getNumPoints(); ++i ){ 
      outC[0*out_feedbackcontrol.getNumPoints() + i] = out_feedbackcontrol.getTime(i); 
      for( int j=0; j<out_feedbackcontrol.getNumValues(); ++j ){ 
        outC[(1+j)*out_feedbackcontrol.getNumPoints() + i] = out_feedbackcontrol(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"CONTROLS",OutC );
    mxArray *OutP = NULL;
    double  *outP = NULL;
    OutP = mxCreateDoubleMatrix( out_feedbackparameter.getNumPoints(),1+out_feedbackparameter.getNumValues(),mxREAL ); 
    outP = mxGetPr( OutP );
    for( int i=0; i<out_feedbackparameter.getNumPoints(); ++i ){ 
      outP[0*out_feedbackparameter.getNumPoints() + i] = out_feedbackparameter.getTime(i); 
      for( int j=0; j<out_feedbackparameter.getNumValues(); ++j ){ 
        outP[(1+j)*out_feedbackparameter.getNumPoints() + i] = out_feedbackparameter(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"PARAMETERS",OutP );
    mxArray *OutZ = NULL;
    double  *outZ = NULL;
    OutZ = mxCreateDoubleMatrix( out_algstates.getNumPoints(),1+out_algstates.getNumValues(),mxREAL ); 
    outZ = mxGetPr( OutZ );
    for( int i=0; i<out_algstates.getNumPoints(); ++i ){ 
      outZ[0*out_algstates.getNumPoints() + i] = out_algstates.getTime(i); 
      for( int j=0; j<out_algstates.getNumValues(); ++j ){ 
        outZ[(1+j)*out_algstates.getNumPoints() + i] = out_algstates(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"ALGEBRAICSTATES",OutZ );
    mxArray *OutConv = NULL;
    if ( returnvalue == SUCCESSFUL_RETURN ) { OutConv = mxCreateDoubleScalar( 1 ); }else{ OutConv = mxCreateDoubleScalar( 0 ); } 
    mxSetField( plhs[0],0,"CONVERGENCE_ACHIEVED",OutConv );


    clearAllStaticCounters( ); 
 
} 

