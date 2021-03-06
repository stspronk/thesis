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
    DifferentialState phi;
    DifferentialState theta;
    DifferentialState psi;
    DifferentialState THRUST;
    Control delta_phi;
    Control delta_theta;
    Control delta_psi;
    Control delta_THRUST;
    Parameter Time; 
    DifferentialEquation acadodata_f1(0, Time);
    acadodata_f1 << dot(x) == (cos(psi)*v_x-sin(psi)*v_y);
    acadodata_f1 << dot(y) == (cos(psi)*v_y+sin(psi)*v_x);
    acadodata_f1 << dot(z) == v_z;
    acadodata_f1 << dot(v_x) == (1/5.00000000000000000000e-01*THRUST*sin(theta)-5.00000000000000000000e-01*v_x);
    acadodata_f1 << dot(v_y) == ((-sin(phi))/5.00000000000000000000e-01*THRUST-5.00000000000000000000e-01*v_y);
    acadodata_f1 << dot(v_z) == ((4.90500000000000024869e+00+THRUST*cos(phi)*cos(theta))/5.00000000000000000000e-01-5.00000000000000000000e-01*v_z);
    acadodata_f1 << dot(phi) == delta_phi;
    acadodata_f1 << dot(theta) == delta_theta;
    acadodata_f1 << dot(psi) == delta_psi;
    acadodata_f1 << dot(THRUST) == delta_THRUST;

    OCP ocp1(0, Time, 50);
    ocp1.minimizeMayerTerm(Time);
    ocp1.subjectTo(acadodata_f1);
    ocp1.subjectTo(AT_START, x == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, y == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, z == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, v_x == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, v_y == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, v_z == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, phi == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, theta == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, psi == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_START, THRUST == (-4.91600000000000036948e+00));
    ocp1.subjectTo(AT_END, x == 3.00000000000000000000e+00);
    ocp1.subjectTo(AT_END, y == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_END, z == (-1.50000000000000000000e+00));
    ocp1.subjectTo(AT_END, v_x == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_END, v_y == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_END, v_z == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_END, phi == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_END, theta == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_END, psi == 0.00000000000000000000e+00);
    ocp1.subjectTo(AT_END, THRUST == (-4.91600000000000036948e+00));
    ocp1.subjectTo((-1.60000000000000008882e+00) <= z <= 0.00000000000000000000e+00);
    ocp1.subjectTo((-4.36332312998582383390e-01) <= phi <= 4.36332312998582383390e-01);
    ocp1.subjectTo((-4.36332312998582383390e-01) <= theta <= 4.36332312998582383390e-01);
    ocp1.subjectTo((-4.36332312998582383390e-01) <= psi <= 4.36332312998582383390e-01);
    ocp1.subjectTo((-2.30000000000000000000e+01) <= delta_phi <= 2.30000000000000000000e+01);
    ocp1.subjectTo((-2.30000000000000000000e+01) <= delta_theta <= 2.30000000000000000000e+01);
    ocp1.subjectTo((-2.30000000000000000000e+01) <= delta_psi <= 2.30000000000000000000e+01);
    ocp1.subjectTo((-4.00000000000000000000e+01) <= delta_THRUST <= 4.00000000000000000000e+01);
    ocp1.subjectTo((-1.96209999999999986642e+01) <= THRUST <= 0.00000000000000000000e+00);
    ocp1.subjectTo(0.00000000000000000000e+00 <= Time <= 5.00000000000000000000e+00);


    OptimizationAlgorithm algo1(ocp1);
    algo1.set( KKT_TOLERANCE, 1.000000E-10 );
    algo1.set( MAX_NUM_ITERATIONS, 400 );
    returnValue returnvalue = algo1.solve();

    VariablesGrid out_states; 
    VariablesGrid out_parameters; 
    VariablesGrid out_controls; 
    VariablesGrid out_disturbances; 
    VariablesGrid out_algstates; 
    algo1.getDifferentialStates(out_states);
    algo1.getControls(out_controls);
    algo1.getParameters(out_parameters);
    const char* outputFieldNames[] = {"STATES", "CONTROLS", "PARAMETERS", "DISTURBANCES", "ALGEBRAICSTATES", "CONVERGENCE_ACHIEVED"}; 
    plhs[0] = mxCreateStructMatrix( 1,1,6,outputFieldNames ); 
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
    OutC = mxCreateDoubleMatrix( out_controls.getNumPoints(),1+out_controls.getNumValues(),mxREAL ); 
    outC = mxGetPr( OutC );
    for( int i=0; i<out_controls.getNumPoints(); ++i ){ 
      outC[0*out_controls.getNumPoints() + i] = out_controls.getTime(i); 
      for( int j=0; j<out_controls.getNumValues(); ++j ){ 
        outC[(1+j)*out_controls.getNumPoints() + i] = out_controls(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"CONTROLS",OutC );
    mxArray *OutP = NULL;
    double  *outP = NULL;
    OutP = mxCreateDoubleMatrix( out_parameters.getNumPoints(),1+out_parameters.getNumValues(),mxREAL ); 
    outP = mxGetPr( OutP );
    for( int i=0; i<out_parameters.getNumPoints(); ++i ){ 
      outP[0*out_parameters.getNumPoints() + i] = out_parameters.getTime(i); 
      for( int j=0; j<out_parameters.getNumValues(); ++j ){ 
        outP[(1+j)*out_parameters.getNumPoints() + i] = out_parameters(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"PARAMETERS",OutP );
    mxArray *OutW = NULL;
    double  *outW = NULL;
    OutW = mxCreateDoubleMatrix( out_disturbances.getNumPoints(),1+out_disturbances.getNumValues(),mxREAL ); 
    outW = mxGetPr( OutW );
    for( int i=0; i<out_disturbances.getNumPoints(); ++i ){ 
      outW[0*out_disturbances.getNumPoints() + i] = out_disturbances.getTime(i); 
      for( int j=0; j<out_disturbances.getNumValues(); ++j ){ 
        outW[(1+j)*out_disturbances.getNumPoints() + i] = out_disturbances(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"DISTURBANCES",OutW );
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

