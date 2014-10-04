/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *  This program includes the functions relevant for hexagonal template placement.
 *
 *  If you edit this code, add the date and your name (also your e-mail address)
 *  to the following list and describe the difference from the current version.
 *
 *  Creation Date and Author: 
 *  2014-06-26 ; First version by Koh Ueno (ueno@vega.ess.sci.osaka-u.ac.jp)
 *
 */

#include <kagali/KGLStdlib.h>
#include <kagali/nrutil.h>
#include <kagali/KGLParameters.h>


void KGLReprodHexagonGrids( //begin{proto}
    int    imode,    /**< (1 or 2) choose one of the two eigen vectors */
    double mm,       /**< minimal match */
    double lambda1,  /**< eigen value 1 of the metric */
    double lambda2,  /**< eigen value 2 of the metric */
    double theta,    /**< inclination angle of ellipse */
    double tau0,     /**< mass parameter which is used to define the metric */
    double tau3,     /**< mass parameter which is used to define the metric */
    double *tau0d,   /**< (out) mass parameter grids which are newly generated */
    double *tau3d    /**< (out) mass parameter grids which are newly generated */
    ) //end{proto}
{
  double cost,sint;
  double cost2,sint2;
  double dX1,dX2;
  double dtau0,dtau3;
  cost=cos(theta);
  sint=sin(theta);
  cost2=cost*cost;
  sint2=sint*sint;
  dX1=sqrt((1-mm)/lambda1);
  dX2=sqrt((1-mm)/lambda2);
  
  if(imode==1){
      dtau0=sint*dX2;
      dtau3=cost*dX2;
      tau0d[1]=tau0+sqrt(3.)*dtau0;
      tau3d[1]=tau3+sqrt(3.)*dtau3;
      tau0d[2]=tau0-sqrt(3.)*dtau0;
      tau3d[2]=tau3-sqrt(3.)*dtau3;
      
      dtau0=dX1*cost*sqrt(3.)/2.+dX2*sint/2.;
      dtau3=-dX1*sint*sqrt(3.)/2.+dX2*cost/2.;
      tau0d[3]=tau0+sqrt(3.)*dtau0;
      tau3d[3]=tau3+sqrt(3.)*dtau3;
      tau0d[4]=tau0-sqrt(3.)*dtau0;
      tau3d[4]=tau3-sqrt(3.)*dtau3;
      
      dtau0=dX1*cost*sqrt(3.)/2.-dX2*sint/2.;
      dtau3=-dX1*sint*sqrt(3.)/2.-dX2*cost/2.;
      tau0d[5]=tau0+sqrt(3.)*dtau0;
      tau3d[5]=tau3+sqrt(3.)*dtau3;
      tau0d[6]=tau0-sqrt(3.)*dtau0;
      tau3d[6]=tau3-sqrt(3.)*dtau3;
  
  }else if(imode==2){
      dtau0=cost*dX1;
      dtau3=-sint*dX1;
      tau0d[1]=tau0+sqrt(3.)*dtau0;
      tau3d[1]=tau3+sqrt(3.)*dtau3;
      tau0d[2]=tau0-sqrt(3.)*dtau0;
      tau3d[2]=tau3-sqrt(3.)*dtau3;
      
      dtau0=dX1*cost/2.+dX2*sint*sqrt(3.)/2.;
      dtau3=-dX1*sint/2.+dX2*cost*sqrt(3.)/2.;
      tau0d[3]=tau0+sqrt(3.)*dtau0;
      tau3d[3]=tau3+sqrt(3.)*dtau3;
      tau0d[4]=tau0-sqrt(3.)*dtau0;
      tau3d[4]=tau3-sqrt(3.)*dtau3;
      
      dtau0=-dX1*cost/2.+dX2*sint*sqrt(3.)/2.;
      dtau3=dX1*sint/2.+dX2*cost*sqrt(3.)/2.;
      tau0d[5]=tau0+sqrt(3.)*dtau0;
      tau3d[5]=tau3+sqrt(3.)*dtau3;
      tau0d[6]=tau0-sqrt(3.)*dtau0;
      tau3d[6]=tau3-sqrt(3.)*dtau3;    

  }else{
    fprintf(stderr,"HexagonGrids: imode should be 1 or 2.");
    exit(EXIT_FAILURE);
  }

  return;
}


void KGLReprodHexagonVertices( //begin{proto}
    int    imode,    /**< (1 or 2) choose one of the two eigen vectors */
    double mm,       /**< minimal match */
    double lambda1,  /**< eigen value 1 of the metric */
    double lambda2,  /**< eigen value 2 of the metric */
    double theta,    /**< inclination angle of ellipse */
    double tau0,     /**< mass parameter which is used to define the metric */
    double tau3,     /**< mass parameter which is used to define the metric */
    double *tau0v,   /**< (out) Hexagon vertices around a given (tau0,tau3) */
    double *tau3v    /**< (out) Hexagon vertices around a given (tau0,tau3) */
    ) //end{proto}
{
  double cost,sint;
  double cost2,sint2;
  double dX1,dX2;
  double dtau0,dtau3;
  cost=cos(theta);
  sint=sin(theta);
  cost2=cost*cost;
  sint2=sint*sint;
  dX1=sqrt((1-mm)/lambda1);
  dX2=sqrt((1-mm)/lambda2);
  
  if(imode==1){
    dtau0=dX1*cost/2.+dX2*sint*sqrt(3.)/2.;
    dtau3=-dX1*sint/2.+dX2*cost*sqrt(3.)/2.;
    tau0v[1]=tau0+dtau0;
    tau3v[1]=tau3+dtau3;
    
    dtau0=-dX1*cost/2.+dX2*sint*sqrt(3.)/2.;
    dtau3=dX1*sint/2.+dX2*cost*sqrt(3.)/2.;
    tau0v[2]=tau0+dtau0;
    tau3v[2]=tau3+dtau3;
    
    dtau0=dX1*cost/2.-dX2*sint*sqrt(3.)/2.;
    dtau3=-dX1*sint/2.-dX2*cost*sqrt(3.)/2.;
    tau0v[3]=tau0+dtau0;
    tau3v[3]=tau3+dtau3;
    
    dtau0=-dX1*cost/2.-dX2*sint*sqrt(3.)/2.;
    dtau3=dX1*sint/2.-dX2*cost*sqrt(3.)/2.;
    tau0v[4]=tau0+dtau0;
    tau3v[4]=tau3+dtau3;
    
    dtau0=dX1*cost;
    dtau3=-dX1*sint;
    tau0v[5]=tau0+dtau0;
    tau3v[5]=tau3+dtau3;
    tau0v[6]=tau0-dtau0;
    tau3v[6]=tau3-dtau3;
    
  }else if(imode==2){
    dtau0=dX1*cost*sqrt(3.)/2.+dX2*sint/2.;
    dtau3=-dX1*sint*sqrt(3.)/2.+dX2*cost/2.;
    tau0v[1]=tau0+dtau0;
    tau3v[1]=tau3+dtau3;
    
    dtau0=-dX1*cost*sqrt(3.)/2.+dX2*sint/2.;
    dtau3=dX1*sint*sqrt(3.)/2.+dX2*cost/2.;
    tau0v[2]=tau0+dtau0;
    tau3v[2]=tau3+dtau3;
    
    dtau0=dX1*cost*sqrt(3.)/2.-dX2*sint/2.;
    dtau3=-dX1*sint*sqrt(3.)/2.-dX2*cost/2.;
    tau0v[3]=tau0+dtau0;
    tau3v[3]=tau3+dtau3;
    
    dtau0=-dX1*cost*sqrt(3.)/2.-dX2*sint/2.;
    dtau3=dX1*sint*sqrt(3.)/2.-dX2*cost/2.;
    tau0v[4]=tau0+dtau0;
    tau3v[4]=tau3+dtau3;
    
    dtau0=dX2*sint;
    dtau3=dX2*cost;
    tau0v[5]=tau0+dtau0;
    tau3v[5]=tau3+dtau3;
    tau0v[6]=tau0-dtau0;
    tau3v[6]=tau3-dtau3;
    
  }else{
    fprintf(stderr,"HexagonVertices: imode should be 1 or 2.");
    exit(EXIT_FAILURE);
  }
  
  return;
}
