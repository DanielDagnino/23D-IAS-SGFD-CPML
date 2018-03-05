/*************************************************************************/
// Authors of the sunmo modification: Daniel Dagnino.
// Authors of the sunmo: See bellow.
/*************************************************************************/
// This modification allows to make a snmo for a difference between to velocities V_{NMO}
// It is designed to be used in a Fortran main program.

// tr input : The shotgather calculated with the velocity V_0
//    output: The shotgather that would correspond to the velocity V_p
// In other words, if V_p > V_0 => The TT arrive earlier (shotgather compressed).
//                    V_p < V_0 => The TT arrive later (shotgather expanded).
/*************************************************************************/

/*************************************************************************/
/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUNMO: $Revision: 1.31 $ ; $Date: 2013/03/06 20:35:27 $    */

/* Credits:
 *  SEP: Shuki Ronen, Chuck Sword
 *  CWP: Shuki Ronen, Jack, Dave Hale, Bjoern Rommel
 *      Modified: 08/08/98 - Carlos E. Theodoro - option for lateral offset
 *      Modified: 07/11/02 - Sang-yong Suh -
 *    added "upward" option to handle decreasing velocity function.
 *      CWP: Sept 2010: John Stockwell
 *    1. replaced Carlos Theodoro's fix 
 *    2. added  the instruction in the selfdoc to use suazimuth to set 
 *        offset so that it accounts for lateral offset. 
 *        3. removed  Bjoren Rommel's anisotropy stuff. sunmo_a is the 
 *           version with the anisotropy parameters left in.
 *        4. note that scalel does not scale the offset field in
 *           the segy standard.
 * Technical Reference:
 *  The Common Depth Point Stack
 *  William A. Schneider
 *  Proc. IEEE, v. 72, n. 10, p. 1238-1254
 *  1984
 *
 * Trace header fields accessed: ns, dt, delrt, offset, cdp, scalel
 */

#include <stdio.h>
#include <string.h>

#include "./su.h"
#include "./segy.h"
// #include "./cwp.h"

// #pragma GCC diagnostic ignored "-Wwrite-strings"

// // define MARK which you can place throughout your code to find what line 
// // it is crashing on
// #ifndef MARK
// #define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
// #endif

char *sdoc[] = {" ",NULL};

/*********************** self documentation ******************************/
void nmo_dag( void *ptr, int nt, int ntr, void *pvoffset, void *ptnmo1, void *ptnmo2, void *pvnmo1, void *pvnmo2, int ntnmo, float dt, int invert ){
  float *ovvt1=NULL;  /* array[nt] of sloth for a particular trace */
  float *ovvt2=NULL;  /* array[nt] of sloth for a particular trace */
  float tn;    /* NMO time (time after NMO correction) */
  float *qtn=NULL;  /* NMO-corrected trace q(tn) */
  float *ttn=NULL;  /* time t(tn) for NMO */
  float *atn=NULL;  /* amplitude a(tn) for NMO */
  float *qt=NULL;    /* inverse NMO-corrected trace q(t) */
  float *tnt=NULL;  /* time tn(t) for inverse NMO */
  float *at=NULL;    /* amplitude a(t) for inverse NMO */
  float temp;    /* temporary float */
  float tsq;    /* temporary float */
  float offset;    /* value of offset honoring scalel */
  
  int it, itr, itmute;
  float lmute, smute, osmute;
  float v, tsq1, tsq2;
  float **tr_aux;
  
  double *tr;
  float *tnmo1, *tnmo2, *vnmo1, *vnmo2, *voffset;
  
//   FILE *fileID;
  
//   //----------------------------------------------------------//
//   fileID = fopen("hola.txt","w");
//   fprintf(fileID,"C HOLA 0\n");  
//   MARK
  
  tr = (double*) ptr;
  voffset = (float*) pvoffset;
  tnmo1 = (float*) ptnmo1;
  vnmo1 = (float*) pvnmo1;
  tnmo2 = (float*) ptnmo2;
  vnmo2 = (float*) pvnmo2;
  
  //----------------------------------------------------------//
//   MARK
  
  tr_aux = ealloc2float(nt,ntr);
  
  for (it=0;it<nt;it++){
  for (itr=0;itr<ntr;itr++){
    tr_aux[itr][it] = (float) tr[itr*nt+it];
    }
    }
  
  //----------------------------------------------------------//
//   MARK
  
  /* allocate workspace */
  ovvt1 = ealloc1float(nt);
  ovvt2 = ealloc1float(nt);
  
  for (it=0,tn=0.; it<nt; ++it,tn+=dt){
    intlin(ntnmo,tnmo1,vnmo1,vnmo1[0],vnmo1[ntnmo-1],1,&tn,&v);
    ovvt1[it] = 1./(v*v);
    }
  
  for (it=0,tn=0.; it<nt; ++it,tn+=dt){
    intlin(ntnmo,tnmo2,vnmo2,vnmo2[0],vnmo2[ntnmo-1],1,&tn,&v);
    ovvt2[it] = 1./(v*v);
    }
  
  //----------------------------------------------------------//
//   MARK
  
  lmute = 25.;
  smute = 1.5;
  
  /* allocate workspace */
  ttn = ealloc1float(nt);
  atn = ealloc1float(nt);
  qtn = ealloc1float(nt);
  tnt = ealloc1float(nt);
  at  = ealloc1float(nt);
  qt  = ealloc1float(nt);
  
  for (itr=0; itr<ntr; ++itr){
    
    /* compute time t(tn) (normalized) */
    offset = voffset[itr];
    temp = ((float) offset*offset)/(dt*dt);
    for (it=0,tn=0.; it<nt; ++it,tn+=1.){
      tsq1 = temp*ovvt1[it];
      tsq2 = temp*ovvt2[it];
      ttn[it] = it - ( sqrt( tn*tn + tsq2 ) - sqrt( tn*tn + tsq1 ) );
      }
    
    /* compute inverse of stretch factor a(tn) */
    atn[0] = ttn[1]-ttn[0];
    for (it=1; it<nt; ++it) atn[it] = ttn[it]-ttn[it-1];
    
    /* determine index of first sample to survive mute */
    osmute = 1./smute;
    for (it=0; it<nt-1 && atn[it]<osmute; ++it);
    itmute = it;
    
    //----------------------------------------------------------//
    /* if inverse NMO will be performed */
    if (invert) {
      
      /* compute tn(t) from t(tn) */
//   yxtoxy    compute a regularly sampled function x(y) from a regularly
//             sampled, monotonically increasing function y(x)
//   void yxtoxy( int nx, float dx, float fx, float y[], 
//        int ny, float dy, float fy, float xylo, float xyhi, float x[]);
      
      yxtoxy( nt-itmute, 1., itmute, &ttn[itmute], nt-itmute, 1., itmute, -nt, nt, &tnt[itmute] );
      
      /* adjust mute time */
      itmute = 1. + ttn[itmute];
      itmute = MIN(nt-2,itmute);
      
      /* compute a(t) */
      for (it=itmute+1; it<nt; ++it) at[it] = tnt[it]-tnt[it-1];
      at[itmute] = at[itmute+1];
      
    }
    
    /* if forward (not inverse) nmo */
    if (!invert) {
      
      /* do nmo via 8-point sinc interpolation */
      ints8r(nt,1.,0.,tr_aux[itr],0.,0.,nt-itmute,&ttn[itmute],&qtn[itmute]);
      
      /* apply mute */
      for (it=0; it<itmute; ++it) qtn[it] = 0.;
      
      /* apply linear ramp */
      for (it=itmute; it<itmute+lmute && it<nt; ++it) qtn[it] *= (float)(it-itmute+1)/(float)lmute;
      
      /* if specified, scale by the NMO stretch factor */
      for (it=0; it<nt; ++it) qtn[it] *= atn[it];
      
      /* copy NMO corrected trace to output trace */
      memcpy( (void *) tr_aux[itr], (const void *) qtn, nt*sizeof(float) );
      
    /* else inverse nmo */
    } else {
      
      /* do inverse nmo via 8-point sinc interpolation */
      ints8r(nt,1.,0.,tr_aux[itr],0.,0.,nt-itmute,&tnt[itmute],&qt[itmute]);
      
      /* apply mute */
      for (it=0; it<itmute; ++it) qt[it] = 0.;
      
      /* if specified, undo NMO stretch factor scaling */
      for (it=itmute; it<nt; ++it) qt[it] *= at[it];
      
      /* copy inverse NMO corrected trace to output trace */
      memcpy( (void *) tr_aux[itr], (const void *) qt, nt*sizeof(float) );
      
    }
  
  }
  
  for (itr=0; itr<ntr; ++itr){
  for (it=0; it<nt; ++it){
    tr[itr*nt+it] = (double) tr_aux[itr][it];
    }
    }
  
//   fclose(fileID);
  
  //----------------------------------------------------------//
  //----------------------------------------------------------//
  //----------------------------------------------------------//
//   MARK
  free2float(tr_aux);
  
//   MARK
  free1float(ovvt1);
//   MARK
  free1float(ovvt2);
  
//   MARK
  free1float(ttn);
//   MARK
  free1float(atn);
//   MARK
  free1float(qtn);
//   MARK
  free1float(tnt);
//   MARK
  free1float(at);
//   MARK
  free1float(qt);
  
//   MARK
  
  return;
  }


