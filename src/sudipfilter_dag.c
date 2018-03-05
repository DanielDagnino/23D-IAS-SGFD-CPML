/*************************************************************************/
// Authors of the sunmo modification: Daniel Dagnino.
// Authors of the sunmo: See bellow.
/*************************************************************************/
// 
/*************************************************************************/

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUDIPFILT: $Revision: 1.20 $ ; $Date: 2011/11/12 00:09:00 $    */

#include "./su.h"
#include "./segy.h"
#include "./header.h"

#include <signal.h>

// #pragma GCC diagnostic ignored "-Wwrite-strings"

// // define MARK which you can place throughout your code to find what line 
// // it is crashing on
// #ifndef MARK
// #define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
// #endif

// char *sdoc[] = {" ",NULL};

/* Credits:
 *
 *  CWP: Dave (algorithm--originally called slopef)
 *       Jack (reformatting for SU)
 *
 * Trace header fields accessed: ns, dt, d2
 */

/*********************** self documentation ******************************/

/* prototypes */
void slopefilter (int nslopes, float *slopes, float *amps, float bias,
  int nt, float dt, int ntr, float dx, float **tr_aux);

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
void dipfilter_dag( void *ptr, float dx, float dt, int nt, int ntr, int nslopes, void *pslopes, void *pamps ){
  double *tr;
//   float *slopes;  /* slopes at which amplitudes are specified */
//   float *amps;  /* amplitudes corresponding to slopes   */
  float **tr_aux;
  int it, itr;
  float bias=0.000667;
  
  //----------------------------------------------------------//
  tr = (double*) ptr;
//   slopes = (float*) pslopes;
//   amps = (float*) pamps;
  
  //----------------------------------------------------------//
  tr_aux = ealloc2float(nt,ntr);
  
  // Save tr in auxiliar array.
  for (it=0;it<nt;it++)
    for (itr=0;itr<ntr;itr++) tr_aux[itr][it] = (float) tr[it*ntr+itr];
  
  //----------------------------------------------------------//
  /* Apply slope filter */
  slopefilter(nslopes,pslopes,pamps,bias,nt,dt,ntr,dx,tr_aux);
  
  //----------------------------------------------------------//
  /* Output filtered traces */
  for (itr=0; itr<ntr; ++itr)
    for (it=0; it<nt; ++it) tr[it*ntr+itr] = (float) tr_aux[itr][it];
  
  //----------------------------------------------------------//
  free2float(tr_aux);
  
  return;
}



/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
void slopefilter (int nslopes, float *slopes, float *amps, float bias,
  int nt, float dt, int nx, float dx, float **tr_aux){
/******************************************************************************
apply slope filter in frequency-wavenumber domain
*******************************************************************************
Input:
nslopes   number of slopes (and amplitudes) specified
slopes    slopes at which amplitudes are specified (see notes below)
amps    amplitudes corresponding to slopes (see notes below)
bias    linear moveout slope before and after filtering
nt    number of time samples
dt    time sampling interval
nx    number of traces
dx    trace space (spatial sampling interval)
tr_aux   file pointer to data to be filtered

Output:
tr_aux   file pointer to filtered data
*******************************************************************************
Notes:
Linear interpolation and constant extrapolation are used to
determine amplitudes for slopes that are not specified.
******************************************************************************/
  int ntfft;    /* nt after padding for FFT */
  int nxfft;    /* nx after padding for FFT */
  float sfft;   /* scale factor for FFT */
  int nw;     /* number of frequencies */
  float dw;   /* frequency sampling interval */
  float fw;   /* first frequency */
  int nk;     /* number of wavenumbers */
  float dk;   /* wavenumber sampling interval */
  float w,k;    /* frequency and wavenumber */
  int it,ix,iw,ik;  /* sample indices */
  float slope,amp;  /* slope and amplitude for particular w,k */
  complex **cpfft;  /* complex FFT workspace */
  float **pfft;   /* float FFT workspace */
  float phase;    /* phase shift for bias */
  complex cshift;   /* complex phase shifter for bias */
  
  /* determine lengths and scale factors for prime-factor FFTs */
  ntfft = npfar(nt);
  nxfft = npfa(nx);
  sfft = 1.0/(ntfft*nxfft);
  
  /* determine frequency and wavenumber sampling */
  nw = ntfft/2+1;
  dw = 2.0*PI/(((float) ntfft)*dt);
  fw = 0.000001*dw; /* non-zero to avoid divide by zero w */
  nk = nxfft;
  dk = 2.0*PI/(((float) nxfft)*dx);
  
  /* allocate real and complex workspace for FFTs */
  cpfft = alloc2complex(nw,nk);
  pfft = alloc2float(ntfft,nxfft);
  
  
  
  
  
  /* DAGNINO start: To avoid non-initialized values and condition jumps it is necessary to zero out. */
  for (ix=0; ix<nk; ix++)
    for (it=0; it<nw; it++){
      cpfft[ix][it].r = 0.;
      cpfft[ix][it].i = 0.;
    }
  
  for (ix=0; ix<nxfft; ix++)
    for (it=0; it<ntfft; it++) pfft[ix][it] = 0.;
  /* DAGNINO end. */
  
  
  
  
  
  /* copy data from input to FFT array and pad with zeros */
  for (ix=0; ix<nx; ix++) {
    for (it=0; it<nt; it++) pfft[ix][it] = tr_aux[ix][it];
    for (it=nt; it<ntfft; it++) pfft[ix][it] = 0.0;
    }
  for (ix=nx; ix<nxfft; ix++)
    for (it=0; it<ntfft; it++) pfft[ix][it] = 0.0;
  
  /* Fourier transform t to w */
  pfa2rc(1,1,ntfft,nx,pfft[0],cpfft[0]);
  
  /* do linear moveout bias via phase shift */
  for (ix=0; ix<nx; ix++) {
    for (iw=0,w=0.0; iw<nw; iw++,w+=dw) {
      phase = -((float) ix)*dx*w*bias;
      cshift = cmplx(cos(phase),sin(phase));
      cpfft[ix][iw] = cmul(cpfft[ix][iw],cshift);
    }
  }
  
  /* Fourier transform x to k */
  pfa2cc(-1,2,nw,nxfft,cpfft[0]);
  
  /* loop over wavenumbers */
  for (ik=0; ik<nk; ik++) {
  
    /* determine wavenumber */
    k = (ik<=nk/2) ? ((float) ik)*dk : ((float) (ik-nk))*dk;
    
    /* loop over frequencies */
    for (iw=0,w=fw; iw<nw; iw++,w+=dw) {
    
      /* determine biased slope */
      slope = k/w+bias;
      
      /* linearly interpolate to find amplitude */
      intlin(nslopes,slopes,amps,amps[0],amps[nslopes-1],1,&slope,&amp);
      
      /* include fft scaling */
      amp *= sfft;
      
      /* filter real and imaginary parts */
      cpfft[ik][iw].r *= amp;
      cpfft[ik][iw].i *= amp;
    }
  }
  
  /* Fourier transform k to x */
  pfa2cc(1,2,nw,nxfft,cpfft[0]);
  
  /* undo linear moveout bias via phase shift */
  for (ix=0; ix<nx; ix++) {
    for (iw=0,w=0.0; iw<nw; iw++,w+=dw) {
      phase = ((float) ix)*dx*w*bias;
      cshift = cmplx(cos(phase),sin(phase));
      cpfft[ix][iw] = cmul(cpfft[ix][iw],cshift);
    }
  }
  
  /* Fourier transform w to t */
  pfa2cr(-1,1,ntfft,nx,cpfft[0],pfft[0]);
  
  /* copy filtered data from FFT array to output */
  for (ix=0; ix<nx; ix++)
    for (it=0; it<nt; it++) tr_aux[ix][it] = pfft[ix][it];
  
  /* free workspace */
  free2complex(cpfft);
  free2float(pfft);
}



