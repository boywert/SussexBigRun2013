#ifndef KERNEL_H
#define KERNEL_H

#ifdef QUINTIC_KERNEL

#if !defined(TWODIMS) && !defined(ONEDIM)
#define  NORM 2187.0/(40.0*M_PI)	/*!< For 3D-normalized kernel */
#else
#ifndef  ONEDIM
#define  NORM 15309.0/(478.0*M_PI)	/*!< For 2D-normalized kernel */
#else
#define  NORM 243.0/40.0        	/*!< For 1D-normalized kernel */
#endif
#endif

#else /* of QUINTIC_KERNEL, so here comes cubic kernel */

#if !defined(TWODIMS) && !defined(ONEDIM)
#define  NORM 8.0/M_PI                   /*!< For 3D-normalized kernel */
#else
#ifndef  ONEDIM
#define  NORM 40.0/(7.0*M_PI)	         /*!< For 2D-normalized kernel */
#else
#define  NORM  4.0/3.0        	         /*!< For 1D-normalized kernel */
#endif
#endif

#endif /* QUINTIC_KERNEL */

static inline void kernel_hinv(double h, double *hinv, double *hinv3, double *hinv4)
{
  *hinv = 1.0 / h;
#ifndef  TWODIMS
#ifndef  ONEDIM
  *hinv3 = *hinv * *hinv * *hinv;
#else
  *hinv3 = *hinv;
#endif
#else
  *hinv3 = *hinv * *hinv / boxSize_Z;
#endif

  *hinv4 = *hinv3 * *hinv;
} 

/* Attention: Here we assume that kernel is only called with range 0..1 for u  
   as done in hydra or density !! 
   Call with mode 0 to calculate dwk and wk
   Call with mode -1 to calculate only wk
   Call with mode +1 to calculate only dwk */

static inline void kernel_main(double u, double hinv3, double hinv4, double *wk, double *dwk, int mode)
{
#ifdef QUINTIC_KERNEL
  double t1 = (1.0 - u);
  double t2 = t1 * t1;
  double t4 = t2 * t2;

  if(mode >= 0) *dwk = -5.0 * t4;
  if(mode <= 0) *wk = t4 * t1;

  if (u < 2.0/3.0)
    {
      t1 = (2.0/3.0 - u);
      t2 = t1 * t1;
      t4 = t2 * t2;
      if(mode >= 0) *dwk += 30.0 * t4;
      if(mode <= 0) *wk -= 6.0 * t4 * t1;
    }
  if (u < 1.0/3.0)
    {
      t1 = (1.0/3.0 - u);
      t2 = t1 * t1;
      t4 = t2 * t2;
      if(mode >= 0) *dwk -= 75.0 * t4;
      if(mode <= 0) *wk += 15.0 * t4 * t1;
    }

#else /* of QUINTIC_KERNEL, so here comes cubic kernel */

  if(u < 0.5)
    {
      if(mode >= 0) *dwk = u * (18.0 * u - 12.0);
      if(mode <= 0) *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
  else
    {
      double t1 = (1.0 - u);
      double t2 = t1 * t1;
      if(mode >= 0) *dwk = -6.0 * t2;
      if(mode <= 0) *wk = 2.0 * t2 * t1;
    }
#endif /* QUINTIC_KERNEL */

  if(mode >= 0) *dwk *= NORM * hinv4;
  if(mode <= 0) *wk *= NORM * hinv3;
}

#endif
