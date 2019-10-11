      integer, parameter ::totim=600*120
      INTEGER NZ,NZM1,NZP1,NDIM,NX,NY,NPTS
      PARAMETER(NZ = 39, NZM1 = NZ-1, NZP1 = NZ+1 )
      PARAMETER(NDIM = 1, NX = 1, NY = 1, NPTS = NX*NY)
      INTEGER NVEL,NSCLR,NVP1,NSP1,NSB
      PARAMETER(NVEL =  2 , NSCLR= 2 , NVP1 = NVEL+1,NSP1=NSCLR+1)
      PARAMETER(NSB  =  1 ) !NSCLR-2
c      INTEGER NRECS
c      PARAMETER(NRECS = 1440)
      INTEGER itermax
      REAL hmixtolfrac
      PARAMETER(hmixtolfrac = 1.7          , itermax = 200 )
c     temporary grid
      INTEGER NGRID,NZL,NZU,NZDIVmax,NZtmax,NZP1tmax,igridmax
      PARAMETER(NGRID   = NZ )
      PARAMETER(NZL     =  1 , NZU  = 2    , NZDIVmax = 8 )
      PARAMETER(NZtmax  = NZ +(NZL+NZU)*(NZDIVmax-1), NZP1tmax=NZtmax+1)
      PARAMETER(igridmax= 5)
c fluxes and forcing
      INTEGER NSFLXS,NJDT,NSFLXSM1,NSFLXSP2,NDHARM
      PARAMETER(NSFLXS  =  9 , NJDT = 1    , NSFLXSM1 = NSFLXS-1 ,
     +                                       NSFLXSP2 = NSFLXS+2 )
      PARAMETER(NDHARM  =  5) 
c ocean advection
c      INTEGER maxmodeadv
c      PARAMETER(maxmodeadv = 6)	
c richardson mixing
      INTEGER MR,MRP1
      PARAMETER(MR      =100 , MRP1 = MR+1)
c rad/conv model
*      INTEGER NPLEV,NPSAVE
*      PARAMETER(NPLEV   = 38 )
*      PARAMETER(NPSAVE  = 10 )
c output buffer
*      INTEGER NDOUT,NBUFF
*      PARAMETER(NDOUT   = 10+NZP1)
*      PARAMETER(NBUFF   = NZP1*(NVp1+NSP1) + NZP1*(NVEL+NSCLR)
*     +                  + NDOUT + 3*NZ + 5*NSFLXS )
c    +                  + 2*NPLEV + NPSAVE ) ! to store rad/conv output
c
c  Parameters for regional coupling
      INTEGER NX_GLOBE,NY_GLOBE
      PARAMETER(NX_GLOBE=1, NY_GLOBE=1)
      INTEGER*8 NPTS_GLOBE
      PARAMETER (NPTS_GLOBE= NX_GLOBE*NY_GLOBE)
c-----------------------------------------------------------------------
c Note in this version:
c    -albedo for ocean is set to 0.06, which is used for QSW from fcomp
c        or fread when rad/conv is not running:
c        albocn=0.06                         (in subroutine init cnsts)
c-----------------------------------------------------------------------
c
c main (permanent) grid
c
c     NZ    : number of layers      
c     NZP1  : number of grid points 
c     NDIM  & NX & NY : dimension of the model(not used in this version)
c     NVEL  : number of velocity components, i.e. 2 for U and V
c     NSCLR : number of scalars, i.e. T, S, and additional scalars
c     NSB   : number of biological scalars, used in "biocommon.inc"
c     hmixtolfrac : convergence tolerance for hmix(new)-hmix(old) 
c             iteration in ocnstep: fraction of layer thickness hm(kmix)
c     itermax : maximum number of hmix iterations (on main or temporary
c             grids.
c     NRECS : number of time records in forcing flux file ! bxp 20-Aug-2015
c
c temporary grid
c
c     NGRID : number of grids = permanent grid + temporary grids,
c             if only p-grid used: NGRID = 1, if t-grid used: NGRID = NZ
c     NZL   & NZU : refinement interval: it is defined on permanent-grid
c             from (kmix+NZU) to (kmix-NZL)
c     NZDIVmax: maximum number of fine(temporary) grid intervals
c             per permanent grid interval. It is needed to dimension the
c             necessary arrays; the number to be used is being read in 
c             as NZDIV.
c     NZtmax: maximum number of layers on temporary grid,
c             the number to be used in this run is read in as NZT.
c     igridmax : maximum number of new grid refinements in ocntgrid
c
c fluxes and forcing
c 
c     NSFLXS: number of fluxes: sflux(NSFLXS,5,0:NJDT)
c     NJDT  : number of older flux values used to extrapolate new fluxes
c     NDHARM: maximum number of harmonics used to specify forcing
c
c ocean advection
c     maxmodeadv: maximum number of different modes for advection
c richardson mixing
c
c     MR    : dimension of "F of Ri" in ri_mix routine
c
c rad/conv model
c
c     NPLEV : number of levels in atmospheric model.
c             (Note: NPLEV=plev necessary in "rad.par")
c     NPSAVE: number of extra atmospheric variables to be saved
c
c dimension of output array
c
c     NDOUT : dimension of output array "dout", so that a number of 
c             scalar parameters, as well as one extra profile 
c             (on layer or grid) can be stored for diagnostic purposes.
c     NBUFF : dimension of data array "buffer"
