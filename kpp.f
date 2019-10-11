      subroutine kppmix(
     $                  zgrid , hwide , Shsq, dVsq, 
     $                  ustar , Bo    , Bosol , alphaDT , betaDS ,
     $                  dbloc , Ritop , Coriol, jwtype,
     $                  visc  , difs  , dift  , ghats , hbl , kbl) 

c.......................................................................
c
c     Main driver subroutine for kpp vertical mixing scheme and 
c     interface to greater ocean model
c
c     written by : bill large,   june,  6, 1994
c     modified by: jan morzel,   june, 30, 1994
c                  bill large, august, 11, 1994
c                  bill large, november 1994,   for 1d code
c
c.......................................................................
c     
        include 'parameter.inc'
      parameter (km = NZ, kmp1 = nzp1, imt = NX*NY) !kmp1 =km+1
      parameter (mdiff = 3)  ! number of diffusivities for local arrays
c
c input
      real zgrid(kmp1)   ! vertical grid (<= 0)            (m)
      real hwide(kmp1)   ! layer thicknesses               (m)
      real Shsq(kmp1)    ! (local velocity shear)^2       (m/s)^2
      real dVsq(kmp1)    ! (velocity shear re sfc)^2      (m/s)^2
      real ustar        ! surface friction velocity       (m/s)
      real Bo           ! surface turbulent buoy. forcing (m^2/s^3)
      real Bosol        ! radiative buoyancy forcing      (m^2/s^3)
      real alphaDT(kmp1) ! alpha * DT  across interfaces
      real betaDS(kmp1)  ! beta  * DS  across interfaces
      real dbloc(km)     ! local delta buoyancy            (m/s^2)
      real Ritop(km)     ! numerator of bulk Richardson Number (m/s)^2
c          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
      real Coriol            ! Coriolis parameter              (s^{-1})
      integer jwtype    ! Jerlov water type               (1 -- 5)
c
c output
      real visc(0:kmp1)  ! vertical viscosity coefficient  (m^2/s)
      real difs(0:kmp1)  ! vertical scalar diffusivity     (m^2/s)
      real dift(0:kmp1)  ! vertical temperature diffusivity(m^2/s)
      real ghats(km)     ! nonlocal transport              (s/m^2)
      real hbl          ! boundary layer depth (m)
      integer kbl       ! index of first grid level below hbl
c
c local
      real bfsfc        ! surface buoyancy forcing        (m^2/s^3)
      real ws           ! momentum velocity scale
      real wm           ! scalar   velocity scale 
      real caseA        ! = 1 in case A; =0 in case B
      real stable       ! = 1 in stable forcing; =0 in unstable
      real dkm1(mdiff)   ! boundary layer difs at kbl-1 level
      real gat1(mdiff)   ! shape function at sigma=1
      real dat1(mdiff)   ! derivative of shape function at sigma=1
      real blmc(km,mdiff)! boundary layer mixing coefficients
      real sigma        ! normalized depth (d / hbl)
      real Rib(2)        ! bulk Richardson number


c zero the mixing coefficients 
        do 310 ki=0,km
        visc(ki) = 0.0
        difs(ki) = 0.0
        dift(ki) = 0.0
 310  continue

c compute RI and IW interior diffusivities everywhere
!      IF(LRI) THEN
      call ri_iwmix ( km  , kmp1, imt   ,
     $                Shsq,  dbloc , zgrid ,
     $                visc, difs, dift  )
!      ENDIF

c fill the bottom kmp1 coefficients for blmix
         visc (kmp1) = visc(km)
         difs (kmp1) = difs(km)
         dift (kmp1) = dift(km)
 
c  compute boundary layer mixing coefficients ??
!      IF(LKPP) THEN
c diagnose the new boundary layer depth

      call  bldepth (km   , kmp1 , imt   , zgrid, hwide, dVsq,   
     $               dbloc, Ritop, ustar , Bo   , Bosol, Coriol, jwtype,
     $               hbl  , bfsfc, stable, caseA, kbl  ,
     $               Rib  , sigma, wm    , ws   , ocdepth)

c compute boundary layer diffusivities
      call blmix   (km   , kmp1 , imt , mdiff , zgrid, hwide , 
     $              ustar, bfsfc, hbl , stable, caseA, 
     $              visc , difs , dift, kbl   , 
     $              gat1 , dat1 , dkm1, blmc  , ghats, 
     $              sigma, wm   , ws  )
 
c enhance diffusivity at interface kbl - 1
      call enhance (km   , kmp1 , imt , mdiff , dkm1  , visc ,
     $              difs , dift , hbl , kbl   , zgrid , caseA,
     $              blmc ,ghats )

c combine interior and boundary layer coefficients and nonlocal term
      do ki= 1,km
            if(ki.lt.kbl) then
               visc(ki)=blmc(ki,1)
               difs(ki)=blmc(ki,2)
               dift(ki)=blmc(ki,3)
            else
               ghats(ki)=0.
            endif
      end do

c
c     NPK 25/9/08.  Trap for negative values of diffusivities.
c     If negative, set to a background value of 1E-05.
c
      DO ki= 1,km
            IF (dift(ki) .LT. 0) dift(ki)=1E-05
            IF (difs(ki) .LT. 0) difs(ki)=1E-05
            IF (visc(ki) .LT. 0) visc(ki)=1E-05
      end do
!     ENDIF              ! of LKPP

      return
      end

c ********************************************************************

      subroutine  bldepth (
     $          km   , kmp1 , imt   , zgrid, hwide, dVsq , 
     $          dbloc, Ritop, ustar , Bo   , Bosol, Coriol, jwtype,
     $          hbl  , bfsfc, stable, caseA, kbl  ,
     $          Rib  , sigma, wm    , ws   ,ocdepth)
c
c     the oceanic planetray boundary layer depth, hbl, is determined as
c     the shallowest depth where the bulk richardson number is
c     equal to the critical value, Ricr.
c
c     bulk richardson numbers are evaluated by computing velocity and
c     buoyancy differences between values at zgrid(kl) < 0 and surface
c     reference values.
c     in this configuration, the reference values are equal to the
c     values in the surface layer.  
c     when using a very fine vertical grid, these values should be 
c     computed as the vertical average of velocity and buoyancy from 
c     the surface down to epsilon*zgrid(kl).
c
c     when the bulk richardson number at k exceeds Ricr, hbl is
c     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
c
c     The water column and the surface forcing are diagnosed for 
c     stable/ustable forcing conditions, and where hbl is relative 
c     to grid points (caseA), so that conditional branches can be 
c     avoided in later subroutines.
c
c
c  model  
      include 'times.com'
      integer km,kmp1      ! number of vertical levels
      integer imt          ! number of horizontal grid points
      real zgrid(kmp1) ! vertical grid (<= 0)              (m)
      real hwide(kmp1) ! layer thicknesses                 (m)
c
c  input
      real dVsq(kmp1)  ! (velocity shear re sfc)^2      (m/s)^2
      real dbloc(km)   ! local delta buoyancy              (m/s^2)
      real Ritop(km)   ! numerator of bulk Richardson Number (m/s)^2
c          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
      real ustar      ! surface friction velocity         (m/s)
      real Bo         ! surface turbulent buoyancy forcing(m^2/s^3)
      real Bosol      ! radiative buoyancy forcing        (m^2/s^3)
      real Coriol     ! Coriolis parameter                (1/s)
      integer jwtype  ! Jerlov water type                 (1 to 5)
      real ocdepth

c  output
      real hbl        ! boundary layer depth              (m)
      real bfsfc      ! Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)
      real stable     ! =1 in stable forcing; =0 unstable
      real caseA      ! =1 in case A, =0 in case B 
      integer kbl     ! index of first grid level below hbl 
c 
c  local
      real Rib(2)      ! Bulk Richardson number
      real sigma      ! normalized depth (d/hbl)
      real wm,ws ! turbulent velocity scales         (m/s)
      real dmo(2)      ! Monin-Obukhov Depth
      real hek        ! Ekman depth
      logical LEK, LMO     ! flags for MO and Ekman depth checks

      LOGICAL L_INITFLAG
      common /initflag/L_INITFLAG

      save epsln,Ricr,epsilon,cekman,cmonob,cs,cv,vonk,hbf
c
      data epsln           /  1.e-16 /
      data DelVmin         /  .005   /
      data Ricr            /  0.30   /
      data epsilon         /  0.1    /
      data cekman          /  0.7    /
      data cmonob          /  1.0    /
      data cs              / 98.96   /
      data cv              /  1.6    /
      data vonk            /  0.4    /
      data hbf             /  1.0    /

c       Set MO and Ekman depth flags
      LEK = .true. 
      LMO = .true. 
!      jwtype=1
      ocdepth=-1000.
 
c find bulk Richardson number at every grid level find implied hri
c find Monin-Obukvov depth at every grid level find L
c Compare hri, L and hek to give hbl and kbl
c
c note: the reference depth is -epsilon/2.*zgrid(k), but the reference
c       u,v,t,s values are simply the surface layer values,
c       and not the averaged values from 0 to 2*ref.depth,
c       which is necessary for very fine grids(top layer < 2m thickness)
c note: max values when Ricr never satisfied are
c                                      kbl(i)=km and hbl(i) -zgrid(km)
c       min values are                 kbl(i)=2      hbl(i) -zgrid(1)
 
      Vtc =  cv * sqrt(0.2/cs/epsilon) / vonk**2 / Ricr
  
c     indices for array Rib(i,k), the bulk Richardson number.
      ka = 1
      ku = 2

c     initialize hbl and kbl to bottomed out values
         Rib(ka) = 0.0
         dmo(ka) = -zgrid(kmp1)
         kbl    = km
         hbl    = -zgrid(km)
         hek =  cekman * ustar / (abs(Coriol) + epsln)

      do 200 kl = 2,km

c compute bfsfc = sw fraction at hbf * zgrid
c To optimize the code choose the "swfrac_opt" ?
         call swfrac_opt(hbf,zgrid(kl),kl,jwtype,bfsfc)
           IF(kbl.ge.km) THEN
 
c           use caseA as temporary array for next call to wscale
            caseA  = -zgrid(kl)
 
c           compute bfsfc= Bo + radiative contribution down to hbf * hbl
c Bosol or B0sol
            bfsfc  = Bo 
     $                  + Bosol * (1. - bfsfc)
            stable = 0.5 + SIGN( 0.5, bfsfc+epsln )
            sigma  = stable * 1. + (1.-stable) * epsilon

           ENDIF

c        compute velocity scales at sigma, for hbl= caseA = -zgrid(kl)
         call wscale( sigma, caseA, ustar, bfsfc,   wm, ws)

           IF(kbl.ge.km) THEN
 
c         compute the turbulent shear contribution to Rib
            bvsq =0.5*
     $            ( dbloc(kl-1) / (zgrid(kl-1)-zgrid(kl))+
     $              dbloc(kl) / (zgrid(kl)-zgrid(kl+1)) )
            Vtsq = - zgrid(kl) * ws * sqrt(abs(bvsq)) * Vtc
c         compute bulk Richardson number at new level, dunder
            Rib(ku) = Ritop(kl) / (dVsq(kl)+Vtsq+epsln)
            Rib(ku) = MAX( Rib(ku), Rib(ka) + epsln)
c         linear interpolate to find hbl where Rib = Ricr
            hri   = -zgrid(kl-1) + (zgrid(kl-1)-zgrid(kl)) *
     $                  (Ricr - Rib(ka)) / (Rib(ku)-Rib(ka))

c         compute the Monin Obukov length scale 
c            fmonob    = stable(i) * LMO
            fmonob    = stable * 1.0
            dmo(ku) = cmonob * ustar * ustar * ustar 
     >                   / vonk / (abs(bfsfc) + epsln)
            dmo(ku) = fmonob * dmo(ku) - (1.-fmonob) *zgrid(kmp1) 
              if(dmo(ku).le.(-zgrid(kl))) then
               hmonob =(dmo(ku)-dmo(ka))/(zgrid(kl-1)-zgrid(kl))
               hmonob =(dmo(ku)+hmonob*zgrid(kl)) / (1.-hmonob)
              else
               hmonob = -zgrid(kmp1)
              endif
 
c     compute the Ekman depth
c     fekman  =  stable(i) * LEK
              fekman  =  stable * 1.0
              hekman  = fekman * hek - (1.-fekman) * zgrid(kmp1)
              
c     compute boundary layer depth
              hmin  = MIN(hri, hmonob,  hekman, -ocdepth)
!	write(*,*)hmin,hri, hmonob, hekman
              if(hmin .lt. -zgrid(kl) ) then
                 hbl = hmin
                 kbl = kl
              endif
!        write(*,*)hmin,hri,hmonob,hekman,ocdepth      
           ENDIF                !kbl
        ksave = ka
        ka    = ku
        ku    = ksave
        
 200  continue
      
      call swfrac(-1.0,hbl,jwtype,bfsfc)
      
      bfsfc  = Bo + Bosol * (1. - bfsfc)
      stable = 0.5 + SIGN( 0.5, bfsfc )
      bfsfc  = bfsfc + stable * epsln !ensures bfsfc never=0
c determine caseA and caseB
         caseA  = 0.5 + 
     $       SIGN( 0.5,-zgrid(kbl) -0.5*hwide(kbl) -hbl)

      return
      end

c *********************************************************************

      subroutine wscale( sigma, hbl, ustar, bfsfc,
     $                  wm , ws   )
c
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c
c
c lookup table
      parameter ( ni = 890,              ! number of values for zehat
     $            nj = 48)              ! number of values for ustar

      real wmt(0:ni+1,0:nj+1)           ! lookup table for wm
      real wst(0:ni+1,0:nj+1)           ! lookup table for ws
      real deltaz                       ! delta zehat in table
      real deltau                       ! delta ustar in table
      real zmin,zmax                    ! zehat limits for table
      real umin,umax                    ! ustar limits for table
      logical firstf
      save wmt,wst,deltaz,deltau,zmin,zmax,umin,umax,firstf
c
      data zmin,zmax  / -4.e-7, 0.0   / ! m3/s3
      data umin,umax  /  0.   , .04   / ! m/s
      data firstf     / .true.        /

c  model
      integer imt          ! number of horizontal grid points

c  input
      real sigma      ! normalized depth (d/hbl)
      real hbl        ! boundary layer depth (m)
      real ustar      ! surface friction velocity         (m/s)
      real bfsfc      ! total surface buoyancy flux       (m^2/s^3)

c  output
      real wm,ws ! turbulent velocity scales at sigma

c local
      real zehat           ! = zeta *  ustar**3
      real zeta            ! = stability parameter d/L

      save epsln,c1,am,cm,c2,zetam,as,cs,c3,zetas,vonk

      data epsln           /   1.0e-20/
      data c1              /   5.0   /
      data am,cm,c2,zetam  /   1.257 ,  8.380 , 16.0 , - 0.2  /
      data as,cs,c3,zetas  / -28.86  , 98.96  , 16.0 , - 1.0  /
      data vonk            /   0.40  /
c
c construct the wm and ws lookup tables
c
      if(firstf) then

         deltaz = (zmax-zmin)/(ni+1) 
         deltau = (umax-umin)/(nj+1)

         do 100 i=0,ni+1
            zehat = deltaz*(i) + zmin
            do 90 j=0,nj+1
               usta = deltau*(j) + umin
               zeta = zehat/(usta**3+epsln)

               if(zehat.ge.0.) then
                  wmt(i,j) = vonk*usta/(1.+c1*zeta)
                  wst(i,j) = wmt(i,j)
               else
                  if(zeta.gt.zetam) then
                     wmt(i,j) = vonk* usta * (1.-c2*zeta)**(1./4.)
                  else
                     wmt(i,j) = vonk* (am*usta**3 - cm*zehat)**(1./3.)
                  endif
                  if(zeta.gt.zetas) then
                     wst(i,j) = vonk* usta * (1.-c3*zeta)**(1./2.)
                  else
                     wst(i,j) = vonk* (as*usta**3 - cs*zehat)**(1./3.)
                  endif
               endif   
 90         continue   
 100     continue
         firstf=.false.
      endif       
 
c use lookup table for zehat < zmax  ONLY;  otherwise use stable formulae
         zehat = vonk * sigma * hbl * bfsfc

         IF (zehat .le. zmax) THEN

         zdiff  = zehat-zmin
         iz = int( zdiff/deltaz )
         iz = min( iz , ni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff/deltau)
         ju = min( ju , nj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff/deltaz - float(iz)
         ufrac = udiff/deltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was

         ELSE

         ucube = ustar**3
         wm = vonk * ustar * ucube / (ucube + c1 * zehat)
         ws = wm

         ENDIF   
      return
      end

c **********************************************************************
 
      subroutine ri_iwmix ( km  , kmp1, imt   ,
     $                      Shsq, dbloc , zgrid ,
     $                      visc, difs, dift  )
c
c     compute interior viscosity diffusivity coefficients due to
c     shear instability (dependent on a local richardson number)
c     and due to background internal wave activity.
c

c  input

      real Shsq(kmp1)    ! (local velocity shear)^2          (m/s)^2
      real dbloc(km)     ! local delta buoyancy              (m/s^2)
      real zgrid(kmp1)   ! vertical grid (<= 0)              (m)
      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
 
c output
      real visc(0:kmp1)  ! vertical viscosivity coefficient  (m^2/s)
      real difs(0:kmp1)  ! vertical scalar diffusivity       (m^2/s)
      real dift(0:kmp1)  ! vertical temperature diffusivity  (m^2/s)
 
c local variables
      real Rig,Rigg          ! local richardson number
      real fri               ! function of Rig
 
      save epsln,Riinfty,Ricon,difm0,difs0,difmcon,difscon,
     &            difmiw,difsiw,c1

      data  epsln   / 1.e-16 /   ! a small number          
      data  Riinfty /  0.8     /  ! LMD default was = 0.7
      data  Ricon   / -0.2    /  ! note: exp was repl by multiplication
      data  difm0   / 0.005  /  ! max visc due to shear instability
      data  difs0   / 0.005  /  ! max diff ..  .. ..    ..
!     data  difmiw  / 0.0001  /  ! background/internal waves visc(m^2/s)
      data  difmiw  / 0.0005  /  ! background/internal waves visc(m^2/s) DD change
!      data  difsiw  / 0.00001 /  ! ..         ..       ..    diff(m^2/s)
      data  difsiw  / 0.00005 /  ! ..         ..       ..    diff(m^2/s) DD change
      data  difmcon / 0.0000   /  ! max visc for convection  (m^2/s)
      data  difscon / 0.0000   /  ! max diff for convection  (m^2/s)
      data  c1/ 1.0/
      data  c0/ 0.0/
      data  mRi/ 1 /               ! number of vertical smoothing passes
c
c     compute interior gradient Ri at all interfaces, except surface

c-----------------------------------------------------------------------
c     compute interior gradient Ri at all interfaces ki=1,km, (not surface)
c       use visc(ki=1,km) as temporary storage to be smoothed
c       use dift(ki=1,km) as temporary storage of unsmoothed Ri
c       use difs(ki=1,km) as dummy in smoothing call
 
      do ki = 1, km
            Rig  = dbloc(ki) * (zgrid(ki)-zgrid(ki+1)) /
     $                          ( Shsq(ki) + epsln)
            dift(ki) = Rig          
	    visc(ki) = dift(ki)                        
      end do
 
c-----------------------------------------------------------------------
c     vertically smooth Ri mRi times
      do mr = 1,mRi
        call z121(kmp1,c0,Riinfty,visc,difs)
      enddo

c-----------------------------------------------------------------------
c                           after smoothing loop
      do ki = 1, km
c  evaluate f of unsmooth Ri (fri) for convection        store in fcon
c  evaluate f of   smooth Ri (fri) for shear instability store in fri
 
            Rigg  = AMAX1( dift(ki) , Ricon )
            ratio = AMIN1( (Ricon-Rigg)/Ricon , c1 )
            fcon  = (c1 - ratio*ratio)
            fcon  = fcon * fcon * fcon
 
            Rigg  = AMAX1( visc(ki) , c0 )
            ratio = AMIN1( Rigg/Riinfty , c1 )
            fri   = (c1 - ratio*ratio)
            fri   = fri * fri * fri

c ----------------------------------------------------------------------
c            evaluate diffusivities and viscosity
c    mixing due to internal waves, and shear and static instability
 
          visc(ki) = (difmiw + fcon * difmcon + fri * difm0)
          difs(ki) = (difsiw + fcon * difscon + fri * difs0)
          dift(ki) = difs(ki)

      end do
c ------------------------------------------------------------------------
c         set surface values to 0.0
 
                 visc(0)    = c0
                 dift(0)    = c0
                 difs(0)    = c0
      return
      end

c *********************************************************************
      Subroutine z121 (kmp1,vlo,vhi,V,w)

c    Apply 121 smoothing in k to 2-d array V(i,k=1,km)
c     top (0) value is used as a dummy
c     bottom (kmp1) value is set to input value from above.
 
c  input

      real V(0:kmp1)  ! 2-D array to be smoothed in kmp1 direction
      real w(0:kmp1)  ! 2-D array of internal weights to be computed
 
      km  = kmp1 - 1
 
      w(0)    =   0.0         
      w(kmp1) =   0.0             
      V(0)    =   0.0      
      V(kmp1) =   0.0

      do k=1,km  
         if((V(k).lt.vlo).or.(V(k).gt.vhi)) then
              w(k) = 0.0
         else 
              w(k) = 1.0
         endif
      enddo
 
      do k=1,km
          tmp = V(k)
          V(k) = w(k-1)*V(0)+2.*V(k)+w(k+1)*V(k+1)
          wait   = w(k-1) + 2.0 + w(k+1)
          V(k) = V(k) / wait             
          V(0) = tmp
      enddo

      return
      end
c *********************************************************************
c *********************************************************************

      subroutine blmix 
     $             (km   , kmp1 , imt , mdiff , zgrid, hwide , 
     $              ustar, bfsfc, hbl , stable, caseA, 
     $              visc , difs , dift, kbl   , 
     $              gat1 , dat1 , dkm1, blmc  , ghats, 
     $              sigma, wm   , ws  )
c     mixing coefficients within boundary layer depend on surface
c     forcing and the magnitude and gradient of interior mixing below
c     the boundary layer ("matching").

CAUTION if mixing bottoms out at hbl = -zgrid(km) THEN
c   fictitous layer kmp1 is needed with small but finite width (eg. 1.e-10)
c model

      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
      integer mdiff          ! number of viscosities + diffusivities
      real zgrid(kmp1)   ! vertical grid (<=0)               (m)
      real hwide(kmp1)   ! layer thicknesses                 (m)
c
c input
      real ustar        ! surface friction velocity         (m/s)
      real bfsfc        ! surface buoyancy forcing        (m^2/s^3)
      real hbl          ! boundary layer depth              (m)
      real stable       ! = 1 in stable forcing
      real caseA        ! = 1 in case A
      real visc(0:kmp1)  ! vertical viscosity coefficient    (m^2/s)
      real difs(0:kmp1)  ! vertical scalar diffusivity       (m^2/s)
      real dift(0:kmp1)  ! vertical temperature diffusivity  (m^2/s)
      integer kbl       ! index of first grid level below hbl
c
c output
      real gat1(mdiff)
      real dat1(mdiff)
      real dkm1(mdiff)   ! boundary layer difs at kbl-1 level
      real blmc(km,mdiff)! boundary layer mixing coefficients(m^2/s)
      real ghats(km)     ! nonlocal scalar transport
c
c  local
      real sigma        ! normalized depth (d / hbl)
      real ws, wm  ! turbulent velocity scales         (m/s)
 
      save epsln,epsilon,c1,am,cm,c2,zetam,as,cs,c3,zetas,
     $     cstar,grav,vonk

      data epsln             /   1.e-20 /
      data epsilon           /   0.1    /
      data c1                /   5.0    /
      data am,cm,c2,zetam    /   1.257  ,  8.380, 16.0, - 0.2 / !7-24-92
      data as,cs,c3,zetas    / -28.86   , 98.96 , 16.0, - 1.0 /
      data cstar             /    5.    /
!      data cstar             /    1.    /  ! BP testing
      data grav              /   9.816  /
      data vonk              /   0.4    /
c
      cg = cstar * vonk * (cs * vonk * epsilon)**(1./3.)
       
c compute velocity scales at hbl
      sigma = stable * 1.0 + (1.-stable) * epsilon
 
      call wscale( sigma, hbl, ustar, bfsfc,   wm, ws)
      kn    = ifix(caseA+epsln) *(kbl -1) +
     $            (1-ifix(caseA+epsln)) * kbl
 
c find the interior viscosities and derivatives at hbl(i) 
      delhat = 0.5*hwide(kn) - zgrid(kn) - hbl
      R      = 1.0 - delhat / hwide(kn)
      dvdzup = (visc(kn-1) - visc(kn)) / hwide(kn) 
      dvdzdn = (visc(kn)   - visc(kn+1)) / hwide(kn+1)
      viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $                        R  * (dvdzdn + abs(dvdzdn)) )
 
      dvdzup = (difs(kn-1) - difs(kn)) / hwide(kn) 
      dvdzdn = (difs(kn)   - difs(kn+1)) / hwide(kn+1)
      difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $                        R  * (dvdzdn + abs(dvdzdn)) )
 
      dvdzup = (dift(kn-1) - dift(kn)) / hwide(kn) 
      dvdzdn = (dift(kn)   - dift(kn+1)) / hwide(kn+1)
      diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $                        R  * (dvdzdn + abs(dvdzdn)) )
c
      visch  = visc(kn) + viscp * delhat
      difsh  = difs(kn) + difsp * delhat
      difth  = dift(kn) + diftp * delhat
c
      f1 = stable * c1 * bfsfc / (ustar**4+epsln) 
      gat1(1) = visch / hbl / (wm+epsln)
      dat1(1) = -viscp / (wm+epsln) + f1 * visch
      dat1(1) = min(dat1(1),0.) 

      gat1(2) = difsh  / hbl / (ws+epsln)
      dat1(2) = -difsp / (ws+epsln) + f1 * difsh 
      dat1(2) = min(dat1(2),0.) 
  
      gat1(3) = difth /  hbl / (ws+epsln)
      dat1(3) = -diftp / (ws+epsln) + f1 * difth 
      dat1(3) = min(dat1(3),0.) 

      do ki = 1,km       
c     compute turbulent velocity scales on the interfaces
            sig     = (-zgrid(ki) + 0.5 * hwide(ki)) / hbl
            sigma= stable*sig + (1.-stable)*AMIN1(sig,epsilon)
             call wscale( sigma, hbl, ustar, bfsfc,   wm,  ws)
c
c     compute the dimensionless shape functions at the interfaces
c
            sig = (-zgrid(ki) + 0.5 * hwide(ki)) / hbl
            a1 = sig - 2.
            a2 = 3.-2.*sig
            a3 = sig - 1.
c
            Gm = a1 + a2 * gat1(1) + a3 * dat1(1) 
            Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
            Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
c
c     compute boundary layer diffusivities at the interfaces
c
            blmc(ki,1) = hbl * wm * sig * (1. + sig * Gm)
            blmc(ki,2) = hbl * ws * sig * (1. + sig * Gs)
            blmc(ki,3) = hbl * ws * sig * (1. + sig * Gt)
c
c     nonlocal transport term = ghats * <ws>o
            ghats(ki) = (1.-stable) * cg / (ws * hbl +epsln)

      end do
c find diffusivities at kbl-1 grid level 
         sig      =  -zgrid(kbl-1)  / hbl
         sigma =stable * sig + (1.-stable) * AMIN1(sig,epsilon)
c
         call wscale( sigma, hbl, ustar, bfsfc,   wm, ws)
c
         sig = -zgrid(kbl-1) / hbl
         a1= sig - 2.
         a2 = 3.-2.*sig
         a3 = sig - 1.
         Gm = a1 + a2 * gat1(1) + a3 * dat1(1)
         Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
         Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
         dkm1(1) = hbl * wm * sig * (1. + sig * Gm)
         dkm1(2) = hbl * ws * sig * (1. + sig * Gs)
         dkm1(3) = hbl * ws * sig * (1. + sig * Gt)

      return
      end

c ******************************************************************

      subroutine enhance (km   , kmp1  , imt   , mdiff , dkm1  , visc ,
     &                    difs , dift  , hbl   , kbl   , zgrid , caseA,
     &                    blmc , ghats )
c
c enhance the diffusivity at the kbl-.5 interface
c
c input

      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
      integer mdiff          ! number of viscosities + diffusivities
      integer kbl       ! grid above hbl
      real hbl          ! boundary layer depth             (m)
      real dkm1(mdiff)   ! bl diffusivity at kbl-1 grid level
      real zgrid(kmp1)   ! vertical grid (<= 0)             (m)
      real visc(0:kmp1)  ! enhanced viscosity               (m^2/s) 
      real difs(0:kmp1)  ! enhanced thermal diffusivity     (m^2/s)
      real dift(0:kmp1)  ! enhanced scalar  diffusivity     (m^2/s)
      real caseA        ! = 1 in caseA, = 0 in case B
 
c input/output
      real ghats(km)     ! nonlocal transport               (s/m**2)
c                              modified ghats at kbl(i)-1 interface
c output
      real blmc(km,mdiff)! enhanced bound. layer mixing coeff.
c
c local
      real delta             ! fraction hbl lies beteen zgrid neighbors
c
      do ki=1,km-1

          if(ki .eq. (kbl - 1) ) then

            delta = (hbl+zgrid(ki)) / (zgrid(ki)-zgrid(ki+1))

            dkmp5 = caseA * visc(ki) + (1.-caseA) * blmc(ki,1)
            dstar = (1.-delta)**2 * dkm1(1) + delta**2 * dkmp5      
            blmc(ki,1) = (1.-delta) * visc(ki) + delta * dstar

            dkmp5 = caseA * difs(ki) + (1.-caseA) * blmc(ki,2)
            dstar = (1.-delta)**2 * dkm1(2) + delta**2 * dkmp5    
            blmc(ki,2) = (1.-delta) * difs(ki) + delta * dstar

            dkmp5 = caseA * dift(ki) + (1.-caseA) * blmc(ki,3)
            dstar = (1.-delta)**2 * dkm1(3) + delta**2 * dkmp5     
            blmc(ki,3) = (1.-delta) * dift(ki) + delta * dstar
            
            ghats(ki) = (1.-caseA) * ghats(ki)

          endif
      end do

      return
      end

c***********************************************************************
