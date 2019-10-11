      SUBROUTINE fluxes(taux,tauy,solar,nsolar,pme)

      IMPLICIT NONE
      
      include 'parameter.inc'
      include 'flx_sfc.com'
      include 'ocn_paras.com'
      include 'flx_paras.com'
      include 'timocn.com'
      include 'times.com'

      REAL taux,tauy,solar,nsolar, pme

!      call forcing(taux,tauy,solar,nsolar,pme)
            IF ((taux .EQ. 0.0) .AND. (tauy .EQ. 0.0)) THEN
               taux=1.e-10
            ENDIF
            sflux(1,5,0)=taux
            sflux(2,5,0)=tauy
            sflux(3,5,0)=solar             
            sflux(4,5,0)=nsolar
            sflux(5,5,0)=0.0 ! Melting of sea-ice = 0.0               
            sflux(6,5,0)=pme/100000 ! assuming rain = P-E          
            call ntflx
      
      RETURN
      END

************************************************************************
      SUBROUTINE ntflx

      IMPLICIT NONE

      include 'parameter.inc'
      include 'times.com'
      include 'flx_sfc.com'
      include 'flx_profs.com'
      include 'vert_pgrid.com'
      include 'ocn_paras.com'

      INTEGER k
      REAL SWDK
      REAL SWDK_OPT(0:NZ)
      EXTERNAL SWDK

      COMMON /SWDK_SAVE/SWDK_OPT

      IF (ntime .le. 1) THEN
      DO k=0,NZ
         swdk_opt(k)=swdk(-dm(k))
      ENDDO
      ENDIF
      DO k=0,NZ
         IF (ntime .ge. 1) THEN 
            wXNT(k,1)=-sflux(3,5,0)*swdk_opt(k)/(rho(0)*CP(0))
         ENDIF
      ENDDO

      RETURN
      END

*******************************************************************

      REAL FUNCTION SWDK(z)
      include 'parameter.inc'
      include 'proc_pars.com'

      parameter(max=5)
      real Rfac(max),a1(max),a2(max)
c         types =  I       IA      IB      II      III
c             j =  1       2       3       4       5
      data Rfac /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1   /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2   / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /

      j = jerlov

      SWDK =        Rfac(j)  * dexp(dble(z/a1(j))) 
     >       + (1.0-Rfac(j)) * dexp(dble(z/a2(j)))

      return
      end
**************************************************************
