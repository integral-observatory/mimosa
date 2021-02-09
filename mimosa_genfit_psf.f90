
!************************************************
! PSF FIT MODULES : DECLARATIONS AND INTERFACES
!************************************************


!****************************
MODULE IBIS_PSF_AUX
!****************************
INTERFACE

SUBROUTINE IBIS_PSF_MINIMIZE(Status)
!------------------------------------
INTEGER  :: Status
END SUBROUTINE IBIS_PSF_MINIMIZE

SUBROUTINE IBIS_PSF_HESSIAN(status)
!------------------------------------
IMPLICIT NONE
INTEGER :: Status
END SUBROUTINE IBIS_PSF_HESSIAN

END INTERFACE

!****************************
END MODULE IBIS_PSF_AUX
!****************************







!*********************
MODULE FIT_PSF_MODULE
!********************

INTERFACE
SUBROUTINE PSF_FIT(FixedPerm,inCatNumber,SKYCLE,VAR,FILTER, xc,yc,xs,ys,&
                   fitflux,locerr,chi2red,Status)
!-------------------------------------------------------------------
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
REAL(kind=4) ,dimension(:,:),pointer :: SKYCLE,VAR
LOGICAL      ,dimension(:,:),pointer :: Filter
REAL(kind=8) :: xc,yc  ! fit centre
REAL(kind=8) :: xs,ys  ! source fit position
REAL(kind=8) :: fitflux,locerr,chi2red
INTEGER :: FixedPerm,inCatNumber,Status
END SUBROUTINE PSF_FIT
END INTERFACE
!*********************
END MODULE FIT_PSF_MODULE
!********************



!*********************
MODULE PSF_FIT_SUB_MODULE
!********************

INTERFACE

!PSF FIT SUBROUTINES
SUBROUTINE IBIS_PSF_INITMINI(FixedPerm,IncatNumber,SKY,VAR,sources_to_fit,xc,yc,xs,ys,print)
!-------------------------------------------------------------------
  IMPLICIT NONE
  ! Parameter
  REAL      (kind=4), dimension(:,:), pointer :: SKY,VAR
  INTEGER :: FixedPerm,IncatNumber,sources_to_fit
  REAL      (kind=8),dimension(sources_to_fit) :: xs,ys
  REAL      (kind=8) :: xc,yc
  INTEGER :: print 

END SUBROUTINE IBIS_PSF_INITMINI

SUBROUTINE IBIS_PSF_FIT_MINIMAP(sources_to_fit,xs,ys,fitflux,locerr,Status)
!-------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: sources_to_fit
  REAL      (kind=8),dimension(sources_to_fit) :: xs,ys
  REAL      (kind=8),dimension(sources_to_fit) :: fitflux,locerr
  INTEGER :: Status 

END SUBROUTINE IBIS_PSF_FIT_MINIMAP


END INTERFACE
!*************************
END MODULE PSF_FIT_SUB_MODULE
!*************************








! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! PSF FIT PROCEDURES CODE
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!======================================================================
SUBROUTINE PSF_FIT(FixedPerm,inCatNumber,SKYCLE,VAR,FILTER,xc,yc,xs,ys,&
                   fitflux,locerr,chi2red,Status)
!======================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FIT_PAR_MODULE
USE PSF_FIT_SUB_MODULE
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
REAL(kind=4) ,dimension(:,:),pointer :: SKYCLE,VAR
LOGICAL      ,dimension(:,:),pointer :: Filter

REAL(kind=8) :: xc,yc  ! fit centre
REAL(kind=8) :: xs,ys  ! source fit position
REAL(kind=8) :: fitflux,locerr,chi2red
INTEGER :: FixedPerm,inCatNumber,Status,result
!LOCAL VARIABLES
LOGICAL :: jestpeak
REAL(kind=8) ,dimension(:),pointer :: fitfluxtab,locerrtab,xfit,yfit
REAL(kind=8) :: snrmax
INTEGER :: iok,icsky,jcsky,idim,jdim,imax,jmax
CHARACTER(len=20) :: procName = 'FIT_PSF'

Status = 0
result = 0

ALLOCATE(xfit(1),yfit(1),locerrtab(1),fitfluxtab(1),stat = iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif 

fitfluxtab(:) = 0.
locerrtab(:)  =  0.
xfit(1) = xc
yfit(1) = yc
! FITTING to get 1 fine source location :xs,ys
CALL IBIS_PSF_INITMINI(FixedPerm,inCatNumber,SKYCLE,VAR,1,xc,yc,xfit,yfit,DebugMode)  
CALL Ibis_Psf_Fit_MiniMap(1,xfit,yfit,fitfluxtab,locerrtab,Status)
if(Status.ne.0)then
   result = 1
   Status = ISDC_OK
   xs = xc;ys = yc
   fitflux = skycle(nint(xs),nint(ys))-rmean
   locerr = 0.
else
   xs = xfit(1)
   ys = yfit(1)
   fitflux = fitfluxtab(1)
   locerr = locerrtab(1)
endif

if(.not.FILTER(nint(xs),nint(ys)))then ! forbidden zone
   call WAR_MESSAGE(procName,&
        'Fit in forbidden zone - reset',0,Status)
   result = 1
   xs = xc
   ys = yc
endif
chi2red = spchisq

if(result==1)then !bad first fit result
   call resid(xc,yc,jestpeak,imax,jmax,snrmax)
   if(jestpeak)then ! unusual peak in resmap

      idim = size(SkyCle,1)
      jdim = size(SkyCle,2)
      icsky = (idim+1)/2
      jcsky = (jdim+1)/2
  
      !two source fit
      DEALLOCATE(xfit,yfit,locerrtab,fitfluxtab)
      ALLOCATE(xfit(2),yfit(2),locerrtab(2),fitfluxtab(2),stat = iok)
      IF (iok /= 0) then
         call MESSAGE(procName,'Allocation problem',AllocError,Status)
         return
      endif
      fitfluxtab(:) = 0.
      locerrtab(:)  =  0.
      xfit(1) = xc
      yfit(1) = yc
     ! second_peak
      xfit(2) = imax
      yfit(2) = jmax
      CALL IBIS_PSF_INITMINI(FixedPerm,inCatNumber,SKYCLE,VAR,2,xc,yc,xfit,yfit,DebugMode)  
      CALL Ibis_Psf_Fit_MiniMap(2,xfit,yfit,fitfluxtab,locerrtab,Status)
      if(Status.eq.0)then
         chi2red = spchisq
         xs = xfit(1); ys = yfit(1)
         locerr = locerrtab(1); fitflux = fitfluxtab(1)
      else
         call WAR_MESSAGE(procName,&
              ' Failure of  PSF fit procedure ',1,Status)
      endif 
      
   endif ! unusual peak in resmap
endif !res==0

DEALLOCATE(xfit,yfit,locerrtab,fitfluxtab)


CONTAINS

!.........................................
subroutine resid(xc,yc,jest,imax,jmax,snrmax)
!.........................................
USE FIT_PAR_MODULE
IMPLICIT NONE
Logical :: jest
INTEGER   :: imax,jmax
REAL(kind=8)    ::xc,yc
REAL(kind=8)    ::mean,absmean,dev,snr,snrmax
INTEGER   ::  i,j,idim,jdim

idim = imapmax-imapmin+1
jdim = jmapmax-jmapmin+1

mean = sum(resmap)/idim/jdim
dev  = 0.
do i=imapmin,imapmax
do j=jmapmin,jmapmax
   dev = dev+ (resmap(i,j)-mean)**2
enddo
enddo
dev = sqrt(dev/idim/jdim)
absmean = sum(abs(resmap(imapmin:imapmax,jmapmin:jmapmax)))/idim/jdim
jest = .false.
snrmax = 0.
do i=imapmin,imapmax
do j=jmapmin,jmapmax
   snr = (resmap(i,j)-mean)/dev
   if(snr  > 6.)then
      if(snr > snrmax)then
         snrmax = snr
         imax = i
         jmax = j
      endif
      jest = .true.
   endif
enddo
enddo
if(jest)then
  i = idim/2
  j = jdim/2
  imax = xc-i+imax-1
  jmax = yc-j+jmax-1
endif
!.........................................
end subroutine resid
!.........................................


!======================
END SUBROUTINE PSF_FIT
!======================


!#####################################
! SUBROUTINES FROM PSF_FIT_SUB_MODULE
!#####################################


!|| --------------------------------------------------------------------------
!|| SUB IBIS_PSF_INITMINI                                                  1.2
!||
!||   SUB IBIS_PSF_INITMINI(SKY,VAR,xc,yc,xs,yc)
!||
!||   Initialise values for minimization and prepare minicarte
!||   Values/minicarte are stored in variables of the module
!||     IBIS_PSF_FIT_PAR
!||
!||   Input:
!||    SKY         R4  Array Poi     Sky image
!||    VAR         R4  Array Poi     Variance image
!||    xc,yc       R4                Center of minimap in sky coord.
!||    xs,ys       R4                Input source position in sky coord.
!||
!||   History:
!||   1.2 - 21/09/00 AG : C, dbStaSpRe, dbw defined before, change name
!||                       from INITMINI2
!||   1.1 - 31/08/00 AG : C, include iDetTyp to allow PICSIT psf fit,
!||                       explicit type of R8 variab/values
!||
!|| --------------------------------------------------------------------------
!======================================================
SUBROUTINE IBIS_PSF_INITMINI(FixedPerm,IncatNumber,SKY,VAR,&
     sources_to_fit,xc,yc,xs,ys,print)
!======================================================
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
  USE FIT_PAR_MODULE
  USE COMMON_PAR
  IMPLICIT NONE

  ! Parameter
  REAL      (kind=4), dimension(:,:), pointer :: SKY,VAR
  INTEGER :: FixedPerm,IncatNumber,sources_to_fit
  REAL      (kind=8),dimension(sources_to_fit) :: xs,ys
  REAL      (kind=8) :: xc,yc
  INTEGER :: print
!     LOCAL VARIABLES
  INTEGER :: isky,jsky,iexfov,isx,jsx,jfitdim
  INTEGER :: ns,ncomp,i,j,l,k,ki,lj,ifgx,ifgy,ieqspre
  REAL      (kind=8) :: xcsky,ycsky,varia
  REAL      (kind=8) :: dbStaSpRe

  method = 1
  f04asfprint = print
  isky=SIZE(SKY,1)
  jsky=SIZE(SKY,2)
  xcsky=REAL(isky+1)/2.
  ycsky=REAL(jsky+1)/2.

  !     SOME PARAMETERS and initializations
  ! Double Precision STArting SPatial REsolution 
  ! (~ 0 for perfect binning, ~ 0.4 for ISGRI, ~ 0.2 for Picsit)
  dbStaSpRe = 0.01D00
  IF (iDetTyp > 0) dbStaSpRe = 0.430D00 / DBLE(iDetTyp)
!!      dbw=2.435D00/DBLE(iDetTyp)
  IF (dbw == 0.D00) dbw=2.435D00
  
  fitmap=0.
  error=0.

  ! FINE DEC PSF
  iexfov=1
  idec=-1

!     INPUTS
  ifitdim = mapdim
  jfitdim = mapdim
  if(ifitdim.gt.mapdim)then
!!$     print*,' Fit region too big :',&
!!$          ifitdim,' Size set to max :',mapdim
     ifitdim=mapdim
  endif
 if(jfitdim.gt.mapdim)then
!!$     print*,' Fit region too big :',&
!!$          jfitdim,' Size set to max :',mapdim
     jfitdim=mapdim
  endif 

!  Values of variance and error at minimap center (not at source position !)
  isx=NINT(xc)
  jsx=NINT(yc)
  if(iexfov.eq.1)then
     varia=VAR(isx,jsx)
     error=SNGL(DSQRT(DBLE(varia)))
  endif

!  WRITE IMAGE SECTOR IN THE ARRAY FITMAP
  ifd2=ifitdim/2
  jfd2=jfitdim/2
  ifgx=isx-ifd2
  ifgy=jsx-jfd2
  k=ifgx-1
  l=ifgy-1
  do i=1,ifitdim
     ki=k+i
     do j=1,jfitdim
        lj=l+j
        fitmap(i,j)=SKY(ki,lj)
        varmap(i,j)=VAR(ki,lj)
     enddo
  enddo

  ! COORDINATES X0 Y0
  x00=REAL(ifgx-1)
  y00=REAL(ifgy-1)

  ! PIXEL COORDINATES OF THE TELESCOPE AXIS :
  if(idec.eq.-2)then
     xta=xcsky-0.5
     yta=ycsky-0.5
  elseif(idec.eq.-1)then
     xta=xcsky
     yta=ycsky
  endif
!!           print*,' Source',isx,jsx,SKY(isx,jsx),error,x00,y00,ifgx,&
!!                  (k+ifitdim),ifgy,(l+ifitdim)

!  DEFINITION OF TOTAL SPSF
  nsouf = sources_to_fit


  PSFP(:,:) = 0.
!  Enter Bkg Level (cts/pix)
  PSFP(1,1) = 0.D0
  PSFP(1,2) = 1.D0

  ieqspre=0
  do ns=1,nsouf
     ncomp=(ns-1)*3+1
! Enter photon index (> 0)
     spin(ns) = 0.
! Enter Starting Y T.A. Pos.
     PSFP(ncomp+2,1) = xs(ns)
     PSFP(ncomp+2,1) = PSFP(ncomp+2,1)-x00
     PSFP(ncomp+2,2) = 1.D0
! Enter Starting Z T.A. Pos.
     PSFP(ncomp+3,1) = ys(ns)
     PSFP(ncomp+3,1) = PSFP(ncomp+3,1)-y00
     PSFP(ncomp+3,2) = 1.D0
        if((ns==1).and.(FixedPerm==1))then
        ! only first sou pos can be fixed
        ! and permission to fix
        ! not for model
        if(FitPosFixed)then
           PSFP(ncomp+2:ncomp+3,2) = 0.0
           ! all source positions fixed
        else
           if(inCatNumber > 0)then
              !catalog source
              if(IncatFixed(InCatNumber))then
                 PSFP(ncomp+2:ncomp+3,2) = 0.0
          
              endif
           endif
        endif

         
     endif  ! only first sou pos can be fixed
     
! Enter Starting Source Intensity
     PSFP(ncomp+1,1) = 0.D0
     PSFP(ncomp+1,2) = 1.D0
     if(spin(ns).eq.0)ieqspre=1
  enddo

  if(ieqspre.eq.1)then
! Starting Resolution (sigma in pix)  ==> ',$)
     PSFP((nsouf*3+2),1) = dbStaSpRe
     PSFP((nsouf*3+2),2) = 1
  endif

! Confidence level for error estimation
  confl= 0.683
  np   = 1
!=================================
 END SUBROUTINE IBIS_PSF_INITMINI
!=================================

!|| --------------------------------------------------------------------------
!|| SUB IBIS_PSF_FIT_MINIMAP                                               1.1
!||  
!||   SUB IBIS_PSF_FIT_MINIMAP
!||
!||   Description:
!||     Fit of a minimap 
!|| 
!||   Input:
!||     Input and output in variable of module IBIS_PSF_FIT_PAR
!||
!||   Dependencies:
!||     Variable must be deifned for example by a call of INITMINI
!||  
!||   Sub-History :
!||   1.1 - 21/09/00 AG : C, iPsfTest
!||   1.0 ibis_decon1  02/09/99-AG : C, correct mask/reg error, iPsfTest
!||
!|| --------------------------------------------------------------------------
!=========================================================================
SUBROUTINE IBIS_PSF_FIT_MINIMAP(sources_to_fit,xs,ys,fitflux,locerr,Status)
!========================================================================
 USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
 USE FIT_PAR_MODULE
  USE IBIS_PSF_AUX
  IMPLICIT NONE
  INTEGER :: sources_to_fit
  REAL      (kind=8),dimension(sources_to_fit) :: xs,ys
  REAL      (kind=8),dimension(sources_to_fit) :: fitflux,locerr
  INTEGER :: Status
  ! Local variables
  REAL      (kind=8) :: fsumsq, bestchi
  INTEGER :: i,j,ncomp,ns,ieqspre

  Status = 0
! CHECK STARTING POINT
  IF (iPsfTest == 10) THEN       
     ieqspre=0
     print*,''
     print*,' Starting parameters and freedom :'
     ncomp=1
     write(6,1020) ncomp, PSFP(1,1), INT(PSFP(1,2))
     do ns=1,nsouf
        ncomp=(ns-1)*3+1
        do i=1,3
           write(6,1020)(ncomp+i),PSFP(ncomp+i,1),INT(PSFP(ncomp+i,2))
        enddo
        if(spin(ns).eq.0)ieqspre=1
     enddo
     ncomp=nsouf*3+2
     if(ieqspre.eq.1)then
        write(6,1020) ncomp,PSFP(ncomp,1),INT(PSFP(ncomp,2))
     endif
1020 FORMAT(2X,' Parameter N.',I3,' :',E15.6,'  Free (1=y 0=n) :',I3)
  ENDIF

!  FILL THE DOUBLE PRECISION MINIMAP ARRAY
  imapmin = 1
  imapmax = ifitdim
  jmapmin = 1
  jmapmax = ifitdim

!  MINIMIZATION
  CALL IBIS_PSF_MINIMIZE(Status)
  If(Status.ne.0)return

  
  fsumsq=spchisq*dof   
  bestchi=fsumsq

!  Model in psfmap
  do j=1,ifitdim
     do i=1,ifitdim
        psfmap(i,j)=(fitmap(i,j)-resmap(i,j)*error)
     enddo
  enddo

!  ERRORS IN THE ASSUMPTION OF QUADRATIC APPROXIMATION AND USING THE PSF
!  FOR DELTA-DECODING 
  CALL IBIS_PSF_HESSIAN(Status)
  If(Status.ne.0)return
  !  RESULT OF FIT
  IF (iPsfTest == 10) THEN  

     print*,' Best fit parameters :'
     ncomp=1
     write(6,1020) ncomp,PSFP(1,1),INT(PSFP(1,2))
     do ns=1,nsouf
        ncomp=(ns-1)*3+1
        do i=1,3
           write(6,1020)(ncomp+i),PSFP(ncomp+i,1),INT(PSFP(ncomp+i,2))
              
        enddo
     enddo
     ncomp=nsouf*3+2
     if(ieqspre.eq.1)then
        write(6,1020) ncomp,PSFP(ncomp,1),INT(PSFP(ncomp,2))
     endif
     fsumsq=spchisq*dof   
     print*,' Fit par. (iter,ndof,nparm,fsumsq,rchisq):',&
          iter,INT(dof),(nparm+lparm),fsumsq,spchisq

     !CALL IBIS_PSF_OUTPUT

  ENDIF

  if(RawFitFluxVerif)then
     if(PSFP(2,1) <   maxval(fitmap(ifd2-2:ifd2+2,jfd2-2:ifd2+2)))then
        ! filled flux too low
        status = 15
        return
     endif
  endif

 do ns=1,nsouf
    
    xs(ns)=SNGL(PSFP(3,ns))+x00
    ys(ns)=SNGL(PSFP(4,ns))+y00
    fitflux(ns)=SNGL(PSFP(2,ns))
    locerr(ns)=SQRT(SNGL(PSFE(3,ns))**2 + SNGL(PSFE(4,ns))**2)
    
enddo
!==================================
END SUBROUTINE IBIS_PSF_FIT_MINIMAP
!==================================


!#############################################
! PSF FIT, PROCEDURES FROM IBIS_PSF_AUX MODULE
!#############################################
!|| --------------------------------------------------------------------------
!|| SUB IBIS_PSF_MINIMIZE                                                  1.1
!||
!||   SUB IBIS_PSF_MINIMIZE
!||
!||   Perform minimization
!||   Input /output stored in variables of the module
!||     IBIS_PSF_FIT_PAR
!||
!||   Input:
!||
!||
!||   History:
!||   1.1 - 21/09/00 AG : C
!||                    
!||
!|| --------------------------------------------------------------------------

!==================================
SUBROUTINE IBIS_PSF_MINIMIZE(Status)
!==================================
  USE FIT_PAR_MODULE
  IMPLICIT NONE
  INTEGER :: Status
  ! Local parameters and arrays
  INTEGER ,parameter :: liw=mnnlpar+2, lw=MAX0((liw-2)*(liw-3)/2+12*liw,13)
  REAL      (kind=8), parameter :: bulim=15.D00, bllim=0.05D00
  INTEGER :: iw(liw)
  REAL      (kind=8) :: xx(mnnlpar),ww(lw)

  ! Local variables
  INTEGER :: ns,ncomp,i,ifail,ndof,ieqspre,idiffuse,ibound
  REAL      (kind=8) :: fsumsq, rchisq

  Status = 0
!  Bound type for non-linear parameters: individual bounds given in bl & bu
  ibound=0

!  Number of free linear & non linear parameters and set xx
  nparm=0
  lparm=0
  if(PSFP(1,2).eq.1)lparm=lparm+1
  ieqspre=0
  do ns=1,nsouf
     ncomp=1+(ns-1)*3
     if(spin(ns).eq.0)ieqspre=1
!  Intensity
     if(PSFP(ncomp+1,2).eq.1)lparm=lparm+1
!  Non linear parameters: position
     do i=2,3
        if(PSFP(ncomp+i,2).eq.1)THEN
           nparm=nparm+1
           xx(nparm)=PSFP(ncomp+i,1)
           bl(nparm)=DBLE(imapmin)
           bu(nparm)=DBLE(imapmax)
        endif
     enddo
  enddo
!  Non linear parameters: spatial resolution
  if(ieqspre.eq.1)THEN
     ncomp=3*nsouf+2
     if(PSFP(ncomp,2).eq.1)THEN
        nparm=nparm+1
        xx(nparm)=PSFP(ncomp,1)
        bl(nparm)=bllim
        bu(nparm)=bulim
     endif
  endif
!  Linear param: Diffuse component
  ncomp=3*nsouf+3
  IF (PSFP(ncomp,2).EQ.1) THEN
     lparm=lparm+1
     idiffuse=1
  ENDIF

  iter = 0
  fsumsq = 0.D0
 

!  Degrees of freedom
  ndof = (imapmax-imapmin+1)*(jmapmax-jmapmin+1)-nparm-lparm

!  NAG Routine for minimization: call FUNCT1 for model
  f04asfresult= 0
  ifail = -1
  call X04aaf(1,-1)
  call X04ABF(1,-1)
  CALL E04JAF(nparm,ibound,bl,bu,xx,fsumsq,iw,liw,ww,lw,ifail)
!!      if(ifail.ne.0)print*,'Minimization Terminated with IFAIL =',ifail
  f04asfresult = ifail
  if(f04asfresult.ne.0)then
     Status = 1
     return
  endif
  CALL FUNCT1(nparm,xx,fsumsq)
  if(f04asfresult.ne.0)then
     Status = 1
     return
  endif
  rchisq = fsumsq / DBLE(ndof)
 
!c      print*,' Best fit values (non-lin):'
!c      print*,xx
!c      print*,' Fit parameters (iter,ndof,nparm,fsumsq,rchisq): '
!c      print*,iter,ndof,(nparm+lparm),fsumsq,rchisq
!c      print*,''
 
  dof= REAL(ndof)
  spchisq= SNGL(rchisq)
!===============================
END SUBROUTINE IBIS_PSF_MINIMIZE
!================================

      
!|| --------------------------------------------------------------------------
!|| SUB IBIS_PSF_HESSIAN                                                   1.0
!||
!||   SUB IBIS_PSF_HESSIAN
!||
!||   Description:
!||     Computation of the Hessian of the chisquare to determine
!||     errors in the quadratic approximation
!||     It is assumed that H=Jt*J where Jt is the transposte of J and J the
!||     Jacobian of the residuals 
!||     Functional form of PSF used is the delta decoding PSF
!|| 
!||   Sub-History :
!||     1.0 - 21/09/00 AG : C,
!||     0.0   02/09/99 - AG : iw
!||
!|| --------------------------------------------------------------------------

!===================================
SUBROUTINE IBIS_PSF_HESSIAN(status)
!===================================

USE FIT_PAR_MODULE
USE PSF_SUB
IMPLICIT NONE
INTEGER :: Status

! Local parameters and arrays
INTEGER, parameter :: nmax=mnpar,ip=nmax,m=mapdim*mapdim
REAL      (kind=8) :: DJB(nmax,m),HINV(nmax,ip),H(nmax,ip)
!!      DIMENSION WSP(ip)
REAL      (kind=8) :: PSFX(mapdim),PSFY(mapdim),DERX(mapdim),DERY(mapdim)
REAL      (kind=8) :: DERSX(mapdim),DERSY(mapdim)

! Local variables
INTEGER :: ns,ncomp,i,j,k,ifail
INTEGER   (kind=8) :: inc,nnc,n1,n2,ich
REAL      (kind=8) :: x,y,dxc,dyc,sig,dbwei

INTEGER :: idim,imin1,ntpar,nsigp,l0,l,kmax,iw
REAL      (kind=8) :: derror,cur
REAL      (kind=8) :: G01FCF,G01ECF
REAL      (kind=8) :: deltachi,prob
REAL      (kind=8) :: ddiffuse(mapdim,mapdim)
     
Status = 0
! Test variable
iw=0

! Error 
derror=DBLE(error)
if(IABS(idec).eq.1)derror=DSQRT(derror*derror/(dbw*dbw))

ddiffuse = 0.

idim=imapmax-imapmin+1
imin1=imapmin-1
ntpar=nparm+lparm
if(PSFP(3*nsouf+2,2).EQ.1)nsigp=ntpar

DO j=jmapmin,jmapmax
   l0=(j-jmapmin)*idim-imin1
   DO i=imapmin,imapmax
      l=l0+i
      do k=1,ntpar
         DJB(k,l)=0.
      enddo
   enddo
enddo

inc=0
DO nnc=1,(nsouf+1)
   ns=nnc-1
   ncomp=IABS(2+(ns-1)*3)

! CALCUL of Jacobian DJB
   if(ns.eq.0)then
      IF(PSFP(1,2).EQ.1)THEN
         inc=inc+1
         DO j=jmapmin,jmapmax
            l0=(j-jmapmin)*idim-imin1
            DO i=imapmin,imapmax
               l=l0+i
               DJB(1,l)= 1./derror
            ENDDO
         ENDDO
      ENDIF
   elseif(ns.ge.1)then
      dxc=PSFP(ncomp+1,1)
      dyc=PSFP(ncomp+2,1)
      if(spin(ns).ne.0)then
         n1=i1
         n2=i2
      else
         n1=1
         n2=1
         sig=PSFP((nsouf*3+2),1)
         dbwei=1.D00
      endif
      DO ICH=n1,n2
         if(spin(ns).ne.0)then
            sig=DBLE(SIGSI(ich))
            dbwei=DBLE(weig(ich,ns))
         endif
         DO i = imapmin,imapmax
            x=DBLE(i)
            psfx(i) = psfd(x,dxc,sig,dbw)
            derx(i) = derpsf(x,dxc,sig,dbw)
            if(spin(ns).eq.0)dersx(i) = dersig(x,dxc,sig,dbw)
         ENDDO
         DO j = jmapmin,jmapmax
            y=DBLE(j)
            psfy(j) = psfd(y,dyc,sig,dbw)
            dery(j) = derpsf(y,dyc,sig,dbw)
            if(spin(ns).eq.0)dersy(j) = dersig(y,dyc,sig,dbw)
         ENDDO
         DO j = jmapmin,jmapmax
            l0=(j-jmapmin)*idim-imin1
            DO i = imapmin,imapmax
               l=l0+i
               k=inc
               IF(PSFP(ncomp,2).EQ.1)THEN
                  k=k+1
                  DJB(k,l) = DJB(k,l) + DBWEI * &
                       psfx(i)*psfy(j)/derror
               ENDIF
               IF(PSFP(ncomp+1,2).EQ.1)THEN
                  k=k+1
                  DJB(k,l) = DJB(k,l) + DBWEI * &
                       PSFP(ncomp,1)*derx(i)*psfy(j)/derror
               ENDIF
               IF(PSFP(ncomp+2,2).EQ.1)THEN
                  k=k+1
                  DJB(k,l) = DJB(k,l) + DBWEI * &
                      PSFP(ncomp,1)*psfx(i)*dery(j)/derror
               ENDIF
               IF(PSFP(nsouf*3+2,2).EQ.1)THEN
                  if(spin(ns).eq.0)then
                     DJB(nsigp,l) = DJB(nsigp,l) + PSFP(ncomp,1)* &
                          (psfx(i)*dersy(j)+dersx(i)*psfy(j))/derror
                  endif
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF(PSFP(ncomp,2).EQ.1)inc=inc+1
      IF(PSFP(ncomp+1,2).EQ.1)inc=inc+1
      IF(PSFP(ncomp+2,2).EQ.1)inc=inc+1
    endif
ENDDO

IF(PSFP(nsouf*3+2,2).EQ.1)inc=inc+1
!  Contribution Diffuse
ncomp=3*nsouf+3
IF(PSFP(ncomp,2).EQ.1)THEN
   inc=inc+1
   DO j=jmapmin,jmapmax
      l0=(j-jmapmin)*idim-imin1
      DO i=imapmin,imapmax
         l=l0+i
         DJB(inc,l)=DJB(inc,l)+ddiffuse(i,j)/derror
      ENDDO
   ENDDO
ENDIF


! CHECK OF THE CALCUL OF DERIVATIVES
!if(inc.ne.ntpar)print*,' Number free parameter wrong !'
!C         l0=(j-jmapmin)*idim-imin1
!C         DO i=imapmin,imapmax
!C            l=l0+i
!C            write(6,270) l,(DJB(k,l),k=1,inc)
!C         ENDDO
!C      ENDDO
!C 270  FORMAT(I5,10D13.6,/)

! Computing Hessian matrix
H=1.D00
kmax = (imapmax-imapmin+1) * (jmapmax-jmapmin+1)
DO j = 1,inc
   DO i = 1,j
      cur = 0.D0
      DO k = 1,kmax
         cur = cur + DJB(i,k) * DJB(j,k)
      ENDDO
      H(i,j) = cur
   ENDDO
ENDDO
!C      CALL F01CLF(H,DJB,DJB,nmax,ip,m,ifail)
 201  FORMAT(10D13.6,/)
      
!print*,''
!print*,' Hessian matrix : '
!do i=1,inc
!   write(*,201)(H(i,j),j=1,inc)
!enddo
      
     
! Inversion Hessian Matrix
ifail = 0 ! cholesky factorisation 
CALL F07FDF('U',inc,H,nmax,ifail)
Status = ifail
if(Status.ne.0)then
   !print *,' not pos-def matrix'
   return
endif

HINV= 0.D00
DO i = 1,inc
   HINV(i,i) = 1.D0
ENDDO
CALL F07FEF('U',inc,inc,H,nmax,HINV,nmax,ifail)
Status = ifail
if(Status.ne.0)then
   !print *,'F07FEF problem '
   return
endif
!C      CALL F01AAF(H,nmax,inc,HINV,nmax,WSP,ifail)
IF (iw == 1) THEN
   print*,''
   print*,' Curvature matrix : '
   do i=1,inc
      write(6,201)(HINV(i,j),j=1,inc)
   enddo

   print*,''
   print*,' Single parameter errors (68.3%) for free parameters: '
   print*,(DSQRT(HINV(i,i)),i=1,inc)
   print*,''
ENDIF
       
!      prob=1-SNGL(G01ECF('U',DREAL(deltachi),DBLE(np),ifail))
!      deltachi=G01CCF(DBLE(confl),np,ifail)
      deltachi=SNGL(G01FCF(DBLE(confl),DBLE(np),ifail))
      prob= SNGL(G01ECF('L',DBLE(deltachi),DBLE(np),ifail))

      IF (iw == 1) THEN
         print*,'' 
         print 2,confl,np,deltachi,prob

2        format(1X,' Deltachi at', &
         F6.3,' conf. level for ',I3,' param = ',F8.5,' (Prob=',F6.3,')')
      ENDIF

!  PUT errors in PSFE
      inc=0
      DO nnc=1,(nsouf+1)
         ns=nnc-1
         ncomp=IABS(2+(ns-1)*3)
         if(ns.eq.0)then
           IF(PSFP(1,2).EQ.1)THEN
             PSFE(1,1)=DSQRT(HINV(1,1))*DSQRT(DBLE(deltachi))
             inc=inc+1
           ENDIF
         elseif(ns.ge.1)then
           do i=1,3
             IF(PSFP(ncomp+i-1,2).EQ.1)THEN
                inc=inc+1
                PSFE(ncomp+i-1,1)=DSQRT(HINV(inc,inc))*DSQRT(DBLE(deltachi))
             ENDIF
           enddo
         endif
      enddo
      IF(PSFP(nsouf*3+2,2).EQ.1)THEN
         inc=inc+1
         PSFE(nsouf*3+2,1)=DSQRT(HINV(inc,inc))*DSQRT(DBLE(deltachi))
      ENDIF
      IF(PSFP(nsouf*3+3,2).EQ.1)THEN
         inc=inc+1
         PSFE(nsouf*3+3,1)=DSQRT(HINV(inc,inc))*DSQRT(DBLE(deltachi))
ENDIF
!===============================
END SUBROUTINE IBIS_PSF_HESSIAN
!===============================

!##############################################
! PSF FIT; PROCEDURES FROM PSF_FIT_AUX1 MODULE
!###############################################



