
!*********************
MODULE OVER_FIT_SUB_MODULE
!********************

INTERFACE

!GAUSSIAN FIT SUBROUTINES

SUBROUTINE OVER_FIT(FixedPerm,InCatNumber,FitDataType,anglestat,SKY,VAR,Filter,xc,yc,xs,ys,fitflux,locerr,&
     fluxerr,Status)
!------------------------------------------------------------------
IMPLICIT NONE
REAL(kind=4), dimension(:,:), pointer :: SKY,VAR
LOGICAL      ,dimension(:,:),pointer :: Filter
REAL(kind=8) :: xs,ys,xc,yc,fitflux,locerr,fluxerr
INTEGER :: FixedPerm,InCatNumber,FitDataType,anglestat,Status
END SUBROUTINE OVER_FIT

END INTERFACE
!*************************
END MODULE OVER_FIT_SUB_MODULE
!*************************




!***************************************************
! OVER FIT MODULES : DECLARATIONS AND INTERFACES
!***************************************************



!*********************
MODULE OVER_FIT_AUX
!*********************
INTERFACE

SUBROUTINE OVER_MINIMIZE(Status)
!-------------------------------------
IMPLICIT NONE
INTEGER :: Status
END SUBROUTINE OVER_MINIMIZE



SUBROUTINE OVER_HESSIAN(status)
!---------------------------------
IMPLICIT NONE
INTEGER :: Status
END SUBROUTINE OVER_HESSIAN


END INTERFACE
!*********************
END MODULE OVER_FIT_AUX
!*********************

!*********************
MODULE OVER_FIT_SUB
!*********************
INTERFACE

SUBROUTINE OVER_FIT_INIT(FixedPerm,InCatNumber,SKY,VAR,sources_to_fit,&
           xc,yc,xs,ys,fitdec)
!------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL(kind=4), dimension(:,:), pointer :: SKY,VAR
INTEGER :: FixedPerm,InCatNumber,sources_to_fit
REAL (kind=8),dimension(sources_to_fit) :: xs,ys
REAL(kind=8) :: xc,yc
INTEGER,dimension(4) :: fitdec
END SUBROUTINE OVER_FIT_INIT

SUBROUTINE OVER_FIT_MAP(sources_to_fit,xs,ys,fitflux,locerr,fluxerrtab,Status)
!-----------------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: sources_to_fit
REAL      (kind=8),dimension(sources_to_fit) :: xs,ys
REAL      (kind=8),dimension(sources_to_fit) :: fitflux,locerr,fluxerrtab
INTEGER :: Status
END SUBROUTINE OVER_FIT_MAP

END INTERFACE
!*********************
END MODULE OVER_FIT_SUB
!*********************



!#################################
! MODULE OVER_FIT_AUX SUBROUTINES
!#################################

!==================================
SUBROUTINE OVER_MINIMIZE(Status)
!==================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FIT_PAR_MODULE
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: Status
!LOCAL VARIABLES
!NAG definitions
INTEGER ,parameter :: liw=MaxNonLinPar+2, &
                       lw=MAX0((liw-2)*(liw-3)/2+12*liw,13)
INTEGER :: iw(liw)
REAL      (kind=8) :: ww(lw)
! non-linear params bounds
! angle,sigmax,sigmay
REAL      (kind=8), parameter :: Bulim(3)=(/360.0d0,2.D00,2.D00/), &
                                 Bllim(3)=(/-360.0d0,0.01D00,0.01D00/)
INTEGER :: ipar,ns,ndof,ibound,ifail,run,i,j
REAL      (kind=8) :: actNLPar(MaxNonLinPar)
REAL      (kind=8) :: fsumsq,redchisq
Logical :: good
CHARACTER(len=20) :: procname

Status = 0
procname = 'OVER_MINIMIZE'
!  Bound type for non-linear parameters: 
!     individual bounds given in Bl & bu
ibound=0

! non linear parameters
   nparm=0 
   ! angle,sigmax,sigmay
   do ipar=2,4
      if(fitperm(ipar) == 1)then
         nparm = nparm+1
         actNLPar(nparm) = FitVal(ipar)
         LowBounds(nparm)=Bllim(ipar-1)
         UpperBounds(nparm)=bulim(ipar-1)
      endif
   enddo

   ! sources ( x,y)
   do ns=1,nsouf
      ipar = 4+(ns-1)*3+2
      if(fitperm(ipar) == 1)then
         nparm = nparm+1
         actNLPar(nparm) = FitVal(ipar)
         LowBounds(nparm)=DBLE(imapmin)
         UpperBounds(nparm)=DBLE(imapmax)
      endif
      ipar = ipar+1
      if(fitperm(ipar) == 1)then
         nparm = nparm+1
         actNLPar(nparm) = FitVal(ipar)
         LowBounds(nparm)=DBLE(imapmin)
         UpperBounds(nparm)=DBLE(imapmax)
      endif
   enddo

!linear parameters
   lparm=0 
   !background
   if(fitperm(1) == 1)lparm = lparm+1
    ! sources intensity
   do ns=1,nsouf
      ipar = 4+(ns-1)*3+1
      if(fitperm(ipar) == 1)then
         lparm = lparm+1
      endif
   enddo
  iter = 0
  fsumsq = 0.D0
 

!  Degrees of freedom
  ndof = (imapmax-imapmin+1)*(jmapmax-jmapmin+1)-nparm-lparm



  !  NAG Routine for minimization: call FUNCT1 for model

  good = .true.
  f04asfresult= 0
  ifail = -1
  call X04aaf(1,-1)
  call X04ABF(1,-1)
  CALL E04JAF(nparm,ibound,LowBounds,UpperBounds,&
        actNLPar,fsumsq,iw,liw,ww,lw,ifail)
  
  f04asfresult = ifail
  if((ifail.ge.5).and.(ifail.le.8)) then
     ! doubts about minimum
     ifail = -1
     call X04aaf(1,-1)
      call X04ABF(1,-1)
      CALL E04JAF(nparm,ibound,LowBounds,UpperBounds,&
           actNLPar,fsumsq,iw,liw,ww,lw,ifail)
      f04asfresult = ifail
   endif
   if((ifail.ge.5).and.(ifail.le.8)) then
      call war_message(procname,&
           ' Some doubt about E04JAF',0,Status)
      f04asfresult=0
   endif

   if(f04asfresult.ne.0)then
      Status = 1
      return
   endif


   CALL FUNCT1(nparm,actNLPar,fsumsq)
   if(f04asfresult.ne.0)then
      Status = 1
      return
   endif

psfmap(:,:) = 0.
do j=1,ifitdim
   do i=1,ifitdim
      if(filmap(i,j))then
         psfmap(i,j)=resmap(i,j)/errormap(i,j)
      endif
   enddo
enddo

! verification of resid at the peak position
!if(abs(psfmap(ifd2,jfd2))>9.) Status = 20




!reduced chi2
redchisq = fsumsq / DBLE(ndof)
dof= REAL(ndof)
spchisq= SNGL(redchisq)
!==================================
END SUBROUTINE OVER_MINIMIZE
!==================================

!===================================
SUBROUTINE OVER_HESSIAN(status)
!===================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FIT_PAR_MODULE
USE OVER_INTPROC
IMPLICIT NONE
INTEGER :: Status

! Local variaBles

REAL      (kind=8) :: GRADIENT(MaxPar,mapdim,mapdim),HesRawNorm
REAL      (kind=8) :: Hessian(MaxPar,MaxPar),InvHessian(MaxPar,MaxPar)
INTEGER :: ns,i,j,k1,k2,k,ifail,npar,ipar,nc
REAL      (kind=8) :: x,y,cur
INTEGER :: iw
REAL      (kind=8) :: G01FCF,G01ECF
REAL      (kind=8) :: deltachi,prob
REAL      (kind=8) :: a,b,sina,cosa,sigmax,sigmay,hsou
Logical            :: HesNorm
character(len=20) :: procname   
  
Status = ISDC_OK
procname = 'OVER_HESSIAN'

sina = sin(fitval(2))
cosa = cos(fitval(2))
sigmax = FitVal(3)
sigmay = FitVal(4)

GRADIENT(:,:,:) = 0.

npar = 0
! bckg derivative
if(Fitperm(1).eq.1)then
   npar = npar+1
   DO j=jmapmin,jmapmax
   DO i=imapmin,imapmax
      if(filmap(i,j))GRADIENT(1,i,j)=1./errormap(i,j) ! bckg derivative
         do k=2,Maxpar
            GRADIENT(k,i,j)=0.0d0
         enddo
   enddo
   enddo
else
 DO j=jmapmin,jmapmax
   DO i=imapmin,imapmax
      do k=1,Maxpar
         GRADIENT(k,i,j)=0.0d0
      enddo
   enddo
   enddo
endif

! angle derivatives
if(Fitperm(2).eq.1)then
   npar = npar+1
   DO ns=1,nsouf
      ipar = 4+(ns-1)*3+1
      if(FitPerm(ipar) == 1)then ! active (fitted) source
         a = FitVal(ipar+1)
         b = FitVal(ipar+2)
         hsou = FitVal(ipar)
         DO j=jmapmin,jmapmax
         DO i=imapmin,imapmax
            x = dBle(i)
            y = dBle(j)
            if(filmap(i,j))&
            GRADIENT(npar,i,j)= GRADIENT(npar,i,j)+&
                 FitVal(ipar)*Der_Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y,5)&
                 /errormap(i,j)
            
         enddo
         enddo
      endif
   enddo
endif


! sigmax derivatives
if(Fitperm(3).eq.1)then
   npar = npar+1
   DO ns=1,nsouf
      ipar = 4+(ns-1)*3+1
      if(FitPerm(ipar) == 1)then ! active (fitted) source
         a = FitVal(ipar+1)
         b = FitVal(ipar+2)
         hsou = FitVal(ipar)
         DO j=jmapmin,jmapmax
         DO i=imapmin,imapmax
            x = dBle(i)
            y = dBle(j)
            if(filmap(i,j))&
            GRADIENT(npar,i,j)= GRADIENT(npar,i,j)+&
                 FitVal(ipar)*Der_Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y,3)&
                 /errormap(i,j)
         enddo
         enddo
      endif
   enddo
endif


! sigmay derivatives
if(Fitperm(4).eq.1)then
   npar = npar+1
   DO ns=1,nsouf
      ipar = 4+(ns-1)*3+1
      if(FitPerm(ipar) == 1)then ! active (fitted) source
         a = FitVal(ipar+1)
         b = FitVal(ipar+2)
         hsou = FitVal(ipar)
         DO j=jmapmin,jmapmax
         DO i=imapmin,imapmax
            x = dBle(i)
            y = dBle(j)
            if(filmap(i,j))&
            GRADIENT(npar,i,j)= GRADIENT(npar,i,j)+&
                 FitVal(ipar)*Der_Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y,4)&
                 /errormap(i,j)
         enddo
         enddo
      endif
   enddo
endif

! sources : intensity, posx, posy deriv
DO ns=1,nsouf
   ipar = 4+(ns-1)*3+1
   if(FitPerm(ipar) == 1)then ! active (fitted) source
      a = FitVal(ipar+1)
      b = FitVal(ipar+2)
      hsou = FitVal(ipar)
      DO j=jmapmin,jmapmax
      DO i=imapmin,imapmax
         x = dBle(i)
         y = dBle(j)
         ! int deriv
          if(filmap(i,j))then
             GRADIENT(npar+1,i,j)= GRADIENT(npar+1,i,j)+&
              Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y)&
                 /errormap(i,j)
             ! posx deriv
             GRADIENT(npar+2,i,j)= GRADIENT(npar+2,i,j)+&
              FitVal(ipar)*Der_Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y,1)&
                 /errormap(i,j)
             !posy deriv
             GRADIENT(npar+3,i,j)= GRADIENT(npar+3,i,j)+&
              FitVal(ipar)*Der_Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y,2)&
                 /errormap(i,j)
          endif
      enddo
      enddo
      npar = npar+3
   endif
enddo
      


! Computing Hessian matrix
Hessian(:,:)=1.D00
do k1=1,npar
do k2=1,npar
   cur = 0.D0
  DO j=jmapmin,jmapmax
  DO i=imapmin,imapmax
     cur = cur + GRADIENT(k1,i,j) * GRADIENT(k2,i,j)
  enddo
  enddo
  Hessian(k1,k2) = cur
ENDDO
enddo

!!$InvHessian = Hessian
!!$call F02AAF(InvHessian,MaxPar,npar,EigenVal,we,ifail)
!!$print *,' Eigenvalues min :',minval(Eigenval(1:npar))


HesNorm = .true.
do k1=1,npar
   HesRawNorm = sum(abs(hessian(k1,1:k1)))
   if(HesRawNorm < 1.e-15)then
       HesNorm = .false.
   endif
enddo
!SPR 3943
if(.not.HesNorm)then
   str250 = ' Too small Hessian array raw norm '
   call war_message(procname,str250,0,status)
   Status = 7
   return
endif

! Inversion of Hessian Matrix
! First: Cholesky factorisation 
!        result in Hessian
ifail = 0 
CALL F07FDF('U',npar,Hessian,MaxPar,ifail)

if(ifail.ne.0)then
   write(str250,'("F07FDF fail :",I6)')ifail
   call war_message(procname,str250,0,status)
   Status = ifail
   return
endif

InvHessian= 0.D00
DO i = 1,npar
   InvHessian(i,i) = 1.D0
ENDDO
!inversion  
ifail = 0 
CALL F07FEF('U',npar,npar,Hessian,MaxPar,InvHessian,MaxPar,ifail)

if(ifail.ne.0)then
   write(str250,'("F07FEF fail :",I6)')ifail
   call war_message(procname,str250,0,status)
   Status = ifail
   return
endif

       
ifail = 0 
deltachi=SNGL(G01FCF(DBLE(confl),1.0d0,ifail))

deltachi=0.83323169

!!problem with isdcmath 
!!$write(str250,*)'1 G01ECF(L,DBLE(deltachi),1.0d0,confl,deltachi,ifail : ',confl,deltachi,ifail
!!$call message(procname,str250,zeroerror,status)
!!$if(ifail.ne.0)then
!!$   write(str250,'("G01FCF fail :",I6)')ifail
!!$   call war_message(procname,str250,0,status)
!!$   Status = ifail
!!$  return
!!$endif
!!$
!!$prob= SNGL(G01ECF('L',DBLE(deltachi),1.0d0,ifail))
!!$write(str250,*)'2 G01ECF(L,DBLE(deltachi),1.0d0,deltachi,prob,ifail : ',deltachi,prob,ifail
!!$call message(procname,str250,zeroerror,status)
!!$if(ifail.ne.0)then
!!$  
!!$   write(str250,'("G01ECF fail :",I6)')ifail
!!$   call war_message(procname,str250,0,status)
!!$   Status = ifail
!!$  return
!!$endif
     
ipar=0
do i=1,4
   if(FitPerm(i).eq.1)then
      ipar = ipar+1
      FitErr(i)  =DSQRT(INVHESSIAN(ipar,ipar))*DSQRT(DBLE(deltachi))
   endif
enddo

do ns=1,nsouf
   nc = 4+(ns-1)*3+1
   if(FitPerm(nc).eq.1)then ! active source
      ipar = ipar+1
      FitErr(nc) = DSQRT(INVHESSIAN(ipar,ipar))*DSQRT(DBLE(deltachi))
      nc = nc+1
      ipar = ipar+1
      FitErr(nc) = DSQRT(INVHESSIAN(ipar,ipar))*DSQRT(DBLE(deltachi))
      nc = nc+1
      ipar = ipar+1
      FitErr(nc) = DSQRT(INVHESSIAN(ipar,ipar))*DSQRT(DBLE(deltachi))
   endif
enddo
!===============================
END SUBROUTINE OVER_HESSIAN
!===============================



!################################
!MODULE OVER_FIT_SUB SUBROUTINES 
!################################


!=============================================
SUBROUTINE OVER_FIT_INIT(FixedPerm,InCatNumber,SKY,VAR,sources_to_fit,&
           xc,yc,xs,ys,fitdec)
!============================================ 
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FITS_DECLARATIONS
USE FIT_PAR_MODULE
USE COMMON_PAR
USE OVER_FIT_AUX
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL(kind=4), dimension(:,:), pointer :: SKY,VAR
INTEGER :: FixedPerm,InCatNumber,sources_to_fit
REAL (kind=8),dimension(sources_to_fit) :: xs,ys
REAL(kind=8) :: xc,yc
!LOCAL VARIABLES
REAL(kind=8) :: xcsky,ycsky,varakt
INTEGER :: isky,jsky,iexfov,isx,jsx,jfitdim,ns,ipar
 INTEGER :: ncomp,i,j,l,k,ki,lj,ifgx,ifgy,ieqspre
INTEGER,dimension(4) :: fitdec



method = 3
 isky=SIZE(SKY,1)
 jsky=SIZE(SKY,2)
 xcsky=REAL(isky+1)/2.
 ycsky=REAL(jsky+1)/2.

 ifitdim=mapdim
 jfitdim=mapdim

 fitmap=0.
 filmap = .true.



!  fit map center in the image
 isx=NINT(xc)
 jsx=NINT(yc)

!  WRITE IMAGE SECTOR IN THE ARRAY FITMAP
 ifd2=ifitdim/2
 jfd2=jfitdim/2

! lower left corner of the fit map in the image
 ifgx=max(2,isx-ifd2)
 ifgy=max(2,jsx-jfd2)

! k,l - image coordiantes
! i,j - fit map coordinates
! fitmap,varmap,errormap  creation

 k=ifgx-1
 l=ifgy-1

 freedom = 0
 do i=1,ifitdim
    ki=k+i
    do j=1,jfitdim
       lj=l+j
!print *,i,j,ifitdim,jfitdim,ki,lj
       fitmap(i,j)=dble(SKY(ki,lj))/multicoeff
       varakt = dble(VAR(ki,lj))
       if(varakt.gt.0.)then
          varmap(i,j)=   varakt/multicoeff**2
          errormap(i,j) = sqrt(varakt) /multicoeff
          freedom = freedom + 1
       else
         filmap(i,j) = .false.
       endif
    enddo
 enddo
!print *,'over_fit_init no var '

 ! COORDINATES X0 Y0 - fit map lower left corner in the image
 x00=dble(ifgx-1)
 y00=dble(ifgy-1)

nsouf = sources_to_fit
fitval(:) = 0

! 1   - bckg       lin
! 2   - angle     nlin
! 3   - x width   nlin
! 4   - y width   nlin
! i   - int.       lin
! i+1 - posx      nlin
! i+2 - posy      nlin
fitperm(1:4) = fitdec(1:4)
fitval(1:4) = initval(1:4)

fitperm(5:4+nsouf*3) = 1



do ns=1,nsouf
   ipar = 4+(ns-1)*3+1
   fitval(ipar) = 0.0d0
   fitval(ipar+1) = xs(ns)-x00
   fitval(ipar+2) = ys(ns)-y00
     if((ns==1).and.(FixedPerm==1))then
        ! only first sou pos can be fixed
        ! and permission to fix ( not map fit)
        ! not for model
        if(FitPosFixed)then
           fitperm(ipar+1:ipar+2) = 0
           ! all source positions fixed
        else
           if(inCatNumber > 0)then
              !catalog source
              if(IncatFixed(InCatNumber))then
                 !position of this source fixed
                 fitperm(ipar+1:ipar+2) = 0
              endif
           endif!catalog source
        endif!FitPosFixed=false

         
     endif  ! only first sou pos can be fixed

enddo

! Confidence level for error estimation
confl= 0.683
!=============================
END SUBROUTINE OVER_FIT_INIT
!=============================



!=========================================================================
SUBROUTINE OVER_FIT_MAP(sources_to_fit,xs,ys,fitflux,locerr,fluxerr,Status)
!========================================================================
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FIT_PAR_MODULE
USE OVER_FIT_AUX

IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: sources_to_fit
REAL      (kind=8),dimension(sources_to_fit) :: xs,ys
REAL      (kind=8),dimension(sources_to_fit) :: fitflux,locerr,fluxerr
INTEGER :: Status
!LOCAL VARIABLES
INTEGER :: i,j,ipar,ns,ic,jc
REAL      (kind=8) :: fsumsq, bestchi,peakval

Status = 0
!  FILL THE DOUBLE PRECISION MINIMAP ARRAY
imapmin = 1
imapmax = ifitdim
jmapmin = 1
jmapmax = ifitdim
ic = ifitdim/2
jc = ifitdim/2

!  MINIMIZATION
CALL OVER_MINIMIZE(Status)




If(Status.ne.0)return

fsumsq=spchisq*dof   
bestchi=fsumsq



 ! too big gaussian width
   if( FitVal(3)> widthTol)then
      status = 4
      return
   endif
   if( FitVal(4)> widthTol)then
      status = 4
      return
   endif



! verification of 1st sourcefit  quality only
do ns=1,1
    ipar = 4+(ns-1)*3+2
    peakval =maxval(fitmap(ifd2-2:ifd2+2,jfd2-2:ifd2+2))

    if(RawFitFluxVerif)then
       if (FitVal(ipar-1) < peakval)then
          ! fited flux too low
          status = 15
          return
       endif
    endif
 
    if(peakval > 0.)then
       if(abs(FitVal(ipar-1)-peakval)/peakval > 1.5 )then
          ! peak too high
          status = 5
          return
       endif
    else
       status = 6
       return
    endif
    ! too big position error

   if(abs(xs(ns)-x00-FitVal(ipar)) > posTol)then
      status = 3
      return
   endif

   if(abs(ys(ns)-y00- Fitval(ipar+1))> posTol)then
      status = 3
      return
   endif

enddo


!  ERRORS IN THE ASSUMPTION OF QUADRATIC APPROXIMATION AND USING gaussian

CALL OVER_HESSIAN(Status)
If(Status.ne.0)return

! VERIFICATION

!too big fit error
 if( FitErr(3)> ErrTol)then
      status = 7
      return
   endif
   if( FitErr(4)> ErrTol)then
      status = 7
      return
   endif

do ns=1,nsouf
    ipar = 4+(ns-1)*3+2

    xs(ns)=SNGL(FitVal(ipar))+x00
    ys(ns)=SNGL(Fitval(ipar+1))+y00
    fitflux(ns)=SNGL(FitVal(ipar-1))
    !15.05.05 corrected
    locerr(ns)=SQRT(SNGL(FitErr(3)**2 + FitErr(4))**2)
    fluxerr(ns)=SNGL(FitErr(ipar-1))
 enddo



!==================================
END SUBROUTINE OVER_FIT_MAP
!==================================



!================================================================
SUBROUTINE OVER_FIT(FixedPerm,InCatNumber,FitDataType,anglestat,SKY,VAR,Filter,xc,yc,xs,ys,fitflux,locerr,&
     fluxerr,Status)
!================================================================== 
USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FIT_PAR_MODULE
USE FITS_DECLARATIONS
USE COMMON_PAR
USE OVER_FIT_SUB
USE GENERAL_SUB_MODULE ,only: ResidPeak
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
REAL(kind=4), dimension(:,:), pointer :: SKY,VAR
LOGICAL      ,dimension(:,:),pointer :: Filter
REAL(kind=8) :: xs,ys,xc,yc,fitflux,locerr,fluxerr
INTEGER :: FixedPerm,InCatNumber,FitDataType,anglestat,Status
! FitDataType = 0  data
! FitDataType = 1  model
!LOCAL VARIABLES
REAL(kind=8) ,dimension(:),pointer :: fitfluxtab,locerrtab,&
     fluxerrtab,xfit,yfit
REAL(kind=8) ,dimension(:),pointer :: x_deb_old,y_deb_old,x1,y1
INTEGER :: result,iloc(2),iok,idim,jdim,icsky,jcsky
INTEGER :: id,jd,im,jm,imap,jmap
INTEGER :: result1,nsou
REAL(kind=8) :: chi2red,chi2red1,peak,p1,dist,d1
INTEGER,dimension(4) :: fitdec
CHARACTER(len=20) :: procname
!FixedPerm 
! 0 - no position fixed permision ( for model fitting)
! 1 - position fixed permitted
procname = 'OVER_FIT'
Status = ISDC_OK
result = 0
result1 = 0
id = 5
im = mapdim/2
ALLOCATE(xfit(1),yfit(1),locerrtab(1),fitfluxtab(1),&
     fluxerrtab(1),stat = iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif 

fitfluxtab(:) = 0.0d0
locerrtab(:)  =  0.0d0
fluxerrtab(:) = 0.0d0
xfit(1) = xc
yfit(1) = yc

multicoeff = 1.

fitdec(:) = 1 ! b,alpha,sigmax,sigmay to be fitted

if(anglestat==0)fitdec(2) = 0
initval(1:4) = (/0.0d0,0.0d0,1.034d0,1.034d0/)

! ATTENTION : OTHER PROCEDURES OF FIT
! IF THE GAUSS_ORG.F90 - OVER_FIT_INIT AND OVER_FIT_MAP


CALL OVER_FIT_INIT(FixedPerm,IncatNumber,SKY,VAR,1,xc,yc,xfit,yfit,fitdec)    
CALL OVER_FIT_MAP(1,xfit,yfit,fitfluxtab,locerrtab,fluxerrtab,Status)
! Status : 
!         7 - too big error on width  
!        15 - fit flux too low 
if(Status.ne.ISDC_OK)then
   If((Status == 3).or.(Status==4).or.(Status==8).or.(Status==20))then
      ! 1st fit not accepted 
      ! 3 : fitted peak > 2 pixel dist. 
      ! 4 : gaussian width > 3 pixels 
      ! 5 : peak too high
      ! 8 : not converged
      ! 20 : residuals too big in the peak zone
      result=0
      result1 = 1
      xs = xc
      ys = yc
      fitflux = sky(nint(xs),nint(ys))/multicoeff
      locerr = nanf
      fluxerr = nanf
   else
   
      result = 1
      xs = xc;ys = yc
      fitflux = sky(nint(xs),nint(ys))
      locerr = nanf
      fluxerr = nanf
   endif
  
else
   xs = xfit(1)
   ys = yfit(1)
   fitflux = fitfluxtab(1)
   locerr = locerrtab(1)
   fluxerr = fluxerrtab(1)
endif

if(.not.FILTER(nint(xs),nint(ys)))then ! forbidden zone
   call WAR_MESSAGE(procName,&
        'Fit in forbidden zone - reset',0,Status)
   Status = 0
   result = 1
   xs = xc
   ys = yc
endif

chi2red = spchisq
chi2red1 = chi2red+1

if(FitDataType==0)then 
   ! model fitting - only one source model allowed
   return
endif


nsou = 1
ALLOCATE(x_deb_old(1),y_deb_old(1),&
           stat = iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif
x_deb_old(1)=xc
y_deb_old(1)=yc

if((result==0).and.(result1==1))then 
   ! first fit converged but not accepted
   Status = 0
   exclumap = .true.
   imap = mapdim/2+1
   jmap =imap
   exclumap(imap,jmap)=.false.

   call ResidPeak(im,jm,peak,dist)
   

   do while ((peak >10.).and.(chi2red >1.3).and.(chi2red < chi2red1)&
        .and.(Status .ge.ISDC_OK).and.(nsou < MaxSouNumber) )
 
      Status = ISDC_OK
      chi2red1 = chi2red
      Status  = 0
     
     
     ALLOCATE(x1(nsou),y1(nsou),stat = iok)
      IF (iok /= 0) then
         call MESSAGE(procName,'Allocation problem',AllocError,Status)
         return
      endif
      x1 = x_deb_old
      y1 = y_deb_old

      nsou = nsou+1
      
      DEALLOCATE(xfit,yfit,locerrtab,fitfluxtab,fluxerrtab,x_deb_old,y_deb_old)
      ALLOCATE(xfit(nsou),yfit(nsou),locerrtab(nsou),&
           fitfluxtab(nsou),x_deb_old(nsou),y_deb_old(nsou),&
           fluxerrtab(nsou),stat = iok)
      IF (iok /= 0) then
         call MESSAGE(procName,'Allocation problem',AllocError,Status)
         return
      endif
      x_deb_old(1:nsou-1) = x1
      y_deb_old(1:nsou-1) = y1
      x_deb_old(nsou) = im-imap+xc
      y_deb_old(nsou) = jm-jmap+yc
      deallocate(x1,y1)
      fitfluxtab(:) = 0.
      fluxerrtab(:) = 0.
      locerrtab(:)  =  0.
      xfit = x_deb_old
      yfit = y_deb_old
     
      Status = ISDC_OK

     
      CALL OVER_FIT_INIT(FixedPerm,inCatNumber,SKY,VAR,nsou,xc,yc,xfit,yfit,fitdec)    
      CALL OVER_FIT_MAP(nsou,xfit,yfit,fitfluxtab,locerrtab,fluxerrtab,Status) 
     
     
     
      if(Status.eq.0)then
         !fit accepted
         if((spchisq >0.8).and.(spchisq < chi2red1))then
            !better fit
            xs = xfit(1); ys = yfit(1)
            locerr = locerrtab(1)
            fitflux = fitfluxtab(1)
            fluxerr = fluxerrtab(1)
            result1 = 0
            chi2red = spchisq
         endif
      else
         chi2red1 = chi2red1+1
         
      endif
     
      exclumap(im,jm) = .false.

      call ResidPeak(im,jm,peak,dist)

  enddo



  if((Status.eq.ISDC_OK).and.(result1 ==0))then
     !one of fits  accepted
     chi2red = spchisq
     xs = xfit(1); ys = yfit(1)
     locerr = locerrtab(1)
     fitflux = fitfluxtab(1)
     fluxerr = fluxerrtab(1)
  else
     !no fits  accepted
     if(status == ISDC_OK)status = 1
     xs = xc;ys = yc
     fitflux = sky(nint(xs),nint(ys))/multicoeff
     locerr = nanf
     fluxerr = nanf
    
  endif ! no fits  accepted
     
 
endif !res==0 - first fit converged

DEALLOCATE(xfit,yfit,locerrtab,fitfluxtab,fluxerrtab,x_deb_old,y_deb_old)

!print *,' obtained chi2 red :',chi2red

!call FITS_FILE8(1,psfmap,1,nom_out='psfmap.fits')
!==========================
END SUBROUTINE OVER_FIT
!==========================
