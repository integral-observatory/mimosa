


! GENERAL DECLARATIONS
!*****************
MODULE COMMON_PAR
!*****************
IMPLICIT NONE
INTEGER,save :: method
REAL(kind=8) ,save    ::coeffnorm
!********************
END MODULE COMMON_PAR
!********************



!*******************
MODULE FIT_PAR_MODULE
!*******************
INTEGER,save :: funct1no = 0
INTEGER,parameter :: MaxSouNumber=3
INTEGER,parameter :: MaxParPerSou=5
INTEGER,parameter :: MaxNonLinParPerSou = 4
INTEGER,parameter :: MaxLinParPerSou = 1
INTEGER,parameter :: MaxParNum = MaxSouNumber*MaxParPerSou+1
INTEGER,parameter :: MaxNonLinparNum = MaxSouNumber* MaxNonLinParPerSou
INTEGER,parameter :: MaxLinparNum = MaxSouNumber* MaxLinParPerSou+1
!OVER FIT DEFINITIONS - BEGIN
INTEGER,parameter :: MaxPar = 4+MaxSouNumber*3,&
                     MaxNonLinPar = 3+MaxSouNumber*2,&
                     MaxModelNumber = MaxSouNumber+1
REAL(kind=8),dimension(4+3*MaxSouNumber)    :: fitval,fiterr
REAL(kind=8),dimension(4) :: initval
INTEGER,dimension(4+3*MaxSouNumber) :: fitperm
REAL(kind=8),dimension(MaxNonLinPar) :: LowBounds,UpperBounds
!OVER FIT DEFINITIONS - END
REAL(kind=8),parameter :: PosTol = 2. ! max position tol. in GaussFit
REAL(kind=8),parameter :: WidthTol = 3. ! max width  tol in GussFit
REAL(kind=8),parameter :: ErrTol = 3.  ! max error on width  in GussFit
REAL(kind=8),save ::spi
REAL(kind=8),dimension(MaxSouNumber,MaxParPerSou) :: &
     GaussSources,GaussErrors
REAL(kind=8) :: BkgLevel,BkgError
INTEGER, parameter :: mns=10, mnpar=3*mns+3, mnnlpar=2*mns+1
INTEGER, parameter :: mapdim= 40
REAL    (kind=8), save :: profile(mapdim)
REAL      (kind=8) :: MODEL(MaxModelNumber+1,mapdim,mapdim)
REAL    (kind=8), save :: fitmap(mapdim,mapdim), varmap(mapdim,mapdim)
REAL    (kind=8), save :: psfmap(mapdim,mapdim), resmap(mapdim,mapdim)
REAL    (kind=8), save :: errormap(mapdim,mapdim)
LOGICAL         , save :: filmap(mapdim,mapdim),exclumap(mapdim,mapdim)
REAL    (kind=8), save :: bl(mnnlpar),bu(mnnlpar),dbw
REAL    (kind=8), save :: bl_gen(MaxNonLinparNum),bu_gen(MaxNonLinparNum)
REAL    (kind=8), save :: PSFP(mnpar,2),PSFE(mnpar,2)
INTEGER, save :: nsouf,nparm,lparm,freedom
INTEGER, save :: ifitdim,imapmin,imapmax,jmapmin,jmapmax,ifd2,jfd2
REAL    (kind=8), save :: error,x00,y00,xta,yta
REAL    (kind=8), save :: dof,spchisq,confl,multicoeff
REAL    (kind=8), save :: spin(mns),SIGSI(95),weig(95,mns)
INTEGER, save :: i1,i2,idec,np,iter
INTEGER, save :: iDetTyp, iPsfTest
REAL    (kind=8), save :: f04asfresult,f04asfprint,rmean,esd,stdev,xsta,ysta,time
CONTAINS 

!..........................
FUNCTION Gauss(x,mu,sigma)
!...........................
IMPLICIT NONE
REAL(kind=8) :: x,mu,sigma,Gauss

REAL(kind=8) :: coeff,arg

if(sigma.eq.0.)then
  ! print *,' error in sigma '
   Gauss = 0.
   return
endif
coeff = 1.0d0/spi/sigma
arg = -0.50d0*(x-mu)**2/sigma**2
gauss = coeff*exp(arg)

!..................
END FUNCTION Gauss
!..................

!***********************
END MODULE FIT_PAR_MODULE
!************************


!***********************
MODULE OVER_INTPROC
!***********************

INTERFACE

FUNCTION Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y)
!---------------------------------------------------
IMPLICIT NONE
REAL(kind=8) :: a,b,sina,cosa,sigmax,sigmay,x,y,Over_Gauss
END FUNCTION Over_Gauss

FUNCTION Der_Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y,ktora)
!------------------------------------------------------------
IMPLICIT NONE
REAL(kind=8) :: a,b,sina,cosa,sigmax,sigmay,x,y,Der_Over_Gauss
INTEGER :: ktora
END FUNCTION Der_Over_Gauss

END INTERFACE

!***********************
END MODULE OVER_INTPROC
!***********************


!***************
MODULE PSF_SUB
!***************
INTERFACE

FUNCTION PSFD(x,xc,sigx,w)
!-------------------------
IMPLICIT NONE
REAL      (kind=8) :: PSFD
REAL      (kind=8) :: xc,x,sigx,w
END FUNCTION PSFD

FUNCTION  PSFF(x,xc,sigx,w)
!-------------------------
IMPLICIT NONE
REAL      (kind=8) :: PSFF
REAL      (kind=8) :: xc,x,sigx,w
END FUNCTION  PSFF

FUNCTION DERPSF(x,xc,sigx,w)
!-------------------------
IMPLICIT NONE
REAL      (kind=8) :: DERPSF
REAL      (kind=8) :: xc,x,sigx,w
END FUNCTION DERPSF

FUNCTION DERSIG(x,xc,sigx,w)
!-------------------------
IMPLICIT NONE
REAL      (kind=8) :: DERSIG
REAL      (kind=8) :: xc,x,sigx,w
END FUNCTION DERSIG

FUNCTION TPSF(x,y,ns1,ns2)
!-------------------------
IMPLICIT NONE
REAL      (kind=8) :: TPSF, x,y
INTEGER :: ns1,ns2
END FUNCTION TPSF

END INTERFACE
!***************
END MODULE PSF_SUB
!***************



!***********************
 MODULE FIT_FUNCT1
!***********************

INTERFACE
  
SUBROUTINE FUNCT1(n,xc,fc)
!-------------------------
IMPLICIT NONE
INTEGER  :: n
REAL      (kind=8) :: xc(n),fc
END SUBROUTINE FUNCT1

END INTERFACE

!***********************
END  MODULE FIT_FUNCT1
!***********************

!****************************
MODULE PSF_FIT_AUX1
!****************************
INTERFACE

SUBROUTINE PSF_FUNCT1(n,xc,fc)
!-------------------------
IMPLICIT NONE
INTEGER :: n
REAL   (kind=8) :: xc(n),fc
END SUBROUTINE PSF_FUNCT1

END INTERFACE
!****************************
END MODULE PSF_FIT_AUX1
!****************************

!****************************
MODULE OVER_FIT_AUX1
!****************************
INTERFACE

SUBROUTINE OVER_FUNCT1(n,xc,fc)
!--------------------------
IMPLICIT NONE
INTEGER              :: n
REAL      (kind=8) :: xc(n),fc
END SUBROUTINE OVER_FUNCT1
END INTERFACE
!****************************
END MODULE OVER_FIT_AUX1
!****************************


!#################################
! MODULE OVER_INTPROC SUBROUTINES
!#################################

!===================================================
FUNCTION Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y)
!===================================================
IMPLICIT NONE
REAL(kind=8) :: a,b,sina,cosa,sigmax,sigmay,x,y,Over_Gauss

REAL(kind=8) :: argx,argy

if(sigmax*sigmay.le.0.)then
!!$   print *,' a,b,sina,cosa,sigmax,sigmay,x,y',&
!!$           a,b,sina,cosa,sigmax,sigmay,x,y,&
!!$           '--> error in sigma '
   Over_Gauss = 0.
   return
endif
argx = ((x-a)*cosa-(y-b)*sina)/sigmax
argy = ((x-a)*sina+(y-b)*cosa)/sigmay

Over_Gauss = exp(-0.50d0*(argx**2+argy**2))

!========================
END FUNCTION Over_Gauss
!========================


!=============================================================
FUNCTION Der_Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y,ktora)
!=============================================================
IMPLICIT NONE
REAL(kind=8) :: a,b,sina,cosa,sigmax,sigmay,x,y,Der_Over_Gauss
INTEGER :: ktora
! 1 - a
! 2 - b
! 3 - sigmax
! 4 - sigmay
! 5 - angle
REAL(kind=8) :: argx,argy,func

if(sigmax*sigmay.le.0.)then
!!$   print *,' a,b,sina,cosa,sigmax,sigmay,x,y',&
!!$           a,b,sina,cosa,sigmax,sigmay,x,y,&
!!$           '--> error in sigma '
   Der_Over_Gauss = 0.
   return
endif
!!$argx = (a+x*cosa-y*sina)/sigmax
!!$argy = (b+x*sina+y*cosa)/sigmay
argx = ((x-a)*cosa-(y-b)*sina)/sigmax
argy = ((x-a)*sina+(y-b)*cosa)/sigmay

func = exp(-0.50d0*(argx**2+argy**2))
select case (ktora)
case(1)
!!$   Der_Over_Gauss= -func*argx*(argx-a/sigmax)
   Der_Over_Gauss= func*(argx*cosa/sigmax + argy*sina/sigmay)
case(2)
!!$   Der_Over_Gauss= -func*argy*(argy-b/sigmay)
   Der_Over_Gauss= -func*(argx*sina/sigmax-argy*cosa/sigmay )
case(3)
!!$   Der_Over_Gauss= func*argx**2/sigmax
   Der_Over_Gauss= func*argx**2/sigmax
case(4)
!!$   Der_Over_Gauss= func*argy**2/sigmay
   Der_Over_Gauss= func*argy**2/sigmay
case(5)
!!$   Der_Over_Gauss= -func*&
!!$        (-argx*(x*sina+y*cosa)/sigmax+argy*(x*cosa-y*sina)/sigmay)
   Der_Over_Gauss= func*argx*argy*((sigmay**2-sigmax**2)/sigmax/sigmay)
end select
!=============================
END FUNCTION  Der_Over_Gauss
!=============================

!**************************
MODULE GENERAL_SUB_MODULE
!**************************

INTERFACE


SUBROUTINE ResidPeak(im,jm,peak,dist)
!-----------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: im,jm
REAL      (kind=8) :: peak,dist
END SUBROUTINE ResidPeak

END INTERFACE

!**************************
END MODULE GENERAL_SUB_MODULE
!**************************



!#############################
! SUBROUTINE FUNCT1
!############################

!=========================
SUBROUTINE FUNCT1(n,xc,fc)
!=========================
USE COMMON_PAR
USE OVER_FIT_AUX1
USE PSF_FIT_AUX1

IMPLICIT NONE
INTEGER  :: n
REAL      (kind=8) :: xc(n),fc

select case(method)
case(1)
  call PSF_FUNCT1(n,xc,fc)

case(3)
   call OVER_FUNCT1(n,xc,fc)
end select

!=========================
END SUBROUTINE FUNCT1
!=========================

!===============================
SUBROUTINE PSF_FUNCT1(n,xc,fc)
!===============================
USE FIT_PAR_MODULE
USE PSF_SUB
IMPLICIT NONE

! Local parameters and arrays
INTEGER ,parameter :: ia=mnpar-mnnlpar
INTEGER :: n
REAL      (kind=8) :: xc(n),xx(ia)
REAL      (kind=8) :: dbmap2(mapdim,mapdim)
REAL      (kind=8) :: PHF(mns+1,mapdim,mapdim),AF(ia,ia),BF(ia)
REAL      (kind=8) :: DX(mns,mapdim),DY(mns,mapdim)
REAL      (kind=8) :: W1(ia),W2(ia)

! Local variables
INTEGER :: ns,ncomp,i,j,k,ifail
INTEGER :: iadec,nc,inc,nnc,n1,n2,ich,ic,id,idf
REAL      (kind=8) :: x,y,dxc,dyc,sig,dbwei,da,s1,s2,s,ski2,ski,fc


! Dec type
iadec=IABS(idec)

!  Non linear parameters 
k=0
do ns=1,nsouf
   ncomp=1+(ns-1)*3
   do i=2,3
      if(PSFP(ncomp+i,2).eq.1)THEN
         k=k+1
         PSFP(ncomp+i,1)=xc(k)
      endif
   enddo
enddo
ncomp=3*nsouf+2
if(PSFP(ncomp,2).eq.1)THEN
   k=k+1
   PSFP(ncomp,1)=xc(k)
endif

DO i=imapmin,imapmax
   DO j=jmapmin,jmapmax
      dbmap2(i,j)=DBLE(fitmap(i,j))
   ENDDO
ENDDO
 
 
! COEFFICIENTS POUR LE CALCUL DES MOINDRES CARRES
nc=lparm
do k=1,nc
   do i=imapmin,imapmax
      do j=jmapmin,jmapmax
         PHF(k,i,j) = 0.D00
      enddo
   enddo
enddo
      
inc=0
DO nnc=1,(nsouf+1)
   ns=nnc-1
   ncomp=IABS(2+(ns-1)*3)

!  VALUE OF PSF
   if(ns.ge.1)then
      dxc=PSFP(ncomp+1,1)
      dyc=PSFP(ncomp+2,1)
      IF(PSFP(ncomp,2).EQ.1)inc=inc+1
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
            if(iadec.eq.2)then
               DX(nnc,i) = PSFD(x,dxc,sig,dbw)
            elseif(iadec.eq.1)then
               DX(nnc,i) = PSFF(x,dxc,sig,dbw)
            endif
         ENDDO
         DO j = jmapmin,jmapmax
            y=DBLE(j)
            if(iadec.eq.2)then
               DY(nnc,j) = PSFD(y,dyc,sig,dbw)
            elseif(iadec.eq.1)then
               DY(nnc,j) = PSFF(y,dyc,sig,dbw)
            endif
         ENDDO
         IF(PSFP(NCOMP,2).EQ.1)THEN
            DO i=imapmin,imapmax
               DO j=jmapmin,jmapmax
                  PHF(inc,i,j) = PHF(inc,i,j) + &
                       dbwei*DX(nnc,i)*DY(nnc,j)
               ENDDO
            ENDDO
         ELSE
            da=PSFP(ncomp,1)
            DO i=imapmin,imapmax
               DO j=jmapmin,jmapmax
                  dbmap2(i,j) = dbmap2(i,j) - da* &
                       dbwei*DX(nnc,i)*DY(nnc,j)
               ENDDO
            ENDDO
         ENDIF
    ENDDO
     else
        IF(PSFP(1,2).EQ.1)THEN
           inc=inc+1
           DO i=imapmin,imapmax
              DO j=jmapmin,jmapmax
                 PHF(1,i,j) = 1.0D00
              ENDDO
           ENDDO
        ELSE
           da=PSFP(ncomp,1)
           DO i=imapmin,imapmax
              DO j=jmapmin,jmapmax
                 dbmap2(i,j)=dbmap2(i,j)-da
              ENDDO
           ENDDO
        ENDIF
      endif
ENDDO

!  Initializes ARRAYS AF and BF of EQUAT to 0
DO IC = 1,NC
   BF(IC) = 0.0
   DO ID = 1,NC
      AF(IC,ID) = 0.0
   ENDDO
ENDDO

DO IC = 1,NC
   S1 = 0.0D0
   S2 = 0.0D0
   DO j = jmapmin,jmapmax
      DO i = imapmin,imapmax
         S1 = S1 + (PHF(IC,I,J))**2
         S2 = S2 + PHF(IC,I,J)*dbmap2(i,j)
      ENDDO
   ENDDO
   AF(IC,IC) = AF(IC,IC) + S1
   BF(IC) = BF(IC) + S2
!!!         BF(IC) = BF(IC) - S2       !
ENDDO
NNC = NC -1
DO IC = 1,NNC
   IDF = IC + 1
   DO ID = IDF,NC
      S = 0.0
      DO j = jmapmin,jmapmax
         DO i = imapmin,imapmax
            S = S + PHF(IC,I,J) * PHF(ID,I,J)
         ENDDO
      ENDDO
      AF(IC,ID) = AF(IC,ID) + S
      AF(ID,IC) = AF(ID,IC) + S
   ENDDO
ENDDO


100   FORMAT(1X,6(2X,1PE8.1))
   
! CALCUL OF CONTRIBUTIONS

!  Trough NAG Routine (negative values permitted)
ifail = -1
call X04aaf(1,-1)
call X04ABF(1,-1)
CALL F04ASF(AF,ia,BF,nc,XX,W1,W2,ifail)
f04asfresult = ifail
if(f04asfresult.ne.0)then 
   if(f04asfprint.gt.0)then
      print *,' uncorrect result from F04asf ifail = ',ifail
      print *,' matrix : '
      do i=1,ia
         print *,(af(i,j),j=1,ia)
      enddo
      print *,' b:'
      print *,' bf(1:ia)'
   endif
   return
endif
! CALCUL OF SUM OF SQUARES AND OF KI2
 
iter = iter + 1
SKI2 = 0.0
DO j = jmapmin,jmapmax
   DO i = imapmin,imapmax
      S = 0.0D0
      DO IC = 1,NC
         S = S + PHF(IC,i,j) * XX(IC)
      ENDDO
      SKI = dbmap2(i,j) - S
      SKI2 = SKI2 + SKI**2
      if(error >0)resmap(i,j)=REAL(ski)/error
   ENDDO
ENDDO
if(error > 0.)&
ski2 = ski2 / DBLE(error)**2
fc = ski2
 
! PARAMETERS AT EACH STEP

ic=0
if(PSFP(1,2).eq.1)then
   ic=1
   PSFP(1,1)=xx(ic)
endif
do ns=1,nsouf
   ncomp=2+(ns-1)*3
   if(PSFP(ncomp,2).eq.1)THEN
      ic=ic+1
      PSFP(ncomp,1)=xx(ic)
   endif
enddo
     
!      rchisq = ski2 / DFLOAT(ndof)
!      print*,' Chisquare =',sky2
!      print*,' Parameters :',(xc(k),k=1,nparm),(xx(k),k=1,lparm)

!======================
END SUBROUTINE PSF_FUNCT1
!========================

!===============================
SUBROUTINE OVER_FUNCT1(n,NonLinPars,fc)
!===============================
USE FIT_PAR_MODULE
USE OVER_INTPROC
IMPLICIT NONE

! Local parameters and arrays
INTEGER ,parameter :: LinparNumber=MaxPar-MaxNonLinPar
INTEGER :: n
REAL      (kind=8) :: NonLinPars(n),linPars(LinparNumber)
REAL      (kind=8) :: dbmap2(mapdim,mapdim)
REAL      (kind=8) :: AF(LinparNumber,LinparNumber),BF(LinparNumber)
! rotated gauss models
REAL      (kind=8) :: DXY(MaxModelNumber,mapdim,mapdim)
! NAG arrays
REAL      (kind=8) :: W1(LinparNumber),W2(LinparNumber)

! Local variaBles
INTEGER :: ns,i,j,k,ifail,ipar,npar,nnc
INTEGER :: inc,ic,id,idf
REAL      (kind=8) :: x,y,da,s1,s2,s,ski2,ski,fc
REAL      (kind=8) :: a,b,sina,cosa,sigmax,sigmay




!   Non linear parameters 

npar = 0

! angle,sigmax,sigmay
do ipar=2,4
   if(fitperm(ipar) == 1)then
      npar = npar+1
      fitval(ipar) = NonLinPars(npar)
   endif
enddo

! source positions
do ns=1,nsouf
   ipar = 4+(ns-1)*3+2
   do i=ipar,ipar+1
      if(fitperm(i).eq.1)THEN
         npar = npar+1
         FitVal(i)=NonLinPars(npar)
      endif
   enddo
enddo

! current data - constant if all FitPerm = 1
DO i=imapmin,imapmax
   DO j=jmapmin,jmapmax
      dbmap2(i,j)=DBLE(fitmap(i,j))
   ENDDO
ENDDO
 
 
! COEFFICIENTS POUR LE CALCUL DES MOINDRES CARRES

sina = sin(fitval(2))
cosa = cos(fitval(2))
sigmax = FitVal(3)
sigmay = FitVal(4)

do k=1,lparm
   do i=imapmin,imapmax
      do j=jmapmin,jmapmax
         MODEL(k,i,j) = 0.D00
      enddo
   enddo
enddo

inc = 0 
! background term - 1st model


IF(Fitperm(1).EQ.1)THEN
   inc=inc+1
   DO i=imapmin,imapmax
      DO j=jmapmin,jmapmax
           if(filmap(i,j))then
              MODEL(1,i,j) = 1.0D00
           endif
      ENDDO
   ENDDO
ELSE
   da=FitVal(1)
   DO i=imapmin,imapmax
      DO j=jmapmin,jmapmax
          if(filmap(i,j))then
             dbmap2(i,j)=dbmap2(i,j)-da
          endif
      ENDDO
   ENDDO
ENDIF

      
! source intensities
DO ns=1,nsouf
   ipar = 4+(ns-1)*3+1
   ! model calcul
   a = DBLE(Fitval(ipar+1))
   b = DBLE(FitVal(ipar+2))
   DO i = imapmin,imapmax
      x=DBLE(i)
      DO j = jmapmin,jmapmax
         y=DBLE(j)
         if(filmap(i,j))then
            DXY(ns,i,j) = Over_Gauss(a,b,sina,cosa,sigmax,sigmay,x,y)
         endif
      ENDDO
   enddo
   if(fitperm(ipar) == 1)then !linear  parameter to fit - ns-th model
      inc=inc+1
      DO i=imapmin,imapmax
      DO j=jmapmin,jmapmax
         if(filmap(i,j))then
            MODEL(inc,i,j) =  DXY(ns,i,j)
         endif
      ENDDO
      ENDDO
   else ! linear  parameter fixed
      da=FitVal(ipar)
      DO i=imapmin,imapmax
      DO j=jmapmin,jmapmax
           if(filmap(i,j))then
              dbmap2(i,j) = dbmap2(i,j) - da* &
                   DXY(inc,i,j)
           endif
      ENDDO
      ENDDO
   endif
 ENDDO


 


!  Initializes ARRAYS AF and BF of EQUAT to 0
DO IC = 1,lparm
   BF(IC) = 0.0d0
   DO ID = 1,lparm
      AF(IC,ID) = 0.0d0
   ENDDO
ENDDO

DO IC = 1,lparm
   S1 = 0.0D0
   S2 = 0.0D0
   DO j = jmapmin,jmapmax
      DO i = imapmin,imapmax
         if(filmap(i,j))then
            S1 = S1 + (MODEL(IC,I,J))**2/varmap(i,j)
            S2 = S2 + MODEL(IC,I,J)*dbmap2(i,j)/varmap(i,j)
         endif
      ENDDO
   ENDDO
   AF(IC,IC) = AF(IC,IC) + S1
   BF(IC) = BF(IC) + S2

ENDDO

NNC = lparm -1
DO IC = 1,NNC
   IDF = IC + 1
   DO ID = IDF,lparm
      S = 0.0
      DO j = jmapmin,jmapmax
         DO i = imapmin,imapmax
             if(filmap(i,j))then
                S = S + MODEL(IC,I,J) * MODEL(ID,I,J)/varmap(i,j)
             endif
         ENDDO
      ENDDO
      AF(IC,ID) = AF(IC,ID) + S
      AF(ID,IC) = AF(ID,IC) + S
   ENDDO
ENDDO



   
! CALCUL OF CONTRIBUTIONS

!  Trough NAG Routine (negative values permitted)
ifail = -1
call X04aaf(1,-1)
call X04ABF(1,-1)
CALL F04ASF(AF,LinParNumber,BF,lparm,LINPARS,W1,W2,ifail)
f04asfresult = ifail
if(f04asfresult.ne.0)then 
   if(f04asfprint.gt.0)then
      print *,' uncorrect result from F04asf ifail = ',ifail
      print *,' matrix : '
      do i=1,LinparNumber
         print *,(af(i,j),j=1,LinparNumber)
      enddo
      print *,' b:'
      print *,' bf(1:LinparNumber)'
   endif
   return
endif
! CALCUL OF SUM OF SQUARES AND OF KI2
 
iter = iter + 1
SKI2 = 0.0
DO j = jmapmin,jmapmax
   DO i = imapmin,imapmax
      S = 0.0D0
      DO IC = 1,lparm
         S = S + MODEL(IC,i,j) * LINPARS(IC)
      ENDDO
      if(filmap(i,j))then
         SKI = dbmap2(i,j) - S
         SKI2 = SKI2 + SKI**2/varmap(i,j)
         resmap(i,j)=ski
      endif
   ENDDO
ENDDO

fc = ski2
 
! PARAMETERS AT EACH STEP

ic=0
if(FitPerm(1).eq.1)then
   ic=1
   FitVal(1)=linPars(1)
endif
do ns=1,nsouf
    ipar = 4+(ns-1)*3+1
  
   if(FitPerm(ipar).eq.1)THEN
      ic=ic+1
      Fitval(ipar)=linPars(ic)
   endif
enddo
     

!print '("over_funct1 ",2(F15.6,1x),2x,f20.10)',fitval(6:7),fc
!======================
END SUBROUTINE OVER_FUNCT1
!========================

!=============================================
SUBROUTINE ResidPeak(im,jm,peak,dist)
!============================================ 
USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FIT_PAR_MODULE
USE COMMON_PAR
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: im,jm
REAL      (kind=8) :: peak,dist
!LOCAL VARIABLES
INTEGER :: id,jd,imap,jmap
REAL      (kind=8) :: p1,d1

imap = mapdim/2+1
jmap =imap
 peak = 0.
 im = 0
 jm = 0
 dist = 100.
 do id=1,mapdim
    do jd=1,mapdim
      if( exclumap(id,jd))then
         p1 = psfmap(id,jd)
         d1 = sqrt((real(id-imap))**2+(real(jd-jmap))**2)
         if((p1 > 7.).and.(p1 > peak))then
            im = id
            jm = jd
            peak = p1
            dist = d1
         endif
      endif
   enddo
 enddo

!===================
END SUBROUTINE ResidPeak
!===================


!#########################################
! PSF FIT ; PROCEDURES FROM PSF_SUB MODULE
!#########################################
!|| -----------------------------------------------------------------------
!||                           PSF  FUNCTIONS
!|| -----------------------------------------------------------------------
!|| 
!||
!|| -----------------------------------------------------------------------
!||     FUNCTION PSFD(x,xc,sigx,w)
!||  
!||  Delta decoding PSF
!|| 
!==========================
FUNCTION PSFD(x,xc,sigx,w)
!==========================
IMPLICIT NONE

      REAL      (kind=8) :: S15ADF
      REAL      (kind=8) :: PSFD,SQRT2,SQRT2SI,ROOT2PI,ROOT2PIW,ROOTPI
      REAL      (kind=8) :: xc,x,sigx,w,u0,u1,u2,u3
      REAL      (kind=8) :: erfcu0, erfcu1, erfcu2, erfcu3
      REAL      (kind=8) :: derfc12, derfc01, derfc23
      INTEGER :: ifail
            
      SQRT2    = DSQRT(2.0d0)
      SQRT2SI  = SQRT2*sigx
      ROOT2PI  = 2.506628275D0
      ROOT2PIW = ROOT2PI*w
      ROOTPI   = 1.772453851D00
      
! Respect to the formula notations:  u0=u2  u1=u+ u2=u- u3=u1
       u0=(x-xc-(w+1)/2)/SQRT2SI
       u1=(x-xc-(w-1)/2)/SQRT2SI
       u2=(x-xc+(w-1)/2)/SQRT2SI
       u3=(x-xc+(w+1)/2)/SQRT2SI

      ifail=1
      erfcu0=S15ADF(u0,ifail)
      !if(ifail.eq.1)PRINT *,' ifail =',ifail
      ifail=1
      erfcu1=S15ADF(u1,ifail)
      !if(ifail.eq.1)PRINT *,' ifail =',ifail
      ifail=1
      erfcu2=S15ADF(u2,ifail)
     ! if(ifail.eq.1)PRINT *,' ifail =',ifail
      ifail=1
      erfcu3=S15ADF(u3,ifail)
      !if(ifail.eq.1)PRINT *,' ifail =',ifail

      derfc12=erfcu1-erfcu2
      derfc01=erfcu0-erfcu1
      derfc23=erfcu2-erfcu3

      PSFD = sigx/ROOT2PIW*(DEXP(-u0**2)-DEXP(-u1**2)-DEXP(-u2**2) + &
             DEXP(-u3**2)) + .5D00/w*derfc12 - sigx/(w*SQRT2)*       &
             (u0*derfc01-u3*derfc23)
!==================
END FUNCTION PSFD
!==================


!|| -----------------------------------------------------------------------
!||     FUNCTION  PSFF(x,xc,sigx,w)
!|| 
!||  Fine cross correlation PSF
!|| 
      
!==========================
FUNCTION  PSFF(x,xc,sigx,w)
!==========================
IMPLICIT NONE

      REAL      (kind=8) :: S15ADF
      REAL      (kind=8) :: PSFF,SQRT2,SQRT2SI,ROOT2PI,ROOT2PIW,ROOTPI
      REAL      (kind=8) :: xc,x,sigx,w,umw,uc,upw
      REAL      (kind=8) :: erfcuc, erfcumw, erfcupw, derfc1, derfc2
      INTEGER :: ifail
      
      SQRT2    = DSQRT(2.0d0)
      SQRT2SI  = SQRT2*sigx
      ROOT2PI  = 2.506628275D0
      ROOT2PIW = ROOT2PI*w
      ROOTPI   = 1.772453851D00
      
       umw=(x-xc+w)/SQRT2SI
       uc=(x-xc)/SQRT2SI
       upw=(x-xc-w)/SQRT2SI
       
      ifail=1
      erfcuc=S15ADF(uc,ifail)
     ! if(ifail.eq.1)PRINT *,' ifail =',ifail
      ifail=1
      erfcumw=S15ADF(umw,ifail)
     ! if(ifail.eq.1)PRINT *,' ifail =',ifail
      ifail=1
      erfcupw=S15ADF(upw,ifail)
     ! if(ifail.eq.1)PRINT *,' ifail =',ifail

      derfc1=erfcuc-erfcumw
      derfc2=erfcupw-erfcuc

      PSFF = sigx/(ROOT2PIW) *( ROOTPI*(umw*derfc1 - upw*derfc2) &
             -  2.D00*DEXP(-uc**2) + DEXP(-umw**2) + DEXP(-upw**2) )
!==================
END FUNCTION PSFF
!==================

!|| -------------------------------------------------------------------------
!||     FUNCTION DERPSF(x,xc,sigx,w)
!|| 
!||  Partial derivative of PSFD(x) respect to xc
!|| 
!============================
FUNCTION DERPSF(x,xc,sigx,w)
!============================
IMPLICIT NONE

      REAL      (kind=8) :: S15ADF
      REAL      (kind=8) :: DERPSF,SQRT2,SQRT2SI,ROOT2PI,ROOT2PIW,ROOTPI
      REAL      (kind=8) :: xc,x,sigx,w,u0,u1,u2,u3
      REAL      (kind=8) :: erfcu0, erfcu1, erfcu2, erfcu3
      INTEGER :: ifail
      
      SQRT2    = DSQRT(2.0d0)
      SQRT2SI  = SQRT2*sigx
      ROOT2PI  = 2.506628275D0
      ROOT2PIW = ROOT2PI*w
      ROOTPI   = 1.772453851D00
      
! Respect to the formula notations:  u0=u2  u1=u+ u2=u- u3=u1
       u0=(x-xc-(w+1)/2)/SQRT2SI
       u1=(x-xc-(w-1)/2)/SQRT2SI
       u2=(x-xc+(w-1)/2)/SQRT2SI
       u3=(x-xc+(w+1)/2)/SQRT2SI

      ifail=1
      erfcu0=S15ADF(u0,ifail)
     ! if(ifail.eq.1)PRINT *,' ifail =',ifail
      ifail=1
      erfcu1=S15ADF(u1,ifail)
     ! if(ifail.eq.1)PRINT *,' ifail =',ifail
      ifail=1
      erfcu2=S15ADF(u2,ifail)
     ! if(ifail.eq.1)PRINT *,' ifail =',ifail
      ifail=1
      erfcu3=S15ADF(u3,ifail)
     ! if(ifail.eq.1)PRINT *,' ifail =',ifail

      DERPSF=1./(2*w) *(erfcu0-erfcu1-erfcu2+erfcu3)
!==================
END FUNCTION DERPSF
!==================

!|| -------------------------------------------------------------------------
!||     FUNCTION DERSIG(x,xc,sigx,w)
!||
!||  Partial derivative of PSFD(x) respect to sigx
!|| 
!===========================
FUNCTION DERSIG(x,xc,sigx,w)
!===========================
      
IMPLICIT NONE

      REAL      (kind=8) :: DERSIG,SQRT2,SQRT2SI,ROOT2PI,ROOT2PIW,ROOTPI
      REAL      (kind=8) :: xc,x,sigx,w,u0,u1,u2,u3
      
      
      SQRT2    = DSQRT(2.0d0)
      SQRT2SI  = SQRT2*sigx
      ROOT2PI  = 2.506628275D0
      ROOT2PIW = ROOT2PI*w
      ROOTPI   = 1.772453851D00
      
!   Respect to the formula notations:  u0=u2  u1=u+ u2=u- u3=u1
       u0=(x-xc-(w+1)/2)/SQRT2SI
       u1=(x-xc-(w-1)/2)/SQRT2SI
       u2=(x-xc+(w-1)/2)/SQRT2SI
       u3=(x-xc+(w+1)/2)/SQRT2SI

      DERSIG = 1./(ROOT2PI*w) * &
        (DEXP(-u0**2)-DEXP(-u1**2)-DEXP(-u2**2)+DEXP(-u3**2))
!===================
END FUNCTION DERSIG
!===================

!|| -------------------------------------------------------------------------
!||     FUNCTION TPSF(x,y,ns1,ns2)
!||
!||  Total PSF value at x,y :
!||  sum of SPSF of sources from ns1 to ns2 + background level
!==========================
FUNCTION TPSF(x,y,ns1,ns2)
!==========================
USE FIT_PAR_MODULE
IMPLICIT NONE

      REAL      (kind=8) :: PSFD, PSFF
      REAL      (kind=8) :: TPSF, a,xc,yc,x,y,sig,dbwei,spsfx,spsfy
      INTEGER :: iadec,ns,n1,n2,ns1,ns2,ncomp,ich

      TPSF=0.D00
      iadec=IABS(idec)
      IF(NS1.NE.0)THEN
      DO ns=ns1,ns2
        ncomp=1+(ns-1)*3
        a=PSFP(ncomp+1,1)
        xc=PSFP(ncomp+2,1)
        yc=PSFP(ncomp+3,1)
        if(spin(ns).ne.0)then
          n1=i1
          n2=i2
        else
          n1=1
          n2=1
          sig=PSFP(nsouf*3+2,1)
          dbwei=1.D00
        endif
        DO ich=n1,n2
            if(spin(ns).ne.0)then
              sig=DBLE(SIGSI(ich))
              dbwei=DBLE(weig(ich,ns))
            endif
            if(iadec.eq.2)then
              spsfx = PSFD(x,xc,sig,dbw)
              spsfy = PSFD(y,yc,sig,dbw)
            elseif(iadec.eq.1) then
              spsfx = PSFF(x,xc,sig,dbw)
              spsfy = PSFF(y,yc,sig,dbw)
            endif
            TPSF = TPSF + dbwei*a*spsfx*spsfy
        ENDDO
      ENDDO
      ENDIF
      TPSF = TPSF + PSFP(1,1)
!==================
END FUNCTION TPSF
!==================



