



!******************
MODULE CHI2_FIT_PAR
!******************
  
IMPLICIT NONE

!!$      ! Basic parameters:
      INTEGER, parameter :: MaxNumSou=20, &
                                     MaxNumPar=3*MaxNumSou+3,&
                                     MaxNumNonLinPar=2*MaxNumSou+1
     
      INTEGER          :: iSpaceDim,jSpaceDim
      REAL    (kind=4),dimension(:,:),pointer, save  :: FittingCarte
      REAL    (kind=4),dimension(:,:),pointer, save :: errorCarte,lsvar
      REAL    (kind=4),dimension(:,:),pointer, save :: bkgCarte
      REAL    (kind=4),dimension(:,:),pointer, save :: diffCarte
      REAL    (kind=4),dimension(:,:),pointer, save :: psfCarte,effiCarte
      REAL    (kind=4),dimension(:,:),pointer, save :: resCarte
      REAL    (kind=4),dimension(:,:),pointer, save :: difCarte
      LOGICAL         ,dimension(:,:),pointer, save :: filterCarte
      
      REAL    (kind=8), save :: bl(MaxNumNonLinPar),bu(MaxNumNonLinPar)
      REAL    (kind=8), save :: XTX(MaxNumNonLinPar,MaxNumNonLinPar)
      REAL    (kind=8), save :: FITRESTAB(MaxNumPar,2)
      INTEGER, save :: nsouf,nparm,lparm,tot_mes
      INTEGER, save :: ifitdim,jfitdim,iCartemin,iCartemax,jCartemin,jCartemax
      REAL    (kind=4), save :: error
      REAL    (kind=4), save :: dof,spchisq,confl
      REAL    (kind=4), save :: spin(MaxNumSou),SIGSI(95),weig(95,MaxNumSou)
      INTEGER, save :: i1,i2,idec,np,iter,good_mesures
      INTEGER, save :: iDetTyp, PrintLevel=1,chi2_status
      
      ! Output param stat,
      REAL    (kind=4), save :: esd,stdev,xsta,ysta,time
!*********************
END MODULE CHI2_FIT_PAR
!*********************


!$$$$$$$$$$$$$$$$$$$$$$$$
! MODULES FOR LS FIT   
!$$$$$$$$$$$$$$$$$$$$$$$$

!**************************
MODULE LS_FIT_MODULE
!**************************
INTERFACE



 SUBROUTINE LS_FIT_INIT(DET,DETVAR,Effi,BCKGMODEL,DIFFMODEL,&
                          sou_number,Status)
!-------------------------------------------------
IMPLICIT NONE

REAL(kind=4), dimension(:,:), pointer  :: DET,DETVAR,Effi
REAL(kind=4), dimension(:,:), pointer,optional ::&
                                               BCKGMODEL,DIFFMODEL

INTEGER                     :: sou_number
INTEGER                  :: Status
END SUBROUTINE LS_FIT_INIT

SUBROUTINE LS_FIT_SHD(Det,Detn,LSBckg,LSBckgError,LSSources,&
                                LSSourceErrors,LSRes,Status)
!----------------------------------------------------------

IMPLICIT NONE
REAL(kind=4), dimension(:,:), pointer   :: Det,Detn
REAL(kind=4)       :: LSBckg,LSBckgError
REAL(kind=4), dimension(:), pointer   :: LSSources,LSSourceErrors
INTEGER :: LSRes,Status

END SUBROUTINE LS_FIT_SHD

END INTERFACE

!**************************
END MODULE LS_FIT_MODULE
!**************************


!**************************
 MODULE LS_FIT_INTERNAL
!*************************

INTERFACE

 SUBROUTINE LS_FIT(status)
!--------------------------
IMPLICIT NONE
INTEGER :: Status
END  SUBROUTINE LS_FIT


SUBROUTINE LS_VARIANCE(Status)
!----------------------------
IMPLICIT NONE
INTEGER :: Status
END SUBROUTINE LS_VARIANCE

END INTERFACE

!**************************
END MODULE LS_FIT_INTERNAL
!**************************

!###############################################
! CODE OF EXTERNAL USE PROCEDURES FOR Least squares FIT 
!               - INTERFACES IN LS_FIT_MODULE
!###############################################





!========================================================
 SUBROUTINE LS_FIT_INIT(DET,DETVAR,Effi,BCKGMODEL,DIFFMODEL,&
                          sou_number,Status)
!========================================================

      USE ISDC
      USE MIMOSA_CONTROL_MODULE
      USE MIMOSA_GLOBVAR_MODULE
      USE MIMOSA_USE_MODULE
      USE COMMON_PAR
      USE CHI2_FIT_PAR
      USE LS_FIT_INTERNAL
      IMPLICIT NONE

      ! Parameter
      REAL(kind=4), dimension(:,:), pointer  :: DET,DETVAR,Effi
      REAL(kind=4), dimension(:,:), pointer,optional ::&
                                               BCKGMODEL,DIFFMODEL
      INTEGER                       :: sou_number
      INTEGER :: Status
      ! dimension of source position arrays - in nsouf
!     LOCAL VARIABLES
      REAL(kind=4), dimension(:), pointer :: dectab
      INTEGER :: iexfov,iok
      INTEGER :: ns,ncomp,i,ieqspre
      REAL      (kind=4) :: xcdim,ycdim
      REAL      (kind=8) :: dbStaSpRe
      CHARACTER(len=20)::procName

     
      Status  = ISDC_OK

      method = 4 ! parameter for FUNCT1
      chi2_status = 0
      procName = 'LS_FIT_INIT'
      if(DebugMode.eq.3)&
         call Message(procName,' ',ZeroError,Status)

      nsouf = min(MaxNumSou,sou_number)
      if(nsouf < sou_number)then
         write(str250,'(" Too many sources :",I5," .Maximum ",I4," allowed - no simultaneous fit" )')sou_number,MaxNumSou
         call Message(procname,str250,Zeroerror,Status)
         Status=1
         return
      endif

      i = size(fitdectab,1)
      allocate(dectab(i),stat=iok)

      IF (iok /= 0) then
      call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
      endif

      dectab = FitdecTab
    
      iSpaceDim = size(DET,1)
      jSpaceDim = size(DET,2)
      tot_mes = iSpaceDim*jSpaceDim
      coeffnorm = iSpaceDim*jSpaceDim/duration(1)

      allocate(FittingCarte(iSpaceDim,jSpaceDim),&
               bkgCarte(iSpaceDim,jSpaceDim),LSVar(nsouf+1,nsouf+1),&
               effiCarte(iSpaceDim,jSpaceDim),&
               errorCarte(iSpaceDim,jSpaceDim),&
               diffCarte(iSpaceDim,jSpaceDim),&
               filterCarte(iSpaceDim,jSpaceDim),stat=iok)
      IF (iok /= 0) then
      call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
      endif

      FittingCarte = det

      bkgCarte=0.;LSVar = 0.
      diffCarte=0.
      effiCarte = effi
    
      filterCarte = .true.
      where(DETVAR.eq.0.)
           filterCarte = .false.
      endwhere
     
      filterCarte(65:66,:) = .false.
      filterCarte(:,33:34) = .false.
      filterCarte(:,67:68) = .false.
      filterCarte(:,101:102) = .false.

      errorCarte=1.
      good_mesures = sum(errorCarte,mask=filterCarte)
   
      if(ScwFitMode==1)then ! CHi2 fit
         where(filtercarte)
            errorcarte=sqrt(DETVAR)
         endwhere
      endif

      if(present(BCKGMODEL)) then
         bkgCarte = BCKGMODEL*effiCarte
      else
         bkgCarte = efficarte
      endif

      if(present(DIFFMODEL)) then

         diffCarte = DIFFMODEL*effiCarte
      else
         if(dectab(3*nsouf+3) .ge.0) then
            diffCarte  = effiCarte
      endif
      endif
     

     

      xcdim=REAL(iSpaceDim+1)/2.
      ycdim=REAL(jSpaceDim+1)/2.

!     SOME PARAMETERS and initializations
      ! Double Precision STArting SPatial REsolution 
      ! (~ 0 for perfect binning, ~ 0.4 for ISGRI, ~ 0.2 for Picsit)
      dbStaSpRe = 0.01D00
      IF (iDetTyp > 0) dbStaSpRe = 0.430D00 / DBLE(iDetTyp)
!!      dbw=2.435D00/DBLE(iDetTyp)
 !     IF (dbw == 0.D00) dbw=2.435D00
      ifitdim=iSpaceDim
      jfitdim=jSpaceDim
      xcdim = ifitdim/2.
      ycdim = jfitdim/2.
      
      iexfov=1
      idec=-1


   

!  Bkg Level (cts/pix)
         FITRESTAB(1,1) = 0.D0
         FITRESTAB(1,2) = dectab(1)

        
         do ns=1,nsouf
            ncomp=(ns-1)*3+1
! Enter photon index (> 0)
            spin(ns) = 0.
! Enter Starting Y T.A. Pos.
            FITRESTAB(ncomp+2,1) = FitPosX(ns)
           
            FITRESTAB(ncomp+2,2) = dectab(ncomp+2)
! Enter Starting Z T.A. Pos.
            FITRESTAB(ncomp+3,1) = FitPosY(ns)
          
            FITRESTAB(ncomp+3,2) = dectab(ncomp+3)
! Enter Starting Source Intensity
            FITRESTAB(ncomp+1,1) = 0.D0
            FITRESTAB(ncomp+1,2) = dectab(ncomp+1)
            if(spin(ns).eq.0)ieqspre=1
         enddo

        
! Starting Resolution (sigma in pix)  ==> ',$)
           FITRESTAB((nsouf*3+2),1) = dbStaSpRe
           FITRESTAB((nsouf*3+2),2) = dectab(nsouf*3+2)
       

!diffuse emission
           
           FITRESTAB((nsouf*3+3),1) = 1.
           FITRESTAB((nsouf*3+3),2) = dectab(nsouf*3+3)
         

! Confidence level for error estimation
         confl= 0.683
         np   = 1
         deallocate(dectab)
!==========================
END SUBROUTINE LS_FIT_INIT
!==========================

!==============================================================
SUBROUTINE LS_FIT_SHD(DET,DETN,LSBckg,LSBckgError,LSSources,&
                                LSSourceErrors,LSRes,Status)
!==============================================================

USE ISDC
USE COMMON_PAR
USE CHI2_FIT_PAR
USE LS_FIT_INTERNAL
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE

IMPLICIT NONE
REAL(kind=4), dimension(:,:), pointer   :: Det,Detn
REAL(kind=4)       :: LSBckg,LSBckgError
REAL(kind=4), dimension(:), pointer   :: LSSources,LSSourceErrors
INTEGER :: LSRes,Status

! Local variables

INTEGER :: i,j
CHARACTER(len=20)   :: procname

procname = 'LS_FIT_SHD'
Status = ISDC_OK
LSRes = 0
    


!  FILL THE DOUBLE PRECISION MINIMAP ARRAY
         iCartemin = 1
         iCartemax = ifitdim
         jCartemin = 1
         jCartemax = jfitdim

!  MINIMIZATION
CALL LS_FIT(Status)
if(Status.ne.0)then
  call MESSAGE(procName,' Error in LS fit ',FitError,Status)
  return
endif 
 
!variance
call LS_VARIANCE(Status)

!  RESULT OF FIT
!call FIT_RESULTS(status)

LSBckg = FITRESTAB(1,1)
LSBckgError = sqrt(LSVar(1,1))
do i=1,nsouf
   j = (i-1)*3+2
   LSSources(i) = FITRESTAB(j,1)
   LSSourceErrors(i) =sqrt(LSVar(i+1,i+1))
enddo


      
     
Detn = Det - LSBckg
do i=1,nsouf
   Detn = Detn -LSSources(i)*SourceModel(i,:,:)
enddo

if( associated(bkgCarte))deallocate(bkgCarte)
if( associated(diffCarte))deallocate(diffCarte)
if( associated(LSVar))deallocate(LSVar)
if( associated(errorCarte))deallocate(errorCarte)
if( associated(FittingCarte))deallocate(FittingCarte)
if( associated(filterCarte))deallocate(filterCarte)
if( associated(effiCarte))deallocate(effiCarte)

!==========================
END SUBROUTINE LS_FIT_SHD
!==========================

!#######################################################################
! CODE OF INTERNAL USE SUBROUTINES FOR LS FIT                         
! INTERFACES IN LS_FIT_INTERNAL                                       
!#######################################################################



!=======================
 SUBROUTINE LS_FIT(status)
!=======================
USE ISDC
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
USE COMMON_PAR
USE CHI2_FIT_PAR
USE FITS_DECLARATIONS
IMPLICIT NONE
INTEGER :: Status

! Local variables
INTEGER,parameter :: ia = MaxNumNonLinPar 
REAL(kind=8) :: xx(MaxNumNonLinPar)
REAL(kind=8) :: dbCarte2(iSpaceDim,jSpaceDim)
REAL(kind=8) :: PHF(MaxNumSou+2,iSpaceDim,jSpaceDim),BF(ia)
REAL(kind=8) :: W1(ia),W2(ia)
INTEGER :: ifail,ieqspre,idiffuse
INTEGER :: nc,inc,nnc,ic,id,idf
REAL(kind=8) :: dxc,dyc,da
INTEGER :: ns,ncomp
character(len=20)::procname
    
procname = 'LS_FIT'
Status = 0


dbCarte2(:,:)=DBLE(FittingCarte(:,:))


!  Number of free linear  parameters and set xx
lparm=0
if(FITRESTAB(1,2).eq.1)lparm=lparm+1
ieqspre=0
do ns=1,nsouf
   ncomp=1+(ns-1)*3
   if(spin(ns).eq.0)ieqspre=1
!  Intensity
   if(FITRESTAB(ncomp+1,2).eq.1)lparm=lparm+1

enddo
ncomp=3*nsouf+3
IF (FITRESTAB(ncomp,2).EQ.1) THEN
         lparm=lparm+1
         idiffuse=1
ENDIF


 
! COEFFICIENTS POUR LE CALCUL DES MOINDRES CARRES
nc=lparm
PHF(:,:,:) = 0.D00
      
      
inc=0
DO nnc=1,(nsouf+1)
   ns=nnc-1
   ncomp=IABS(2+(ns-1)*3)
   if(ns.ge.1)then
       dxc=FITRESTAB(ncomp+1,1)
       dyc=FITRESTAB(ncomp+2,1)
       IF(FITRESTAB(ncomp,2).EQ.1)inc=inc+1
       IF(FITRESTAB(NCOMP,2).EQ.1)THEN
          PHF(inc,:,:) = PHF(inc,:,:) + SourceModel(ns,:,:)*effiCarte(:,:)
       ELSEIF(FITRESTAB(NCOMP,2).EQ.0)THEN
           da=FITRESTAB(ncomp,1)
           dbCarte2(:,:) = dbCarte2(:,:) - da*SourceModel(ns,:,:)*effiCarte(:,:)
       ENDIF
         
    else
       IF(FITRESTAB(1,2).EQ.1)THEN
          inc=inc+1
          PHF(1,:,:) = bkgCarte(:,:)
         ! call FITS_FILE(1,bkgCarte,1,'bkg_ls.fits')
          
       ELSEIF(FITRESTAB(1,2).EQ.0)THEN
          da=FITRESTAB(1,1)
          dbCarte2(:,:)=dbCarte2(:,:)-da
       ENDIF
     endif
ENDDO

! Diffuse emission

IF(FITRESTAB(3*nsouf+3,2).EQ.1)THEN
             inc=inc+1
             PHF(inc,:,:) = diffCarte(:,:)
             
ELSEif(FITRESTAB(3*nsouf+3,2).EQ.0)THEN
             da=FITRESTAB(1,1)
             dbCarte2(:,:)=dbCarte2(:,:)-da
            
ENDIF

!  Initializes ARRAYS AF and BF of EQUAT to 0
DO IC = 1,NC
         BF(IC) = 0.0
         DO ID = 1,NC
            XTX(IC,ID) = 0.0
         ENDDO
ENDDO   

DO IC = 1,NC

   where(filtercarte)
      PHF(IC,:,:) = PHF(IC,:,:)/errorcarte(:,:)
   endwhere

   XTX(IC,IC) = sum(PHF(IC,:,:)**2,mask=filterCarte)
   BF(IC) = sum(PHF(IC,:,:)*dbCarte2(:,:)/errorcarte(:,:),mask=filterCarte) 
ENDDO 
  
NNC = NC -1
DO IC = 1,NNC
   IDF = IC + 1
   DO ID = IDF,NC
       XTX(IC,ID) = sum(PHF(IC,:,:) * PHF(ID,:,:),&
                        mask=filterCarte)
       XTX(ID,IC) = XTX(IC,ID)
   ENDDO
ENDDO  

! CALCUL OF CONTRIBUTIONS
   ifail = -1
   call X04aaf(1,-1)
   call X04ABF(1,-1)
! CALCULATES AF*XX=BF
    CALL F04ASF(XTX,ia,BF,nc,XX,W1,W2,ifail)
 if(ifail.ne.0)then
    call WAR_MESSAGE(procname,' Error in F04ASF',0,Status)
    select case(ifail)
    case(1)
       write(str250,*)' LS matrix is singular '
    case(2)
       write(str250,*)' LS matrix too ill-conditioned '
    case(3)
       write(str250,*)' uncorrect input for matrix inversion '
    end select
   call WAR_Message(procName,str250,0,Status)
   
    chi2_Status=FitError
    return
endif

 
! PARAMETERS AT EACH STEP

ic=0
if(FITRESTAB(1,2).eq.1)then
 ic=1
 FITRESTAB(1,1)=xx(ic)
endif
do ns=1,nsouf
         ncomp=2+(ns-1)*3
         if(FITRESTAB(ncomp,2).eq.1)THEN

              ic=ic+1
              FITRESTAB(ncomp,1)=xx(ic)
         endif
enddo


!==========================
END SUBROUTINE LS_FIT
!==========================

!=============================
SUBROUTINE LS_VARIANCE(Status)
!=============================
USE ISDC
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
USE COMMON_PAR
USE CHI2_FIT_PAR
IMPLICIT NONE
INTEGER :: Status
!LOCAL VARIABLES
INTEGER :: n,k,ifail,nnc,ns,ncomp
INTEGER,PARAMETER :: nmax= MaxNumNonLinPar 
REAL(kind=8)              :: A(nmax,nmax),BB(NMAX),C(NMAX),WK1(NMAX),WK2(NMAX),res(iSpaceDim,jSpaceDim),s2
!REAL(kind=8)              ::invert(nmax)
CHARACTER(len=80) :: str180
CHARACTER(len=20) :: procname

Status = 0
procname = 'LS_VARIANCE'


n = nsouf+1

!INVERSION OF RTR

do k=1,n
   a =XTX
   ifail = 1
   bb = 0.0d0
   bb(k) = 1.0d0
   IFAIL = -1
   call X04aaf(1,-1)
   call X04ABF(1,-1)
   call F04ASF(a,nmax,BB,n,C,WK1,WK2,IFAIL)
   if(ifail.ne.0) then
      select case(ifail)
         case(1)
              write(str180,*)' RTR matrix is singular '
         case(2)
              write(str180,*)' RTR matrix too ill-conditioned '
         case(3)
              write(str180,*)' uncorrect input for matrix inversion '
      end select
      str180=str180// ' Variance matrix set to MaxVar'
      call WAR_Message(procName,str180,0,Status)
     LSVar=MaxMaxLvar
      status = 1
      return
  else
    LSVar(1:n,k) = c(1:n)

   
  endif
enddo

!ESTIMATION OF S^2
res(:,:)=DBLE(FittingCarte(:,:))

DO nnc=1,(nsouf+1)
   ns=nnc-1
   ncomp=IABS(2+(ns-1)*3)
   if(ns.ge.1)then
     res(:,:) = res(:,:) - FITRESTAB(ncomp,1)*SourceModel(ns,:,:)*effiCarte(:,:)
   else
     res(:,:) = res(:,:) - FITRESTAB(1,1)
     endif
ENDDO

! Diffuse emission excluded for the moment

where(filtercarte)
   res = res/errorcarte
endwhere


s2 = sum(res*res,mask=filterCarte)/(good_mesures-n)

LSVar = LSVar*s2
!=============================
END SUBROUTINE LS_VARIANCE
!=============================
