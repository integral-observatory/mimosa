
!***********************************
MODULE IBIS_DECON_AUX
!***********************************
INTERFACE
SUBROUTINE DoMotif(motif,is,js,rep,replic)
!------------------------------
    IMPLICIT NONE
    REAL(kind=4),dimension(:,:),pointer :: motif
    INTEGER :: is,js,rep,replic
 END SUBROUTINE DoMotif
 END INTERFACE
!***********************************
END MODULE IBIS_DECON_AUX
!***********************************


!***********************************
MODULE IBIS_DECON_LIB_MOD
!***********************************
  INTERFACE
   
  
   !IBIS_IMAGING_INIT
    SUBROUTINE IBIS_IMAGING_INIT(iDet,iprint)
      IMPLICIT NONE

      ! I/F Parameters
      INTEGER :: iDet
      INTEGER,optional :: iprint
    END  SUBROUTINE IBIS_IMAGING_INIT


    ! IBIS_IMASTA
 
      SUBROUTINE IBIS_IMASTA(IMAGE,WEI,ImaStaPar,Masque,iprint,Status)
      IMPLICIT NONE
      ! Input
      REAL     (kind=4), dimension(:,:),   pointer :: IMAGE,WEI
      LOGICAL  , dimension(:,:), optional, pointer :: Masque
      INTEGER  , optional ::                          iprint
      ! Output
      REAL     (kind=4), dimension(:),     pointer :: ImaStaPar
      INTEGER :: Status
      END SUBROUTINE IBIS_IMASTA
 

  ! IBIS_IMASN

    SUBROUTINE IBIS_IMASN(IMAGE,VAR,WEI,mask,reg,soupar,imastat,sousnr,&
                          iprint,Status)
      IMPLICIT NONE
      ! INTERFACE variables
      REAL     (kind=4), dimension(:,:), pointer :: IMAGE
      REAL     (kind=4), dimension(:,:), optional, pointer :: VAR,WEI,reg
      LOGICAL  , dimension(:,:), optional,  pointer :: mask
      REAL     (kind=4), dimension(:), optional, pointer :: soupar, imastat
      REAL     (kind=4), dimension(:,:), optional, pointer :: sousnr
      INTEGER, optional :: iprint
      INTEGER :: Status
    END SUBROUTINE IBIS_IMASN



 
END INTERFACE

!*****************************
END MODULE  IBIS_DECON_LIB_MOD
!*****************************


!**************************
 MODULE FLUX_APP_MOD
!**************************

INTERFACE
FUNCTION UnitSource(x,y) 
!-----------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL(kind=4)     :: x,y,UnitSource
END FUNCTION UnitSource
END INTERFACE

!*****************************
END MODULE  FLUX_APP_MOD
!**************************


!**************************
 MODULE DeconIma_AUX
!**************************
INTERFACE
SUBROUTINE Deconvolution(WEI,DetEffi,DET,&
                        RWORK,RSKY,CWORK,CWORK2,SKY,&
                        ishif,ifCorr,idim,jdim,status)
!-----------------------------------------------------

IMPLICIT NONE

! GLOBAL variables
INTEGER                                  :: ifCorr,idim,jdim,status
INTEGER, dimension(2):: ishif
REAL(kind=4)   , dimension(:,:), pointer :: WEI,DET,DetEffi
REAL(kind=4)   , dimension(:,:), pointer :: SKY
COMPLEX   (kind=4), dimension(:,:) ,pointer      :: CWORK,CWORK2
REAL      (kind=4), dimension(:,:), pointer      :: RWORK,RSKY
END SUBROUTINE Deconvolution

END INTERFACE
!**************************
END MODULE DeconIma_AUX
!**************************

!**************************
 MODULE DeconIma_MODULE
!**************************

INTERFACE

SUBROUTINE DeconIma(WEI,DET,DetEffi,DETVAR,SKY,VAR,&
     ivar,ifCorr,ifExpo,ifNorm,Status)
!--------------------------------------------------------
IMPLICIT NONE

! GLOBAL variables
REAL(kind=4)   , dimension(:,:), pointer :: WEI,DET,DetEffi,DETVAR
REAL(kind=4)   , dimension(:,:), pointer :: SKY,VAR
INTEGER                          :: ivar,ifCorr,ifExpo,ifNorm
INTEGER                                  ::status
END SUBROUTINE DeconIma

END INTERFACE
!**************************
END  MODULE DeconIma_MODULE
!**************************






!===========================
SUBROUTINE Peaks(SigmaSky,filter)
!===========================
USE ISDC
USE FITS_DECLARATIONS
IMPLICIT NONE
REAL(kind=4) ,dimension(:,:),pointer :: SigmaSky
LOGICAL  ,dimension(:,:),pointer :: Filter


!LOCAL VARIABLES
INTEGER :: ifan,jfan,idim,jdim,i,j,i1,i2,j1,j2,i0,j0,isou,jsou


REAL(kind=4),dimension(3,3) :: modelval
INTEGER,dimension(3) ::modeli,modelj



idim = size(SigmaSky,1)
jdim = size(SigmaSky,2)
i0 = (idim+1)/2
j0 = (jdim+1)/2

do isou=1,idim
do jsou=1,jdim
if(SigmaSky(isou,jsou) > 0)then
   modeli = 0;modelj = 0;modelval = 0.
   ifan = isou

   do while((ifan.gt.0).and.(Filter(ifan,j0)))
     ifan = ifan-129
   end do
   ifan = ifan+129

   modeli(1) = ifan
   i=1
   do while((ifan.le.idim).and.(Filter(ifan,j0)))
     ifan = ifan+129
     i=i+1
     modeli(i) = ifan
   enddo

   jfan = jsou
   do while((jfan.gt.0).and.(Filter(i0,jfan)))
    jfan = jfan-129
   enddo
   jfan = jfan+129

   modelj(1) = jfan
   j=1
   do while((jfan.le.jdim).and.(Filter(i0,jfan)))
     jfan = jfan+129
    j=j+1
     modelj(j) = jfan
   enddo

   do i=1,3
   do j=1,3
    if((modeli(i).gt.0).and.(modelj(j).gt.0).and.(filter(modeli(i),modelj(j))))then
          i1 =max( modeli(i)-1,1)
          i2 = min(modeli(i)+1,idim)
          j1 =max( modelj(j)-1,1)
          j2 = min(modelj(j)+1,jdim)
          modelval(i,j)=maxval(sigmasky(i1:i2,j1:j2),filter(i1:i2,j1:j2))
     endif
  
   enddo
   enddo


   do i=1,3
   do j=1,3
      i1 = modeli(i)
      j1 = modelj(j)
      if(i1*j1.gt.0)then
         if((sigmasky(i1,j1).gt.0).and.(filter(i1,j1)))then
         if((abs(isou-i1).gt.2).and.(abs(jsou-j1).gt.2)) then
!!$          print *,' confusion of positions : ',isou,jsou,sigmasky(isou,jsou)
!!$           print *,'                         ',i1,j1,sigmasky(i1,j1)
         endif
        endif
     endif
  enddo
  enddo

!iloc = maxloc(modelval)
!imax = iloc(1)
!jmax = iloc(2)

endif
enddo
enddo



!====================
END SUBROUTINE Peaks
!====================

!===========================
SUBROUTINE FantomLines(Sky)
!===========================
USE ISDC
USE FITS_DECLARATIONS
IMPLICIT NONE
REAL(kind=4),dimension(:,:),pointer :: Sky

INTEGER :: idim,jdim,i,j
REAL(kind=4) :: mean,dev,err
REAL(kind=4),dimension(:),pointer :: hortab,vertab
idim = size(Sky,1)
jdim = size(Sky,2)
mean = sum(Sky)/idim/jdim
allocate(hortab(jdim),vertab(jdim))
!horizonthal lines
do j=1,jdim
  dev = sum((Sky(:,j)-mean)**2)/idim
  err = sqrt(dev)
  hortab(j)= err
enddo
do i=1,idim
  dev = sum((Sky(i,:)-mean)**2)/jdim
  err = sqrt(dev)
  vertab(i)= err
enddo
open(53,file='hor_err.fits',status='unknown')
write(53,*)hortab
close(53)
open(54,file='ver_err.fits',status='unknown')
write(54,*)vertab
close(54)
deallocate(hortab)
!===========================
END SUBROUTINE FantomLines
!===========================



!========================
FUNCTION UnitSource(x,y) 
!========================
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL(kind=4)     :: x,y,UnitSource
!LOCAL VARIABLES
REAL(kind=4),dimension(4),parameter :: fluxpeak=(/440.98,442.55,442.08,441.55/)
INTEGER  :: ip,jp
REAL(kind=4)     :: dx,dy,dz,fx,fy,fz

ip = nint(x)
jp = nint(y)
dx = abs(ip-x)
dy = abs(jp-y)
dz = sqrt(dble(dx**2+dy**2))
fx = fluxpeak(1)+(fluxpeak(2)-fluxpeak(1))*dx/0.5
fy = fluxpeak(1)+(fluxpeak(3)-fluxpeak(1))*dy/0.5
fz = fluxpeak(1)+(fluxpeak(4)-fluxpeak(1))*dz/0.5
UnitSource = (fx+fy+fz)/3. ! peak approximation of unit source
!========================
END FUNCTION UnitSource
!========================





!|| --------------------------------------------------------------------------
!|| SUB IBIS_IMAGING_INIT                                                  1.2 
!||
!||  SUB IBIS_IMAGING_INIT(iDet)
!||
!||  Defines basic parameters of a coded mask imaging system included in the
!||  mod IBIS_IMAGING_PARAM for the IBIS telescope
!||
!||  Input :
!||  - iDet    I4     Specifies detector layer: 0=external 1=ISGRI,2=PICSIT
!||                   if < 0 => perfect system with integer number of pix per
!||                   per mask element = IABS(iDet)
!||
!||  Output : set the global variables of the IBIS_IMAGING_PARAM as defined 
!||           below
!||  
!||  Error Codes:
!||  - not yet defined
!||
!||  History:
!||  1.2 31/08/00 AG : size detector variable
!||  1.2 07/07/00 AG : idet=0 => externally defined
!||  1.1 24/02/00 AG : redefinition of variable names (...Dim in pixel units)
!|| 

SUBROUTINE IBIS_IMAGING_INIT(iDet,iprint)
      USE MIMOSA_GLOBVAR_MODULE
      USE IBIS_IMAGING_PARAM

      IMPLICIT NONE

      ! I/F Parameters
      INTEGER :: iDet
      INTEGER,optional :: iprint
      ! Local parameters
      ! IBIS Geometry (mm, degrees)


      ! Picsit Parameters
!!$      INTEGER, parameter :: PicsNumPix = 4096
!!$      INTEGER, parameter :: PicsYDim=64, PicsZDim=64
      REAL    (kind=4), parameter :: PicsPixSiz=9.2, PicsPixEffSiz=9.2
      REAL    (kind=4), parameter :: PicsPixThi=30., PicsPixDZSiz=0.
      REAL    (kind=4), parameter :: PicsRotAng=0. 
      REAL    (kind=4), parameter :: PicsMasDis=3110. ! Det-top to Mask-bottom
      
      ! LOCAL VARIABLES
      ! FOV
      INTEGER::iFovDim,jFovDim,iFCFovDim,jFCFovDim,i50FovDim,j50FovDim
      INTEGER :: ireb       ! rebinning for optimum system

     
      ! Defines IBIS Mask dimensions
      IF (idet >= 0)THEN      
         iMaBPDim = 53 ; jMaBPDim = 53
         iMasDim = 95 ; jMasDim = 95
      ENDIF

      ! Detector dependent parameters
      IF (idet == 0)THEN
      ! Externally defined values for IBIS system
         Detector='IBIS/External'

      ELSEIF (idet < 0) THEN
      ! Ideal detector (integer number of pixels per mask element)

         Detector='Perfect System'

         ! Detector values in pixel units
         ireb=IABS(idet)

         iDetDim = iMaBPDim*ireb
         jDetDim = jMaBPDim*ireb
!!$         PixDZDim  = 0.
!!$         PixEffDim = 1.
!!$         PixThiDim = 0.

         ! Mask Values in pixel units
         MaElPixDim = REAL(ireb)

         MaElThiDim = MaElThi / IsgrPixSiz
         MasDetRotAng = 0.
         MasDetDisDim = IsgrHalfMasDis/(MaElSiz/REAL(ireb))

      ELSEIF (idet == 1) THEN
      ! ISGRI Detector Layer

         Detector='IBIS/ISGRI'

         ! Detector values in pixel units
         iDetDim = 130
         jDetDim = 134
         iActiveDetDim = 128
         jActiveDetDim = 128
!!$         PixDZDim  = IsgrPixDZSiz / IsgrPixSiz
!!$         PixDzHalfDim  = PixDZDim/2.
!!$         PixEffDim = IsgrPixEffSiz / IsgrPixSiz
!!$         PixThiDim = IsgrPixThi / IsgrPixSiz

         ! Mask Values in pixel units
         MaElPixDim = MaElSiz / IsgrPixSiz
     !    MaElPixDim = 129./53.

         MaElThiDim = MaElThi / IsgrPixSiz
         MasDetRotAng = MasRotAng + IsgrRotAng
         MasDetDisDim = IsgrHalfMasDis / IsgrPixSiz

      ELSEIF (idet == 2) THEN
      ! PICSIT Detector Layer

         Detector='IBIS/PICSIT'

         ! Detector values in pixel units
         iDetDim = 65
         jDetDim = 67
!!$         PixDZDim  = PicsPixDZSiz / PicsPixSiz
!!$         PixEffDim = PicsPixEffSiz / PicsPixSiz
!!$         PixThiDim = PicsPixThi / PicsPixSiz

         ! Mask Values in pixel units
         MaElPixDim = MaElSiz / PicsPixSiz
         MaElPixDim = 64.5/53.

         MaElThiDim = MaElThi / PicsPixSiz
         MasDetRotAng = MasRotAng + PicsRotAng
         MasDetDisDim = PicsMasDis / PicsPixSiz

      ENDIF
      
      ! Derived Imaging System parameters
      PixAngDim = atan(1./MasDetDisDim)*180./PI

      MasPixDimI=REAL(iMasDim)*MaElPixDim
      MasPixDimJ=REAL(jMasDim)*MaElPixDim
      iMasPixDim=INT(MasPixDimI)+1
      jMasPixDim=INT(MasPixDimJ)+1

    

     
END SUBROUTINE IBIS_IMAGING_INIT


!|| --------------------------------------------------------------------------
!|| SUB IBIS_IMASTA                                                        1.1
!||
!||   SUB IMASTAT(IMAGE,WEI,ImaStaPar,Masque,iprint,Status)
!||
!||   Description:
!||     Simple and weighted statistics of an image sector of the input 
!||     image array IMAGE.
!||  
!||   Parameters
!||     IMAGE          REAL    POINT     IN    input image 
!||     WEI          REAL    POINT     IN    input weighting array (0.-1.)
!||                                          where wei=0 => filter=f
!||     Masque       LOGICAL POINT OPT IN    input filter for stat calculation
!||                                    OUT   filter used (IN & WEI>0)
!||     iprint       INT           OPT IN    print stat (=0 no, 1=yes)
!||
!||     ImaStaPar    REAL    POINT     OUT   Vector of statistical parameters
!||                                          including sequencially:
!||
!||   -  num         =  number of pixels
!||   -  tot         =  total number of counts
!||   -  rmean       =  mean value
!||   -  stdev       =  standard deviation
!||   -  imn,jmn     =  position of min value
!||   -  rmin        =  minimum value
!||   -  imx,jmx     =  max position
!||   -  rmax        =  max value
!||
!||  History:
!||  1.1 09/03/00 AG : rectangular region, C, output param
!||  1.0 24/02/00 AG : from ibis_osm_isgr, C, iprint
!||    
!||
!|| --------------------------------------------------------------------------

      SUBROUTINE IBIS_IMASTA(IMAGE,WEI,ImaStaPar,Masque,iprint,Status)

      IMPLICIT NONE

      ! INTERFACE 
      ! Input
      REAL     (kind=4), dimension(:,:),   pointer :: IMAGE,WEI
      LOGICAL  , dimension(:,:), optional, pointer :: Masque
      INTEGER  , optional :: iprint
     
      ! Output
      REAL     (kind=4), dimension(:),     pointer :: ImaStaPar
      INTEGER :: Status
      ! Local
      LOGICAL  , dimension(:,:), pointer :: Filter
      INTEGER                  :: npix,nregpix,iloc(2),imasque
      REAL     (kind=4)                  :: tot,tot2,rmean,var,stdev,sterr
      INTEGER                  :: i,j,idim,jdim,ii,if,ji,jf
      INTEGER                  :: imn,jmn,imx,jmx,iw,iok
      REAL     (kind=4)                  :: rmin,rmax,efficiency

      Status = 0

      imasque=0
      idim=SIZE(IMAGE,1)
      jdim=SIZE(IMAGE,2)
       
       ALLOCATE(Filter(1:idim,1:jdim),stat=iok)
       IF (iok /= 0)then
         if(present(iprint))print *,'Allocation problem for filter'
         status = 1
         return
       endif
      filter=.true.
      ! Allocation and Default values for filter
      IF (PRESENT(masque)) THEN
         IF(.not.ASSOCIATED(masque)) THEN
           ALLOCATE(masque(1:idim,1:jdim),stat=iok)
           IF (iok /= 0)then
              if(present(iprint))print *,'Allocation problem for filter'
              status = 1
             return
           endif
          Masque=Filter
          imasque=1
         ELSE
           Filter=Masque
         ENDIF
      ENDIF
      ! Value of iw from iprint and default
      IF (PRESENT(iprint)) THEN
         iw=iprint
      ELSE
         iw=0
      ENDIF

      ! Number of pixels in filter of analysis
      nregpix=COUNT(FILTER)

      ! Use only non-noisy pixel with weight larger than 0
      WHERE(WEI <= 0.) Filter = .false.

      ! Compute statistics for image without weigth
      npix=0 ; tot=0. ; tot2=0. ; var=0. ; stdev=0.
      ! Region of analysis
      npix=COUNT(filter)
      ii=idim ; if=1 ; ji=jdim ; jf=1
      do i=1,idim
         do j=1,jdim
            IF(filter(i,j)) THEN
              ii=MIN(ii,i)
              if=MAX(if,i)
              ji=MIN(ji,j)
              jf=MAX(jf,j)
            ENDIF
         enddo
      enddo

      tot=SUM(IMAGE,filter)
      tot2=SUM((IMAGE*IMAGE),filter)
      var = 0.
      stdev = 0.

      if(npix > 0) then
         IF (npix > 1) var=(tot2-(tot*tot)/REAL(npix))/REAL(npix-1)
         IF (var >= 0.) stdev=SQRT(var)
         rmean=tot/REAL(npix)
         sterr=stdev/SQRT(REAL(npix))
      else
         Status = 2
         return
      endif

      rmax=MAXVAL(IMAGE,filter)
      iloc=MAXLOC(IMAGE,filter)
      imx=iloc(1) ; jmx=iloc(2)
      rmin=MINVAL(IMAGE,filter)
      iloc=MINLOC(IMAGE,filter)
      imn=iloc(1) ; jmn=iloc(2)

      ImaStaPar(1) = REAL(npix)
      ImaStaPar(2) = rmean
      ImaStaPar(3) = stdev
      ImaStaPar(4) = rmax
      ImaStaPar(5) = rmin

      Efficiency = SUM(WEI,filter)/REAL(npix)

      IF(iw == 2)THEN
        PRINT *,''
        PRINT *,' Image Statistics:'
        PRINT *,''
        PRINT *,' Rect. Region    ',ii,if,ji,jf
        PRINT *,' Region Dim.     ',(if-ii+1),(jf-ji+1),nregpix
        PRINT *,' Total pixels %R ',npix,(REAL(npix)/REAL(nregpix)*100.)
        PRINT *,' Total cts, Eff. ',tot,efficiency
        PRINT *,' Mean intensity  ',rmean
        PRINT *,' Standard dev/err',stdev,sterr
        PRINT *,' Max pos / val   ',imx,jmx,' ',rmax
        PRINT *,' Min pos / val   ',imn,jmn,' ',rmin
        PRINT *,''
      ENDIF

      ! Compute statistiques for image with weigth if WEI /= 1.
      IF(Efficiency /= 1. ) THEN
        tot=0. ; tot2=0. ; var=0. ; stdev=0.
        tot=SUM(IMAGE*WEI,filter)
        tot2=SUM((IMAGE*IMAGE*WEI*WEI),filter)
        if(npix > 0)then
           IF (npix > 1) var=(tot2-(tot*tot)/REAL(npix))/REAL(npix-1)
           IF (var >= 0.) stdev=SQRT(var)
           rmean=tot/REAL(npix)
           sterr=stdev/SQRT(REAL(npix))
        else
           status=2
           return
        endif
        rmax=MAXVAL(IMAGE*WEI,filter)
        iloc=MAXLOC(IMAGE*WEI,filter)
        imx=iloc(1) ; jmx=iloc(2)
        rmin=MINVAL(IMAGE*WEI,filter)
        iloc=MINLOC(IMAGE*WEI,filter)
        imn=iloc(1) ; jmn=iloc(2)

        ImaStaPar(6) = REAL(npix)
        ImaStaPar(7) = rmean
        ImaStaPar(8) = stdev
        ImaStaPar(9) = rmax
        ImaStaPar(10) = rmin

      ENDIF

      IF(imasque == 1) masque=filter
      deallocate(filter)
      END SUBROUTINE IBIS_IMASTA

!|| --------------------------------------------------------------------------
!|| SUB IBIS_IMASN                                                         2.1
!||
!||  SUB IBIS_IMASN(IMA,VAR,WEI,mask,reg,soupar,imastat,sousnr,iprint,Status)
!||
!||  IMAge statistics and source Signal to Noise (SN) estimates
!||  in a rectangular image sector, for given number of sources and where
!||  pixels whose positions correspond to .false. values of logical array
!||  mask are not considered in the analysis.
!|| 
!||  Input :
!||  - IMAGE       R4  Array Poi      imagege array to be analysed 
!||  - VAR       R4  Array Poi Opt  image of variance (dim=IMA def=0.)
!||  - WEI       R4  Array Poi Opt  image of weights for image (dim=IMA def=1)
!||  - mask      L   Array Poi Opt  logical filter for IMAGE (dim=IMA def=T)
!||  - reg       R4  Array Poi Opt  limits in i and j  of the image 
!||                                 sector to analyse (def=IMA lim,)
!||  - soupar    R4  10 Vec    Opt  parameters for the stat. analysis
!||     soupar(1) =  number of source to search (def=0)
!||     soupar(2) =  source radius in pix. (def=0.)
!|| 
!||  - sousnr    R4  Array Poi Opt  if sousnr(1:3,) /= 0. contains in input 
!||                                 positions of sources to search for,
!||                                 other parameters ignored
!||                                 (meaning indiv. values: see below)
!||  - iprint    I4            Opt  print stat (0/1/2=no/yshort/ycomp, def=0)
!||
!||    
!||  Output :
!||  - imastat   R4  10 Vec    Opt  Parameters for image general statistics:
!||     imastat(1)  = num         :  number of pixels
!||     imastat(2)  = rtot        :  total number of counts
!||     imastat(3)  = spmean      :  mean value
!||     imastat(4)  = spsdev      :  standard deviation
!||     imastat(5&6)= imn,jmn     :  min position
!||     imastat(7)  = rmin        :  minimum value
!||     imastat(8&9)= imx,jmx     :  max position
!||     imastat(10) = rmax        :  max value
!||  - sousnr    R4  array Poi Opt  Result on sources (10 x nsou) 
!||                                 If snr =/ 0. for ns1 values in input -> 
!||        
!||     Y Z         = source position respect to IMA center in pix
!||     flux        = source image value 
!||     error       = source value error (=SQRT(Var(Y,Z)) or StDev if no VAR)
!||     
!||                 = source total counts 
!||
!||  Parameters :
!||  -  
!||
!||  Test parameters :
!||

!||  Sub-History :
!||    2.1  11/02/00 - AG : correct spvar,iprint
!||    ibis_decon2      06/09/99-AG : C, imastat
!||    ibis_decon1      02/09/99-AG : C, correct mask/reg error,
!||    ibis_decon_lib0  13/08/99 AG : 

     SUBROUTINE IBIS_IMASN(IMAGE,VAR,WEI,mask,reg,soupar,imastat,sousnr,&
                           iprint,Status)

      IMPLICIT NONE

      ! INTERFACE variables
      REAL     (kind=4), dimension(:,:), pointer :: IMAGE
      REAL     (kind=4), dimension(:,:), optional, pointer :: VAR,WEI,reg
      LOGICAL  , dimension(:,:), optional, pointer :: mask
      REAL     (kind=4), dimension(:), optional, pointer :: soupar, imastat
      REAL     (kind=4), dimension(:,:), optional, pointer :: sousnr
      INTEGER  , optional :: iprint
      INTEGER :: Status
      ! Local variables
      LOGICAL  , dimension(:,:), pointer :: filter
      REAL     (kind=4), dimension(:,:), pointer :: snr
      INTEGER :: ii,if,ji,jf,imn,jmn,imx,jmx,iw
      REAL     (kind=4) :: rtot,spmean,spvar,spsdev,exsdev,val
      

      ! Local
       INTEGER,parameter::soumax = 20
      INTEGER :: idim,jdim,ic,jc,is,js,i,j,iloc(2),ili,ilf,jli,jlf
      INTEGER :: ns,nsou,nsousnr,nrpix,ns1,npix,iok
      REAL     (kind=4) :: rmin,rmax,rtot2,spsterr
      REAL     (kind=4) :: rad
    

      Status = 0
      ! Initializations
      iw = 0
    
      

      ! Dimensions
      idim=SIZE(IMAGE,1)
      jdim=SIZE(IMAGE,2)
      ! Telescope axis in image center
      ic=(idim+1)/2
      jc=(jdim+1)/2

      ! Allocation and Default values for filter
      ALLOCATE(filter(1:idim,1:jdim),stat=iok)
      IF (iok /= 0) then
         Status = 1
         if(present(iprint))print *,'Allocation problem for filter'
         return
      endif
      filter=.false.
      ! Default values for region
      IF (.not.PRESENT(reg))THEN
         ii=1 ; if=idim ; ji=1 ; jf=jdim
      ELSE
         ii=INT(reg(1,1)) ; if=INT(reg(2,1)) ; ji=INT(reg(2,1))
         jf=INT(reg(2,2))
      ENDIF
      ! Set filter to true in the interested region
      filter(ii:if,ji:jf)=.true.

      IF (PRESENT(mask)) THEN
         IF(.not.ASSOCIATED(mask)) THEN
         ! problem to know if it will always give ...
           ALLOCATE(mask(1:idim,1:jdim),stat=iok) 
           IF (iok /= 0) then
             Status = 1
             if(present(iprint)) print *,'Allocation problem for mask'
             return
           endif
           mask=filter
         ELSE
           filter=mask
         ENDIF
      ENDIF

      ! Default values for source parameters
      IF (.not.PRESENT(soupar))THEN
         ns=0
         rad=0.
      ELSE
         ns=INT(soupar(1))
         ! ns - number of peaks to be retrieved
         if(ns > souMax)ns = soumax
         rad=soupar(2)
      ENDIF

      ! Default values for source snr

      nsousnr=10  ! Y,Z,flux,err,snr,ctsondet - number of info per source
      !ns - number of source peaks to be found
      ! both outpout source peak array sousnr and work array snr
      ! will have the same size ns= number of peaks asked in soupar(1)
      IF (ns >0 )THEN
         ALLOCATE(snr(1:nsousnr,1:ns),stat=iok)
         snr(:,:)=0.
         IF (PRESENT(sousnr))THEN
         
            ! input array of sources exist
            IF(.not.ASSOCIATED(sousnr)) THEN
               ! output peak list will be put to wanted size
              ALLOCATE(sousnr(1:nsousnr,1:ns),stat=iok) 
              IF (iok /= 0) then
                Status = 1
                if(present(iprint))print *, 'Allocation problem '
                return
             endif         
              sousnr(:,:)=0.         
            ELSE
              ! find number of sources defined in sousnr
              ns1=SIZE(sousnr,2) 
              ! ns1 - number of sources as defined in input array sousnr
              ns1=MIN(ns1,ns)
              ! ns - max number of sources to be searched for
              snr(:,1:ns1)=sousnr(:,1:ns1)
              DEALLOCATE(sousnr)
              ALLOCATE(sousnr(1:nsousnr,1:ns),stat=iok)
              IF (iok /= 0) then
                 Status = 1
                if(present(iprint))print*, 'Allocation problem '
                return
             endif
            ENDIF
         ENDIF
      ENDIF ! sources to be searched for ( ns > 0)


      IF (PRESENT(VAR)) THEN
         exsdev=VAR(ic,jc)
         if (exsdev >= 0. ) exsdev=SQRT(exsdev)
      ELSE
         exsdev=0.
      ENDIF


      

      DO nsou=0,ns
         rmax = -1.e10
         rmin = 1.e10
         npix = 0
         rtot = 0.
         rtot2 = 0.
         do i=1,idim
            do j=1,jdim
               
               if(filter(i,j))then
                  npix = npix+1
                  val = image(i,j)
                  rtot = rtot+val
                  rtot2 = rtot2+val**2
                  if (rmax < val)then
                     rmax = val
                     imx = i
                     jmx = j
                     
                  endif
                  if(rmin > val)then
                     imn = i
                     jmn = j
                     rmin = val
                    
                  endif
               endif
            enddo
         enddo
       
         if(npix > 0) then
            IF (npix > 1) then
               spvar=(rtot2-(rtot*rtot)/REAL(npix))/REAL(npix-1)
            else 
               spvar = 0.
            endif
            
            IF (spvar >= 0.) spsdev=SQRT(spvar)
            spmean=rtot/REAL(npix)
            !      spsterr=spmean/SQRT(REAL(npix))
         else
            if(nsou == 0)then
               Status = 2
               if(present(iprint))print *,'npix = 0 in IBIS_IMASN'
               return
            else
               exit
            endif
         endif
        

         IF (nsou > 0) THEN
            ! Case sousnr(:,nsou)=(flux) is  0. -> get max position
            IF(snr(1,nsou)==0. .and. snr(2,nsou)==0 .and. snr(3,nsou)==0.)THEN
               ! maximum found
               is=imx 
               js=jmx
               snr(1,nsou)=is-ic
               snr(2,nsou)=js-jc             !
            ELSE
            ! Case sousnr(3,nsou)=(flux) is /= 0. ->
               is=NINT(snr(1,nsou))+ic ;  js=NINT(sousnr(2,nsou))+jc
            ENDIF
            snr(3,nsou)=IMAGE(is,js)
            if (PRESENT(VAR)) snr(4,nsou)=SQRT(VAR(is,js))

            ! Set filter to .false. in a rad-wide square around source
            ili=MAX((is-INT(rad)),ii) ; ilf=MIN((is+INT(rad)),if)
            jli=MAX((js-INT(rad)),ji) ; jlf=MIN((js+INT(rad)),jf)
            do i=ili,ilf
            do j=jli,jlf
              
               filter(i,j)=.false.
            enddo
            enddo
          
         ENDIF
         IF (nsou==0 .or. nsou==ns) THEN
            Nrpix=((if-ii+1)*(jf-ji+1))
          
         ENDIF
      ENDDO

      ! Compute flux with final mean, S/N
      ! & Print result if required (iw=2) 
      
      do i=1,ns
         is=snr(1,i)+ic
         js=snr(2,i)+jc
         if (snr(4,i) /= 0.) snr(5,i)=snr(3,i)/snr(4,i)
      enddo

      


      ! Set mask = filter (filter in output)
      if (PRESENT(mask)) mask=filter

      ! source peak list  - if was done
    
      if(ns > 0 ) then 
         if (PRESENT(sousnr)) then
            
            sousnr=snr
         endif
      endif

      IF (PRESENT(imastat)) THEN
         imastat(1) = REAL(npix)
         imastat(2) = rtot
         imastat(3) = spmean
         imastat(4) = spsdev
         imastat(5) = REAL(imx) ; imastat(6) = jmx ;  imastat(7) = rmax
         imastat(8) = REAL(imn) ; imastat(9) = jmn ;  imastat(10) = rmin
      ENDIF

    if(associated(filter)) deallocate(filter)
    if(associated(snr)) deallocate(snr)
    END SUBROUTINE IBIS_IMASN

