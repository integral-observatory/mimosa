

!****************************
MODULE IbisMapSearch_DEFS
!****************************
INTEGER,parameter :: auxDim=10
INTEGER  :: MinCatNumber
INTEGER :: nclea,itest,MapFitres
REAL(kind=4) :: radius,siglim1,siglim2,radidentif,xs,ys,xc,yc
INTEGER :: isky,jsky,icsky,jcsky,nSouToSearch,nSouInCat,is,js,irad,ivar
REAL(kind=4), dimension(:,:), pointer :: tempSnr
INTEGER,parameter:: tempSnrSize = 11
LOGICAL, dimension(:,:), pointer :: FILTER,FILTER2,Filter0
REAL(kind=4), dimension(:,:), pointer :: SKYMOD,BORDERS
REAL(kind=4), dimension(:,:), pointer ::snr
REAL(kind=4),dimension(:,:),pointer ::FanMap
REAL(kind=4), dimension(:), pointer :: dista
INTEGER, dimension(:), pointer :: index
INTEGER, dimension(:), pointer :: sourceList
REAL(kind=4), dimension(:), pointer :: soupar,imastat,catSouFlux
!****************************
END MODULE IbisMapSearch_DEFS
!****************************



!****************************
MODULE SEARCH_CLEAN_AUX_MODULE
!****************************

INTERFACE

SUBROUTINE BETTER_MEAN(scwmapflag,fanflag,idim,jdim,sky,filter0,&
                       FanMap,mean,status)
!----------------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER  :: scwmapflag,fanflag
INTEGER  :: idim,jdim,status
REAL(kind=4),dimension(:,:),pointer :: sky
LOGICAL,dimension(:,:),pointer ::FILTER0
REAL(kind=4),dimension(:,:),pointer :: FanMap
REAL(kind=4) :: mean
END SUBROUTINE BETTER_MEAN



END INTERFACE

!*********************************
END MODULE SEARCH_CLEAN_AUX_MODULE
!*********************************

!**************************
 MODULE SOURCE_IDENT1_MOD
!**************************

INTERFACE

 

   SUBROUTINE SourceGhosts(SigmaSky,Filter,xs,ys,&
             modeli,modelj,modelval,fancount,IfSource) 
   !--------------------------------------------------------------
   IMPLICIT NONE
   REAL(kind=4) ,dimension(:,:),pointer :: SigmaSky
   LOGICAL  ,dimension(:,:),pointer :: Filter
   REAL(kind=4):: xs,ys
   REAL(kind=4),dimension(3,3) :: modelval
   INTEGER,dimension(3) ::modeli,modelj
    INTEGER :: fancount
   LOGICAL        :: IfSource
   END SUBROUTINE SourceGhosts


    SUBROUTINE Identification(xs,ys,radidentif,dista,&
           MinCatNumber,NIdenSou)
    !------------------------------------------------
    IMPLICIT NONE
    ! INPUT/OUTPUT VARIABLES
    REAL(kind=4)                           :: xs,ys,radidentif
    ! xs,ys - source fine position
    ! radidentif - identification radius in pixels
    REAL(kind=4), dimension(:), pointer    :: dista
    ! dista - 0 or distance from all identified catalogue sources
    INTEGER                        :: MinCatNumber,NIdenSou
    ! MinCatNumber - position in catalogue of the nearest source
    ! NIdenSou - number of sources identified with this position
    END SUBROUTINE Identification

    SUBROUTINE  MapFlux(xs,ys,RefFrameNumber,MapPoint,count)
    !------------------------------------------------------------
      IMPLICIT NONE
      !INPUT/OUTPUT VARIABLES
      REAL(kind=4)     :: xs,ys ! position with respect to the image centre
     INTEGER  :: RefFrameNumber
      REAL(kind=4), dimension(:,:), pointer ::MapPoint
      REAL(kind=4)     :: count
     
    END SUBROUTINE  MapFlux

    SUBROUTINE IsScwSource(eband,x,y,scwno)
    !--------------------------------------
      IMPLICIT NONE

      !INPUT/OUTPUT VARIABLES
      INTEGER  :: eband,scwno
      REAL(kind=4):: x,y

   END SUBROUTINE IsScwSource

 END INTERFACE

!****************************
END MODULE SOURCE_IDENT1_MOD
!****************************

!****************************
MODULE SOURCE_IDENT2_MOD
!****************************

INTERFACE

SUBROUTINE SourceFantoms(ncatsou,icsky,jcsky,xs,ys,modeli,modelj,modelval,&
                        IdentNum,PosError,CatNumber,Result,Status)
!----------------------------------------------------------------
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE SOURCE_IDENT1_MOD
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER  :: ncatsou,icsky,jcsky
REAL(kind=4)     :: xs,ys
REAL(kind=4),dimension(3,3)  :: modelval
INTEGER,dimension(3) :: modeli,modelj
INTEGER  :: IdentNum,CatNumber,Result
REAL(kind=4)     :: PosError 
INTEGER  :: Status
END SUBROUTINE SourceFantoms


END INTERFACE

!****************************
END MODULE SOURCE_IDENT2_MOD
!****************************



!******************************
MODULE MAP_SEARCH_AUX_MODULE
!******************************
INTERFACE



SUBROUTINE MapAllocInitPrints(Sky,SourceQualityFlag,Status)
!---------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: Status
INTEGER, dimension(:), pointer :: SourceQualityFlag
REAL(kind=4),dimension(:,:),pointer :: SKY
END SUBROUTINE MapAllocInitPrints


SUBROUTINE MapPeakIdent(OutBand,MapExpo,CatSource,Peak,IdentSou,QQche,Status)
!----------------------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
Logical         :: QQche,CatSource,peak
REAL(kind=4), dimension(:,:), pointer :: MapExpo
INTEGER :: OutBand,IdentSou,Status
END SUBROUTINE MapPeakIdent

SUBROUTINE MapResults(RefFrameNumber,isou,IdentSou,flux,fitflux,locerr,&
     fluxerr,theorFluxErr,sigma,mapexpo,Status)
!----------------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: RefFrameNumber,isou,IdentSou,Status
REAL(kind=4)    :: flux,fitflux,locerr,fluxerr,sigma,theorFluxErr
REAL(kind=4), dimension(:,:), pointer :: mapexpo
END SUBROUTINE MapResults


END INTERFACE
!******************************
END MODULE MAP_SEARCH_AUX_MODULE
!******************************




!###############################################
!# SUBROUTINES CODE FROM SEARCH_CLEAN_AUX_MODULE
!###############################################

!=============================================================
SUBROUTINE BETTER_MEAN(scwmapflag,fanflag,idim,jdim,sky,filter0,&
                       FanMap,mean,status)
!==============================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER  :: scwmapflag,fanflag
INTEGER  :: idim,jdim,status
REAL(kind=4),dimension(:,:),pointer :: sky
LOGICAL,dimension(:,:),pointer ::FILTER0
REAL(kind=4),dimension(:,:),pointer :: FanMap
REAL(kind=4) :: mean
!LOCAL VARIABLES
REAL(kind=4),parameter :: eps = 1.e-4
REAL(kind=4) :: dev,depx,depy
logical :: jest
integer :: i,j,ile,ns,nsou,iok
integer :: ic,jc,is,js,i1,i2,j1,j2,i0,j0
logical,dimension(:,:),pointer :: filter
character(len=20) :: procname

status = ISDC_OK
procname = 'BETTER_MEAN'

ALLOCATE(filter(idim,jdim),stat=iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif 
filter = filter0


if(FanFlag==0)then !FanMap creation
    if(scwmapflag==0)then
       depx = Xdisi
       depy = Xdisj
    else
       depx = Mapdisi
       depy = Mapdisj 
    endif
   FanMap = 0
   ic = (idim+1)/2
   jc = (jdim+1)/2
   ! source and fantom positions exclusion
   nsou = size(InScwCat,2)
   do ns = 1,nsou
      is = nint(InScwCat(1,ns)+ic+depx)
      js = nint(InScwCat(2,ns)+jc+depy) 

      do i=-129,0,129
         do j=-129,0,129
            i0 = is+i
            j0 = js+j
            i1 = max(1,i0-3)
            i2 = min(idim,i0+3)
            j1 = max(1,j0-3)
            j2 = min(jdim,j0+3)
            i1 = min(max(1,i1),idim)
            i2 = min(max(1,i2),idim)
            j1 = min(max(1,j1),idim)
            j2 = min(max(1,j2),idim)
            filter(i1:i2,j1:j2) = .false.
           
            if((i == 0).and.(j==0))then ! source position
               FanMap(i1:i2,j1:j2) = FanMap(i1:i2,j1:j2)+10
            else
               FanMap(i1:i2,j1:j2) = FanMap(i1:i2,j1:j2)+1
            endif
         enddo
      enddo
enddo
else
   where(fanmap > 0)filter = .false.
endif



jest = .true.
do while(jest) 
   mean = sum(sky,filter)/idim/jdim
   dev = 0.
   do i=1,idim
      do j=1,jdim
         dev= dev+(sky(i,j)-mean)**2
      enddo
   enddo
   dev = sqrt(dev/idim/jdim)
!!$   IF(DebugMode.gt.1)then
!!$      print *,' BETTER_MEAN mean dev active : ',mean,dev,&
!!$           real(count(filter))/real(idim)/real(jdim)
!!$   endif
if(dev .lt.eps)return
   jest = .false.
ile = 0
   do i=1,idim
      do j=1,jdim
         if(filter(i,j))then
            if(abs(sky(i,j)-mean)/dev > 5.5)then
               filter(i,j) = .false.
               jest = .true.
               ile = ile+1
            endif
         endif
      enddo
   enddo
enddo

IF(ASSOCIATED(filter))DEALLOCATE (filter)
!============================
END SUBROUTINE BETTER_MEAN
!============================






!##################################################
!   SUBROUTINES CODE FROM  MODULE SOURCE_IDENT1_MOD
!##################################################

!===========================================================================
SUBROUTINE SourceGhosts(SigmaSky,Filter,xs,ys,modeli,modelj,modelval,fancount,&
                      IfSource) 
!=============================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
IMPLICIT NONE
REAL(kind=4) ,dimension(:,:),pointer :: SigmaSky
LOGICAL  ,dimension(:,:),pointer :: Filter
REAL(kind=4):: xs,ys
REAL(kind=4),dimension(3,3) :: modelval
INTEGER,dimension(3) ::modeli,modelj
INTEGER :: fancount
LOGICAL        :: IfSource
!LOCAL VARIABLES
INTEGER :: isou,jsou,ifan,jfan,idim,jdim,&
                   i,j,i1,i2,j1,j2,i0,j0,aktm,&
                   totfan,iloc(2),imax,jmax

Logical :: wpolu

modeli = 0;modelj = 0;modelval = 0.
isou = nint(xs)
jsou = nint(ys)

!fanhei = 0.65

totfan = 0
idim = size(SigmaSky,1)
jdim = size(SigmaSky,2)
i0 = (idim+1)/2
j0 = (jdim+1)/2

fancount = 0


ifan = isou

wpolu = .false.
! szukanie najmniejszego fantomu
do while((ifan.gt.0).and.(Filter(ifan,j0)))
  ifan = ifan-129
  wpolu = .true.
end do
if(wpolu)ifan = ifan+129

modeli(1) = ifan
i=2
ifan = ifan+129
do while((ifan.le.idim).and.(Filter(ifan,j0)).and.(i < 4))
  modeli(i) = ifan
  ifan = ifan+129
  i=i+1
enddo

jfan = jsou
wpolu = .false.
! szukanie najmniejszego fantomu
do while((jfan.gt.0).and.(Filter(i0,jfan)))
   jfan = jfan-129
   wpolu = .true.
enddo
if(wpolu)jfan = jfan+129

modelj(1) = jfan
j=2
jfan = jfan+129
do while((jfan.le.jdim).and.(Filter(i0,jfan)).and.(j < 4))
   modelj(j) = jfan
   jfan = jfan+129
   j=j+1
enddo


do i=1,3
do j=1,3
 aktm = 0.
 if((modeli(i).gt.0).and.(modelj(j).gt.0).and.(filter(modeli(i),modelj(j))))then
       i1 =max( modeli(i)-GhostPeakWidth,1)
       i2 = min(modeli(i)+GhostPeakWidth,idim)
       j1 =max( modelj(j)-GhostPeakWidth,1)
       j2 = min(modelj(j)+GhostPeakWidth,jdim)
       
       modelval(i,j)=maxval(sigmasky(i1:i2,j1:j2),filter(i1:i2,j1:j2))
       aktm =maxval(sigmasky(i1:i2,j1:j2),filter(i1:i2,j1:j2))
  endif
! print *,i,j,modeli(i),modelj(j),modelval(i,j),aktm
        
enddo
enddo

fancount = 0
totfan = 0
do i=1,3
do j=1,3
   if((modeli(i).gt.0).and.(modelj(j).gt.0).and.(filter(modeli(i),modelj(j))))then
      totfan = totfan+1
      if((modelval(i,j).gt.4.5))fancount = fancount+1
   endif
enddo
enddo

iloc = maxloc(modelval)
imax = iloc(1)
jmax = iloc(2)

!print *,' totfan : ',totfan,' valids : ',fancount

ifsource = .false.
if(fancount.ge.4)ifSource = .true.
!if(.not.ifsource)then
!  print *,' no '
!endif
select case(totfan)
case(9,8,7,6,5)
  
 !  print *,'can be FCFOV source'
   select case(fancount)
   case(5,6,7,8,9)
   !   print *,'at least 5 ghosts found - FCFOV source' 
      if((imax.ne.2).or.(jmax.ne.2))then
     !    print *,' problem : max fantom : ',imax,jmax
      endif
      if(sqrt(real((isou-modeli(2))**2+(jsou-modelj(2))**2)).lt.3)then
        ! print *,' source pos OK'
      else
       ! print *,' incorrect source pos, should be :',modeli(2),modelj(2)
     endif
  case(4)
     ! print *,' probable source position :',modeli(imax),modelj(jmax)
  case(3,2,1,0)
   ! print *,' not enough of fantoms '
end select


case(4)
  ! print *,' PCFOV source - possibility of confusion'
   if(fancount.lt.3)then
   !  print *,' not enough of fantoms '
   else
  !   print *,' PCFOV source at  : ',modeli(imax),modelj(jmax)
   endif
case default
   ! print *,' incorrect number of fantoms - not a source'
endselect

!==========================
END SUBROUTINE SourceGhosts
!==========================
!===============================================================
SUBROUTINE Identification(xs,ys,radidentif,dista,&
           MinCatNumber,NIdenSou)
!================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
REAL(kind=4)                           :: xs,ys,radidentif
! xs,ys - source fine position
! radidentif - identification radius in pixels
REAL(kind=4), dimension(:), pointer    :: dista
! dista - 0 or distance from all identified catalogue sources
INTEGER                        :: MinCatNumber,NIdenSou
! MinCatNumber - position in catalogue of the nearest source
! NIdenSou - number of sources identified with this position

!LOCAL VARIABLES
INTEGER      :: ncatsou,icatsou 
REAL(kind=4) :: minDist

nidensou=0
minCatNumber = 0


minDist = radidentif 



!idim = size(dista)
ncatsou = size(InScwCat,2) ! number of sources in the input catalogue
if(ncatsou >0)then ! not empty input catalogue
   dista(:) = SQRT((InScwCat(1,:)-xs)**2 +(InScwCat(2,:)-ys)**2)

   ! search for the best candidate in the input
   !source catalogue

   DO icatsou=1,ncatsou
      IF ( dista(icatsou) < radidentif ) THEN
         if( dista(icatsou) < mindist)then
            minDist=dista(icatsou)
            minCatNumber = icatsou
         endif
         nidensou=nidensou+1
      else
         dista(icatsou) = 0.
      endif
   ENDDO
endif! not empty input catalogue
!=============================
END SUBROUTINE Identification
!=============================


!===============================================================
 SUBROUTINE  MapFlux(xs,ys,RefFrameNumber,MapPoint,&
              count)
!=============================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FLUX_APP_MOD
USE ATTI_DECLARATIONS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL(kind=4)     :: xs,ys ! position with respect to the image centre
INTEGER  :: RefFrameNumber
REAL(kind=4), dimension(:,:), pointer ::MapPoint
REAL(kind=4)     :: count

!LOCAL VARIABLES
INTEGER  :: frame,impsize,jmpsize,i1,i2,j1,j2,ip,jp,k
INTEGER  :: n1,n2,t
REAL(kind=4)     :: x,y,alpha,delta,xp,yp,pci,pcj,mean,flux
INTEGER  :: Status
CHARACTER(len=20):: procName

procName = 'MapFlux'
Status = ISDC_OK

flux = 0.

impsize = size(MapPoint,1)
jmpsize = size(MapPoint,2)
pci=REAL(impsize+1)/2.
pcj=REAL(jmpsize+1)/2.
mean = sum(MapPoint)/impsize/jmpsize
! from map to sphere
call atti_conv_yzad(xs,ys,alpha,delta,RefFrameNumber,ProjType,MagnifFactor,Status)
if(Status.ne.0)then
   call WAR_MESSAGE(procName,&
        'Attitude singularity - cannot approximate the flux',0,Status)
  
   return
endif


k=0
do frame= 1,RefFrameNumber-1
   
   ! from sphere to individual image
   call atti_conv_adyz(alpha,delta,x,y,frame,1,ImageMagnifFactor,Status)
    if(Status.ne.0)then
      call WAR_MESSAGE(procName,&
     'Attitude singularity - cannot approximate the flux',0,Status)
      return
   endif
   flux = flux+UnitSource(x-Xdisi,y-Xdisj)
   k=k+1
enddo
flux = flux/real(k)
!from sphere to point map image
call atti_conv_adyz(alpha,delta,xp,yp,RefFrameNumber,ProjType,PointMagnifFactor,Status)
if(Status.ne.0)then
   call WAR_MESSAGE(procName,&
        'Attitude singularity - cannot approximate the flux',0,Status)
   return
endif
ip = nint(xp+pci+Xdisi)
jp=nint(yp+pcj+Xdisj)
i1 = max(1,ip-1);i2 = min(impsize,ip+1)
j1 = max(1,jp-1);j2 = min(jmpsize,jp+1)
t=0.
do n1=i1,i2
do n2=j1,j2
   if(Mappoint(n1,n2) >0.)t=t+1
enddo
enddo
if(t> 3)then
   count = maxval(MapPoint(i1:i2,j1:j2))
   count = count/flux
else
   count = 0.
endif


!=======================
END  SUBROUTINE  MapFlux
!=======================


!=======================================
SUBROUTINE IsScwSource(eband,x,y,scwno)
!========================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DECLARATIONS
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER  :: eband,scwno
REAL(kind=4):: x,y
!LOCAL VARIABLES
INTEGER  :: nsou,scw ,ns,status
REAL(kind=4):: dist,distrad
LOGICAL :: jest
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'IsScwSource'
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)

distrad = 2.
scwno = 0


scw = 0
jest = .false.
do while((scw  < ScwNumber) .and. (.not.jest))
   scw = scw+1
   nsou = AllSouList(scw,eband)%SouNumber
   ns = 0
   do while((ns  < nsou) .and. (.not.jest))
      ns = ns+1
      dist = sqrt((x-AllSouList(scw,eband)%xmap(ns))**2 + &
           (y-AllSouList(scw,eband)%ymap(ns))**2 )
      if(dist < distrad)then
         jest = .true.
         scwno = scw
      endif
   enddo
enddo
!=========================
END SUBROUTINE IsScwSource
!=========================



!###############################################
!SUBROUTINES CODE FROM MODULE SOURCE_IDENT2_MOD
!###############################################

!=========================================================================
SUBROUTINE SourceFantoms(ncatsou,icsky,jcsky,xs,ys,modeli,modelj,modelval,&
                        IdentNum,PosErr,CatNumber,Result,Status)
!=========================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE SOURCE_IDENT1_MOD
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER  :: ncatsou,icsky,jcsky
REAL(kind=4)     :: xs,ys
REAL(kind=4),dimension(3,3)  :: modelval
INTEGER,dimension(3) :: modeli,modelj
INTEGER  :: IdentNum,CatNumber,Result
REAL(kind=4)     :: PosErr 
INTEGER  :: Status

!LOCAL VARIABLES
REAL(kind=4) :: radidentif,siglim1,souPosErr,fanPosErr,maxFanValue,xsf,ysf
REAL(kind=4), dimension(:), pointer :: dista
INTEGER :: iok, MinCatNumberSou,IdentSou,identFantoms,maxFani,maxFanj,if,jf
INTEGER ::MinCatNumberFan,NIdenFan
INTEGER,dimension(3,3) ::fantomSource
character(len=20)::procName

Status=ISDC_OK
procName = 'Sourcefantoms'
Result = 0
radidentif=FineTol*CleanParTab(6)
siglim1=CleanParTab(4)
ALLOCATE(dista(1:ncatsou),stat=iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif 
dista=0. 
souposerr = 0.
! Source Identification
call Identification(xs-icsky-Xdisi,ys-jcsky-Xdisj,radidentif,dista,&
                    MinCatNumberSou,IdentSou)
if(IdentSou>0)souPosErr =dista( MinCatNumberSou)
!print *,' Sourcefantoms xs,ys dista :',xs,ys,dista( MinCatNumberSou)
! Eventual Ghost Identification
identFantoms = 0
if(FantomAnalysis==1)then
   fantomSource(:,:) = 0
   identFantoms = 0
   maxFanValue = 0.
   maxFani = 0
   maxFanj = 0
   do if=1,3
   do jf=1,3
      if((modeli(if).ne.0).and.(modelj(jf).ne.0).and.&
         (modelval(if,jf).gt.siglim1)) then
         xsf = modeli(if)
         ysf = modelj(jf)
         if((abs(xs-xsf).gt.100).or.(abs(ys-ysf).gt.100))then
             ! Fantom  Identification
              call Identification(xsf-icsky-Xdisi,ysf-icsky-Xdisi,radidentif,dista,&
                    MinCatNumberFan,NIdenFan) 
              if(NIdenFan.gt.0)then ! catalogue source on fantom position
                 fantomSource(if,jf) = MinCatNumberFan
                 identFantoms = identFantoms+1
                  if(modelval(if,jf).gt.maxFanValue)then
                     maxFanValue = modelval(if,jf)
                     maxFani = if
                     maxFanj = jf
                     fanPosErr = dista(MinCatNumberFan)
                  endif
               endif ! catalogue source on fantom position
          endif
       endif
    enddo
    enddo
endif !Fantom Analysis

Result = -1
IdentNum=IdentSou+IdentFantoms
if(IdentSou .gt.0)then ! max peak identified
   Result = 0
   posErr = souPosErr
   CatNumber =  MinCatNumberSou
   if(identFantoms.gt.0)then
      Result = 1
      call MESSAGE(procName,&
           'Possible source at a ghost position:',ZeroError,Status)
      do if=1,3
      do jf=1,3
         if(fantomSource(if,jf).gt.0)then
            write(str250,'(2i4,3x," cat. number: ",i3)',err=460)modeli(if),&
                  modelj(jf),Int(InScwCat(10,fantomSource(if,jf)))
            goto 461
460         str250 = errorstr//' 460'
461         call MESSAGE(procName,str250,ZeroError,Status)
         endif
      enddo
      enddo
   endif ! idenFantoms.gt.0
else   ! max peak NOT identified
     if(identFantoms.gt.0)then !identified source at ghost position
         Result = 2
         write(str250,&
              '("Main peak ",2(f5.1,2x)," not identified but ")',err=462)xs,ys
         goto 463
462      str250 = errorstr//' 462'
463      call MESSAGE(procName,str250,ZeroError,Status)
         if = modeli(maxFani)
         jf = modelj(maxFanj)
         posErr = fanPosErr
         CatNumber = FantomSource(MaxFani,MaxFanj)
         write(str250,&
                    '(" catalogue source ",i3,&
                   & " identified at the ghost position:",2i5)',err=464)&
                     fantomSource(maxFani,maxFanj),if,jf
         goto 465
464      str250 = errorstr//' 464'
465      call MESSAGE(procName,str250,ZeroError,Status)
         xs = if
         ys = jf  
     endif !identified source at ghost position
 endif  ! max peak NOT identified

deallocate(dista)
!============================
END SUBROUTINE SourceFantoms
!============================



!#######################################################
! CODE OF SUBROUTINES FROM  MODULE MAP_SEARCH_AUX_MODULE
!########################################################

!========================================================
SUBROUTINE MapAllocInitPrints(Sky,SourceQualityFlag,Status)
!=========================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IbisMapSearch_DEFS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: Status
INTEGER, dimension(:), pointer :: SourceQualityFlag
REAL(kind=4),dimension(:,:),pointer :: SKY
!LOCAL VARIABLES
INTEGER :: iloc(2),iok,i
REAL(kind=4) :: val
character(len=20) :: procname
procname = 'MapAllocInitPrints'
Status = ISDC_OK

val    = CleanParTab(2)   
nclea  = NINT(val)
! source radius for excluding
radius = CleanParTab(3)/MagnifFactor+1

siglim1=CleanParTab(4)
siglim2=CleanParTab(5)
! source identification radius
radidentif=CleanParTab(6)/MagnifFactor
val    = CleanParTab(10)
itest  =  NINT(val)

!25.08 to make coherent with Ibis SkyClean search
irad=INT(radius) !+1

! Sky dim, center (note integer center pixel for even images)
isky=SIZE(SKY,1) ; jsky=SIZE(SKY,2)
icsky=(isky+1)/2 ; jcsky=(jsky+1)/2


!number of sources in the input catalogue
nSouInCat=SIZE(InScwCat,2)


!number of sources to search for from .par file
nSouToSearch=IABS(nclea)


if(DebugMode.ge.1)then
   str250 = 'Sources in Fov:'
  call message(procname,str250,ZeroError,Status)
endif
do i=1,nSouInCat
   is = nint(IOSPHERE*InScwCat(1,i)+(isky+1)/2)
   js = nint(InScwCat(2,i)+(jsky+1)/2)
   if((is-1>0).and.(is+1<isky+1).and.(js-1>0).and.(js+1<jsky+1))then
      iloc = maxloc(Sky(is-1:is+1,js-1:js+1))
      if(DebugMode.ge.1)then
         write(str250,'(I2,3x,I3," Ra Dec : ",&
              &F6.1,2x,F6.1,"  Y Z : ",F6.1,2x,F5.1,&
              &" dist : ",f5.1," peak : ",F8.1,&
              &" (",2(I3,1x),")")',err=330)&
              i,int(InScwCat(10,i)),InScwCat(4,i),InScwCat(5,i),&
              InScwCat(1,i)+icsky+Xdisi,InScwCat(2,i)+jcsky+Xdisj,&
              sqrt(InScwCat(1,i)**2+InScwCat(2,i)**2)&
              ,maxval(Sky(is-1:is+1,js-1:js+1)),&
              iloc(1)+is-2,iloc(2)+js-2
         goto 331
330      str250 = ' '
331      call message(procname,str250,ZeroError,Status)
      endif
   endif
enddo

!----------------------------
!ALLOCATIONS/INITIALISATIONS
!----------------------------
if(associated(SourceQualityFlag))deallocate(SourceQualityFlag)

ALLOCATE(dista(1:nSouInCat),index(1:nSouInCat),catSouFlux(Sourcenumber),&
         BORDERS(1:isky,1:jsky),sourceList(nSouToSearch),&
         FILTER(1:isky,1:jsky),FILTER2(1:isky,1:jsky),&
         FILTER0(1:isky,1:jsky),&
         SKYMOD(1:isky,1:jsky),FanMap(isky,jsky),&
         soupar(1:auxDim),imastat(1:auxDim),snr(1:auxDim,nSouToSearch),&
         tempSnr(1:nSouToSearch,tempsnrsize),&
         SourceQualityFlag(nSouToSearch),stat=iok)
! SourceQualityFlag - should not be de-allocated at the end of
! IbisSkyClean
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif 

FanMap(:,:) = 0.
sourceList(:) = 0
SourceQualityFlag(:) = 4
catSouFlux(:) = 0.
BORDERS(:,:) = 0.
dista=0. ; index=(/ (i, i=1,nSouInCat) /)
 FILTER(:,:)=.true. ; FILTER2(:,:)=.true.
FILTER0(:,:)=.true.  
SKYMOD(:,:)= 0.
soupar(:)=0. ; imastat(:)=0. ; snr(:,:)=0.
tempSnr(:,:) = 0.


soupar(1)=nSouToSearch      ! number of sources to search for 
soupar(2)=radius    ! source radius

ivar=0
!======================================
END SUBROUTINE MapAllocInitPrints
!======================================

!======================================================================
SUBROUTINE MapPeakIdent(OutBand,MapExpo,CatSource,Peak,IdentSou,QQche,Status)
!=======================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IbisMapSearch_DEFS
USE SOURCE_IDENT1_MOD
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
Logical         :: CatSource,peak,QQche
REAL(kind=4), dimension(:,:), pointer :: MapExpo
INTEGER :: OutBand,IdentSou,Status
!LOCAL VARIABLES
INTEGER :: soulim,sousize,scwno,i,j
REAL(kind=4) :: xd,yd
Logical :: okSource,condition
character(len=20)::procName

Status=ISDC_OK
procName = 'MapPeakIdent'


peak = .false.
QQche=.false.

sousize = size(snr,2)


do soulim=1,sousize
        
   xd = snr(1,soulim) 
   yd = snr(2,soulim)
   xc=xd+REAL(icsky)
   yc=yd+REAL(jcsky)
   xs=xc 
   ys=yc  
   is=NINT(xs) 
   js=NINT(ys)
   ! if(Filter(is,js).and.(SKYMOD(is,js) >  0.))  QQche=.true.
   !27/07/2015 negative peak admitted
    if (Filter(is,js))  QQche=.true.

 enddo

if(.not.QQche)then
   return
endif


okSource = .false.
soulim = 0

!changed in 3.6.5 from .true.
if(MapFluxReNorm)then 
   CatSource = .true.
else
   CatSource = .false.
endif

!preference for  catalogue sources
do while((.not.okSource).and.(soulim.lt.sousize).and.(.not.oksource))
   soulim = soulim+1
   xc=snr(1,soulim)+REAL(icsky)
   yc=snr(2,soulim)+REAL(jcsky)
   xs=xc ; ys=yc ; is=NINT(xs) ; js=NINT(ys)
   !   peak = Filter(is,js).and.(SKYMOD(is,js) > siglim1).and.(snr(3,soulim) >0.)
   !27/07/2015 negative peak admitted

   peak = Filter(is,js).and.(SKYMOD(is,js) > siglim1)!   .and.(snr(3,soulim) >0.)

 identsou=0
!!$      call Identification(xs-icsky-Mapdisi,ys-jcsky-Mapdisj,&
!!$           radidentif,dista,MinCatNumber,IdentSou)
      if(peak)then
         !Source position identified 
         okSource = .true.
        
!!$      else
!!$         okSource = .false.
!!$         identSou = 0
!!$         MinCatNumber = 0
      endif
              
   enddo


!!$   if(.not.okSource)then
!!$      ! no catalogue source near peaks positions
!!$
!!$      catSource = .false.
!!$
!!$      soulim = 0
!!$      peak = .false.
!!$      do while((.not.peak).and.(soulim.lt.sousize))
!!$         soulim = soulim+1
!!$          xd = snr(1,soulim) ; yd = snr(2,soulim)
!!$         xc=xd+REAL(icsky)
!!$         yc=yd+REAL(jcsky)
!!$         xs=xc ; ys=yc ; is=NINT(xs) ; js=NINT(ys)
!!$     
!!$         call IsScwSource(Outband,xd,yd,scwno)
!!$         ! source already found in Scw or several scw contributed 
!!$         ! at this position
!!$         condition = (scwno > 0).or.( MapExpo(is,js) > 1.5*minval(Duration))
!!$         peak = (SKYMOD(is,js) > siglim2).and.condition .and.(filter(is,js))
!!$      enddo
!!$   endif 

   if(.not.peak)then
      !no peak found
      
      do soulim=1,sousize
        
          xd = snr(1,soulim) 
          yd = snr(2,soulim)
          xc=xd+REAL(icsky)
          yc=yd+REAL(jcsky)
          xs=xc 
          ys=yc  
          is=NINT(xs) 
          js=NINT(ys)
          ! Set filter to hide the source in next loop
          irad=INT(soupar(2))
          do i=is-irad,is+irad
             do j=js-irad,js+irad
                FILTER(i,j) = .false.
             enddo
          enddo
         
      enddo
   endif
!===========================
END SUBROUTINE MapPeakIdent
!===========================

!========================================================
SUBROUTINE MapResults(RefFrameNumber,isou,IdentSou,flux,fitflux,locerr,&
     fluxerr,theorFluxErr,sigma,mapexpo,Status)
!=========================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_IMAGING_PARAM
USE IbisMapSearch_DEFS
USE ATTI_DEFS
USE ATTI_Internal
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: RefFrameNumber,isou,IdentSou,Status
REAL(kind=4)    :: flux,fitflux,locerr,fluxerr,sigma,theorFluxErr
REAL(kind=4), dimension(:,:), pointer :: mapexpo
!LOCAL VARIABLES
INTEGER :: isn
REAL(kind=4)    :: d,alpha,delta
REAL(kind=8):: a1,d1,d2,a2,dd,l,b
character(len=20) :: procname
procname = 'MapResults'
Status = ISDC_OK

! flux - total flux
tempSnr(isou,2) = IdentSou
tempSnr(isou,3) = xs-icsky-Mapdisi ! coord in the telescope frame
tempSnr(isou,4) = ys-jcsky-Mapdisj

tempSnr(isou,5) = fitflux ! counts/sec

tempSnr(isou,6) = fitflux*totaltime/IActiveDetDim/JActivedetDim/2.!counts/pixel
tempSnr(isou,7) = fitflux
tempSnr(isou,8) = locerr
tempSnr(isou,10)= sigma
tempSnr(isou,11)= theorFluxErr /totaltime    

if(MapFitRes==0)then
   if(sigma >0.)then
      ! fix G.B. problem 9.04
      tempSnr(isou,8)= 22.1*sigma**(-0.95)+0.16
    
      tempSnr(isou,8)=tempSnr(isou,8)/60./PixAngDim
      
      ! estimation of loc err from article       
   endif
endif


! tempsnr contents :
!  1 
!  2 > 0 if catalogue source 
!  3,4, distance in pixel units from the detector centre ( X telescope axis)
!  5  counts/sec  ( possibly modified by spectral fit)
!  6  counts/pixel    "
!  7  fitted flux
!  8 location error
!  9
! 10  significance
! 11 flux error               


IF(IdentSou > 0) THEN !identified source
   isn= Int(InScwCat(10,minCatNumber))
   sourceList(isou) = isn
   
catSouFlux(isn)=catSouFlux(isn)+tempSnr(isou,5)
 
 call Message(procName,'.......................................................',&
        ZeroError,Status)
   
   write(str250,'(" ::  source no. ",I3," signif. :",&
        &F6.1," InCat no ",I3,2x,A20)',err=472)&
        isou,sigma,isn,IncatName(isn)
   goto 473
472 str250 =errorstr//' 472'
473 call Message(procName,str250,ZeroError,Status)
   d = sqrt((real(is-icsky))**2+(real(js-jcsky))**2)
   write(str250,'("    dist. cts/pix  cts/sec Y Z :",&
        &F8.1,2X,2(F7.2,2X),2(F6.1,2X))',err=474)d,&
        tempSnr(isou,6), tempSnr(isou,5),xs,ys  
   goto 475
474 str250 =errorstr//' 474'
475 call Message(procName,str250,ZeroError,Status)


!!$   call ATTI_CONV_YZAD(tempSnr(isou,3),tempSnr(isou,4),&
!!$        alpha,delta,RefFrameNumber,1,MagnifFactor,status)

   call PixSkyProj(tempSnr(isou,3),tempSnr(isou,4),alpha,delta,Status)
   if(Status.ne.0)then
      call  WAR_Message(procName,&
           ' Cannot calculate fitted ra dec - set to 0. ',0,Status)
      alpha = 0.
      delta = 0.
   endif

   a1 = InCatAlpha(isn)*deg_to_rad
   d1 = InCatDelta(isn)*deg_to_rad
   a2 = alpha*deg_to_rad
   d2 = delta*deg_to_rad
   dd=sla_dsep(a1,d1,a2,d2)*rad_to_deg
   if(FitOffset==1)then 
      tempsnr(isou,8)=dd
   endif
   dd=dd*60.*60.

   write(str250,'(    "  Ra Dec cat. and found: ",2(F8.2,2X),2(F10.4,2X)," pix.err. :",F5.2," ang.err :",F7.4,"arc sec")',err=476)&
        InCatAlpha(isn),InCatDelta(isn),alpha,delta,dista(minCatNumber),dd
   goto 477
476 str250 =errorstr//' 476'
477 call Message(procName,str250,ZeroError,Status)

 call D_CONV_AD_LB(dble(InCatAlpha(isn)),dble(InCatDelta(isn)),l,b)
 write(str250,'(a,f7.2,1x,f7.2)')' Galactic coord :',real(l),real(b)
 call Message(procName,str250,ZeroError,Status)

   tempSnr(isou,1) = minCatNumber 
   tempSnr(isou,9) = dista(minCatNumber)
else ! not identified source
   call Message(procName,'........................................',&
        ZeroError,Status)
     
   write(str250,'(" :: source no. ",I3," signif. :",F6.1,&
        &"UNIDENTIFIED ")',err=478)isou,sigma
   goto 479
478 str250 =errorstr//' 478'
479 call Message(procName,str250,ZeroError,Status)
   d = sqrt((real(is-icsky))**2+(real(js-jcsky))**2)
   write(str250,'("    dist. cts/pix  cts/sec Y Z :",&
        &F5.1,2X,2(F7.2,2X),2(F7.2,2X))',err=480)d,&
        tempSnr(isou,6), tempSnr(isou,5),xs,ys 
   goto 481
480 str250 =errorstr//' 480'
481 call Message(procName,str250,ZeroError,Status)
   call ATTI_CONV_YZAD(tempSnr(isou,3),tempSnr(isou,4),&
        alpha,delta,RefFrameNumber,1,MagnifFactor,status)
   if(Status.ne.0)then
      call  WAR_Message(procName,&
           ' Cannot calculate fitted ra dec - set to 0. ',0,Status)
      alpha = 0.
      delta = 0.
   endif
   write(str250,'(" Ra Dec :  found ",2(F10.4,2X) )',err=482)alpha,delta
   goto 483
482 str250 =errorstr//' 482'
483 call Message(procName,str250,ZeroError,Status)
   tempSnr(isou,1) = 0.
   tempSnr(isou,9) = 0.
endif! not identified source

!==========================
END SUBROUTINE MapResults
!==========================

