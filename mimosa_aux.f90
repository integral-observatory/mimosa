


!*************************************
MODULE MIMOSA_AUX_MODULE
!*************************************

INTERFACE

SUBROUTINE ReadOffCorrBand(nBand, IdxPtr,&
     BinNumber, Bins, Arr,Status)
!--------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: nBand, IdxPtr
INTEGER ::  BinNumber,Status
REAL(kind=8),dimension(:,:),pointer::Bins
REAL(kind=4 ),dimension(:,:),pointer :: Arr
END SUBROUTINE ReadOffCorrBand

SUBROUTINE MapFovSources(iMapSize,jMapSize,IEnergyBand,&
                         iDetDim,jDetDim,Status)
!-------------------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)::iMapSize,jMapSize,IEnergyBand,&
                          iDetDim,jDetDim
!OUTPUT VARIABLES
INTEGER,Intent(out) ::Status
END SUBROUTINE MapFovSources

SUBROUTINE MapCoordinates(iMapSize,jMapSize,xymaplist,Status)
!-------------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)::iMapSize,jMapSize
!OUTPUT VARIABLES
INTEGER,Intent(out) ::Status
real(kind=4),dimension(:,:),pointer :: xymaplist
END SUBROUTINE MapCoordinates



SUBROUTINE ReadRealFitsTab(FileName,Tab,Status)
!---------------------------------------------
IMPLICIT NONE
!OUTPUT VARIABLES
CHARACTER(len=*) :: FileName
REAL(kind=4),dimension(:,:),pointer :: Tab
INTEGER :: Status
END SUBROUTINE ReadRealFitsTab


SUBROUTINE  PixSkyProj(y,z,alpha,delta,Status)
!---------------------------------------------------
IMPLICIT NONE

!INPUT/ OUTPUT VARIABLES
REAL(KIND=4)::y,z
REAL(KIND=4)::alpha,delta
INTEGER,Intent(out) ::Status
END SUBROUTINE  PixSkyProj

END INTERFACE
!*************************************
END MODULE MIMOSA_AUX_MODULE
!*************************************




!*************************************
MODULE LOCAL_MODULE
!*************************************
INTERFACE
SUBROUTINE ReadReal4Array(ArrPtr,iDim,jDim,Arr,Status)
!-----------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in) :: ArrPtr,iDim,jDim
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: Arr
INTEGER            :: Status
END SUBROUTINE ReadReal4Array

SUBROUTINE  SkyPixProj(alpha,delta,y,z,Status)
!----------------------------------------------------
IMPLICIT NONE

!INPUT/ OUTPUT VARIABLES
REAL(KIND=4)::y,z
REAL(KIND=4)::alpha,delta
INTEGER,Intent(out) ::Status
END SUBROUTINE SkyPixProj



END INTERFACE

!*************************************
END MODULE LOCAL_MODULE
!*************************************



!##########################
!EXTERNAL SUBROUTINES CODE 
!##########################

!================================================================
SUBROUTINE  SkyPixProj(alpha,delta,y,z,Status)
!================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DECLARATIONS
USE ATTI_DEFS
USE ATTI_INTERNAL
IMPLICIT NONE

!INPUT/ OUTPUT VARIABLES
REAL(KIND=4)::y,z
REAL(KIND=4)::alpha,delta
INTEGER,Intent(out) ::Status

!LOCAL VARIABLES

REAL(KIND=8):: pixel_ang,l,b,yn,zn,ad,dd
CHARACTER(len=20)::procName

procName ='SkyPixProj'
Status = ISDC_OK

pixel_ang=rad_to_deg*atan(pixel_ratio)
ad=dble(alpha)
dd=dble(delta)
if(EquaGal==0)then
   !EQUAT
   if(ProjType==1)then
      !TAN
      CALL ATTI_CONV_ADYZ(alpha,delta,y,z,point_num,ProjType,MagnifFactor,Status)
   else
      !CAR
      Call D_CAR_ADYZ(ad,dd,yn,zn,pointing_table(point_num,1),pointing_table(point_num,2),pixel_ang) 
      y = real(yn)
      z = real(zn)
   endif
else
   !GAL
   call D_CONV_AD_LB(ad,dd,l,b)
   Call D_CAR_ADYZ(l,b,yn,zn,lx(point_num),bx(point_num),pixel_ang) 
   y = real(yn)
   z = real(zn)
endif

!==========================
END SUBROUTINE SkyPixProj
!==========================


!================================================================
SUBROUTINE  PixSkyProj(y,z,alpha,delta,Status)
!================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DECLARATIONS
USE ATTI_DEFS
USE ATTI_INTERNAL
IMPLICIT NONE

!INPUT/ OUTPUT VARIABLES
REAL(KIND=4)::y,z
REAL(KIND=4)::alpha,delta
INTEGER,Intent(out) ::Status

!LOCAL VARIABLES

REAL(KIND=8):: pixel_ang,l,b,yn,zn,ad,dd
CHARACTER(len=20)::procName

procName ='PixSkyProj'
Status = ISDC_OK

pixel_ang=rad_to_deg*atan(pixel_ratio)

yn=dble(y)
zn=dble(z)

if(EquaGal==0)then
   !EQUAT
   if(ProjType==1)then
      !TAN
     call ATTI_CONV_YZAD(y,z,alpha,delta,point_num,1,MagnifFactor,status)
   else
      !CAR
      Call D_CAR_YZAD(yn,zn,ad,dd,pointing_table(point_num,1),pointing_table(point_num,2),pixel_ang)
     alpha = real(ad)
     delta = real(dd)
   endif
else
   !GAL
   Call D_CAR_YZAD(yn,zn,l,b,lx(point_num),bx(point_num),pixel_ang)
   call D_CONV_LB_AD(l,b,ad,dd)
   alpha = real(ad)
   delta = real(dd)
endif

!==========================
END SUBROUTINE PixSkyProj
!==========================


!================================================================
SUBROUTINE MapFovSources(iMapSize,jMapSize,IEnergyBand,&
                         iDetDim,jDetDim,Status)
!================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DECLARATIONS
USE ATTI_DEFS
USE ATTI_INTERNAL
use  LOCAL_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)::iMapSize,jMapSize,iEnergyBand,&
                           iDetDim,jDetDim
!OUTPUT VARIABLES
INTEGER,Intent(out) ::Status


!LOCAL VARIABLES
INTEGER::i,ns,k,i1,i2,iok,scw,sn,badsou
REAL(KIND=4)::alpha,delta,y,z,flux,di,dj,ymap,zmap
REAL(KIND=8):: pixel_ang,l,b,yn,zn
CHARACTER(len=20)::procName
procName ='MapFovSources'
Status = ISDC_OK

if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)

 pixel_ang=rad_to_deg*atan(pixel_ratio)

!half map length
di =  real(iMapSize)/2.
dj = real(jMapSize)/2.         
ns=0
InSourceList = 0.
! SourceNumber - number of sources in input catalogue           
badsou = 0
DO i=1,SourceNumber
 
   alpha=InCatAlpha(i) ; delta=InCatDelta(i)
   call SkyPixProj(alpha,delta,y,z,Status)
   if(Status.ne.0)then
      badsou = badsou+1
      write(str250,*,err=115)&
           ' Coordinates calculation problem for the source ',i  
      goto 116
115   str250 = 'Coordinates calculation problem'&
           //errorstr//' 115'
116   call WAR_MESSAGE(procName,str250,0,Status)
      Status = 0
      y = -iMapSize-10.
      z = -jMapsize-10.
   endif
  

   ! y,z, in telescope frame
   !transformation to the image frame
   ymap = y+MapDisi
   zmap = z+MapDisj
   IF((ABS(ymap).LE.di).AND.(ABS(zmap).LE.dj)) THEN
      ns=ns+1
      InSourceList(1,ns)  = y ! in telescope frame
      InSourceList(2,ns)  = z
      InSourceList(4,ns)  = alpha
      InSourceList(5,ns)  = delta
      InSourceList(10,ns) = i
     
   ENDIF
ENDDO

if(badsou >0)then
      if(MapOnlyMode.ne.3)then 

         call war_message(procname,&
              'Some sources too far from map centre',Zeroerror,status)
      endif
endif 
!number of sources in map FOV
ScwSourceNumber=ns
MapSourceNumber=ns
! output energy band contains input zones between i1 and i2
i = 1
do while((i.le.InEnergyNumber).and.&
         (InToOutBandNums(i).lt.iEnergyBand).and.&
          (InToOutBandNums(i).ne.0))
 i = i+1
enddo
i1 = i
do while((i.le.InEnergyNumber).and.&
         (InToOutBandNums(i).lt.iEnergyBand+1).and.&
          (InToOutBandNums(i).ne.0))
 i = i+1
enddo
i2 = i-1

ns = 0
do i=1,ScwSourceNumber
    ! k -  number of source in the input catalogue
    ! i -  number of source in InSourceList
    k =INT( InSourceList(10,i))
    flux = sum(InCatFlux(k,i1:i2))
    if(FluxMode.ge.1)flux = PredefFlux
    if(SimulMode==0)then
       ns = ns+1
       InSourceList(9,i) = 1.
    else ! Simulations mode
       IF(flux.gt.0)then
          ns = ns+1
          InSourceList(3,i) = flux*FluxNorm
          !total flux
          if(SimulMode==1)then
             InSourceList(6,i) = InSourceList(3,i)*REAL(iDetDim*jDetDim)
          else
             InSourceList(6,i) = InSourceList(3,i)
          endif
          InSourceList(9,i) = 1.
       else
          InSourceList(9,i) = 0.
       endif
    endif ! Simulations mode
enddo
EnergySourceNumber = ns

IF(ASSOCIATED(InScwCat)) DEALLOCATE (InScwCat)
ALLOCATE(InScwCat(1:10,EnergySourceNumber),stat=iok)


IF (iok /= 0) then
      call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
endif
InScwCat(:,:) = 0.

if(.not.associated(AllSouMapFlux))then
   ! allocation for first output energy band
   allocate(AllSouMapFlux(SourceNumber,OutEnergyNumber),&
            AllSouMapSnr(SourceNumber,OutEnergyNumber),stat=iok)
   IF (iok /= 0) then
      call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
   endif
   AllSouMapFlux(:,:) = 0.
   AllSouMapSnr(:,:) = 0.
endif
ns=0
do i=1,ScwSourceNumber
   if(InSourceList(9,i).gt.0)then
       ns = ns+1
       InScwCat(1:10,ns) = InSourceList(1:10,i)

    endif
enddo

!CALCULATION OF FINE MAP COORDIANTES OF ALL SOURCES FOUND
!IN THE GIVEN ENERGY BAND
do scw=1,ScwNumber
   ! loop on scw
  sn = AllSouList(scw,iEnergyBand)%sounumber
  do i=1,sn
     ! loop on sources
     IF (CoorType >= 1) THEN
        alpha=AllSouList(scw,iEnergyBand)%alpha(i) 
        delta=AllSouList(scw,iEnergyBand)%delta(i)
        CALL ATTI_CONV_ADYZ&
             (alpha,delta,y,z,point_num,ProjType,MagnifFactor,Status)
        if(Status.ne.0)then
           write(str250,*,err=117)&
                ' Coordinates calculation problem for the source (2)',i 
           goto 118
117        str250=' Coordinates calculation problem '&
                    //errorstr//' 117'
118        call WAR_MESSAGE(procName,str250,0,Status)
           y = -iMapSize-10.
           z = -jMapsize-10.
        endif
    
     ELSEIF (CoorType >= 0) THEN
        y=1500. ;z=1500.
     ENDIF
     AllSouList(scw,iEnergyBand)%xmap(i) = y
     AllSouList(scw,iEnergyBand)%ymap(i) = z
  enddo
enddo
!==========================
END SUBROUTINE MapFovSources
!==========================



!==========================================
SUBROUTINE ReadRealFitsTab(FileName,Tab,Status)
!==========================================

USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!OUTPUT VARIABLES
CHARACTER(len=*) :: FileName
REAL(kind=4),dimension(:,:),pointer :: Tab
INTEGER :: Status
!LOCAL VARIABLES
INTEGER::ptr
INTEGER,dimension(1:DAL_MAX_ARRAY_DIMENSION)::axes
INTEGER :: type,numAxes,numValues,idim,jdim
INTEGER,dimension(2)::startVal,endVal
CHARACTER(len=20)::procName


Status = ISDC_OK
procName = 'ReadRealFitsTab'
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)

idim = size(Tab,1)
jdim = size(tab,2)

Status = dal_object_open(FileName,ptr,Status)
if(Status.ne.ISDC_OK)then
  call MESSAGE(procName,'cannot open file'//FileName,&
       OpenError,Status)
  return
endif

Status = dal_array_get_struct(ptr,type,numAxes,axes,&
            Status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,' dal_array_get_struct problem ',&
       IsdcExitStatus,Status)
   return
endif
if((type.ne.DAL_FLOAT).or.(numAxes.ne.2).or.&
   (axes(1).ne.idim).or.(axes(2).ne.jdim))then
       call MESSAGE(procName,'  array invalid ',&
       InvalArrayError,Status)
   return
endif
!reading
startVal = 1
endVal(1) = axes(1)
endVal(2) = axes(2)
numValues = axes(1)*axes(2)
Status = dal_array_get_section(ptr,numAxes,startVal,endVal,&
         DAL_Float,numValues,addrof(Tab(1,1)),status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,' dal_array_get_section problem ',&
         IsdcExitStatus,Status)
   return
endif

Status = dal_object_close(ptr,DAL_SAVE,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
  call MESSAGE(procName,'cannot close file'//FileName,&
      IsdcExitStatus,Status)
  return
endif
!===========================
END SUBROUTINE ReadRealFitsTab
!===========================






!================================================================
SUBROUTINE MapCoordinates(iMapSize,jMapSize,xymaplist,Status)
!================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DECLARATIONS
USE ATTI_DEFS
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)::iMapSize,jMapSize
!OUTPUT VARIABLES
INTEGER,Intent(out) ::Status
real(kind=4),dimension(:,:),pointer :: xymaplist

!LOCAL VARIABLES
INTEGER::i,ns,k,i1,i2,iok,scw,sn
REAL(KIND=4)::alpha,delta,y,z,flux,di,dj,ymap,zmap
CHARACTER(len=20)::procName
procName ='MapCoordinates'
Status = ISDC_OK

if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)
!half map length
di =  real(iMapSize)/2.
dj = real(jMapSize)/2.         
ns=0

if (associated(xymaplist))deallocate(xymaplist)
allocate(xymaplist(SourceNumber,2))

! SourceNumber - number of sources in input catalogue           
DO i=1,SourceNumber
  
     alpha=InCatAlpha(i) ; delta=InCatDelta(i)
     CALL ATTI_CONV_ADYZ(alpha,delta,y,z,point_num,ProjType,MagnifFactor,Status)
    if(Status.ne.0)then
       Status = 0
        y = -iMapSize-10.
        z = -jMapsize-10.
       print *,procname,' bad source ',i
     endif
    
  
     
     ! y,z, in telescope frame
     !transformation to the image frame
     ymap = y+MapDisi
     zmap = z+MapDisj
     IF((ABS(ymap).LE.di).AND.(ABS(zmap).LE.dj)) THEN
  
       ns=ns+1
       xymaplist(i,1)  = y + di+mapdisi
       xymaplist(i,2)  = z + dj+mapdisj
  endif    

ENDDO


!==========================
END SUBROUTINE MapCoordinates
!==========================

!===================================================================
SUBROUTINE ReadOffCorrBand(nBand, IdxPtr,&
     BinNumber, Bins, Arr,Status)
!=====================================================================

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE LOCAL_MODULE

IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: nBand, IdxPtr
INTEGER ::  BinNumber,Status

REAL(kind=8),dimension(:,:),pointer::Bins
REAL(kind=4 ),dimension(:,:),pointer :: Arr
!LOCAL Variables

REAL(kind=8) :: e1,e2,knorm,kc,maxnorm,minenergy,maxenergy
INTEGER ::iDim,jDim, n,nFound,ptr,iok
LOGICAL :: jest
character(len=20) ::row
CHARACTER(len=10) :: str
CHARACTER(len=20)  :: procName


procName = 'ReadrOffCorrBand'
Status = ISDC_OK

e1 = OutEnergyBands(nBand,1)
e2 = OutEnergyBands(nBand,2)

idim = size(Arr,1)
jdim = size(Arr,2)

Arr(:,:) = 1.

minenergy = minval(Bins)
maxenergy = maxval(bins)



if((e1 <  minenergy).or.(e2 <  minenergy).or.&
     (e1 > maxenergy).or.(e2 >maxenergy))then
   call MESSAGE(procName,' Cannot find OffCorr  bin ',OffCorrError,Status)
   return
endif

jest = .false.
n=1

do while((n .le. BinNumber) .and.(.not.jest))
   if ((e1 .ne. Bins(n,1)).or.( e2 .ne. Bins(n,2)))then
      n=n+1
   else
      jest = .true.
   endif
enddo



if (.not.jest)then
    call MESSAGE(procName,' Uncorrect OffCorr  maps',OffCorrError,Status)
    return
endif

write(row,'(I1.1)',err=135)n
goto 136
135 call message(procname,' cannot write number in row',&
         StringWriteError,Status)
return
136 row = '#row=='//row
status = dal3gen_index_find_member(IdxPtr,OffCorrStrucName,row,nFound,&
     ptr,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' dal3gen_index_find_member problem for'//OffCorrStrucName,&
        IsdcExitStatus,Status)
   return
endif
  
call ReadReal4Array(ptr,iDim,jDim,Arr,Status)
if(Status.ne.ISDC_OK)return
  


!===========================
END SUBROUTINE ReadOffCorrBand
!===========================



!===================================================
SUBROUTINE ReadReal4Array(ArrPtr,iDim,jDim,Arr,Status)
!===================================================

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in) :: ArrPtr,iDim,jDim
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: Arr
INTEGER            :: Status

!LOCAL VARIABLES
INTEGER,dimension(1:DAL_MAX_ARRAY_DIMENSION)::axes
INTEGER :: type,numAxes,numValues
INTEGER,dimension(2)::startVal,endVal
CHARACTER(len=20) :: procName


procName = 'ReadReal4Array'
Status = ISDC_OK
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)
 
! array size verification
Status = dal_array_get_struct(ArrPtr,type,numAxes,axes,&
            Status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,' dal_array_get_struct problem ',&
       IsdcExitStatus,Status)
   return
endif

if((numAxes.ne.2).or.&
   (axes(1).ne.idim).or.(axes(2).ne.jdim))then
       call MESSAGE(procName,' array invalid ',&
       InvalArrayError,Status)
   return
endif
 
!reading
startVal = 1
endVal(1) = axes(1)
endVal(2) = axes(2)
numValues = axes(1)*axes(2)

type = DAL_FLOAT
Status = dal_array_get_section(ArrPtr,numAxes,startVal,endVal,&
         type,numValues,addrof(Arr(1,1)),status)

if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,' dal_array_get_section problem ',&
         IsdcExitStatus,Status)
   return
endif

if(DebugMode.eq.3)then
   write(str250,'(" read real array  ofdim ",i5,2x,i5)',&
        err = 117)axes(1),axes(2)
   goto 118
117 str250 = errorstr//' : 117'
118 call Message(procName,str250,ZeroError,Status)
endif
if(DebugMode.eq.3)&
  call Message(procName,' end',ZeroError,Status)
!===========================
END SUBROUTINE ReadReal4Array
!===========================


