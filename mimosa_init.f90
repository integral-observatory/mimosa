

!*********************************
MODULE MIMOSA_INIT_MODULE
!*********************************

INTERFACE

   SUBROUTINE InitWork(Status)
!-------------------------------------------
IMPLICIT NONE
INTEGER :: Status
END SUBROUTINE InitWork

SUBROUTINE ReadParFile(Status)
!-------------------------------
IMPLICIT NONE

!OUTPUT VARIABLES
INTEGER   :: Status

END SUBROUTINE ReadParFile

SUBROUTINE ReadOffCorrIndex(IdxPtr,BinNumber,&
     Bins,Status)
!----------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL(kind=8),dimension(:,:),pointer::Bins
INTEGER :: IdxPtr, BinNumber,Status
END SUBROUTINE ReadOffCorrIndex


SUBROUTINE InitArrays(iSky,jSky,iEdge,jEdge,&
            Shd,ShdVar,ShdEffi,&
            GroupShd,GroupShdVar,GroupShdEffi,&
            Wei,Sky,SkyCle,SkyResid,SkySignif,SkyVar,&
            Status)
!-----------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)::iSky,jSky

!OUTPUT VARIABLES
INTEGER,Intent(out)::Status,iEdge,jEdge
REAL(kind=4),dimension(:,:),pointer :: Shd,ShdVar,ShdEffi
REAL(kind=4),dimension(:,:),pointer ::GroupShd,GroupShdVar,GroupShdEffi
REAL(kind=4),dimension(:,:),pointer :: Sky,SkyCle,Wei
REAL(kind=4),dimension(:,:),pointer :: SkyResid,SkySignif,SkyVar
END SUBROUTINE InitArrays

SUBROUTINE CleanStaInit(RaX,DecX,RaZ,DecZ,AttiNumPoi,nclean,Status)
!-----------------------------------------------------------------
IMPLICIT NONE
!INPUT VARIABLES
INTEGER,Intent(in)           ::AttiNumPoi,nclean
REAL(kind=8), dimension(:), pointer ::RaX,DecX,RaZ,DecZ
!OUTPUT VARIABLES
INTEGER,Intent(out)          ::Status
END SUBROUTINE CleanStaInit

SUBROUTINE AttiInitialise(MosaNumber,raX,decX,raZ,decZ,posAngle,&
                          AttiNumPoi,&
                          iMapSize,jMapSize,Status)
!------------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
REAL(kind=8), dimension(:), pointer ::RaX,DecX,RaZ,DecZ,posAngle
INTEGER,Intent(in)  :: MosaNumber
!OUTPUT VARIABLES
INTEGER,Intent(out)  :: AttiNumPoi
INTEGER,Intent(out)  :: iMapSize,jMapSize
INTEGER,Intent(out)  :: Status
END SUBROUTINE AttiInitialise


SUBROUTINE DESCRIPTION(mosanumber,Rax,Decx,imapsize,jmapsize)
!---------------------
IMPLICIT NONE
INTEGER :: mosanumber,imapsize,jmapsize
REAL(KIND=8),DIMENSION(:),POINTER::raX,decX
END SUBROUTINE DESCRIPTION
END INTERFACE

!*********************************
END MODULE MIMOSA_INIT_MODULE
!*********************************

!***************
MODULE INIT_PROC
!***************

INTERFACE

SUBROUTINE CalcMapSize(MosaNumber,ip,mf,iMapSize,jMapSize,iprint,Status)
!-------------------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER,Intent(in)  :: MosaNumber,ip,iprint
REAL(kind=4),Intent(in)     ::  mf
INTEGER,Intent(out) :: iMapSize,jMapSize,Status
END SUBROUTINE CalcMapSize

 SUBROUTINE CalcRefFrame(AttiNumPoi,RaX,DecX,RaZ,DecZ)
 !--------------------------------------------------------
 IMPLICIT NONE
 INTEGER,Intent(in) ::AttiNumPoi
 REAL(kind=8),dimension(:),pointer:: RaX,DecX,RaZ,DecZ
END SUBROUTINE CalcRefFrame

 SUBROUTINE CalcLBRefFrame(AttiNumPoi,RaX,DecX,RaZ,DecZ)
 !--------------------------------------------------------
 IMPLICIT NONE
 INTEGER,Intent(in) ::AttiNumPoi
 REAL(kind=8),dimension(:),pointer:: RaX,DecX,RaZ,DecZ
END SUBROUTINE CalcLBRefFrame
END INTERFACE


!*******************
END MODULE INIT_PROC
!*******************

!***************
MODULE READ_POM
!***************

INTERFACE
SUBROUTINE ReadPar(parName,par,Status)
!---------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER   :: Status
CHARACTER(LEN=30)  :: parName
CHARACTER(LEN=*) :: par
END SUBROUTINE ReadPar

END INTERFACE


!*******************
END MODULE READ_POM
!*******************

!***************
MODULE LOCAL_INIT
!***************
INTERFACE
SUBROUTINE  CalcMapSizeOneLine(num1,ipr,imaCi,imaCj,minI,maxI,minJ,maxJ,&
     pixel1,pixel2,pixel_ang1,pixel_ang2,&
     arefimage,drefimage,lrefimage,brefimage,&
     aref,dref,lref,bref,ri,&
     rot_mat1,rot_mat2,rot_mat3,rot_mat4,Status)
!-------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER  :: num1,ipr,imaCi,imaCj
INTEGER   :: minI,maxI,minJ,maxJ
REAL(kind=8)      :: pixel1,pixel2,pixel_ang1,pixel_ang2
REAL(kind=8)      :: aref,dref,lref,bref,ri
REAL(kind=8)      :: arefimage,drefimage,lrefimage,brefimage
Real(kind=8),Dimension(3,3)::rot_mat1,rot_mat2,rot_mat3,rot_mat4
INTEGER::Status
END SUBROUTINE  CalcMapSizeOneLine

SUBROUTINE CalcMapSizeOneRow(num1,ipr,imaCi,imaCj,minI,maxI,minJ,maxJ,&
     pixel1,pixel2,pixel_ang1,pixel_ang2,&
     arefimage,drefimage,lrefimage,brefimage,aref,dref,lref,bref,rj,&
     rot_mat1,rot_mat2,rot_mat3,rot_mat4,Status)
!-----------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER  :: num1,ipr,imaCi,imaCj
INTEGER   :: minI,maxI,minJ,maxJ
REAL(kind=8)      :: pixel1,pixel2,pixel_ang1,pixel_ang2
REAL(kind=8)      :: aref,dref,lref,bref,rj
REAL(kind=8)      :: arefimage,drefimage,lrefimage,brefimage
Real(kind=8),Dimension(3,3)::rot_mat1,rot_mat2,rot_mat3,rot_mat4
INTEGER::Status
END SUBROUTINE CalcMapSizeOneRow

SUBROUTINE ReadCovr(CovarTab,Covar1Tab,Status)
!.---------------------------------------------  
IMPLICIT NONE
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::CovarTab,Covar1Tab
INTEGER :: Status
END SUBROUTINE  ReadCovr

END INTERFACE
!***************
END MODULE LOCAL_INIT
!***************


!###################################
! MODULE INIT_PROC   SUBROUTINES CODE 
!####################################

!===================================================================
SUBROUTINE CalcMapSizeOneLine(num1,ipr,imaCi,imaCj,minI,maxI,minJ,maxJ,&
     pixel1,pixel2,pixel_ang1,pixel_ang2,&
     arefimage,drefimage,lrefimage,brefimage,aref,dref,lref,bref,ri,&
     rot_mat1,rot_mat2,rot_mat3,rot_mat4,Status)
!===================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER  :: num1,ipr,imaCi,imaCj
INTEGER   :: minI,maxI,minJ,maxJ
REAL(kind=8)      :: pixel1,pixel2,pixel_ang1,pixel_ang2
REAL(kind=8)      :: arefimage,drefimage,lrefimage,brefimage
REAL(kind=8)      :: aref,dref,lref,bref,ri
Real(kind=8),Dimension(3,3)::rot_mat1,rot_mat2,rot_mat3,rot_mat4
INTEGER::Status

!LOCAL VARIABLES

INTEGER   ::num2
Real(KIND=8)      ::rj,xn,yn,pixj
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'CalcMapSizeOneLine'

pixj = RealMapSize(num1,2,1)+0.5
do while (pixj.lt.RealMapSize(num1,2,2)+1.0)
   rj=pixj-imaCj
   call image_to_sky_projection(InputProjType(num1),ProjType,pixel1,pixel2,&
              pixel_ang1,pixel_ang2,arefimage,drefimage,lrefimage,brefimage,&
              aref,dref,lref,bref,&
              rot_mat1,rot_mat2,rot_mat3,rot_mat4,ri,rj,xn,yn,Status)
   if(Status.ne.0)then
      if(ipr==1)then
         write(str250,'(" Mosa no. ",I3,&
              &" too far from map centre ")')num1
         goto 111
110      str250 = 'input mosa too far from map centre'
111      call WAR_Message(procName,str250,0,Status)
         str250 = '- will be only partally projected onto final mosaicks'
         call WAR_Message(procName,str250,0,Status)
        write(str250,*)' Add. infos : arefimage,drefimage,aref,dref, : ',arefimage,drefimage,aref,dref
        call Message(procName,str250,0,Status)
         ipr = 0 
      endif
      Status = 0 
   else
      if(xn.lt.minI)minI=xn
      if(xn.gt.maxI)maxI=xn
      if(yn.lt.minJ)minJ=yn
      if(yn.gt.maxJ)maxJ=yn
!!$      print '(6f10.2)',ri,rj,aproj,dproj,xn,yn
   endif
   pixj = pixj+1. 
enddo


!=================================
END SUBROUTINE CalcMapSizeOneLine
!=================================

!===================================================================
SUBROUTINE CalcMapSizeOneRow(num1,ipr,imaCi,imaCj,minI,maxI,minJ,maxJ,&
     pixel1,pixel2,pixel_ang1,pixel_ang2,&
     arefimage,drefimage,lrefimage,brefimage,aref,dref,lref,bref,rj,&
     rot_mat1,rot_mat2,rot_mat3,rot_mat4,Status)
!===================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER  :: num1,ipr,imaCi,imaCj
INTEGER   :: minI,maxI,minJ,maxJ
REAL(kind=8)      :: pixel1,pixel2,pixel_ang1,pixel_ang2
REAL(kind=8)      :: aref,dref,lref,bref,rj
REAL(kind=8)      :: arefimage,drefimage,lrefimage,brefimage
Real(kind=8),Dimension(3,3)::rot_mat1,rot_mat2,rot_mat3,rot_mat4
INTEGER::Status

!LOCAL VARIABLES

INTEGER   :: num2
Real(KIND=8)      ::ri,xn,yn,pixi
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'CalcMapSizeOneRow'

pixi = RealMapSize(num1,1,1)-0.5
do while (pixi.lt.RealMapSize(num1,1,2)+1.0)
   ri=pixi-imaCi
    call image_to_sky_projection(InputProjType(num1),ProjType,&
              pixel1,pixel2,&
              pixel_ang1,pixel_ang2,arefimage,drefimage,lrefimage,brefimage,&
              aref,dref,lref,bref,&
              rot_mat1,rot_mat2,rot_mat3,rot_mat4,ri,rj,xn,yn,Status)

   if(Status.ne.0)then
      if(ipr==1)then
         write(str250,'(" Mosa no. ",I3,&
              &" too far from map centre ")')num1
         goto 111
110      str250 = 'input mosa too far from map centre'
111      call WAR_Message(procName,str250,0,Status)
         str250 = '- will be only partally projected onto final mosaicks'
         call WAR_Message(procName,str250,0,Status)
         ipr = 0 
      endif
      Status = 0 
   else
      if(xn.lt.minI)minI=xn
      if(xn.gt.maxI)maxI=xn
      if(yn.lt.minJ)minJ=yn
      if(yn.gt.maxJ)maxJ=yn
   endif
   pixi = pixi+1. 
enddo
!=================================
END SUBROUTINE CalcMapSizeOneRow
!=================================

!===================================================================
SUBROUTINE CalcMapSize(MosaNumber,ip,mf,iMapSize,jMapSize,iprint,Status)
!===================================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER,Intent(in)  :: MosaNumber,ip,iprint
REAL(kind=4),Intent(in)     ::  mf
INTEGER,Intent(out) :: iMapSize,jMapSize,Status
!LOCAL VARIABLES
INTEGER   :: imaCi,imaCj,toobig,cutmap,imc,jmc,im1,im2,jm1,jm2
INTEGER   :: i,j,i1,j1,iskysize,jskysize
REAL(kind=8)      :: pixel1,pixel2,pixel_ang1,pixel_ang2
INTEGER   :: minI,maxI,minJ,maxJ,ipr
INTEGER   :: num1,num2
Real(kind=8),Dimension(3,3)::rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4
Real(KIND=8)      :: ri,rj,xn,yn,pixi,pixj,aref,dref,lref,bref
Real(KIND=8)      :: arefimage,drefimage,lrefimage,brefimage
Real(KIND=8)      :: al,del,alpha1,delta1,alpha2,delta2,d,l1 ,l2,b1,b2
REAL(kind=4):: y,z,ytab(4),ztab(4)
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'CalcMapSize'
if(DebugMode.eq.3)&
     call Message(procName,' ',ZeroError,Status)                 


pixel1 = pixel_ratio
pixel2 = pixel_ratio*mf

pixel_ang1 = pixel_ang_size
pixel_ang2 = rad_to_deg*atan(pixel2)
!initialisation of map size
minI = 0;maxI = 0;minJ = 0;maxJ = 0
!last pointing - OG  map centre
num2 =MosaNumber+1
!rot_mat2 = point_rot_mat(num2,1:3,1:3)
do i=1,3
do j=1,3
   rot_mat2(i,j) = point_rot_mat(num2,i,j)
enddo
enddo

rt2(1:3,1:3) = Transpose(rot_mat2(1:3,1:3))
   
aref = POINTING_TABLE(num2,1)
dref = pointing_table(num2,2)

write(str250,*)' num2,aref,dref : ',num2,aref,dref
call message(procname,str250,0,status)

lref = LX(num2)
bref = BX(num2)


ipr = iprint
call message(procname,' Final mosaics size calculating ...',Zeroerror,status)
!verification for all image positions
do num1=1,MosaNumber

   arefimage = POINTING_TABLE(num1,1)
   drefimage = pointing_table(num1,2)
   lrefimage = LX(num1)
   brefimage = BX(num1)
  
   imaCi = defmapsize(num1,1)/2.
   imaCj =defmapsize(num1,2)/2.

    
   do i=1,3
   do j=1,3
      rot_mat1(i,j) = point_rot_mat(num1,i,j)
   enddo
   enddo
   rt(1:3,1:3) = Transpose(rot_mat1(1:3,1:3))
   Call sla_dmxm(rot_mat2,rt,rot_mat3)
   Call sla_dmxm(rot_mat1,rt2,rot_mat4)
       
   pixi = RealMapSize(num1,1,1)- 0.5
   ri=pixi-imaCi
   call CalcMapSizeOneLine(num1,ipr,imaCi,imaCj,&
        minI,maxI,minJ,maxJ,pixel1,pixel2 ,pixel_ang1,pixel_ang2,&
        arefimage,drefimage,lrefimage,brefimage,&
        aref,dref,lref,bref,ri,rot_mat1,rot_mat2,&
        rot_mat3,rot_mat4,Status)
  
   pixi = RealMapSize(num1,1,2)+ 0.5
   ri=pixi-imaCi
   call CalcMapSizeOneLine(num1,ipr,imaCi,imaCj,&
        minI,maxI,minJ,maxJ,pixel1,pixel2 ,pixel_ang1,pixel_ang2,&
         arefimage,drefimage,lrefimage,brefimage,&
        aref,dref,lref,bref,ri,rot_mat1,rot_mat2,&
        rot_mat3,rot_mat4,Status)
  
   pixj = RealMapSize(num1,2,1)- 0.5
   rj=pixj-imaCj
   call CalcMapSizeOneRow(num1,ipr,imaCi,imaCj,minI,maxI,minJ,maxJ,&
     pixel1,pixel2,pixel_ang1,pixel_ang2,&
     arefimage,drefimage,lrefimage,brefimage,aref,dref,lref,bref,rj,&
     rot_mat1,rot_mat2,rot_mat3,rot_mat4,Status)

   pixj = RealMapSize(num1,2,2)+ 0.5
   rj=pixj-imaCj
   call CalcMapSizeOneRow(num1,ipr,imaCi,imaCj,minI,maxI,minJ,maxJ,&
        pixel1,pixel2,pixel_ang1,pixel_ang2,&
        arefimage,drefimage,lrefimage,brefimage,aref,dref,lref,bref,rj,&
        rot_mat1,rot_mat2,rot_mat3,rot_mat4,Status)
enddo ! end of pointing loop

iMapSize = max(abs(minI),maxI)*2+2
JMapSize = max(abs(minJ),maxJ)*2+2

if(2*(iMapSize/2).eq.iMapSize)iMapSize=iMapSize+1
if(2*(jMapSize/2).eq.jMapSize)jMapSize=jMapSize+1
 
!!$iMapSize=890
!!$jmapsize = 1234
!!$print *,' hardcoded map size'
!!$   
toobig = 0
im1 = iMapSize
jm1 = jMapSize
if(iMapSize > MaxMapISize)then
       iMapSize=MaxMapISize
       toobig = 1
endif
if(jMapSize > MaxMapJSize)then
       jMapSize=MaxMapJSize
       toobig = 1
endif
if((toobig==1).and.(MapParDec(3)==0))then
   write(str250,*,err=101)&
        ' Whole map (',im1,'x',jm1,') is too big and will be cut '
   goto 102
101 str250 = ' Whole map is too big and will be cut '
102 if(ip==1)call message(procname,str250,ZeroError,Status)
endif
cutmap=0
if((MapParDec(3)==1).and.(MapParval(3)< 180))then
   !map radius contolled by user
   cutmap=1
     do i=1,4
        select case(i)
        case(1)
           i1 = 1
           j1 = 1
        case(2)
           i1 = -1
           j1 = 1
        case(3)
           i1 = -1
           j1 = -1
        case(4)
           i1 = 1
           j1 = -1
        end select
        al  = aref+ i1*MapParVal(3)
        del = dref+ j1*MapParVal(3)
        if (al <0)al=al+360.0d0
         al = mod(al,360.0d0)
         if(del > deltalimit)del = deltalimit
         if(del < -deltalimit)del = -deltalimit
!!$        call ATTI_CONV_ADYZ(al,del,ytab(i),ztab(i),num2,projType,mf,status)
        call SkyPixProj(real(al),real(del),ytab(i),ztab(i),Status)
        if(status .ne.0) then
           cutmap =0
           status = 0
           exit
        endif

     enddo
     if(cutmap==1)then
        im1 = iMapSize
        jm1 = jMapsize
        iMapSize = int(sum(abs(ytab))/4.)*2
        jMapsize = int(sum(abs(Ztab))/4.)*2
        if(iMapSize > MaxMapISize)then
           iMapSize=MaxMapISize
        endif
        if(jMapSize > MaxMapJSize)then
           jMapSize=MaxMapJSize
        endif
        write(str250,*,err=104)' User map radius = ',MapParVal(3),' ==> map will be cut to ',iMapSize,'x',jMapSize,' pixels'
        goto 105
104     str250 = ' User map radius set ==> map will be cut '
105     if(ip==1)call message(procname,str250,ZeroError,Status)
     endif

endif !map radius contolled by user
!============================
END SUBROUTINE CalcMapSize
!============================

!=====================================================
SUBROUTINE CalcRefFrame(AttiNumPoi,RaX,DecX,RaZ,DecZ)
!=====================================================
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER,Intent(in) ::AttiNumPoi
REAL(kind=8),dimension(:),pointer:: RaX,DecX,RaZ,DecZ
!LOCAL VARIABLES
Integer::i,scwnum,decyd,k,status,iactive
REAL(kind=8)::ra,dec,minra,mindec,maxra,maxdec,val,sumrax,sumdex,valdex,a,d
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'CalcRefFrame'

scwnum=AttiNumPoi-1

minra = 1000.
maxra = -1000.

mindec =  1000.
maxdec = -1000.
sumrax = 0.
sumdex = 0.
iactive = 0


write(str250,*)' scwnum, rax(scwnum-1),rax(scwnum),decx(scwnum-1),decx(scwnum) :',&
     scwnum, rax(scwnum-1),rax(scwnum),decx(scwnum-1),decx(scwnum)
call message(procname,str250,0,status)

do i=1,scwnum
!!$    write(str250,*)' i ra dec : ',i,rax(i),decx(i)
!!$    call message(procname,str250,0,status)
  
     val = RaX(i)
     valdex = DecX(i)
    
     
     iactive = iactive +1
     sumrax = sumrax+val
     sumdex = sumdex + valdex
     if(val < minra)minra = val
     if(val > maxra)maxra = val

     if(valdex < mindec)mindec = valdex
     if(valdex > maxdec)maxdec = valdex

enddo

dec = sumdex/iactive
!!$write(str250,*)' calculated mean dec : ',dec
!!$call message(procname,str250,0,status)
!!$write(str250,*)'maxra, minra : ',maxra, minra
!!$call message(procname,str250,0,status)
if(abs(maxra-minra) .le. 180.)then
    ra = sumrax/iactive
else
   ra = 0.
   do k=1,scwnum
    
         val = rax(k)
         if(val < 180)then
            ra = ra+val
         else
            ra = ra+val-360.
         endif
    
   enddo
   ra = ra/iactive
endif
!!$write(str250,*)' calculated mean ra : ',ra
!!$call message(procname,str250,0,status)

if(ra < 0.)ra = ra+360.0d0
!!$write(str250,*)'Projtype, MapParDec(1:2): ',Projtype,MapParDec(1:2)
!!$call message(procname,str250,0,status)

if(ProjType==1)then
   if(sum(MapParDec(1:2)) < 2)then
      RaX(AttiNumPoi) = ra
      DecX(AttiNumPoi) = dec
   else
      !user map centre
      decyd = 1

    

      if(abs(maxra-minra) .le. 180.)then
         !normal case
         if((MapParVal(1) < minra-MaxDiffRadius).or.&
              (MapParVal(1) > maxra+MaxDiffRadius))then
            decyd = 0
         endif
      
      else
        
         !near 0.
         minra = 400.
         maxra = -400.
         do k=1,scwnum
        
            val = rax(k)
            if(val < 180)val = val+360.
            if(val < minra)minra = val
            if(val > maxra)maxra = val
       
         enddo
         maxra = maxra-360.
         minra = minra-360.
      val = MapParVal(1)
      if(val > 180.)val = val-360.

      if((val < minra-MaxDiffRadius).or.&
         (val > maxra+MaxDiffRadius))then
         decyd = 0
      endif
   endif!near 0.

   if((MapParVal(2) < mindec-MaxDiffRadius).or.&
      (MapParVal(2) > maxdec+MaxDiffRadius))then
      decyd = 0 
   endif

   write(str250,*) ' ra,dec : ',ra,dec
   call message(procname,str250,ZeroError,Status)

   if(decyd==0)then
      str250 = ' User given map centre too far from barycentre ==> wil be re-set'
      call message(procname,str250,ZeroError,Status)
      MapParDec(1:2) = 0
     
      RaX(AttiNumPoi) = ra
      DecX(AttiNumPoi) = dec
   else
      !user map centre accepted
      RaX(AttiNumPoi)= MapParVal(1) 
      DecX(AttiNumPoi) =MapParVal(2)
   endif
 endif !user map centre  


else
   !CAR projection
   if(EquaGal==0)then
      !RA/DEC__CAR
      DecX(AttiNumPoi) = 0.
      if( MapParDec(1)==1)then
         RaX(AttiNumPoi) =MapParVal(1)
      else
         RaX(AttiNumPoi) =ra
      endif
   else
      !LB_CAR
      !will be calculated in CalcLBRefFrame
   endif
endif!CAR projection


!===============================
  END SUBROUTINE CalcRefFrame
!===============================

!=====================================================
SUBROUTINE CalcLBRefFrame(AttiNumPoi,RaX,DecX,RaZ,DecZ)
!=====================================================
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER,Intent(in) ::AttiNumPoi
REAL(kind=8),dimension(:),pointer:: RaX,DecX,RaZ,DecZ
!LOCAL VARIABLES
Integer::i,scwnum,decyd,k,status,iactive
REAL(kind=8)::ra,dec,minra,mindec,maxra,maxdec,val,sumrax,sumdex,valdex,l,b,a,d
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'CalcLBRefFrame'

scwnum=AttiNumPoi-1

do i=1,scwnum
   val = RaX(i)
   valdex = DecX(i)
   call D_CONV_AD_LB(val,valdex,l,b)
   lx(i) = l
   bx(i) = b
enddo

minra = 1000.
maxra = -1000.

mindec =  1000.
maxdec = -1000.
sumrax = 0.
sumdex = 0.
iactive = 0

do i=1,scwnum
 
     val = LX(i)
     valdex =BX(i)
    iactive = iactive +1
     sumrax = sumrax+val
     sumdex = sumdex + valdex
     if(val < minra)minra = val
     if(val > maxra)maxra = val

     if(valdex < mindec)mindec = valdex
     if(valdex > maxdec)maxdec = valdex

enddo

dec = sumdex/iactive

if(abs(maxra-minra) .le. 180.)then
    ra = sumrax/iactive
else
   ra = 0.
   do k=1,scwnum
    
         val = rax(k)
         if(val < 180)then
            ra = ra+val
         else
            ra = ra+val-360.
         endif
    
   enddo
   ra = ra/iactive
endif
if(ra < 0.)ra = ra+360.0d0

l=ra
b=dec

if(ProjType==1)then
   !TAN projection
   if(sum(MapParDec(1:2)) < 2)then
      LX(AttiNumPoi) = l
      BX(AttiNumPoi) = b
   else
      call D_CONV_AD_LB(RaX(AttiNumPoi),DecX(AttiNumPoi),l,b)
      LX(AttiNumPoi)= l
      BX(AttiNumPoi) =b
   endif!user map centre  
 else
    !CAR projection
    if(EquaGal==0)then
      !RA/DEC__CAR
       call D_CONV_AD_LB(RaX(AttiNumPoi),DecX(AttiNumPoi),l,b)
       LX(AttiNumPoi)= l
      BX(AttiNumPoi) =b
   else
       !LB_CAR
      if( MapParDec(1)==1)then
          LX(AttiNumPoi) = MapParVal(1)
          BX(AttiNumPoi) = 0.
          call D_CONV_LB_AD(LX(AttiNumPoi),BX(AttiNumPoi),ra,dec)
          RaX(AttiNumPoi) = ra
          DecX(AttiNumPoi) = dec
       else
          LX(AttiNumPoi) =l
          BX(AttiNumPoi) = b
          call D_CONV_LB_AD(l,b,ra,dec)
          RaX(AttiNumPoi) = ra
          DecX(AttiNumPoi) = dec
       endif
   endif!CAR projection LB_CAR
 endif 




!===============================
  END SUBROUTINE CalcLBRefFrame
!===============================

!##########################
!EXTERNAL SUBROUTINES CODE 
!##########################


!==========================================================
 SUBROUTINE InitWork(Status)
!==========================================================


USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_DECON_LIB_MOD
USE IBIS_IMAGING_PARAM
USE FITS_DECLARATIONS
USE ATTI_DEFS
USE LOCAL_INIT

IMPLICIT NONE

INTEGER :: Status


!LOCAL  VARIABLES

CHARACTER (len=80) :: bpName159,decodName,pcName
INTEGER    ::iok,imadim,jmadim,i,j
INTEGER    ::bpi,bpj,decodi,decodj,pci,pcj
real(kind=8),dimension(-5:5)::u
CHARACTER(len=20)::procName


covpattern(1,1,1,2) = 0.75
covpattern(1,2,1,1) = 0.75
covpattern(2,1,2,2) = 0.75
covpattern(2,2,2,1) = 0.75

covpattern(1,1,2,1) = 0.75
covpattern(1,2,2,2) = 0.75
covpattern(2,1,1,1) = 0.75
covpattern(2,2,1,2) = 0.75

covpattern(1,1,2,2) = 0.75**2
covpattern(1,2,2,1) = 0.75**2
covpattern(2,1,1,2) = 0.75**2
covpattern(2,2,1,1) = 0.75**2

u(0)=1.
u(1) = 0.731111
u(2) = 0.285715
u(3) = 0.0596829
u(4) = 0.00666397
u(5) = 0.000397725
do  i=-5,-1 
   u(i) = u(-i)
enddo
do i=-5,5
do j=-5,5
   ucoeff(i,j) = u(i)*u(j)
enddo
enddo
message_level = 0


SourceCat = 0
ScwDurationMode = 0
AttiMode = 0

GhostPeakWidth = 1
ScwDurationMode = 0 ! scw duration taken from array



  ! Equagal=1

   OneScwMode=0


    DebugMode = 0
    
    ! images created as to seen from inside    
     IOSphere = -1

    ! telescope sight axis displacement with respect to the image centre
    XDisi = 0.5                 
    XDisj = 0.5    
    MapDisi = 0.0                 
    MapDisj = 0.0    
    ! attitude from OG
    AttiMode = 0
    
    ! projection type in the final sky map
    !ProjType = 0
    ! rebinning type 
    if(PixSpread==0)then
       RebinType = 5 !if no variance weighting ==>0 
    else
       RebinType = 7
    endif

  

    ! ratio image/map/point map  pixel /detector pixel size
    ImageMagnifFactor = 1.
    MagnifFactor = 1.0
    PointMagnifFactor=1.0

    ! Scw duration from OG
    if(CalibMode==0)ScwDurationMode = 0
    ! save random generator routines state
    GRAIN_INIT=1
    ! newsequence
    RepetableGeneratorMode = 1
    ! normal treatment
    EmptyMode = 0
    !no prints from atti_lib
    message_level=0

    if(cleanMode == -1) then
       MapMode = 1
    else
       MapMode = 0
    endif
    
  

    ! 0 map from cleaned images
    ! 1 map from raw images
    ! 2 map from resid images

    ! can be next modified by the hidden file
   
    !  Scw trated in individual analysis

   
    OneBandMode = 0
    !print *,procname,' only one band treated'

    BandStart=1
    BandStop=2


    MapOnlyMode = 0
    select case(DoPart2)
    case(0)
        MapOnlyMode = 1
    case(1)
        MapOnlyMode = 0
    case(2)
        MapOnlyMode = 3
    end select

    ! 0 -  Scw+Map
    ! 1 - Scw only
    ! 2 - Map only
    ! 3 - Scw (only shd reading) + Map

    

    !source searching mode - read from .par file
    ! SearchMode = 1
    ! number of sources to be searched - read from .par file
    ! nToClean=1
    ! analysis of source ghosts
    FantomAnalysis = 1
    ! minimum peak height in sigma for catalogue source

    ! PeakMinHeight and NonIdentPeakMinHeight read from .par file 
    !PeakMinHeight = 6.
    ! minimum peak height in sigma for non catalogue source
    !NonIdentPeakMinHeight = 7


    ! source identification radius 
    SouSearchRadius    = 1.5
   
    ! precision increase for fine position source increase
    FineTol=1.   
    ! image simulated  - read from .par file
    !  SimulMode = 1
    ! input OG source catalog if 0
    SourceCat = 0
    soux = -112.0
    souy = 94.5
   
    predefbkg = 0.
    predefflux = 100.
    
    FitType = 1
    ! 0 PSF fit
    ! 1 gauss


    ! flux predefined
    FluxMode  = 1
    !in simulated image dead zones are filling 
    !with detector mean
    DeadZonesFilling = 1

    ! source pattern drawn
    SourcePattern= 0
    ! FFT arrays created
      FFTArrays =0


    !image normalized 
    ! 0 - old formula
    ! 1 - new formula
    !-1 - no normalisation
     AreaMode = 1

   
    !detector type - ISGRI
    DetType=1 
    !deconvolution type
 
  
 call ReadCovr(CovarTab,Covar1Tab,Status)
  if(Status.ne.ISDC_OK)return


If (AttiMode.ne.-1)then
  Interface_type=1
else 
  Interface_type=0
endif


END SUBROUTINE InitWork
!======================================


!..............................................
SUBROUTINE ReadCovr(CovarTab,Covar1Tab,Status)
!...............................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
USE DAL3GEN_F90_API  
IMPLICIT NONE
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::CovarTab,Covar1Tab
INTEGER :: Status
!LOCAL VARIABLES
INTEGER,dimension(1:DAL_MAX_ARRAY_DIMENSION)::axes
INTEGER :: type,numAxes,numValues,iok
INTEGER,dimension(2)::startVal,endVal
INTEGER::idxPtr,ptr,mem,nFound
CHARACTER(len=200):: selection
CHARACTER(len=250) :: str250
CHARACTER(len=20)::procName

Status = ISDC_OK



procName = 'ReadCovr'
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)


allocate(covartab(400,400),covar1tab(400,400),stat=iok)
    if(iok.ne.0)then
      call MESSAGE(procName,'Allocation problem ',AllocError,Status)
      return
endif
covartab(:,:) = 0.
covar1tab(:,:) = 0.

Status = dal_object_open(CovrFile,idxPtr,Status)
if(Status.ne.ISDC_OK)then
  call MESSAGE(procName,'cannot open covariance file',&
       OpenError,Status)
  return
endif

do mem = 1,2

   select case(mem)
   case(1)
      selection = "COVRTYPE=='CARTESIAN'"
   case(2)
      selection = "COVRTYPE=='DIAGONAL'"
   end select

   status = dal3gen_index_find_member(IdxPtr,&
            CovStrucName,selection,nFound,ptr,Status)
  
   if(Status.ne.ISDC_OK) then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
           'dal_object_find_element problem for'//CovStrucName,&
           IsdcExitStatus,Status)
      return
   endif
   if(nFound.ne.1)then
      write(str250,*)' Uncorrect memebr num. :',nfound,' found for ',&
          CovStrucName,selection
       call MESSAGE(procName,str250,&
           IsdcprocError,Status)
      return
   endif
   Status = dal_array_get_struct(ptr,type,numAxes,axes,&
            Status)

   if(Status.ne.ISDC_OK) then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
           ' dal_array_get_struct problem for '//CovStrucName,&
           IsdcExitStatus,Status)
      return
   endif

   if((type.ne.DAL_FLOAT).or.(numAxes.ne.2).or.&
        (axes(1).ne.400).or.(axes(2).ne.400))then
      call MESSAGE(procName,'Covariance  array invalid ',&
           InvalArrayError,Status)
      return
   endif

   !reading
   startVal = 1
   endVal(1) = axes(1)
   endVal(2) = axes(2)
   numValues = axes(1)*axes(2)
   select case(mem)
   case(1)
      Status = dal_array_get_section(ptr,numAxes,startVal,endVal,&
           type,numValues,addrof(CovarTab(1,1)),status)
   case(2)
      Status = dal_array_get_section(ptr,numAxes,startVal,endVal,&
           type,numValues,addrof(Covar1Tab(1,1)),status)
   end select
   if(Status.ne.ISDC_OK) then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
           ' dal_array_get_section problem for '//CovStrucName,&
           IsdcExitStatus,Status)
      return
   endif


enddo

Status = dal_object_close(idxptr,DAL_SAVE,status)
if(Status.ne.ISDC_OK) then
   call WAR_MESSAGE(procName,' cannot close covariance file ',&
       0,Status)
  
endif


END SUBROUTINE ReadCovr


!====================================
SUBROUTINE ReadPar(parName,par,Status)
!====================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER   :: Status
CHARACTER(LEN=30)  :: parName
CHARACTER(LEN=*) :: par
!LOCAL VARIABLES
CHARACTER(LEN=PIL_LINESIZE) :: str
CHARACTER(LEN=20)  :: procName

Status = ISDC_OK
procName = 'ReadPar'

par = ''
str = ' ' 
Status = PilGetString(parName(1:len_trim(parName)), str)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,'PilGetString problem for'//parName,&
        IsdcExitStatus,Status)
   return
endif
par = trim(str)
!======================
END SUBROUTINE ReadPar
!======================

!=========================================================
SUBROUTINE ReadParFile(Status)
!=========================================================

USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE READ_POM
IMPLICIT NONE


!OUTPUT VARIABLES
INTEGER   :: Status


!LOCAL VARIABLES
INTEGER            :: iok,i,l,clpar,dolnum
Real(kind=4),dimension(:),pointer::itab
CHARACTER(LEN=30)  :: parName
CHARACTER(LEN=PIL_LINESIZE) :: str
CHARACTER(LEN=80) :: str80
CHARACTER(LEN=20)  :: procName

Status = ISDC_OK
procName = 'ReadParFile'
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)

JestOffCorrIndex = .false.
offcorrfilename = 'NONE'
parName ='corrDol'
str = ' ' 
Status = PilGetString(parName(1:len_trim(parName)), str)
if(Status.ne.ISDC_OK)then
   ! off-axis correction file parameter does not exist
   call WAR_MESSAGE(procName,&
        'off-axis correction parameter : corrDol cannot be read',&
        ZeroError,Status)
   JestOffCorrIndex = .false.
 
else
   JestOffCorrIndex = .true.
   offcorrfilename = str
   offcorrfile = offcorrfilename
   l=index(offcorrfilename,'(')
   if(l.gt.0) then
      offcorrfile = offcorrfilename(1:l-1)//'[1]'
    endif
    l=index(str,'[')
   if(l.gt.0) then
      offcorrfilename = str(1:l-1)
    endif
endif

parName = 'covrMod'
   call ReadPar(parName,CovrFile,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
    return 
   endif

parName = 'inCat'
   call ReadPar(parName,InCatFile,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
    return 
   endif
call message(procname,incatfile,0,status)
parName = 'outCat'
call ReadPar(parName,OutCatFile,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   return
endif

   parName = 'outMosIma'
   call ReadPar(parName,OutMosImaFile,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
     return
   endif


   parName = 'outMosRes'
   call ReadPar(parName,OutMosResFile,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      return  
   endif




!number of sources to search for

parName = 'ToSearch'
Status = PilGetInt(parName, nToClean)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,'PilGetInt problem for'//parName,&
        IsdcExitStatus,Status)
   return
endif



PixSpread=1
!!$!PixSpread
!!$parName = 'PixSpread'
!!$Status = PilGetInt(parName, PixSpread)
!!$if(Status.ne.ISDC_OK)then
!!$      Status = ISDC_OK     
!!$endif
!!$if(PixSpread .ne. 0) then
!!$     PixSpread=1
!!$endif

ReCall=0
parName = 'ReCall'
Status = PilGetInt(parName, ReCall)
if(Status.ne.ISDC_OK)then
      Status = ISDC_OK  
      call message(procname,'No mosaics recall',   0,Status)
      ReCall=0
endif
if((ReCall < 0) .or.(ReCall > 1)) then
      call message(procname,'Uncorrect value of recall - will be set to 0',   0,Status)
      ReCall=0
endif
if(Recall==1)then
   call message(procname,'Re-call of mosaicking - some data should be already written onto the disk',   0,Status)
endif

DoPart2=1
ParName = 'DoPart2'
Status = PilGetInt(parName,DoPart2)
if(Status.ne.ISDC_OK)then
      Status = ISDC_OK  
      DoPart2=1
endif
if((DoPart2 < 0).or.(DoPart2 > 1))then
   call message(procname,'Uncorrect value of DoPart2 - will be set to 1',   0,Status)
   DoPart2=1   
endif

FitOffset = 0
ParName = 'FitOffset'
Status = PilGetInt(parName,FitOffset)
if(Status.ne.ISDC_OK)then
      Status = ISDC_OK  
      FitOffset=0
endif
if(( FitOffset < 0).or.(FitOffset  > 1))then
   call message(procname,'Uncorrect value of FitOffset - will be set to 0',   0,Status)
   FitOffset=0  
endif

if (FitOffset==1) then
   call message(procname,'Fit offset in deg will be put into fin_rd_error',0,status)
endif


EnergyBand=0
parName = 'EnergyBand'
Status = PilGetInt(parName,EnergyBand)
if(Status.ne.ISDC_OK)then
      Status = ISDC_OK  
      call message(procname,'All energy bands will be treated',   0,Status)
endif
if(Energyband < 0) then
      call message(procname,'Uncorrect value of EnergyBand -All energy bands will be treated ',   0,Status)
      EnergyBand=0
endif



parName = 'UserProj'
Status = PilGetInt(parName, UserProj)
if(Status.ne.ISDC_OK)then
      Status = ISDC_OK     
endif
if((UserProj < 0).or.(UserProj > 2)) then
    UserProj =0
     str250 = 'IUncorrect value of UserProj ==> will be set to 0'
      call WAR_MESSAGE(procName,str250,0,Status)
endif

select case(UserProj)
case(0)
   EquaGal=1
   ProjType=0
case(1)
   EquaGal=0
   ProjType=0
case(2)
   EquaGal=0
   ProjType=1
end select
   

MapParDec(:) = 0
MapParVal(:) = 0.


 parName = 'MapCentre1'
   Status = pilgetreal4(parName,MapParVal(1))
   if(Status.ne.ISDC_OK)then
      Status = ISDC_OK  
      MapParDec(1) = 0
      MapParVal(1) = 0. 
   else
      MapParDec(1) = 1
   endif

if(ProjType==1)then
  
   parName = 'MapCentre2'
   Status = pilgetreal4(parName,MapParVal(2))
   if(Status.ne.ISDC_OK)then
      Status = ISDC_OK  
      MapParDec(2) = 0
      MapParVal(2) = 0. 
   else
      MapParDec(2) = 1
   endif 
else
   MapParDec(1:2)=1
endif
parName = 'MapSize'
 Status = pilgetreal4(parName,MapParVal(3))
   if(Status.ne.ISDC_OK)then
      Status = ISDC_OK  
      MapParDec(3) = 0
      MapParVal(3) = 0. 
   else
      MapParDec(3) = 1
   endif




!map size verification
if(ProjType==1)then
   i = MapParDec(1)+MapParDec(2)
   select case(i)
   case(0,2)
   case(1)
      str250 = 'Incoherence in mosaicked map centre definition ==> will be calculated '
      call WAR_MESSAGE(procName,str250,0,Status)
      MapParDec(1:2) = 0
   case default
      str250 = 'Error in mosaicked map centre definition ==> will be calculated '
      call WAR_MESSAGE(procName,str250,0,Status)
      MapParDec(1:2) = 0
   end select

endif

if( MapParDec(3)==1)then
   if (MapParVal(3) .le.0)  then
       str250 = &
'Error in mosaicked map size definition ==> will be calculated by MIMOSA'
       call WAR_MESSAGE(procName,str250,0,Status)
       MapParDec(3) = 0
    endif
 endif


!!$!map pixel size verification
!!$if( MapParDec(4)==1)then
!!$   if (MapParVal(4) .le.0)  then
!!$       str250 = &
!!$'Error in mosaicked map pixel size definition ==> will be calculated by MIMOSA'
!!$       call WAR_MESSAGE(procName,str250,0,Status)
!!$       MapParDec(4) = 0
!!$    endif
!!$ endif

parName = 'MinCatSouSnr'
Status = pilgetreal4(parName,PeakMinHeight)
if(Status.ne.ISDC_OK)then
   Status = ISDC_OK  
   PeakMinHeight = DefPeakMinHeight
else
   if(PeakMinHeight < 0.)then
      PeakMinHeight = DefPeakMinHeight
      write(str250,*) ' Uncorrect ',parname,' set to default'
      call WAR_MESSAGE(procName,str250,0,Status)
   endif
endif

parName = 'MinNewSouSnr'
Status = pilgetreal4(parName,NonIdentPeakMinHeight)
if(Status.ne.ISDC_OK)then
   Status = ISDC_OK  
   NonIdentPeakMinHeight = DefNonIdentPeakMinHeight
else
   if(NonIdentPeakMinHeight < 0.)then
      NonIdentPeakMinHeight = DefNonIdentPeakMinHeight
      write(str250,*) ' Uncorrect ',parname,' set to default'
      call WAR_MESSAGE(procName,str250,0,Status)
   endif
endif


!===========================
END SUBROUTINE ReadParFile
!===========================


!============================================================
SUBROUTINE ReadOffCorrIndex(IdxPtr,BinNumber,&
     Bins,Status)
!==============================================================
USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE

IMPLICIT NONE

!INPUT/OUTPUT VARIABLES

REAL(kind=8),dimension(:,:),pointer::Bins
INTEGER :: IdxPtr, BinNumber,Status
!LOCAL VARIABLES
INTEGER :: iok,n
REAL(kind=8),dimension(:),pointer :: tab1,tab2
CHARACTER(len=20) :: procname
procname = 'ReadOffCorrIndex'
Status = ISDC_OK

Status = dal_object_open(OffCorrFile,IdxPtr,Status)
if(Status.ne.ISDC_OK)then
  call war_MESSAGE(procName,'cannot open off-axis correction file '//OffCorrFile,&
       0,Status)
  Status = OffCorrError
  return
endif
Status = dal3gen_index_get_num_members(IdxPtr,OffCorrStrucName,&
          BinNumber,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call war_MESSAGE(procName,&
        'dal3gen_index_get_num_members problem  for '//OffCorrIdxStrucName,&
        0,Status)
   Status = IsdcExitStatus
   return
endif
if(BinNumber==0)then
    call war_MESSAGE(procName,&
        '0 elements in  '//OffCorrIdxStrucName,&
        0,Status)
    Status = OffCorrError
   return
endif

allocate(Bins(BinNumber,2),tab1(BinNumber),tab2(BinNumber),stat=iok)
if(iok.ne.0)then
   call war_Message(procName,'Allocation problem for '//OffCorrIdxStrucName,&
            0,Status)
   Status = AllocError
   return
endif

Status = DAL_TABLE_GET_COL(IdxPtr,'E_MIN',&
        1,DAL_DOUBLE,BinNumber,addrof(tab1(1)),Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
      call war_MESSAGE(procName,&
           ' DAL_TABLE_GET_COL problem for E_MIN in '//OffCorrIdxStrucName ,&
           0,Status)
      Status = IsdcExitStatus
      return
   endif

Status = DAL_TABLE_GET_COL(IdxPtr,'E_MAX',&
        1,DAL_DOUBLE,BinNumber,addrof(tab2(1)),Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
      call war_MESSAGE(procName,&
           ' DAL_TABLE_GET_COL problem for E_MAX in '//OffCorrIdxStrucName  ,&
           0,Status)
      Status = IsdcExitStatus
      return
   endif
do n=1,BinNumber
   if((tab1(n).lt.0.).or.(tab2(n).lt.0.))then
      call war_MESSAGE(procname,&
           'Energy band error 1 in '//OffCorrIdxStrucName,0,Status)
      Status = OffCorrError
      return
   else
      if(tab1(n).ge.tab2(n))then
         call MESSAGE(procname,' Bands error 2 in '//OffCorrIdxStrucName,&
         0,Status)
         Status = OffCorrError
         return
      endif
   endif
enddo
Bins(:,1) = tab1
Bins(:,2) = tab2
deallocate(tab1,tab2)
!==============================
END SUBROUTINE ReadOffCorrIndex
!==============================

!===================================================
SUBROUTINE InitArrays(iSky,jSky,iEdge,jEdge,&
            Shd,ShdVar,ShdEffi,&
            GroupShd,GroupShdVar,GroupShdEffi,&
            Wei,Sky,SkyCle,SkyResid,SkySignif,SkyVar,&
            Status)
!===================================================


USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_IMAGING_PARAM
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)::iSky,jSky

!OUTPUT VARIABLES
INTEGER,Intent(out)::Status,iEdge,jEdge
REAL(kind=4),dimension(:,:),pointer :: Shd,ShdVar,ShdEffi
REAL(kind=4),dimension(:,:),pointer ::GroupShd,GroupShdVar,GroupShdEffi
REAL(kind=4),dimension(:,:),pointer :: Sky,SkyCle,Wei
REAL(kind=4),dimension(:,:),pointer :: SkyResid,SkySignif,SkyVar


!LOCAL  VARIABLES
INTEGER::iok,i,j,m,ii1,ii2,jj1,jj2,i1,i2,j1,j2
INTEGER::iEdge1,jEdge1,isqua,ijedge,ijedge1
REAL      (kind=4) :: skyCentreI,skyCentreJ,xd,yd
CHARACTER(len=20)::procName

Status= 0 
procName = 'InitArrays'
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)

!allocate sky arrays
ALLOCATE(Shd(1:iDetDim,1:jDetDim),DeadZonePattern(1:iDetDim,1:jDetDim),&
     ShdVar(1:iDetDim,1:jDetDim),ScwExpo(iSky,jSky),&
     ShdEffi(1:iDetDim,1:jDetDim),nanfilter(iSky,jSky),&
     GroupShd(1:iDetDim,1:jDetDim),&
     GroupShdVar(1:iDetDim,1:jDetDim),&
     GroupShdEffi(1:iDetDim,1:jDetDim),&
     Wei(1:iDetDim,1:jDetDim),&
     FilterSky(iSky,jSky),CurEffArea(iSky,jSky),&
     Sky(iSky,jSky),ThickCorr(iSky,jSky),OffAxisCorr(iSky,jSky),&
     SkyCle(iSky,jSky),SkyResid(iSky,jSky),&
      DetActivePixels(1:iDetDim,1:jDetDim),&
      SkySignif(iSky,jSky),SkyVar(iSky,jSky),ScwNoSou(iSky,jSky),&
      stat=iok)

if(iok.ne.0)then
   call MESSAGE(procName,'Allocation problem ',AllocError,Status)
   return
endif 



CurEffArea(:,:) = 1.
OffAxisCorr(:,:) = 1.
DeadZonePattern(:,:)=.true.
Shd(:,:)=0.
ShdVar(:,:)=0.
ShdEffi(:,:)=0.
GroupShd(:,:)=0.
GroupShdVar(:,:)=0.
GroupShdEffi(:,:)=0.
Wei(:,:)=1.
FilterSky(:,:) = .true.
Sky(:,:)=0.
SkyCle(:,:)=0.
SkyResid(:,:)=0.
SkySignif(:,:)=0.
SkyVar(:,:)=0.       
ThickCorr(:,:)=0.
DetActivePixels(:,:) = .true.
ScwExpo(:,:)=0.

if(ExpoMap)then
   allocate(OrgEfficiency(1:iDetDim,1:jDetDim),stat=iok)
   if(iok.ne.0)then
      call MESSAGE(procName,'Allocation problem ',AllocError,Status)
      return
   endif
   ScwExpo(:,:)=0.
   OrgEfficiency(:,:) = 0.
endif

DetActivePixels(65:66,:) = .false.
DetActivePixels(:,33:34) = .false.
DetActivePixels(:,67:68) = .false.
DetActivePixels(:,101:102) = .false.



!DEAD ZONES
 If(DeadZonesFilling.lt.2)then
    where(.not.DetActivePixels)
       DeadZonePattern = .false.
       WEI = 0.
    endwhere
endif


  


! Here we include the weighing to take into account 
!the not exact mask b-p det dim 
!(NOT implemented for idet=3 !) because wei is used
! to put DZ to 0. in det so wei  for decon == wei
! of det effic.

If (DetType == 3) THEN
  WEI(1,:) = 0.5
  WEI(iDetDim,:) =0.5
  WEI(:,1:5) = 0.5
  WEI(:,jDetdim-5+1:jDetdim) = 0.5
ENDIF


! Temporary set to 0. of WEI because PH image have last! raw and column = 0.
!         IF (ISim == 0) THEN
!            WEI(idim,:)=0.
!            WEI(:,jdim)=0.
!         ENDIF



! IMasPixDim+IDetDim : FOV size 
iEdge=(iSky-(IMasPixDim+IDetDim))/2 +1  ;iEdge1=iEdge+1
jEdge=(jSky-(JMasPixDim+JDetDim))/2 +1 ;jEdge1=jEdge+1



CodedMaxDist = (IMasPixDim+IDetDim)/2-2

isqua=4  ! good for IBIS/ISGRI

IF (iEdge > 0) THEN

!!$   !old fashion - up to 3.9
!!$   ! Cut lateral bands
!!$  FilterSky(:,1:jEdge)=.false. 
!!$  FilterSky(:,jSky-jEdge:jSky)=.false.
!!$  FilterSky(1:iEdge,:)=.false.
!!$  FilterSky(iSky-iEdge:iSky,:)=.false.
!!$
!!$  ! Cut indefined "squares" 
!!$  i1 = iEdge1
!!$  i2 = iEdge1+isqua
!!$  j1 = jEdge1
!!$  j2 = jEdge1+isqua
!!$  FilterSky(i1:i2,j1:j2)=.false.
!!$
!!$  i1 = iEdge1
!!$  i2 = iEdge1+isqua
!!$  j1 = jSky-jEdge1-isqua
!!$  j2 = jSky-jEdge1
!!$  FilterSky(i1:i2,j1:j2)=.false.
!!$
!!$  i1 = iSky-iEdge1-isqua
!!$  i2 = iSky-iEdge1
!!$  j1 = jSky-jEdge1-isqua
!!$  j2 = jSky-jEdge1
!!$  FilterSky(i1:i2, j1:j2) = .false.
!!$
!!$
!!$  i1 = iSky-iEdge1-isqua
!!$  i2 = iSky-iEdge1
!!$  j1 = jEdge1
!!$  j2 = jEdge1+isqua
!!$  FilterSky(i1:i2,j1:j2)=.false.

! enhanced FOV but not corresponding to old images
! Cut lateral bands
!!$  FilterSky(:,1:jEdge)=.false. 
!!$  FilterSky(:,jSky-jEdge+1:jSky)=.false.
!!$  FilterSky(1:iEdge,:)=.false.
!!$  FilterSky(iSky-iEdge+1:iSky,:)=.false.
!!$
!!$  ! Cut indefined "squares" 
!!$  i1 = iEdge1
!!$  i2 = iEdge1+isqua
!!$  j1 = jEdge1
!!$  j2 = jEdge1+isqua
!!$  FilterSky(i1:i2,j1:j2)=.false.
!!$
!!$  i1 = iEdge1
!!$  i2 = iEdge1+isqua
!!$  j1 = jSky-jEdge1-isqua+1
!!$  j2 = jSky-jEdge1+1
!!$  FilterSky(i1:i2,j1:j2)=.false.
!!$
!!$  i1 = iSky-iEdge1-isqua+1
!!$  i2 = iSky-iEdge1+1
!!$  j1 = jSky-jEdge1-isqua+1
!!$  j2 = jSky-jEdge1+1
!!$  FilterSky(i1:i2, j1:j2) = .false.
!!$
!!$  i1 = iSky-iEdge1-isqua+1
!!$  i2 = iSky-iEdge1+1
!!$  j1 = jEdge1
!!$  j2 = jEdge1+isqua
!!$  FilterSky(i1:i2,j1:j2)=.false.

! correct but decreased FOV 
! horizonthal bands 21 pix cuted
! vertical bands 19 pixel cuted
 ! Cut lateral bands
  FilterSky(:,1:jEdge+1)=.false. 
  FilterSky(:,jSky-jEdge:jSky)=.false.
  FilterSky(1:iEdge+1,:)=.false.
  FilterSky(iSky-iEdge:iSky,:)=.false.

  ! Cut indefined "squares" 
  i1 = iEdge1+1
  i2 = iEdge1+isqua+1
  j1 = jEdge1+1
  j2 = jEdge1+isqua+1
  FilterSky(i1:i2,j1:j2)=.false.

  i1 = iEdge1+1
  i2 = iEdge1+isqua+1
  j1 = jSky-jEdge1-isqua
  j2 = jSky-jEdge1
  FilterSky(i1:i2,j1:j2)=.false.

  i1 = iSky-iEdge1-isqua
  i2 = iSky-iEdge1
  j1 = jSky-jEdge1-isqua
  j2 = jSky-jEdge1
  FilterSky(i1:i2, j1:j2) = .false.

  i1 = iSky-iEdge1-isqua
  i2 = iSky-iEdge1
  j1 = jEdge1+1
  j2 = jEdge1+isqua+1
  FilterSky(i1:i2,j1:j2)=.false.
ENDIF


allocate( AllSouList(ScwNumber,OutEnergyNumber),stat=iok)

if(iok.ne.0)then
   call MESSAGE(procName,'Allocation problem ',AllocError,Status)
   return
endif 

do i=1,ScwNumber
   do j=1,OutEnergyNumber
      AllSouList(i,j)%souNumber = 0   
   enddo
enddo

do m=1,8
   select case(m)
   case(1)
      ii1 = 1; ii2 = 64; jj1 = 1; jj2 = 32
   case(2)
      ii1 = 67; ii2 = 130; jj1 = 1; jj2 = 32
   case(3)
      ii1 = 1; ii2 = 64; jj1 = 35; jj2 = 66
   case(4)
      ii1 = 67; ii2 = 130;  jj1 = 35; jj2 = 66
   case(5)
      ii1 = 1; ii2 = 64; jj1 = 69; jj2 = 100
   case(6)
      ii1 = 67; ii2 = 130; jj1 = 69; jj2 = 100
   case(7)
      ii1 = 1; ii2 = 64; jj1 = 103; jj2 = 134
   case(8)
      ii1 = 67; ii2 = 130; jj1 = 103; jj2 = 134
   end select

   ModulesLimits(m,1) = ii1
   ModulesLimits(m,2) = ii2
   ModulesLimits(m,3) = jj1
   ModulesLimits(m,4) = jj2
enddo


!================================
END SUBROUTINE InitArrays
!================================


!====================================================
SUBROUTINE CleanStaInit(RaX,DecX,RaZ,DecZ,AttiNumPoi,nclean,Status)
!====================================================

USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_IMAGING_PARAM
USE ATTI_DECLARATIONS
USE ATTI_DEFS
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)           ::AttiNumPoi,nclean
REAL(kind=8), dimension(:), pointer ::RaX,DecX,RaZ,DecZ
!OUTPUT VARIABLES
INTEGER,Intent(out)          ::Status


!LOCAL  VARIABLES
INTEGER::i,iok,j,k,l,ib,infixsou
REAL(kind=4),dimension(:), pointer :: distance,pomx,pomy
INTEGER,dimension(:), pointer :: order
REAL (kind=4) :: a,d,yt,zt,flux,val
REAL(kind=8)::x,y,G05CAF
CHARACTER(len=100)::filename,title
CHARACTER(len=20)::procName

procName = 'CleanStaInit'
Status = 0
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)

nanf =  dal3aux_nan()


IOView =  IOSphere ! definition for attitude library

! BKGMOYPAR
ALLOCATE(BkgParTab(1:10),stat=iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif
BkgParTab(:)=0.
BkgParTab(1)=BkgMoyPar*FLOAT(iDetdim*jDetdim)     

! SOURCES
ALLOCATE(InSourceList(10,SourceNumber),stat=iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif
InSourceList(:,:) = 0.
 
! Cleaning param
ALLOCATE(CleanParTab(1:10),stat=iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif
CleanParTab(:) = 0.
CleanParTab(1)= 1.
CleanParTab(2)= REAL(nclean) 
CleanParTab(3)= MaElPixDim  ! source radius = 2.435
CleanParTab(4)= PeakMinHeight  ! Minimum peak height
CleanParTab(5)= NonIdentPeakMinHeight
! radius for source raw and fine identification
CleanParTab(6)= SouSearchRadius
CleanParTab(10)= REAL(DebugMode)    ! test parameter


ALLOCATE(ImaStatPar(1:10),stat=iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif
ImaStatPar(:) = 0.
NumSouDes = 4+2*OutEnergyNumber
ALLOCATE(OutSourceCat(1:SourceNumber,NumSouDes),&
         OutSourceCatFlux(1:SourceNumber,OutEnergyNumber),&
         OutSourceCatFluxErr(1:SourceNumber,OutEnergyNumber),&
         OutSourceCatSig(1:SourceNumber),&
         OutSourceCatLocErr(1:SourceNumber),stat=iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',AllocError,Status)
   return
endif
OutSourceCat(:,:) = 0.
OutSourceCatFlux(:,:) = 0.
OutSourceCatFluxErr(:,:) = 0.
OutSourceCatSig(:) = 0.
OutSourceCatLocErr(:) = 0.
if(associated(distance))deallocate(distance)
if(associated(pomx))deallocate(pomx)
if(associated(pomy))deallocate(pomy)
if(associated(order))deallocate(order)


!====================================================
END SUBROUTINE CleanStaInit
!====================================================

!===========================================================
SUBROUTINE AttiInitialise(MosaNumber,raX,decX,raZ,decZ,PosAngle,&
                          AttiNumPoi,iMapSize,jMapSize, Status)
!===========================================================

USE ISDC
USE DAL3AUX_F90_API
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
USE INIT_PROC
IMPLICIT NONE

!INPUT VARIABLES
REAL(kind=8), dimension(:), pointer ::RaX,DecX,RaZ,DecZ,PosAngle
INTEGER,Intent(in)  ::MosaNumber 
!OUTPUT VARIABLES
INTEGER,Intent(out)  :: AttiNumPoi
INTEGER,Intent(out)  :: iMapSize,jMapSize
INTEGER,Intent(out)  :: Status

!LOCAL VARIABLES
INTEGER  :: iok,i,imax
REAL(kind=8) :: projangle
LOGICAL :: badMapSize 
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'AttiInitialise'
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)

if (abs(IOSphere).ne.1)then
   call WAR_Message(procName,&
    'Uncorrect Y axis orientation, set to the telescope frame ',0,Status)
   IOSPhere = 1
endif


IOView=IOSphere

if(mapParDec(4) == 1) then
   !user defined Map pixel size
   MagnifFactor = mapParVal(4)
   PointMagnifFactor =MagnifFactor
   write(str250,*,err=220)' Map pixel size set by user to ',MagnifFactor*5.1,' arcmin'
   goto 221
220 str250 = ' Map pixel size re-set by user'
221 call message(procname,str250,ZeroError,Status)
endif

  
attiNumPoi = MosaNumber+1


call CalcRefFrame(attiNumPoi,raX,decX,raZ,decZ)


call CalcLBRefFrame(AttiNumPoi,RaX,DecX,RaZ,DecZ)

Call ATTI_INIT(attiNumPoi,raX,decX,raZ,decZ,posAngle,status=Status)
if(err.ne.0)then
   if(err.eq.30)then
      call message(procName,'allocation problem in attitude routine ',&
           AttiError,Status)
      return
   else
      call message(procName,'Bad Attitude of  X,Z axes ',&
           AttiError,Status)
      return
   endif
endif

if(ProjType==0)then
  projangle = pointing_table(attiNumPoi,3)
  if(projangle.ne.0.0)then
     call message(procName,'No rotation angle allowed in the final projection ',&
          AttiError,Status)
     return
  endif
endif

!calculation of sky map 
!it should contain all the Scw images

 CALL CalcMapSize(MosaNumber,1,MagnifFactor,iMapSize,jMapSize,1,Status)
 If(Status.ne.ISDC_OK)return



!================================================
END SUBROUTINE AttiInitialise
!================================================

!===========================
SUBROUTINE DESCRIPTION(mosanumber,Rax,Decx,imapsize,jmapsize)
!===========================

USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
INTEGER :: mosanumber,imapsize,jmapsize
REAL(KIND=8),DIMENSION(:),POINTER::raX,decX
!LOCAL VARIABLES
INTEGER :: status,i
CHARACTER(len=80) :: text
CHARACTER(len=2) :: str
CHARACTER(len=4) :: nstr
CHARACTER(len=10) :: estr1,estr2
!CHARACTER(len=80) :: str80

str = '>>'
status = 0
text = 'CURRENT CONFIGURATION :'
call Message(str,text,ZeroError,Status)
call Message(str,'...............................',ZeroError,Status)

write(str250,*,err=156)'Number of Mosaicks found           :',MosaNumber
goto 157
156 str250 = errorstr//' 156'
157 call MESSAGE(str,str250,ZeroError,Status) 
write(str250,*,err= 158)'Number of Energy bands  :',outEnergyNumber
goto 159
158 str250 = errorstr//' 158'
159 call MESSAGE(str,str250,ZeroError,Status)
do i=1,OutEnergyNumber
  write(estr1,'(F7.1)',err = 160)outEnergyBands(i,1)
  goto 161
160 estr1 = ' '
161  write(estr2,'(F7.1)',err = 162)outEnergyBands(i,2)
  goto 163
162  estr2 = ' '
163 str250 = '   ['//estr1//','//estr2//'] '
  call MESSAGE(str,str250,ZeroError,Status) 
enddo

select case(ProjType)
case(1)
   if(EquaGal==0)then
      text = ' RA/DEC - TAN projection '
   else
      text = ' GLON/GLAT - TAN projection '
   endif
case(0)
  if(EquaGal==0)then
      text = ' RA/DEC - CAR projection '
   else
      text = ' GLON/GLAT - CAR projection '
   endif
end select
call Message(str,text,ZeroError,Status)

if(RebinType==0)then
   text =  ' No flux spread'
else
   text =  ' Flux spread'
endif
call Message(str,text,ZeroError,Status)

call Message(str,' Final mosaics centre :',ZeroError,Status)
if(EquaGal==0)then
   write(text,'(a,f6.2,a,f6.2,a)')'       alpha :',rax(mosanumber+1),' ,delta :',decx(mosanumber+1),'[deg]'
else
    write(text,'(a,f6.2,a,f6.2,a)')'       l :',lx(mosanumber+1),' ,b :',bx(mosanumber+1),'[deg]'
endif
call Message(str,text,ZeroError,Status)
call Message(str,'               size  :',ZeroError,Status)
write(text,'(a,I5,a,I5,a)')'         ',imapsize,' x ',jmapsize,' pixels'
call Message(str,text,ZeroError,Status)

if(FitMode)then
   if(FitType==0)then
      text = ' Fitting procedures on : PSF fit'
   else
      text = ' Fitting procedures on : Gaussian fit'
   endif
else
  text = ' Fitting procedures off'
endif
call Message(str,text,ZeroError,Status) 

if(FitMode)then
  if(FitPosFixed)then
     text = ' All source positions fixed in Scw fitting'
     call Message(str,text,ZeroError,Status) 
  else
     if(InFixSouNumber > 0)then
        if(InFixSouNumber .le. SourceNumber)then
           text = ' Some cat. sources positions will be fixed in Scw fitting'
           call Message(str,text,ZeroError,Status) 
        else
           text = ' All cat. sources positions will be fixed in Scw fitting' 
           call Message(str,text,ZeroError,Status) 
        endif
     endif
  endif


endif


   write(nstr,'(I3)',err=152)abs(nToClean)
   goto 153
152 nstr = ' '
153 str250 = ' Search and fit of  '//nstr//' sources in final mosaics'
call Message(str,str250,ZeroError,Status)




call Message(str,'...............................',ZeroError,Status)

!===========================
END SUBROUTINE DESCRIPTION
!===========================


