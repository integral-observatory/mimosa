

!*************************************
MODULE MIMOSA_BASE3_MODULE
!*************************************
INTERFACE


SUBROUTINE ReplaceElement(GrpPtr,ElemPtr,StrucName,&
                       OutElemFile,OutName,clearFlag,Status)
!-----------------------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER        :: GrpPtr,ElemPtr
CHARACTER(len=*)       :: StrucName,OutElemFile,OutName
INTEGER               :: Status
LOGICAL :: clearFlag
END SUBROUTINE ReplaceElement

SUBROUTINE ReplaceIdx(GrpPtr,ElemPtr,StrucName,StrucNameInt,&
                       OutElemFile,OutName,clearflag,Status)
!-----------------------------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER        :: GrpPtr,ElemPtr
CHARACTER(len=*)       :: StrucName,StrucNameInt,OutElemFile,OutName
LOGICAL :: clearFlag
INTEGER               :: Status
END SUBROUTINE ReplaceIdx

SUBROUTINE EnergyBandInfo(GrpPtr,EBins,CHans,Status)
!----------------------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: EBins
INTEGER,dimension(:,:),pointer ::  Chans
INTEGER         :: GrpPtr,Status
END SUBROUTINE EnergyBandInfo


SUBROUTINE FixSouPosInCat(Status)
!-----------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER         :: Status
END SUBROUTINE FixSouPosInCat

SUBROUTINE ReadAttributesAndAtti(imos,GrpPtr,RaX,DecX,RaZ,DecZ,posAngle,&
     incr,Status)
!---------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER         :: imos,GrpPtr
! pointer to output source catalogue
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle
REAL(kind=8)  ::incr
INTEGER         :: Status
END SUBROUTINE  ReadAttributesAndAtti

SUBROUTINE WriteRecallInfo(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle)
!------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: MosaNumber
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle
END SUBROUTINE   WriteRecallInfo


SUBROUTINE ReadRecallInfo(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle)
!----------------------------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: MosaNumber
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle
END SUBROUTINE ReadRecallInfo

END INTERFACE 

!*************************************
END MODULE MIMOSA_BASE3_MODULE
!*************************************


!*************************************
MODULE MIMOSA_BASE2_MODULE
!*************************************
INTERFACE
SUBROUTINE  WriteAttiHeader(imatype,num,nx,ny,ImaPtr,Status)
!..........................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  :: imatype ! 0 - scw, 1 - map
INTEGER  :: num    !pointing number
INTEGER  :: nx,ny  ! image dimensions
INTEGER  :: ImaPtr !image pointer
!OUTPUT VARIABLES
INTEGER  :: Status
END SUBROUTINE  WriteAttiHeader

SUBROUTINE WriteAttributes(ObjPtr,MemberPtr,scw,OutBand,&
     ImaFlag,ImaType,ImaTypeVal,&
           Status)
!------------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   :: ObjPtr,MemberPtr,scw       ! pointer to index member
INTEGER   :: OutBand         ! out energy band number 
INTEGER           :: ImaFlag
! 0 - groupped shd
! 1 - image
! 2 - carte
character(len=*)  :: imaType,imaTypeVal
!OUTPUT VARIABLES
INTEGER   :: Status
END SUBROUTINE  WriteAttributes

SUBROUTINE VerifBins(nrows,tab1,tab2,tabind,bins,bnum,str,Status)
!---------------------------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER  :: nrows,bnum,iok
CHARACTER(len = *) :: str
REAL(kind=4),dimension(:),pointer :: tab1,tab2,tabind
REAL(kind=4),dimension(:,:),pointer ::bins
INTEGER  :: Status
END SUBROUTINE VerifBins

END INTERFACE

!*************************************
END MODULE MIMOSA_BASE2_MODULE
!*************************************

!*********************************
MODULE MIMOSA_BASE1_MODULE
!*********************************


INTERFACE

SUBROUTINE DECO(scwcleaimage)
!_____________________________________________
IMPLICIT NONE

!INPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: scwcleaimage
END SUBROUTINE deco

subroutine RealSize(scwcleaimage,imin,imax,jmin,jmax)
!------------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
REAL,dimension(:,:),pointer :: scwcleaimage
integer ::imin,imax,jmin,jmax
END SUBROUTINE RealSize

subroutine PutCol(Ptr,nsou,tab,colname,colnr,&
                  SourceList,Status)
!........................................................
IMPLICIT NONE
!INPUT VARIABLES
INTEGER                   :: Ptr,nsou,colnr
REAL(kind=4),dimension(:),pointer :: tab
REAL(kind=4),dimension(:,:),pointer :: SourceList
CHARACTER(len=*)                  :: colname
!OUTPUT VARIABLES
INTEGER             :: Status
end subroutine PutCol

subroutine UpdateOutCat(type,Scw,OutBand,nsou,SourceList,&
                        alpha,delta,flux,fluxerr,snr,locerr, Status)
!........................................
IMPLICIT NONE
!INPUT VARIABLES
INTEGER             ::type,Scw,OutBand, nsou
REAL(kind=4),dimension(:,:),pointer :: SourceList
REAL(kind=4),dimension(:),pointer::  alpha,delta,flux,fluxerr,snr,locerr

!OUTPUT VARIABLES
INTEGER             :: Status
end subroutine UpdateOutCat

SUBROUTINE NewIndexMember&
          (IdxPtr,OutBand,StrucName,MemberPtr,Status)
!........................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER             :: IdxPtr,OutBand
CHARACTER(len=*)            :: StrucName
!OUTPUT VARIABLES
INTEGER             :: MemberPtr
INTEGER             :: Status
END SUBROUTINE NewIndexMember



SUBROUTINE AddIndexMemberArray(ObjPtr,IdxPtr,scw,OutBand,ImaFlag,&
           ArrName,ArrType,ArrTypeVal,Arr,MemberPtr,Status)
!........................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER             :: ObjPtr,IdxPtr,scw,OutBand
!OutBand - number of output energy band
INTEGER                     :: ImaFlag
! if yes then write attitude info into the header
! and some additional keywords
CHARACTER(len=*)            :: ArrName,ArrType,ArrTypeVal
REAL(kind=4),dimension(:,:),pointer :: Arr 
!OUTPUT VARIABLES
INTEGER             :: MemberPtr
INTEGER             :: Status
END SUBROUTINE AddIndexMemberArray

SUBROUTINE CreateClear(ParentPtr,ElemType,ClearFlag,&
                       StrucNameExt,StrucNameInt,FitsName,ElemPtr,&
                       Status)
!........................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
INTEGER  ,Intent(in)   :: ElemType
! if ElemType == 0 then element else index
INTEGER  ,Intent(in)   :: ClearFlag
! 0 - no clearing of existing element ( or index)
! 1 - clears rows(members) of element(index) if exist
character(len=*),Intent(in)   :: StrucNameExt
! for ordinary element its structure name
! for index            its index structure name
character(len=*),Intent(in)   :: StrucNameInt 
! for ordinary element not used
! for index            its member structure name
character(len=*),Intent(in)   :: FitsName  
! name of fits file of the element
!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
END SUBROUTINE CreateClear

SUBROUTINE GetClearElement(ParentPtr,StrucName,ElemPtr,&
                       Status)
!........................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
character(len=*),Intent(in)   :: StrucName 

!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
END SUBROUTINE GetClearElement

SUBROUTINE ClearDetachElement(ParentPtr,StrucName,ElemPtr,&
                       Status)
!........................................................

IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
character(len=*),Intent(in)   :: StrucName 

!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
END SUBROUTINE ClearDetachElement

SUBROUTINE CreateClearIndex(ParentPtr,StrucNameExt,StrucNameInt,&
                            FitsName,ElemPtr,Status)
!........................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
character(len=*),Intent(in)   :: StrucNameExt
! for ordinary element its structure name
! for index            its index structure name
character(len=*),Intent(in)   :: StrucNameInt 
! for ordinary element not used
! for index            its member structure name
character(len=*),Intent(in)   :: FitsName  
! name of fits file of the element
!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
END SUBROUTINE CreateClearIndex

SUBROUTINE GetClearIndex(ParentPtr,StrucNameExt,StrucNameInt,ElemPtr,Status)
!........................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
character(len=*),Intent(in)   :: StrucNameExt
! index structure name
character(len=*),Intent(in)   :: StrucNameInt
! index member structure name
!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
END SUBROUTINE GetClearIndex

SUBROUTINE OpenIndexMember(IdxPtr,&
           MemberName,MemberNum,MemberPtr,&
           nFound,Status)
!...........................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   :: IdxPtr       ! pointer to index
CHARACTER(len=*)  :: MemberName   ! MEMBER NAME
INTEGER   :: MemberNum    ! member row number

!OUTPUT VARIABLES
INTEGER   :: MemberPtr !pointer to member
INTEGER   :: nFound
INTEGER   :: Status
END SUBROUTINE OpenIndexMember


SUBROUTINE ReadArray(ArrPtr,iDim,jDim,Arr,Status)
!................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in) :: ArrPtr,iDim,jDim
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: Arr
INTEGER            :: Status
END SUBROUTINE ReadArray

SUBROUTINE ReadUnknownArray(ArrPtr,iDim,jDim,Arr,Status)
!................................................
IMPLICIT NONE

!INPUT VARIABLES
INTEGER :: ArrPtr,iDim,jDim
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: Arr
INTEGER            :: Status
END SUBROUTINE ReadUnknownArray

SUBROUTINE ReadShdAttributes(MemberNum,MemberPtr,eMinVal,&
            eMaxVal,Status)
!-------------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   :: MemberNum,MemberPtr  !index member no, pointer to
!OUTPUT VARIABLES
REAL(kind=4)      :: eMinVal,EMaxVal
INTEGER   :: Status
END SUBROUTINE ReadShdAttributes

SUBROUTINE ReadImaAttributes(MemberPtr,eMinVal,&
            eMaxVal,Status)
!------------------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   ::MemberPtr       ! pointer to index member
!OUTPUT VARIABLES
REAL(kind=4)      :: eMinVal,EMaxVal
INTEGER   :: Status
END SUBROUTINE  ReadImaAttributes

SUBROUTINE ReadScwAttributes(ScwPtr,Status)
!------------------------------------------
IMPLICIT NONE
!INPUT VARIABLES
INTEGER   :: ScwPtr       ! pointer to scw
INTEGER   :: Status
END SUBROUTINE   ReadScwAttributes

SUBROUTINE ReadIndexInfo(ShdIdxPtr,EBins,Chans,Status)
!----------------------------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER                       :: ShdIdxPtr
REAL(kind=4),dimension(:,:),pointer     :: EBins
INTEGER ,dimension(:,:),pointer :: Chans
INTEGER                       :: Status
END SUBROUTINE ReadIndexInfo


SUBROUTINE CreateElem(OutElemFile,StrucName,ElemPtr,Status)
!---------------------------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER        :: ElemPtr
CHARACTER(len=*)       :: StrucName,OutElemFile
INTEGER               :: Status
END SUBROUTINE CreateElem



END INTERFACE
!*************************************
END MODULE MIMOSA_BASE1_MODULE
!*************************************

!*************************************
MODULE MIMOSA_BASE4_MODULE
!*************************************
INTERFACE

SUBROUTINE Read1Mosa(MosaNumber,EBins,Chans,Status)
!---------------------------------------------------
IMPLICIT NONE 
!input/output variables
INTEGER :: MosaNumber,Status
REAL(kind=4),dimension(:,:),pointer :: EBins
INTEGER,dimension(:,:),pointer ::  Chans
END SUBROUTINE Read1Mosa

SUBROUTINE Read1ReCallMosa(MosaNumber,EBins,Chans,Status)
!----------------------------------------------------------
IMPLICIT NONE 
!input/output variables
INTEGER :: MosaNumber,Status
REAL(kind=4),dimension(:,:),pointer :: EBins
INTEGER,dimension(:,:),pointer ::  Chans
END SUBROUTINE Read1ReCallMosa

SUBROUTINE Read2Mosa(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle,incr,Status)
!-------------------------------------------------------------------
IMPLICIT NONE 

!input/output variables
INTEGER :: MosaNumber,Status
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle
REAL(kind=8)  ::incr
END SUBROUTINE Read2Mosa

END INTERFACE

!*********************************
END MODULE MIMOSA_BASE4_MODULE
!*********************************

!*********************************
MODULE MIMOSA_BASE_MODULE
!*********************************

INTERFACE


SUBROUTINE OgOpen(MosaNumber,MosaImaIdxPtr,MosaSouIdxPtr,&
                  OutCatPtr,RaX,DecX,RaZ,DecZ,posAngle,Status)
!------------------------------------------------------
IMPLICIT NONE 
INTEGER :: MosaNumber
!OUTPUT VARIABLES

! pointer to OG
INTEGER,Intent(out) :: MosaImaIdxPtr
! pointer to mosaicked image list
INTEGER,Intent(out) :: MosaSouIdxPtr
! pointer to mosaicked image source list
INTEGER,Intent(out) :: OutCatPtr
! pointer to output source catalogue
REAL(kind=8), dimension(:), pointer ::RaX,DecX,RaZ,DecZ,posAngle
INTEGER,Intent(out) ::Status

END SUBROUTINE OgOpen




SUBROUTINE WriteImages(type,scw,OutBand,ObjPtr,IdxPtr,CleanImage,&
                        VarImage,ResidImage,SignifImage,&
                        Status)
!-------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)   :: type
! 0 - Scw images
! 1 - Map images
INTEGER,Intent(in)   :: scw
!pointing number
INTEGER,Intent(in)   :: OutBand
!out energy band number
INTEGER   ::ObjPtr,IdxPtr
!pointer to the host object and index of Scw image index
REAL(kind=4),dimension(:,:),pointer  ::CleanImage, VarImage,&
                               ResidImage,SignifImage 

!OUTPUT VARIABLES
INTEGER   :: Status
END SUBROUTINE WriteImages



 SUBROUTINE ReadMosaBandImage(mosa,outBand,iDim,jDim,&
                             MosaImage,MosaVar, MosaExpo,Status)
!----------------------------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES

INTEGER ,Intent(in)  :: mosa,OutBand,iDim,jDim
REAL(kind=4),dimension(:,:),pointer  :: MosaImage,MosaVar,MosaExpo
INTEGER             :: Status
END SUBROUTINE ReadMosaBandImage

SUBROUTINE WriteSourceList(type,Scw,SouNumber,ObjPtr,IdxPtr,OutBand,&
           InpSourceList,Status)
!-------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER             :: type,Scw,SouNumber,ObjPtr,IdxPtr,OutBand
REAL(kind=4),dimension(:,:),pointer :: InpSourceList
!OUTPUT VARIABLES
INTEGER             :: Status
END SUBROUTINE WriteSourceList


SUBROUTINE UpdateIdx(IdxPtr,Status)
!----------------------------------------------------
  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER :: IdxPtr
  !OUTPUT VARIABLES
  INTEGER            :: Status
END SUBROUTINE  UpdateIdx


SUBROUTINE WriteOutCat(OutCatPtr,Status)
!--------------------------------------
IMPLICIT NONE

!INPUT VARIABLES
INTEGER :: OutCatPtr
!OUTPUT VARIABLES
INTEGER            :: Status
END SUBROUTINE WriteOutCat


END INTERFACE

!*********************************
END MODULE MIMOSA_BASE_MODULE
!*********************************




!*************************************
MODULE MIMOSA_CAT_MODULE
!*************************************
CONTAINS

!............................................
SUBROUTINE CopyCatalogue(InCatPtr,OutCatPtr,&
           Status)
!........................................... 
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER ,INTENT(IN) :: InCatPtr,OutCatPtr


!OUTPUT VARIABLES

INTEGER::Status
!LOCAL VARIABLES
INTEGER        :: iok,sou,k,col,ien,i,isn
INTEGER        :: startbin,endbin,numvalues
REAL(kind=4),dimension(:),pointer :: tab
REAL(kind=4)::ene1,ene2
REAL(kind=4),dimension(:,:),pointer :: fluxTab,e1Tab,e2Tab
REAL(kind=8),dimension(:),pointer  ::tabd
REAL(kind=4),dimension(:,:),pointer  ::tab8f,tab32f
LOGICAL ::jest
CHARACTER(len=20) :: colName
CHARACTER(len=16),dimension(:),pointer  :: tab16char
CHARACTER(len=32),dimension(:),pointer  :: tab32char
CHARACTER(len=64),dimension(:),pointer  :: tab64char
CHARACTER(len=100),dimension(:),pointer  :: tab100char
CHARACTER(len=128),dimension(:),pointer  :: tab128char
INTEGER(kind=2),dimension(:),pointer :: tabi
CHARACTER(len=20)                 :: procName

Status = ISDC_OK
procName = 'CopyCatalogue'
if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)

allocate(tab(SourceNumber),&
         tab16char(SourceNumber),&
         tab32char(SourceNumber),&
         tab64char(SourceNumber),&
         tab100char(SourceNumber),&
         tab128char(SourceNumber),&
  tab8f(8,SourceNumber),tab32f(32,SourceNumber),&
  tabd(SourceNumber),fluxTab(256,1),e1Tab(256,1),e2Tab(256,1),&
       tabi(SourceNumber),stat=iok)

if(iok.ne.0)then
   call MESSAGE(procName,'allocation  problem',&
                   AllocError,Status)
      return
endif
tab(:)=0.
tab16char(:)=''
tab32char(:)=''
tab64char(:)=''
tab100char(:)=''
tab128char(:)=''
tab8f(:,:)=0.
tab32f(:,:)=0.
tabd(:)=0.
fluxTab(:,:)=0.
e1Tab(:,:)=0.
e2Tab(:,:)=0.
tabi(:)=0


! ADDING TO THE OUTPUT CATALOQUE INPUT SOURCES ROWS
Status = dal_table_add_rows(OutCatPtr,0,SourceNumber,Status)
 if(Status.ne.ISDC_OK) then
       call MESSAGE(procName,'cannot add rows to output cat.',&
                   OutCatError,Status)
      return
 endif
Status = dal_table_get_num_rows(OutCatPtr,k,Status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
  call MESSAGE(procName,'dal_table_get_num_rows problem.',&
                   IsdcExitStatus,Status)
      return
endif
if(k.ne.SourceNumber)then
   call MESSAGE(procName,'incorrect row number in output cat.',&
                   OutCatError,Status)
      return
endif

 ! RIGHT ASCENSION EXTRACTING AND COPYING
  call  CopyColFloat(InCatPtr,OutCatPtr,&
           'RA_OBJ',InCatAlpha,Status)
if(Status.ne.ISDC_OK)then
        call WAR_MESSAGE(procName,&
          ' Cannot copy RA_OBJ column into output cat',0,Status)
           CopyCatStat =CopyCatStat +1
      endif  
!DECLINATION EXTRACTING AND COPYING
  call  CopyColFloat(InCatPtr,OutCatPtr,&
           'DEC_OBJ',InCatDelta,Status)
if(Status.ne.ISDC_OK)then
        call WAR_MESSAGE(procName,&
          'Cannot copy DEC_OBJ column into output cat',0,Status)
           CopyCatStat =CopyCatStat +1
      endif  


  !FLUX EXTRACTING  AND COPYING
   InCatFlux = 0.0
   startBin = 1
   endBin   = 256
   numValues= 256
   !loop on input sources
   do sou=1,SourceNumber
     !E_MIN
     call CopyColBinFloat(InCatPtr,OutCatPtr,&
          'E_MIN',sou,sou,e1tab,Status)
     if(Status.ne.ISDC_OK)then
        call WAR_MESSAGE(procName,&
          'Cannot copy E_MIN column into output cat',0,Status)
           CopyCatStat =CopyCatStat +1
      endif  
     !E_MAX
     call CopyColBinFloat(InCatPtr,OutCatPtr,&
          'E_MAX',sou,sou,e2tab,Status)
      if(Status.ne.ISDC_OK)then
        call WAR_MESSAGE(procName,&
          'Cannot copy E_MAX column into output cat',0,Status)
           CopyCatStat =CopyCatStat +1
      endif  
     !FLUX
     call CopyColBinFloat(InCatPtr,OutCatPtr,&
          'FLUX',sou,sou,fluxtab,Status)
      if(Status.ne.ISDC_OK)then
        call WAR_MESSAGE(procName,&
          'Cannot copy FLUX column into output cat',0,Status)
           CopyCatStat =CopyCatStat +1
      endif   
     do ien = 1,OutEnergyNumber
        ! calculating of flux in each output energy band
        ene1 =OutEnergyBands(ien,1)
        ene2 = OutEnergyBands(ien,2) 
        i=1
        jest = .false.
        do while((i .le.256).and.(.not. jest))
           if(e1tab(i,1).ge.ene1)then
             if(e2tab(i,1).le.ene2) then
                InCatFlux(sou,ien) = InCatFlux(sou,ien)+fluxTab(i,1)
             else
              jest = .true.
            endif
          endif
          i=i+1   
       enddo

     enddo 

   !FLUX_ERR
     call CopyColBinFloat(InCatPtr,OutCatPtr,&
          'FLUX_ERR',sou,sou,fluxtab,Status)
    if(Status.ne.ISDC_OK)then
     call WAR_MESSAGE(procName,&
       'Cannot copy float bin column into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif   
   
   enddo  !sources loop 

! COPYING OF DAY_ID
Status = DAL_TABLE_GET_COL(InCatPtr,&
         'DAY_ID',1,DAL_DOUBLE,SourceNumber,addrof(tabd(1)),Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot get DAY_ID column',0,Status)
        CopyCatStat =CopyCatStat +1
   endif   
Status = DAL_TABLE_PUT_COL(OutCatPtr,&
         'DAY_ID',1,DAL_DOUBLE,SourceNumber,addrof(tabd(1)),Status)

 if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy DAY_ID into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif   

! COLUMNS OF TYPE FLOAT  with several bins inside
  !SPA_PARS
call CopyColBinFloat(InCatPtr,OutCatPtr,&
          'SPA_PARS',1,SourceNumber,tab8f,Status)
 if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy SPA_PARS  into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif   
  !SPE_PARS
 call CopyColBinFloat(InCatPtr,OutCatPtr,&
          'SPE_PARS',1,SourceNumber,tab32f,Status)
    
 if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy SPE_PARS  into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif
!VAR_PARS
 call CopyColBinFloat(InCatPtr,OutCatPtr,&
          'VAR_PARS',1,SourceNumber,tab8f,Status)
 if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy VAR_PARS  into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif

! COLUMNS OF TYPE FLOAT COPYING
!  ERR_RAD
do col=1,9
select case(col)
case(1)
  colName = 'ERR_RAD'
case(2)
  colName = 'SPI_FLUX_1'
case(3)
  colName = 'SPI_FLUX_2'
case(4)
  colName = 'ISGR_FLUX_1'
case(5)
  colName = 'ISGR_FLUX_2'
case(6)
  colName = 'PICS_FLUX_1'
case(7)
  colName = 'PICS_FLUX_2'
case(8)
  colName = 'JEMX_FLUX_1'
case(9)
  colName = 'JEMX_FLUX_2'
end select

call  CopyColFloat(InCatPtr,OutCatPtr,&
           colname,tab,Status)

if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy float column  into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif
enddo



! COLUMNS OF TYPE CHAR
!  SOURCE_ID

call  CopyColChar(InCatPtr,OutCatPtr,&
           'SOURCE_ID',tab16char,Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy SOURCE_ID  into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif
InCatId = tab16char

!  NAME
call  CopyColChar(InCatPtr,OutCatPtr,&
           'NAME',tab100char,Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy NAME into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif
InCatName = tab100char

! SPA_MODL
call  CopyColChar(InCatPtr,OutCatPtr,&
           'SPA_MODL',tab32char,Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy SPa_MODL into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif

! SPE_MODL
call  CopyColChar(InCatPtr,OutCatPtr,&
           'SPE_MODL',tab128char,Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy SPE_MODL into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif
! VAR_MODL
call  CopyColChar(InCatPtr,OutCatPtr,&
           'VAR_MODL',tab128char,Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy VAR_MODL into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif
! COMMENTS
call  CopyColChar(InCatPtr,OutCatPtr,&
           'COMMENTS',tab128char,Status)
  if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy COMMENTS into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif
   

!COLUMNS of type 1U
do col=1,5
   select case(col)
   case(1)
     colName = 'CLASS'
   case(2)
     colName = 'SPA_NPAR'
   case(3)
     colName = 'SPE_NPAR'
   case(4)
     colName = 'VAR_NPAR'
   case(5)
     colName = 'SEL_FLAG'
   end select
   call  CopyCol1U(InCatPtr,OutCatPtr,&
           colname,tabi,Status)
   if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot copy 1U column into output cat',0,Status)
        CopyCatStat =CopyCatStat +1
   endif
   

enddo

! 27.01.05 SCREW 1609
colName = 'ISGRI_FLAG'
call  GetCol1U(InCatPtr,OutCatPtr,&
           colname,tabi,Status)
if(Status.ne.ISDC_OK)then
   call MESSAGE(procName,&
        'above error is not FATAL : Column ISGRI_FLAG not found in the input catalog ',0,Status)
   call MESSAGE(procName,&
        ' no fixing position of individual sources will be done',0,Status)
   
   CopyCatStat =CopyCatStat +1
else
   call  PutCol1U(InCatPtr,OutCatPtr,&
        colname,tabi,Status)
   if(Status.ne.ISDC_OK)then
      call MESSAGE(procName,&
           'Column ISGRI_FLAG cannot be copied into output catalog ',0,Status)
 
      CopyCatStat =CopyCatStat +1
   endif

   ! creating InCatFixedList if not defined in hidden file
   ! otherwise ignore
   if(.not.associated(InCatFixedList))then
      isn = 0.
      do i=1,SourceNumber
         if(tabi(i)==2)then
            isn = isn+1
         endif
      enddo
      if(isn >0)then
         allocate(InCatFixedList(isn),stat=iok)
         if(iok .ne.0)then
            call war_message(procname,&
  ' cannot allocate InCatFixedList, sources positions will not be fixed',&
                 0,status)
         else
            isn = 0.
            do i=1,SourceNumber
               if(tabi(i)==2)then
                  isn = isn+1
                  InCatFixedList(isn) = i
               endif
            enddo
         endif!InCatFixedList allocated
      endif ! some sources positions will be fixed
   endif
endif
   

deallocate(tabi)
deallocate(tab)
deallocate(tabd)
deallocate(tab8f)
deallocate(tab32f)
deallocate(tab16char)
deallocate(tab32char)
deallocate(tab64char)
deallocate(tab100char)
deallocate(tab128char)
deallocate(e1tab)
deallocate(e2tab)
deallocate(fluxtab)
if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!............................
END SUBROUTINE CopyCatalogue
!............................

!...........................................
SUBROUTINE CopyColBinFloat(InCatPtr,OutCatPtr,&
           ColName,StartBin,EndBin,Tab,Status)
!...........................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER ,INTENT(IN)     :: InCatPtr,OutCatPtr
INTEGER ,INTENT(IN)     :: StartBin,EndBin
CHARACTER(len=* )                :: ColName

REAL(kind=4),dimension(:,:),pointer:: Tab
!OUTPUT VARIABLES
INTEGER :: Status
!LOCAL VARIABLES
INTEGER::NumValues
CHARACTER(len=20)  :: procName
Status = ISDC_OK
procName = 'CopyColBinFloat'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)


 Status = DAL_TABLE_GET_COL_BINS(inCatPtr,&
         ColName,1,DAL_FLOAT,StartBin,EndBin,NumValues,&
          addrof(Tab(1,1)),Status)
 if(Status.ne.ISDC_OK)then
   if(WorkMode==1)then
         call MESSAGE(procName,&
             'Cannot get'//colName//' from input cat.',&
              ZeroError,Status)
      
   else
         return
  endif
endif
 Status = DAL_TABLE_PUT_COL_BINS(outCatPtr,&
         ColName,1,DAL_FLOAT,StartBin,EndBin,NumValues,&
          addrof(Tab(1,1)),Status)
 if(Status.ne.ISDC_OK)then
   if(WorkMode==1)then
         call MESSAGE(procName,&
             'Cannot put'//colName//' into out cat.',&
              ZeroError,Status)
      
   else
         return
  endif
endif


!..........................
END SUBROUTINE CopyColBinFloat
!.........................

!...........................................
SUBROUTINE CopyColFloat(InCatPtr,OutCatPtr,&
           ColName,Tab,Status)
!...........................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER ,INTENT(IN)     :: InCatPtr,OutCatPtr
CHARACTER(len=* ),INTENT(IN)     ::ColName
REAL(kind=4),dimension(:),pointer:: tab
!OUTPUT VARIABLES
INTEGER::Status
!LOCAL VARIABLES
INTEGER::ns
CHARACTER(len=20)  :: procName
Status = ISDC_OK
procName = 'CopyColFloat'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)

Status = DAL_TABLE_GET_COL(InCatPtr,&
         ColName,1,DAL_FLOAT,ns,addrof(tab(1)),Status)
 if(Status.ne.ISDC_OK)then
   if(WorkMode==1)then
         call MESSAGE(procName,&
             'Cannot get'//colName//' from input cat.',&
              ZeroError,Status)
      
   else
         return
  endif
endif

Status = DAL_TABLE_PUT_COL(OutCatPtr,&
         ColName,1,DAL_FLOAT,ns,addrof(tab(1)),Status)
if(Status.ne.ISDC_OK)then
   if(WorkMode==1)then
         call MESSAGE(procName,&
             'Cannot put'//colName//' into out cat.',&
              ZeroError,Status)
      
   else
         return
  endif
endif


!............................
END SUBROUTINE CopyColFLOAT
!............................

!...........................................
SUBROUTINE GetColFloat(InCatPtr,ColName,Tab,Status)
!...........................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER ,INTENT(IN)     :: InCatPtr
CHARACTER(len=* )                :: ColName
REAL(kind=4),dimension(:),pointer:: tab
!OUTPUT VARIABLES
INTEGER::Status
!LOCAL VARIABLES
INTEGER:: SN
CHARACTER(len=20)  :: procName
Status = ISDC_OK
procName = 'GetColFloat'


Status = DAL_TABLE_GET_COL(InCatPtr,&
         ColName,1,DAL_FLOAT,SN,addrof(tab(1)),Status)

if(DebugMode.eq.3)&
  call Message(procName,' end  ',ZeroError,Status)
!............................
END SUBROUTINE GetColFLOAT
!............................




!...........................................
SUBROUTINE CopyColChar(InCatPtr,OutCatPtr,&
           ColName,Tab,Status)
!...........................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER ,INTENT(IN) :: InCatPtr,OutCatPtr
CHARACTER(len=* ),INTENT(IN)     ::ColName
CHARACTER(len=* ),dimension(:),pointer :: tab
!OUTPUT VARIABLES
INTEGER::Status
!LOCAL VARIABLES
INTEGER::n,iok,ns
CHARACTER(len=100),dimension(MaxSouInGenCat) :: taba
CHARACTER(len=20)  :: procName
Status = ISDC_OK
procName = 'CopyColChar'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)



Status = dal_table_get_col_strings&
         (InCatPtr,ColName,1,1,0,&
         ns,taba,Status)

tab(1:ns) = taba(1:ns)

if(Status.ne.ISDC_OK)then
   if(WorkMode==1)then
         call MESSAGE(procName,&
             'Cannot get'//colName//' from input cat.',&
              ZeroError,Status)
      
   else
         return
  endif
endif

Status = dal_table_put_col_strings&
         (outCatPtr,ColName,1,1,0,&
          ns,taba,Status)
if(Status.ne.ISDC_OK)then
   if(WorkMode==1)then
         call MESSAGE(procName,&
             'Cannot put'//colName//' into out cat.',&
              ZeroError,Status)
      
   else
         return
  endif
endif

!............................
END SUBROUTINE CopyColChar
!............................

!.........................................
 SUBROUTINE CopyCol1U(InCatPtr,OutCatPtr,&
           ColName,tabi,Status)
!.........................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER ,INTENT(IN) :: InCatPtr,OutCatPtr
CHARACTER(len=* ),INTENT(IN)     ::ColName
INTEGER(kind=2),dimension(:),pointer :: tabi
!OUTPUT VARIABLES
INTEGER::Status
!LOCAL VARIABLES
INTEGER::nsou
CHARACTER(len=20)  :: procName

Status = ISDC_OK
procName = 'CopyCol1U'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)

nsou = size(tabi)
 Status = dal_table_get_col&
          (InCatPtr,ColName,1,DAL_USHORT,nsou,&
          addrof(tabi(1)),Status)
 if(Status.ne.ISDC_OK)then
    IsdcExitStatus=Status
     call MESSAGE(procName,&
         ' dal_table_get_col problem for '//ColName ,&
             IsdcExitStatus,Status)
        return
 endif

 Status = dal_table_put_col&
          (OutCatPtr,ColName,1,DAL_USHORT,nsou,&
          addrof(tabi(1)),Status)
 if(Status.ne.ISDC_OK)then
    IsdcExitStatus=Status
     call MESSAGE(procName,&
         ' dal_table_put_col problem for '//ColName ,&
             IsdcExitStatus,Status)
        return
 endif
!.........................
END SUBROUTINE CopyCol1U
!.........................


!.........................................
 SUBROUTINE PutCol1U(InCatPtr,OutCatPtr,&
           ColName,tabi,Status)
!.........................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER ,INTENT(IN) :: InCatPtr,OutCatPtr
CHARACTER(len=* ),INTENT(IN)     ::ColName
INTEGER(kind=2),dimension(:),pointer :: tabi
!OUTPUT VARIABLES
INTEGER::Status
!LOCAL VARIABLES
INTEGER::nsou
CHARACTER(len=20)  :: procName

Status = ISDC_OK
procName = 'PutCol1U'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)

nsou = size(tabi)


 Status = dal_table_put_col&
          (OutCatPtr,ColName,1,DAL_USHORT,nsou,&
          addrof(tabi(1)),Status)
 if(Status.ne.ISDC_OK)then
     call MESSAGE(procName,&
         ' dal_table_put_col problem for '//ColName ,&
             0,Status)
     status = 1
        return
 endif
!.........................
END SUBROUTINE PutCol1U
!.........................



!.........................................
 SUBROUTINE GetCol1U(InCatPtr,OutCatPtr,&
           ColName,tabi,Status)
!.........................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER ,INTENT(IN) :: InCatPtr,OutCatPtr
CHARACTER(len=* ),INTENT(IN)     ::ColName
INTEGER(kind=2),dimension(:),pointer :: tabi
!OUTPUT VARIABLES
INTEGER::Status
!LOCAL VARIABLES
INTEGER::nsou
CHARACTER(len=20)  :: procName

Status = ISDC_OK
procName = 'GetCol1U'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)

nsou = size(tabi)
 Status = dal_table_get_col&
          (InCatPtr,ColName,1,DAL_USHORT,nsou,&
          addrof(tabi(1)),Status)
 if(Status.ne.ISDC_OK)then
    IsdcExitStatus=Status
     call MESSAGE(procName,&
         ' dal_table_get_col problem for '//ColName ,&
             IsdcExitStatus,Status)
        return
 endif


!.........................
END SUBROUTINE GetCol1U
!.........................
!*************************************
END MODULE MIMOSA_CAT_MODULE
!*************************************

!*************************
MODULE INTER1
!*************************

INTERFACE

SUBROUTINE PutCatNull(OutCatPtr,n1,n2,Status)
!-----------------------------------
IMPLICIT NONE
!INPUT VARIABLES
INTEGER :: OutCatPtr,n1,n2,Status
END SUBROUTINE PutCatNull


END INTERFACE
!******************
END MODULE INTER1
!******************
!##################################################
! SUBROUTINES OF THE MODULE MIMOSA_BASE_MODULE 
!##################################################


!==========================================================
SUBROUTINE OgOpen(MosaNumber,MosaImaIdxPtr,MosaSouIdxPtr,&
                  OutCatPtr,RaX,DecX,RaZ,DecZ,posAngle,Status)
!===========================================================

USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API
USE DAL3CAT_F90_API
USE DAL3HK_F90_API  !NP
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_IMAGING_PARAM
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
USE MIMOSA_BASE1_MODULE
USE MIMOSA_BASE2_MODULE
USE MIMOSA_BASE3_MODULE
USE MIMOSA_BASE4_MODULE
USE MIMOSA_CAT_MODULE


IMPLICIT NONE 
INTEGER :: MosaNumber

!OUTPUT VARIABLES

! pointer to OG
INTEGER,Intent(out)         :: MosaImaIdxPtr
! pointer to mosaicked image list
INTEGER,Intent(out)         :: MosaSouIdxPtr
! pointer to mosaicked image source list
INTEGER,Intent(out)         :: OutCatPtr
! pointer to output source catalogue
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle
INTEGER,Intent(out)         :: Status

!LOCAL  VARIABLES

CHARACTER (len=128),save  :: InGRPpar  = 'inOG' 
CHARACTER (len=128),save  :: OutGRPpar  = 'outOG'
CHARACTER (len=128),save  :: inDOLpar  =  'inCat '
CHARACTER (len=128),save  :: outDOLpar = 'outCat,outMosIma,outMosRes '
INTEGER     :: inCatPtr,GrpPtr
INTEGER     :: auxcatPtr,ptr,infixsou,localstatus,mosaIdxPtr
INTEGER     :: nFound,num,tb,input_proj
INTEGER     :: type_ds,inst,level,iok,i,j,attiNumber,ls,scwVerifRes
CHARACTER (len=DAL_FILE_NAME_STRING) :: name
CHARACTER (len=30) :: strucname


REAL(kind=4),dimension(:,:),pointer :: MapArray
CHARACTER(100)::bla="&\0"
REAL(kind=4),dimension(:,:),pointer :: EBins
INTEGER,dimension(:,:),pointer ::  Chans
CHARACTER(len=100) :: selection,sel
CHARACTER (len=20)  :: procName,colname
CHARACTER (len=10) :: str10

INTEGER  :: imin,imax,jmin,jmax,ii
INTEGER  :: idim,jdim,type,numAxes
INTEGER,dimension(1:DAL_MAX_ARRAY_DIMENSION)::axes 
REAL(kind=8)      ::val,incr

procName = 'OgOpen'
Status = ISDC_OK

if(DebugMode.eq.3)&
  call Message(procName,'begin ',ZeroError,Status)

call CreateElem(OutCatFile,OutCatStrucName,OutCatPtr,Status)
if(Status.ne.ISDC_OK)return

call CreateElem(OutMosImaFile,MosaImaIdxStrucName,MosaImaIdxPtr,Status)
if(Status.ne.ISDC_OK)return
    

call CreateElem(OutMosResFile,MosaSouIdxStrucName,MosaSouIdxPtr,Status)
if(Status.ne.ISDC_OK)return

open(20,file='mimosa.mosa')
read(20,'(i5)')MosaNumber
close(20)

write(str250,*)'Number of scw to be mosaicked  : ',mosanumber
call message(procname,str250,0,status)

if(MosaNUmber .le.0)then
 call war_MESSAGE(procName,'Empty input : no mosaics ',&
        DataError,Status)
   return
endif 

allocate(MosaNames(MosaNumber),DefMapSize(MosaNumber,2),&
     RealMapSize(MosaNumber,2,2),MosaType(MosaNumber),&
     InputProjType(MosaNumber),ActiveScwArray(MosaNumber),stat=iok)
if(iok.ne.0)then
      call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
endif
ActiveScwArray(:) = 1
RealMapSize = 0
DefMapSize = 0
InputProjType = 0
MosaType = 0


attiNumber =MosaNumber
num=attiNumber+1
Allocate(RaX(num),DecX(num),RaZ(num), DecZ(num),posAngle(num),&
     Duration(1:MosaNumber),LX(num),BX(num),ShdTimes(1:MosaNumber,1:7),stat=iok)

IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',&
        AllocError,Status)
   return
endif

ShdTimes(:,:) = 0.
LX(:) = 0.
BX(:) = 0.
RaX(:)=0.
DecX(:)=0.
RaZ(:)=0.
DecZ(:)=0.
posAngle(:)=0.
Duration(:)=0.
incr= 180.0d0*atan(IsgrPixSiz / IsgrHalfMasDis)/acos(-1.0d0)



write(str250,*)' Reading of :',MosaNumber, ' mosaics'
call message(procname,str250,ZeroError,Status)


if (ReCall==0)then
   call message(procname,' Verifying mosaics list',ZeroError,Status)
   call Read1Mosa(MosaNumber,EBins,Chans,Status)
   if(Status .ne. ISDC_OK)return
else
   call ReadRecallInfo(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle)
   call  Read1ReCallMosa(MosaNumber,EBins,Chans,Status)
   if(Status .ne. ISDC_OK)return
endif

if(EnergyBand > 0) then 
   if(EnergyBand > OutEnergyNumber)then
      write(str250,*)'User chosen energy band number :',EnergyBand,' > ',OutEnergyNumber
     call message(procname,str250,ZeroError,Status)
     call message(procname,' All energy bands will be treated',ZeroError,Status) 
     EnergyBand=0
 else
     write(str250,*)'User chosen energy band  :',EnergyBand, ' will be treated'
     call message(procname,str250,ZeroError,Status)
  endif
else
   write(str250,*)'All energy bands   will be treated'
     call message(procname,str250,ZeroError,Status)
endif
!......................
!INPUT SOURCE CATALOG 
!......................

SourceNumber = 0 ! number of sources in the input catalogue


CoorType=1
FluxNorm=1.
BkgMoyPar = PredefBkg ! meanfull for simulated image only

     
name = InCatStrucName
! IPUT SOURCE CATALOGUE OPENING


call MESSAGE(procName,incatfile,0,status)
if(len_trim(InCatFile) > 0)then
   ! outside input file
   Status = dal_object_open(InCatFile,inCatPtr,Status)
   if(Status.ne.ISDC_OK)then
      call MESSAGE(procName,&
           ' dal_object_open problem for input catalog',&
           DataError,Status)
      return
   endif
else
   if(nToClean > 0)then
       call MESSAGE(procName,&
           ' Empty input catalog but wants to search for sources',&
           DataError,Status)
      return
   endif
endif
write(str250,*)' input catalog analysed '
call message(procname,str250,ZeroError,Status)

!......................
!OUTPUT SOURCE CATALOG 
!......................

 
!NUMBER OF SOURCES
Status = dal_table_get_num_rows(inCatPtr,SourceNumber,Status)   
if((Status.ne.ISDC_OK).and.(nToClean > 0))then
   IsdcExitStatus=Status
   call MESSAGE(procName,'dal_table_get_num_rows  problem'//name,&
        IsdcExitStatus,Status)
   return
endif
Status=0
if(SourceNumber == 0)then
   call war_message(procname,' Empty input source catalogue ',0,status)
else

   if(SourceNumber >MaxSouInGenCat )then
      call message(procname,&
           ' Too many sources in the input catalog',InCatError,Status)
      return
   endif

   Allocate(InCatAlpha(SourceNumber),InCatFixed(SourceNumber),&
        InCatDelta(SourceNumber),InCatId(SourceNumber),&
        InCatName(SourceNumber),&
        InCatFlux(SourceNumber,OutEnergyNumber),&
        stat=iok)
   IF (iok /= 0) then
      call MESSAGE(procName,'Allocation problem',&
           AllocError,Status)
      return
   endif
         
     
         

   InCatAlpha(:)=0.
   InCatDelta(:)=0.
        
   InCatId(:)= ''
   InCatName(:)=''
   InCatFlux(:,:)=0.
        
   ! EXTRACTING OF COORDINATES 
   ! AND WRITING A COPY TO A NEW CATALOGUE
   call CopyCatalogue(inCatPtr,OutCatPtr,Status)
   if(Status.ne.ISDC_OK)return
   status = dal_object_close(inCatPtr,DAL_SAVE,Status)
   if(Status.ne.ISDC_OK)then
      call war_message(procname,&
           ' Closing input catalogue problem',0,Status)
   endif

   call  FixSouPosInCat(Status)

endif ! not empty input catalogue

write(str250,*)' output catalog analysed '
call message(procname,str250,ZeroError,Status)

if(ReCall==0)then
   call message(procname,'Reading attributes and setting real mosaics size ',ZeroError,Status)
   call Read2Mosa(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle,incr,Status)
   if(Status.ne.ISDC_OK)return
endif


! CREATES/CLEARS MOSAICKED IMAGE INDEX

call ReplaceIdx(GrpPtr,MosaImaIdxPtr,&
     MosaImaIdxStrucName,MosaImaStrucName,&
     OutMosImaFile,MosaImaIdxFitsName,.true.,Status)
  
if(Status.ne.ISDC_OK)return



! CREATES/CLEARS MOSAICKED IMAGE SOURCE LIST INDEX

 call ReplaceIdx(GrpPtr,MosaSouIdxPtr,MosaSouIdxStrucName,MosaSouStrucName,&
        OutMosResFile,MosaSouIdxFitsName,.true.,Status)
  
 if(Status.ne.ISDC_OK)return

call  WriteRecallInfo(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle)

!=========================================
END SUBROUTINE OgOpen
!=========================================





!======================================================
SUBROUTINE WriteImages(type,scw,OutBand,ObjPtr,IdxPtr,CleanImage,&
                        VarImage,ResidImage,SignifImage,&
                        Status)
!======================================================
USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API  
USE MIMOSA_GLOBVAR_MODULE  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE1_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER,Intent(in)   :: type
! 0 - Scw images
! 1 - Map images
INTEGER,Intent(in)   :: scw
!pointing number
INTEGER,Intent(in)   :: OutBand
!out energy band number
INTEGER   :: ObjPtr,IdxPtr
!pointer to the index of Scw image index
REAL(kind=4),dimension(:,:),pointer  ::CleanImage, VarImage,&
                               ResidImage,SignifImage 

!OUTPUT VARIABLES
INTEGER   :: Status

!LOCAL VARIABLES
INTEGER   :: ptr,num1,num2,flag

CHARACTER(len=20) ::strucName, procname

procName = 'WriteImages '
Status = ISDC_OK
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)


! mosa index

where(VarImage.eq.0)
   CleanImage = nanf
   VarImage = nanf
   SignifImage = nanf
endwhere
strucName = MosaImaStrucName
flag = 2



Status = dal3gen_index_get_num_members&
         (IdxPtr,strucName,num1,status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
  call MESSAGE (procName,&
       'dal3gen_index_get_num_members problem  ',&
        IsdcExitStatus,Status)
   return
endif

call AddIndexMemberArray(ObjPtr,IdxPtr,scw,OutBand,flag,&
     strucName,"IMATYPE","INTENSITY",CleanImage,&
     ptr,Status)
if(Status.ne.ISDC_OK)return

ptr = 0
call AddIndexMemberArray(ObjPtr,IdxPtr,scw,OutBand,flag,&
     strucName,"IMATYPE","VARIANCE",VarImage,ptr,Status)
if(Status.ne.ISDC_OK)return

ptr = 0
call AddIndexMemberArray(ObjPtr,IdxPtr,scw,OutBand,flag,&
     strucName,"IMATYPE","SIGNIFICANCE",SignifImage,ptr,Status)
if(Status.ne.ISDC_OK)return

ptr = 0


call AddIndexMemberArray(ObjPtr,IdxPtr,scw,OutBand,flag,&
     strucName,"IMATYPE","EXPOSURE",ResidImage,&
     ptr,Status)

if(Status.ne.ISDC_OK)return

ptr = 0
Status = dal3gen_index_get_num_members&
         (IdxPtr,strucName,num2,status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE (procName,&
       'dal3gen_index_get_num_members problem  ',&
        IsdcExitStatus,Status)
   return
endif


!==============================
END SUBROUTINE WriteImages
!==============================

!=============================================================
 SUBROUTINE ReadMosaBandImage(mosa,outBand,iDim,jDim,&
                             MosaImage,MosaVar,&
                            MosaExpo,Status)
!=============================================================
USE ISDC
USE DAL3GEN_F90_API  
USE DAL3AUX_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE1_MODULE
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES

INTEGER ,Intent(in)  :: mosa,OutBand,iDim,jDim
REAL(kind=4),dimension(:,:),pointer  :: MosaImage,MosaVar,MosaExpo

INTEGER             :: Status

!LOCAL VARIABLES
INTEGER   :: GrpPtr,i
REAL(kind=4)      :: eMinVal,eMaxVal
INTEGER   :: scwPtr,imaIdxPtr,imaPtr,ptr,detptr
INTEGER   :: nFound,n,imat,iok
CHARACTER(len=20) ::imaType
CHARACTER(len=4),dimension(5) ::imaTypeTab
CHARACTER(len=15) :: stre1,stre2
CHARACTER(len=10) :: strtol
CHARACTER(len=350) :: selection,detselection
LOGICAL  :: okinten,okvar,oksign,okexpo,oknosou
REAL(kind=4),dimension(:,:),pointer  :: detimage
CHARACTER(len=2000)::detname
CHARACTER(len=20) :: procName
CHARACTER(len=30) :: StrucName
Status = ISDC_OK
procName = 'ReadMosaBandImage'

detmean = 0.

imaTypeTab(1) = 'INTE'
imaTypeTab(2) = 'VARI'
imaTypeTab(3) = 'SIGN'
imaTypeTab(4) = 'EXPO'

ogfile=mosaNames(mosa)
Status = dal_object_open(trim(OgFile)//'[1]',GrpPtr,Status)
if(Status.ne.ISDC_OK)then
   print *,'dal_object_open problem , status :',status
   call MESSAGE(procName,'dal_object_open problem for '//trim(ogfile),&
        DataError,Status)
   return
endif

strtol = '0.001'
    write(stre1,'(g14.7)',err=170)OutEnergyBands(OutBand,1)
    goto 171
170 call message(procname,' cannot write OutEnergyBands(OutBand,1)',&
          StringWriteError,Status)
    return
171 write(stre2,'(g14.7)',err=172)OutEnergyBands(OutBand,2)
    goto 173
172 call message(procname,' cannot write OutEnergyBands(OutBand,2) ',&
          StringWriteError,Status)
      return
173 continue
selection = '( ( near(E_MIN,'//stre1//','//strtol//'))'
selection= selection(1:len_trim(selection))//'&&( near(E_MAX,'//stre2//','//strtol//')) ) &&'
detselection = selection
selection= selection(1:len_trim(selection))//&
     "(IMATYPE=='INTENSITY'.or.IMATYPE=='VARIANCE'.or.IMATYPE=='SIGNIFICANCE' "
selection= selection(1:len_trim(selection))//".or.IMATYPE=='EXPOSURE') "


if(mosaType(mosa)==0)then
   !skyImage

   if(RebinType==17)then
      !reading detector data
      detmean = 0.
      i=scan(ogfile,'/',back=.true.)
      if(i==0)then
         call war_message(procname,&
              'Cannot find detector data directory, covar estim will be used',ZeroError,Status)
      else
         detname=ogfile(1:i)
         detname=trim(detname)//'isgri_cor_shad.fits'
         Status = dal_object_open(trim(detname)//'[1]',detPtr,Status)
         if(Status.ne.ISDC_OK)then
            Status = ISDC_OK
            call war_MESSAGE(procName,&
                 'dal_object_open problem for isgri_cor_shad, covar estim will be used',&
                 ZeroError,Status)
         else
            detselection = trim(detselection)// "(SHD_TYPE=='DETECTOR')"
            strucname=shdstrucname
            Status = dal3gen_index_find_member(detPtr,&
                 StrucName,detselection,nFound,Ptr,Status)
            if(Status.ne.ISDC_OK)then
               Status = ISDC_OK
               call WAR_MESSAGE(procName,&
                    'selection problem isgri_cor_shad, covar estim will be used',ZeroError,Status)
            else
               if(nFound.ne.1)then
                  Status = ISDC_OK
                  call MESSAGE(procName,&
                       'selection problem isgri_cor_shad, covar estim will be used', &
                       ZeroError,Status)
              
               else
                  call ReadUnknownArray(ptr,iDim,jDim,detImage,Status)
                  Status = dal_object_close(ptr,DAL_SAVE,Status)
                  Status = dal_object_close(detptr,DAL_SAVE,Status)
                  if(Status .ne.ISDC_OK)then
                     Status = ISDC_OK
                     call MESSAGE(procName,&
                          'data reading problem isgri_cor_shad, covar estim will be used', &
                          ZeroError,Status)
                  else
                     detmean = detimage(33,33)
                  endif
                  if(associated(detimage))deallocate(detimage)
               endif
            endif
         endif
      endif
   endif !detector data reading
      



   strucname=imastrucname
   Status = dal3gen_index_find_member(grpPtr,&
        StrucName,selection,nFound,imaPtr,Status)
   if(Status==0)then
      if(nFound.lt.3)then
         Status = 1
         write(str250,*)' In SkyImage index ',mosa,' too few members found'
         call MESSAGE(procName,str250, NumMembersError,Status)
         return
      else
         if(nFound==3)then
            write(str250,*)' In SkyImage no. ',mosa,'  no true Exposure ==>ONTIME will be taken'
            call MESSAGE(procName,str250, ZeroError,Status)
         endif
      endif
   else
     IsdcExitStatus=Status
     call MESSAGE(procname,' dal3gen_index_find_member problem for '//selection,&
          IsdcExitStatus,Status)
     return
  endif
else
   strucname=mosaimastrucname
   Status = dal3gen_index_find_member(grpPtr,&
        StrucName,selection,nFound,imaPtr,Status)

   if(Status.ne.0)then
      IsdcExitStatus=Status
      call MESSAGE(procname,' dal3gen_index_find_member problem for '//selection,&
           IsdcExitStatus,Status)
      return
   endif
   if(nFound.lt.4)then
      Status = 1
      write(str250,'(" In Mosa ",I5,"members number in the image index : ",I3," lt 3 ")',err=110 )mosa,nFound
      goto 111
110   str250 = 'members number in the image index  ne 4'
111   call MESSAGE(procName,str250, NumMembersError,Status)
      return
   endif
endif


do n=1,nFound !image selection loop
    Status = dal3gen_index_get_member(imaPtr,MosaImaStrucName,&
                 n,ptr,status)
   if(Status.ne.ISDC_OK) then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
       ' dal3gen_index_get_member problem ',&
       IsdcExitStatus,Status)
   return
   endif
   
   call ReadImaAttributes(ptr,eMinVal,eMaxVal,Status)
   if(Status.ne.ISDC_OK) return
  

   ! ENERGY BAND VERIFICATION
  
   if((eMinval.ne.OutEnergyBands(OutBand,1).or.&
       (eMaxVal.ne.OutEnergyBands(OutBand,2))) )then
     call MESSAGE(procName,&
       ' member selection problem :'//&
       '(eMinval.ne.OutEnergyBands(OutBand,1).or.(eMaxVal.ne.OutEnergyBands(OutBand,2))',&
       MemberSelError,Status)
   return
   endif
   !IMA TYPE verification
  
   imaType = adjustl(shd_type)
   if (imaType(1:4).eq.imaTypeTab(1))then
      imat = 1
   else
      if(imaType(1:4).eq.imaTypeTab(2))then
         imat = 2
      else
         if (imaType(1:4).eq.imaTypeTab(3))then
             imat = 3
         else
            if (imaType(1:4).eq.imaTypeTab(4))then
               imat = 4
            else
               call MESSAGE(procName,' image array type invalid ',&
                    InvalArrayError,Status)
               return
            endif
         endif
      endif
   endif

   !reading

   select case(imat)
     case(1) 
       call ReadUnknownArray(ptr,iDim,jDim,mosaImage,Status)
       if(Status .ne.ISDC_OK)return
        okinten=.true.
        allocate(nanfilter(idim,jdim),stat=status)
        if(Status.ne.ISDC_OK)return

     case(2) 
       call ReadUnknownArray(ptr,iDim,jDim,mosaVar,Status)
       if(Status .ne.ISDC_OK)return
        okvar=.true.
    case(3) 
     
    case(4) 
       call ReadUnknownArray(ptr,iDim,jDim,mosaExpo,Status)
        if(Status .ne.ISDC_OK)return
       okexpo=.true.
   
    end select
Status = dal_object_close(Ptr,DAL_DELETE,status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,'Error closing ima ptr  for '//selection,IsdcExitStatus,Status)
   return
endif
!call deco(scwcleaimage)
if(Status.ne.0)return
enddo ! member selection loop
if((imat==3).and.(mosaType(mosa)==0))then
   if(associated(mosaexpo))deallocate(mosaexpo)
   allocate(mosaexpo(idim,jdim))
   mosaExpo(:,:)=shdtimes(mosa,6)
   okexpo=.true.
endif

call deco(mosaimage)
where(nanfilter==0.)
   mosaImage = 0.
   mosaVar = 0.
   mosaExpo = 0.
endwhere






Status = dal_object_close(imaPtr,DAL_DELETE,status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,'Error closing ima ptr  for '//selection,IsdcExitStatus,Status)
   return
endif

 Status = dal_object_close(grpptr,DAL_SAVE,Status)
 if(Status.ne.ISDC_OK)then
    IsdcExitStatus=Status
    call MESSAGE(procName,&
         ' cannot close Grp ptr  '//trim(ogFile),IsdcExitStatus,Status)
    return
 endif

if(.not.okinten)then
    write(str250,'(" In mosa ",I5," intensity image not found  ")',err=410 )mosa
    goto 411
410 str250 = 'intensity image not found '
411   call MESSAGE(procName,str250 ,NumMembersError,Status)
   return
endif

if(.not.okvar)then
    write(str250,'(" In mosa ",I5," variance image not found  ")',err=220 )mosa
    goto 221
220 str250 = 'variance image not found '
221   call MESSAGE(procName,str250 ,NumMembersError,Status)
   return
endif


if(.not.okexpo)then
   ! Exposure map not found
   if(.not. OneExpoMap)then
      write(str250,'(" In mosa ",I5," Exposure map for a given energy band not found ")' )mosa
      call MESSAGE(procName,str250 ,NumMembersError,Status)
      return
   endif
endif



!===============================
END  SUBROUTINE ReadMosaBandImage
!===============================

!=======================================================
SUBROUTINE WriteSourceList(type,Scw,SouNumber,ObjPtr,IdxPtr,OutBand,&
           InpSourceList,Status)
!======================================================
USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE1_MODULE
USE ATTI_DECLARATIONS
USE IBIS_IMAGING_PARAM
IMPLICIT NONE

!INPUT VARIABLES

INTEGER                        :: type,Scw,SouNumber,ObjPtr,IdxPtr,OutBand
REAL(kind=4),dimension(:,:),pointer    :: InpSourceList
!OUTPUT VARIABLES
INTEGER                        :: Status
!LOCAL VARIABLES
INTEGER,parameter                      :: chardim=1000
REAL(kind=4) :: valr
REAL(kind=8)                           ::val,nand
INTEGER                        :: nsou,ptr,iok,ns,k,size2,col,numer
REAL(kind=4),dimension(:),pointer      ::tab
REAL(kind=4),dimension(:,:),pointer    :: SourceList
INTEGER(kind=2),dimension(:),pointer   :: tabi
CHARACTER(len=16),dimension(chardim)   :: tab16char
CHARACTER(len=100),dimension(chardim)   :: tab100char
REAL(kind=4),dimension(:),pointer      ::  alpha,delta,flux,fluxerr,snr,locerr
CHARACTER(len=100) :: str100
character(len=19) ::strutc
character(len=5) :: strdim
CHARACTER(len=20) :: colName,procName
CHARACTER(len=6)  :: str
CHARACTER(len=1)  ::str1
CHARACTER(len=2)  ::str2
CHARACTER(len=3)  ::str3
CHARACTER(len=4)  ::str4


Status = ISDC_OK
procName = 'WriteSourceList'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)




!creates new index member
!adds energy band info

call  NewIndexMember&
     (IdxPtr,OutBand,MosaSouStrucName,ptr,Status)


if(Status.ne.ISDC_OK)return

   end0="&\0"
   end1="&\0"
   do col=1,7
      !TELAPSE,TSTART,TSTOP,ONTIME,EXPOSURE,DEADC
      val = ogtimes(col)
      colname = OGTimeColNames(col)
     
      Status = dal_attribute_put_real&
              (ptr,colName,val,"&","&",Status)
         if(Status.ne.ISDC_OK) then
            write(str250,*,err=177) &
                 ' dal_attribute_put _real problem for ',colName
            goto 178
177         str250 = ' dal_attribute_put _real problem '&
                    //errorstr//' 177'
178         call WAR_MESSAGE(procName,str250,0,Status)
         endif
      
      end0="&\0"
      end1="&\0"
      if((col==4).or.(col==5))then
         select case (col)
         case(4) ! DATE-OBS
            colName = 'DATE-OBS'
            
         case(5)!DATE-END
            colName = 'DATE-END'
         end select
         val = ogtimes(col)
!!$         Status = dal3gen_convert_ijd2utc(val,0,strUtc,Status)
!!$         if(Status.ne.ISDC_OK) then
!!$            write(str250,*,err=179)'idx dal3gen_convert_ijd2utc problem for ',&
!!$                 colName,' time : ',val
!!$            goto 180
!!$179         str250 = 'idx dal3gen_convert_ijd2utc problem'&
!!$                    //errorstr//' 179'
!!$180         call WAR_MESSAGE(procName,str250,0,Status)
!!$  
!!$         else
!!$            Status = dal_attribute_put_char&
!!$                 (ptr,colName,strutc,"&","&",Status)
!!$            if(Status.ne.ISDC_OK) then
!!$               call WAR_MESSAGE(procName,&
!!$                    ' dal_attribute_put_char problem for '//colName,0,&
!!$                    Status)
!!$            endif
!!$         endif
          end0="&\0"
          end1="&\0"
      endif
   enddo
  


colname = 'BKGPARAM'
 Status = dal_attribute_put_char&
      (ptr,colname,BkgFile,"&","&",Status)
 if(Status.ne.ISDC_OK) then
    call WAR_MESSAGE(procName,&
         ' dal_attribute_put_char problem for '//colName,0,&
         Status)
 endif
      
 end0="&\0"
 end1="&\0"

colname = 'OFFCORR'
if(.not.JestOffCorrIndex)then
   offcorrfilename = 'NONE'
endif
Status = dal_attribute_put_char&
      (ptr,colname,OffCorrFilename,"&","&",Status)
 if(Status.ne.ISDC_OK) then
    call WAR_MESSAGE(procName,&
         ' dal_attribute_put_char problem for '//colName,0,&
         Status)
 endif



!!$if(type==0)then ! image source list
!!$  !COPYING SOME ATTRIBUTES
!!$  Status = dal3gen_attribute_copy(ObjPtr,ptr,&
!!$          'REVOL,SWID,SW_TYPE,SWBOUND,OBTSTART,OBTEND',Status)  
!!$  if(Status.ne.ISDC_OK) then
!!$     call WAR_MESSAGE(procName,&
!!$       'cannot copy keywords to the object '//SouStrucName,&
!!$        0,Status)
!!$      WriteSourceListStat = WriteSourceListStat+1
!!$  endif
!!$else !mapsourcelist
!!$
!!$   Status = dal3gen_attribute_copy(ObjPtr,ptr,'OGID',Status)
!!$   if(Status.ne.ISDC_OK) then
!!$     call WAR_MESSAGE(procName,&
!!$       'cannot copy keywords to the object '//SouStrucName,&
!!$        0,Status)
!!$      WriteSourceListStat = WriteSourceListStat+1
!!$  endif
!!$endif




nsou = SouNumber
if(nsou.gt.chardim)then
    write(strdim,*)chardim
    str250 = ' Number of sources to be written into .res file > '//&
         strdim//'==> only first '//strdim//'willbe written'
    call war_MESSAGE(procName,str250,&
                   0,Status)
    nsou = chardim
endif
if(nsou.eq.0)then
   call WAR_MESSAGE(procName,'no sources found in the current FOV',&
                  0,Status)
   Status = dal3gen_index_update(ptr,IdxPtr,Status)
   if(status == -25109) then
      call WAR_MESSAGE(procName,&
           'Source list index updating : keyword missing ',0,&
           Status)
  
      WriteSourceListStat = WriteSourceListStat +1
   else
      if(Status.ne.ISDC_OK)then
         call WAR_MESSAGE(procName,'Source list index updating error ',&
              1,Status)
         WriteSourceListStat = WriteSourceListStat +1
    
      endif
   endif

   return
endif ! empty source list

if(type==0)then !scw source list allocation
   AllSouList(scw,OutBand)%souNumber = nsou
   allocate(AllSouList(scw,OutBand)%alpha(nsou),&
        AllSouList(scw,OutBand)%delta(nsou),&
        AllSouList(scw,OutBand)%xmap(nsou),&
        AllSouList(scw,OutBand)%ymap(nsou),stat=iok)
   if(iok.ne.0)then
      call MESSAGE(procName,'Allocation problem',&
           AllocError,Status)
      return
   endif
  
   AllSouList(scw,OutBand)%xmap = 0.
   AllSouList(scw,OutBand)%ymap = 0.

endif
size2 = size(InpSourceList,2)
allocate(SourceList(nsou,size2),stat=iok)
if(iok.ne.0)then
   call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
endif
SourceList(:,:)=0.
!source list sorting
! nsou - number of sources found in the current FOV
k=0
do ns = 1,nsou
   if(InpSourceList(ns,1).gt.0)then
     !source from input catalogue
     k=k+1
     SourceList(k,1:size2) = InpSourceList(ns,1:size2)
   endif
enddo
do ns = 1,nsou
   if(InpSourceList(ns,1).eq.0)then
     !new source - placed after input catalogue sources
     k=k+1
     SourceList(k,1:size2) = InpSourceList(ns,1:size2)
   endif
enddo 
! SourceList - sorted list of in-FOV found sources
allocate(tabi(nsou),tab(nsou),&
         alpha(nsou),delta(nsou),flux(nsou),fluxerr(nsou),&
         snr(nsou),locerr(nsou),stat=iok)
if(iok.ne.0)then
   call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
endif
tabi(:)=0
tab(:)=0.
alpha(:)=0.
delta(:)=0.
flux(:)=0.
fluxerr(:)=0.
snr(:) = 0.
locerr(:) = 0.

if(RawFluxEstim)then
   if(type==0)then
      !private option in the scw source list
      Status = dal_table_add_col(ptr,'RAW_FLUX',DAL_FLOAT,1,status)
      if(Status.ne.ISDC_OK) then
         call WAR_MESSAGE(procName,&
              'cannot add RAW_FLUX to source list ',ZeroError,Status)
       
      endif
   endif
endif


Status = dal_table_add_rows(ptr,0,nsou,Status)
if(Status.ne.ISDC_OK) then
   call WAR_MESSAGE(procName,'cannot add rows to source list index.',1,Status)
   call WAR_CONT_MESSAGE(procName,'Source list will not be updated',1,Status)
   WriteSourceListStat = WriteSourceListStat+1
   return
endif

Status = dal_table_get_num_rows(ptr,k,Status)
if(Status.ne.ISDC_OK) then
   call WAR_MESSAGE(procName,'dal_table_get_num_rows problem. ',1,Status)
   call WAR_CONT_MESSAGE(procName,'Source list will not be updated',1,Status)
   WriteSourceListStat = WriteSourceListStat+1
  return
endif
if(k.ne.nsou)then
    call WAR_MESSAGE(procName,'cannot add correct number of rows',1,Status)
    call WAR_CONT_MESSAGE(procName,'Source list will not be updated',1,Status)
    WriteSourceListStat = WriteSourceListStat+1
      return
endif
!SOURCE_ID and !RA_OBJ
numer = 0
do ns=1,nsou
  ! k - number in input catalogue
  k = int(SourceList(ns,1))
  if(k.eq.0) then

      numer = numer+1

      if( numer < 10) then
         write(str1,'(i1)')numer
         !ident
         tab16char(ns) = 'NEW_'//str1
         !name
         tab100char(ns) = 'NEW_'//str1
     
      else
         if(numer < 100)then
            write(str2,'(i2)')numer
            !ident
            tab16char(ns) = 'NEW_'//str2
            !name
            tab100char(ns) = 'NEW_'//str2
           
         else
            if(numer < 1000)then
               write(str3,'(i3)')numer
               !ident
               tab16char(ns) = 'NEW_'//str3
               !name
               tab100char(ns) = 'NEW_'//str3
            
            else
               write(str4,'(i4)')numer
                !ident
               tab16char(ns) = 'NEW_'//str4
               !name
               tab100char(ns) = 'NEW_'//str4
               
            endif
         endif
      endif
    
      tab(ns) =nanf
       
  else
      !alpha
      tab(ns) = InScwCat(4,k)
      !ident
      tab16char(ns) = &
          InCatId(Int(InScwCat(10,k)))
      !name
      tab100char(ns) = &
          InCatName(Int(InScwCat(10,k)))
  endif
enddo 
!SOURCE_ID
Status = dal_table_put_col_strings&
         (ptr,'SOURCE_ID',1,1,nsou,nsou,tab16char,Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot put SOURCE_ID'//&
        ' into source list',0,Status)
  WriteSourceListStat = WriteSourceListStat +1
endif
!NAME
Status = dal_table_put_col_strings&
         (ptr,'NAME',1,1,nsou,nsou,tab100char,Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
       'Cannot put NAME'//&
        ' into source list',0,Status)
  WriteSourceListStat = WriteSourceListStat +1
endif

!RA_OBJ
Status = dal_table_put_col(ptr,'RA_OBJ',1,DAL_FLOAT,nsou,&
         addrof(tab(1)),Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
             ' dal_table_put_col problem for RA_OBJ',0,Status)
   WriteSourceListStat = WriteSourceListStat +1
endif

!DEC column
do ns=1,nsou
  ! k - number in input catalogue
  k = int(SourceList(ns,1))
  if(k.eq.0) then
      tab(ns) = nanf
      tabi(ns) = 1
  else
      tab(ns) = InScwCat(5,k)
      tabi(ns) = 0
  endif
enddo 
!NEW_SOURCE
 Status = dal_table_put_col&
          (ptr,'NEW_SOURCE',1,DAL_USHORT,nsou,&
          addrof(tabi(1)),Status)
     if(Status.ne.ISDC_OK)then
        call WAR_MESSAGE(procName,&
             ' dal_table_put_col problem for NEW_SOURCE' ,0,Status)
        WriteSourceListStat = WriteSourceListStat +1
     endif
!DEC_OBJ
Status = dal_table_put_col(ptr,'DEC_OBJ',1,DAL_FLOAT,nsou,&
         addrof(tab(1)),Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
             ' dal_table_put_col problem for DEC_OBJ',0,Status)
   WriteSourceListStat = WriteSourceListStat +1
endif

!Y_FIN column
call  PutCol(Ptr,nsou,tab,'Y_FIN',3,&
                  SourceList,Status)
if(Status.ne.ISDC_OK)then
    Status=0
    WriteSourceListStat = WriteSourceListStat +1
endif

!Z_FIN column

call  PutCol(Ptr,nsou,tab,'Z_FIN',4,&
                  SourceList,Status)
if(Status.ne.ISDC_OK)then
   Status=0
    WriteSourceListStat = WriteSourceListStat +1
endif

!FIN_YZ_ERR column
call  PutCol(Ptr,nsou,tab,'FIN_YZ_ERR',8,&
                  SourceList,Status)
if(Status.ne.ISDC_OK)then
   Status=0
    WriteSourceListStat = WriteSourceListStat +1
endif 

!FIN_RD_ERR column
call  PutCol(Ptr,nsou,tab,'FIN_RD_ERR',8,&
                  SourceList,Status)
if(Status.ne.ISDC_OK)then
   Status=0
    WriteSourceListStat = WriteSourceListStat +1
endif 


!FLUX column
call  PutCol(Ptr,nsou,tab,'FLUX',5,&
                  SourceList,Status)
if(Status.ne.ISDC_OK)then
    Status=0
     WriteSourceListStat = WriteSourceListStat +1
endif

if(RawFluxEstim)then
if(type==0)then
!private estim
   tab(1:nsou) = RawFluxesTab(1:nsou)
   Status = dal_table_put_col(ptr,'RAW_FLUX',1,DAL_FLOAT,nsou,&
         addrof(tab(1)),Status)
   if(Status.ne.ISDC_OK)then
      call MESSAGE(procName,&
           ' dal_table_put_col problem for RAW_FLUX',&
           ZeroError,Status)
   endif
endif
endif


!FLUX_ERR column
call  PutCol(Ptr,nsou,tab,'FLUX_ERR',11,&
                  SourceList,Status)
if(Status.ne.ISDC_OK)then
    Status=0
     WriteSourceListStat = WriteSourceListStat +1
endif


flux(1:nsou) = SourceList(1:nsou,5)
fluxerr(1:nsou) = SourceList(1:nsou,11)
locerr(1:nsou) = SourceList(1:nsou,8)
!DETSIG column
call  PutCol(Ptr,nsou,tab,'DETSIG',10,&
                  SourceList,Status)
if(Status.ne.ISDC_OK)then
    Status=0
     WriteSourceListStat = WriteSourceListStat +1
endif
snr(1:nsou) = SourceList(1:nsou,10)

!RA_FIN DEC_FIN column
if(type==0)then !Scw coordinates always in tan proj
   do ns=1,nsou
      valr = 1.
      call ATTI_CONV_YZAD(SourceList(ns,3),SourceList(ns,4),&
           alpha(ns),delta(ns),Scw,1,valr,Status)
      if(status.ne.0)then
         call WAR_MESSAGE(procName,&
              ' Problem  in fine alpha,delta calculations ',0,&
              Status)
         WriteSourceListStat = WriteSourceListStat +1
     endif
  enddo
AllSouList(scw,OutBand)%alpha(1:nsou) = alpha(1:nsou)
AllSouList(scw,OutBand)%delta(1:nsou) = delta(1:nsou)
else !map coordinates depends an ProjType
   do ns=1,nsou
     call PixSkyProj(SourceList(ns,3),SourceList(ns,4),alpha(ns),delta(ns),Status)

      if(status.ne.0)then
         call WAR_MESSAGE(procName,&
              ' Problem  in fine alpha,delta calculations ',0,&
              Status)
         WriteSourceListStat = WriteSourceListStat +1
      endif
   enddo
endif

tab = alpha

Status = dal_table_put_col(ptr,'RA_FIN',1,DAL_FLOAT,nsou,&
         addrof(tab(1)),Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
             ' dal_table_put_col problem for RA_FIN',0,Status)
   WriteSourceListStat = WriteSourceListStat +1

endif
tab = delta
Status = dal_table_put_col(ptr,'DEC_FIN',1,DAL_FLOAT,nsou,&
         addrof(tab(1)),Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
             ' dal_table_put_col problem for DEC_FIN',0,Status)
  WriteSourceListStat = WriteSourceListStat +1
endif

Status = dal3gen_index_update(ptr,IdxPtr,Status)
if(status == -25109) then
    call WAR_MESSAGE(procName,'Source list index updating : keyword missing ',0,&
                Status)
  
    WriteSourceListStat = WriteSourceListStat +1
else
  if(Status.ne.ISDC_OK)then
     call WAR_MESSAGE(procName,'Source list index updating error ',1,Status)
     call WAR_CONT_MESSAGE(procName,'Source list maybe corrupted ',1,Status)
     WriteSourceListStat = WriteSourceListStat +1
    
  endif
endif
ptr = 0

!UPDATING OF FINAL CATALOGUE - OutSourceCat
if(SourceCat==0)then
  call UpdateOutCat(type,Scw,OutBand,nsou,SourceList,alpha,delta,flux,&
                    fluxerr,snr,locerr,Status)
  if(Status.ne.ISDC_OK)return
endif

deallocate(tab,tabi,alpha,delta,flux,fluxerr,snr,locerr,SourceList)

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!=====================================
END SUBROUTINE WriteSourceList
!=====================================


!========================================
SUBROUTINE UpdateIdx(IdxPtr,Status)
!========================================
USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER :: IdxPtr
!OUTPUT VARIABLES
INTEGER            :: Status
!LOCAL VARIABLES
REAL(kind=8)      :: val
character(len=19) ::strutc
CHARACTER(len=20) :: procName,colName

 procName = 'UpdateIdx'
 Status = ISDC_OK
 if(DebugMode.eq.3)&
      call Message(procName,' begin',ZeroError,Status)

end0="&\0"
end1="&\0"

colName = 'TSTART'
val = ogtimes(4)
Status = dal_attribute_put_real&
            (IdxPtr,colName,val,"&","&",Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
             ' dal_attribute_put_real problem for '//colName,&
              0,Status)
endif

colName = 'TSTOP'
val = ogtimes(5)
Status = dal_attribute_put_real&
            (IdxPtr,colName,val,"&","&",Status)
if(Status.ne.ISDC_OK)then
   call WAR_MESSAGE(procName,&
             ' dal_attribute_put_real problem for '//colName,&
              0,Status)
endif

end0="&\0"
end1="&\0"
!!$colName = 'DATE-OBS'
!!$val = ogtimes(4)
!!$Status = dal3gen_convert_ijd2utc(val,0,strUtc,Status)
!!$if(Status.ne.ISDC_OK) then
!!$   write(str250,*,err=181)'idx dal3gen_convert_ijd2utc problem for ',&
!!$        colName,' time : ',val
!!$   goto 182
!!$181 str250 = 'idx dal3gen_convert_ijd2utc problem '&
!!$                    //errorstr//' 181'
!!$182 call WAR_MESSAGE(procName,str250,0,Status)
!!$  
!!$else
!!$   Status = dal_attribute_put_char&
!!$            (IdxPtr,colName,strutc,"&","&",Status)
!!$   if(Status.ne.ISDC_OK) then
!!$      call WAR_MESSAGE(procName,&
!!$           ' dal_attribute_put_char problem for '//colName,0,&
!!$            Status)
!!$   endif
!!$endif
end0="&\0"
end1="&\0"




!!$colName = 'DATE-END'
!!$val = ogtimes(5)
!!$      
!!$Status = dal3gen_convert_ijd2utc(val,0,strUtc,Status)
!!$if(Status.ne.ISDC_OK) then
!!$   write(str250,*,err=183)'idx dal3gen_convert_ijd2utc problem for ',&
!!$        colName,' time : ',val
!!$   goto 184
!!$183 str250 = 'idx dal3gen_convert_ijd2utc problem '&
!!$                    //errorstr//' 183'
!!$184 call WAR_MESSAGE(procName,str250,0,Status)
!!$  
!!$else
!!$   Status = dal_attribute_put_char&
!!$            (IdxPtr,colName,strutc,"&","&",Status)
!!$   if(Status.ne.ISDC_OK) then
!!$      call WAR_MESSAGE(procName,&
!!$           ' dal_attribute_put_char problem for '//colName,0,&
!!$            Status)
!!$   endif
!!$endif
!!$end0="&\0"
!!$end1="&\0"


!========================
END SUBROUTINE UpdateIdx
!========================


!========================================
SUBROUTINE WriteOutCat(OutCatPtr,Status)
!========================================
  USE ISDC
  USE MIMOSA_CONTROL_MODULE
  USE MIMOSA_GLOBVAR_MODULE
  USE MIMOSA_USE_MODULE
  USE IBIS_IMAGING_PARAM
  USE INTER1
  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER :: OutCatPtr
  !OUTPUT VARIABLES
  INTEGER            :: Status

  !LOCAL VARIABLES
  INTEGER ::valint
  INTEGER,parameter :: chardim = MaxSouInGenCat
  INTEGER(kind=2),dimension(:),pointer   :: tabi
  REAL(kind=4) :: sumf,simsig
  REAL(kind=4),dimension(:),pointer      :: tabr
  REAL(kind=4),dimension(:),pointer      :: tabr4
  CHARACTER(len=16),dimension(chardim) :: tab16char
  CHARACTER(len=100),dimension(chardim) :: tab100char
  INTEGER   :: nsou,ns,news,iok,col,rows,ib,it,iler,newsou,numer
  INTEGER   :: StartBin,Endbin, NumValues,type,binsize,varlength
  CHARACTER(len=20) :: procName,colName,dimstr
  CHARACTER(len=6)  :: str
  CHARACTER(len=1)  ::str1
  CHARACTER(len=2)  ::str2
  CHARACTER(len=3)  ::str3
  CHARACTER(len=4)  ::str4


  procName = 'WriteOutCat'
  Status = ISDC_OK
  if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)
!!$print *,procname,' allsoumapflux',allsoumapflux
  write(dimstr,'(I5)')chardim

  if(SourceCat.ne.0)then
      call WAR_Message(procName,' No output Catalogue created ',0,Status)
      return
  endif

  if(.not.associated( OutSourceCat))return ! empty output catalogue


 !EBIN_NUM
  valint = OutEnergyNumber
  end0="&\0"
  end1="&\0"
  Status = dal_attribute_put_int&
            (OutCatPtr,'EBIN_NUM',valint,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      call WAR_MESSAGE(procName,&
           ' dal_attribute_put_int problem for EBIN_NUM',0,Status)
   endif

   end0="&\0"
   end1="&\0"



  
  nsou = size(OutSourceCat,1)
  if(nsou.gt.chardim)then
      call War_Message(procName,&
           'Too many sources  ==> only first '//dimstr//' will be written into output catalog',ZeroError,Status)
      nsou = chardim
  endif

  if(nsou ==0)then
      call war_message(procname,' Empty output source catalogue ',0,status)
      return
   endif

  allocate(tabi(nsou),tabr(nsou),stat=iok)
  if(iok.ne.0)then
     call MESSAGE(procName,'allocation problem.',&
          AllocError,Status)
     return
  endif
  tabi(:)=0
  tabr(:)=0.


  news = nsou-SourceNumber
  if(news.gt.0)then
     write(str,'(I5)',err=186)news
     goto 187
186  str250 = ' ** new sources found in OG '&
                    //errorstr//' 186'
187  call MESSAGE(procName,&
          ' '//str//' new sources found in OG ',&
          ZeroError,Status)
     Status = dal_table_add_rows( OutCatPtr,SourceNumber,news,Status)
     if(Status.ne.ISDC_OK) then
        call WAR_MESSAGE(procName,&
             'Output catalogue writing :cannot add rows',1,Status)
         call WAR_CONT_MESSAGE(procName,&
             'Output catalogue will not be created',1,Status)
        WriteOutCatStat = WriteOutCatStat+1
        return
     endif
     call  PutCatNull(OutCatPtr,SourceNumber+1,nsou,Status)
  endif !news.gt.0


  do ns = 1,nsou
     if(OutSourceCat(ns,1).gt.0) then
        !calcul of moyen alpha,delta
        !for the instant no weighting
        OutSourceCat(ns,2) = OutSourceCat(ns,2)/OutSourceCat(ns,1) 
        OutSourceCat(ns,3) = OutSourceCat(ns,3)/OutSourceCat(ns,1) 
        iler = 0
        do ib=1,OutEnergyNumber
           it = OutSourceCat(ns,4+ib)+OutSourceCat(ns,4+OutEnergyNumber+ib)
           iler = iler+it
           if(it > 0)then
              OutSourceCatFlux(ns,ib) = OutSourceCatFlux(ns,ib)/it
              OutSourceCatFluxErr(ns,ib) = OutSourceCatFluxErr(ns,ib)/it
             
           endif
        enddo
        OutSourceCatSig(ns) = sqrt(OutSourceCatSig(ns))
        if(iler >0)OutSourceCatLocErr(ns) = OutSourceCatLocErr(ns)/iler
     endif
  enddo

  do col=1,4
     tabi = 0
     select case(col)
     case(1)
        !SCW_NUMC 
        ! number of times that the source was in Scw FOV
        tabi(1:SourceNumber) = INT(OutSourceCat(1:SourceNumber,4))
        colName = 'SCW_NUM_C'

     case(2)
        !SCW_NUMF 
        ! number of times that the source was found in Scw FOV
        do ns=1,nsou
           tabi(ns) = INT(sum(OutSourceCat(ns,5:5+OutEnergyNumber-1)))
        enddo
        colName = 'SCW_NUM_F'

     case(3)
        !OG_NUM 
        ! number of times that the source was found in OG
        do ns=1,nsou
           tabi(ns) = &
           INT(sum(OutSourceCat(ns,5+OutEnergyNumber:5+2*OutEnergyNumber-1)))
        enddo
        colName = 'OG_NUM'
     case(4)
        !NEW_SOURCE
        tabi(1:SourceNumber) = 0
        if(nsou.gt.SourceNumber)tabi(SourceNumber+1:nsou)=1
        colName = 'NEW_SOURCE'
     end select

   
     Status = dal_table_put_col&
          (OutCatPtr,colName,1,DAL_USHORT,nsou,&
          addrof(tabi(1)),Status)
     if(Status.ne.ISDC_OK)then
        call WAR_MESSAGE(procName,&
             'dal_table_put_col problem for'//colName ,0,Status)
        WriteOutCatStat = WriteOutCatStat +1
     endif

  enddo

!SOURCE_ID
tab16char(1:SourceNumber) = InCatId(1:SourceNumber)

if(nsou.gt.SourceNumber)then

   do newsou=SourceNumber+1,nsou

      numer = newsou-SourceNumber

      if( numer < 10) then
         
         write(str1,'(i1)')numer
         tab16char(newsou) = 'NEW_'//str1
      else
         if(numer < 100)then
            write(str2,'(i2)')numer
            tab16char(newsou) = 'NEW_'//str2
         else
            if(numer < 1000)then
               write(str3,'(i3)')numer
               tab16char(newsou) = 'NEW_'//str3
            else
               write(str4,'(i4)')numer
               tab16char(newsou) = 'NEW_'//str4
            endif
         endif
      endif
   enddo
  
endif

Status = dal_table_put_col_strings&         
      (OutCatPtr,'SOURCE_ID',1,1,nsou,nsou,tab16char,Status)
if(Status.ne.ISDC_OK)then
  call WAR_MESSAGE(procName,&
  'Cannot put SOURCE_ID ',0,Status)
   WriteOutCatStat = WriteOutCatStat +1
endif

!NAME
tab100char(1:SourceNumber) = InCatName(1:SourceNumber)
if(nsou.gt.SourceNumber)then
   do newsou=SourceNumber+1,nsou
       numer = newsou-SourceNumber
      if( numer < 10) then
         write(str1,'(i1)')numer
         tab100char(newsou) = 'NEW_'//str1
      else
         if(numer < 100)then
            write(str2,'(i2)')numer
            tab100char(newsou) = 'NEW_'//str2
         else
            if(numer < 1000)then
               write(str3,'(i3)')numer
               tab100char(newsou) = 'NEW_'//str3
            else
               write(str4,'(i4)')numer
               tab100char(newsou) = 'NEW_'//str4
            endif
         endif
      endif
   enddo
   

endif

Status = dal_table_put_col_strings&
         (OutCatPtr,'NAME',1,1,nsou,nsou,tab100char,Status)
if(Status.ne.ISDC_OK)then
  call WAR_MESSAGE(procName,&
  'Cannot put NAME into source list',0,Status)
   WriteOutCatStat = WriteOutCatStat +1
endif
tabr(:) = 0.

  do col=1,4
     select case(col)
     case(1)
        !ra_fin
        tabr = OutSourceCat(1:nsou,2)
        colName = 'RA_FIN'
     case(2)
        !dec_fin
        tabr = OutSourceCat(1:nsou,3)
        colName = 'DEC_FIN'
    case(3)
        !signif
       if(nsou >SourceNumber)then
          tabr(SourceNumber+1:nsou) =&
               OutSourceCatSig(SourceNumber+1:nsou)
       endif
       do it=1,SourceNumber
          if(OutSourceCatSig(it) >0.)then
             ! source already found
             tabr(it) = OutSourceCatSig(it)
          else
             sumf = sum(AllSouMapFlux(it,:))
             if(sumf > 0. )then
                ! upper limit from Map
                tabr(it) = sqrt(sum(AllSouMapSnr(it,:)**2))
            endif
         endif
      enddo
      colName = 'DETSIG'
     case(4)
         tabr = OutSourceCatLocErr(1:nsou)*PixAngDim
        colname = 'FIN_RD_ERR'
     end select
    
     Status = dal_table_put_col&
          (OutCatPtr,colName,1,DAL_FLOAT,nsou,&
          addrof(tabr(1)),Status)
     if(Status.ne.ISDC_OK)then
        call WAR_MESSAGE(procName,&
             ' dal_table_put_col problem for '//colName,0,Status)
        WriteOutCatStat = WriteOutCatStat +1
     endif
  enddo



 !flux
 colName = 'FLUX'
 Status = dal_table_get_col_struct(outCatPtr,&
          0,NumValues,ColName,type,binsize,varlength,status)
allocate(tabr4(binsize),stat=iok)
  if(iok.ne.0)then
     call MESSAGE(procName,'allocation problem.',&
          AllocError,Status)
     return
  endif


NumValues = binsize

  do ns = 1,nsou
     StartBin =ns
     Endbin = ns
     tabr4(:) = 0.
     if(ns .le. SourceNumber)then 
        !catalogue source
        do ib=1,OutEnergyNumber
           if(OutSourceCatFlux(ns,ib) >0.)then 
              ! source found somewhere
              tabr4(ib) = OutSourceCatFlux(ns,ib)
           else
              ! upper limit from Map
              tabr4(ib) = AllSouMapFlux(ns,ib)
           endif
       ! print *,procname,ib,ns,OutSourceCatFlux(ns,ib),AllSouMapFlux(ns,ib)
        enddo
     else 
        !new source
        tabr4(1:OutEnergyNumber) = OutSourceCatFlux(ns,1:OutEnergyNumber)
     endif
     Status = DAL_TABLE_PUT_COL_BINS(outCatPtr,&
          ColName,1,DAL_FLOAT,StartBin,EndBin,NumValues,&
          addrof(Tabr4(1)),Status)
     if(Status.ne.ISDC_OK)then
        call MESSAGE(procName,&
                'Cannot put'//colName//' into out cat.',&
              ZeroError,Status)
      
     endif
  enddo




  !FLUX_ERR
  colname = 'FLUX_ERR'
  do ns = 1,nsou
     StartBin =ns
     Endbin = ns
     tabr4(:) = 0.
     tabr4(1:OutEnergyNumber) = OutSourceCatFluxErr(ns,1:OutEnergyNumber)
     Status = DAL_TABLE_PUT_COL_BINS(outCatPtr,&
          ColName,1,DAL_FLOAT,StartBin,EndBin,NumValues,&
          addrof(Tabr4(1)),Status)
     if(Status.ne.ISDC_OK)then
        call MESSAGE(procName,&
                'Cannot put'//colName//' into out cat.',&
              ZeroError,Status)
      
     endif
  enddo


  !E_MIN
  colname = 'E_MIN'
  do ns = 1,nsou
     StartBin =ns
     Endbin = ns
     tabr4(:) = 0.
     tabr4(1:OutEnergyNumber) =OutEnergyBands(1:OutEnergyNumber,1)
     Status = DAL_TABLE_PUT_COL_BINS(outCatPtr,&
          ColName,1,DAL_FLOAT,StartBin,EndBin,NumValues,&
          addrof(Tabr4(1)),Status)
     if(Status.ne.ISDC_OK)then
        call MESSAGE(procName,&
                'Cannot put'//colName//' into out cat.',&
              ZeroError,Status)
      
     endif
  enddo


 !E_MAX
  colname = 'E_MAX'
  do ns = 1,nsou
     StartBin =ns
     Endbin = ns
     tabr4(:) = 0.
     tabr4(1:OutEnergyNumber) =OutEnergyBands(1:OutEnergyNumber,2)
     Status = DAL_TABLE_PUT_COL_BINS(outCatPtr,&
          ColName,1,DAL_FLOAT,StartBin,EndBin,NumValues,&
          addrof(Tabr4(1)),Status)
     if(Status.ne.ISDC_OK)then
        call MESSAGE(procName,&
                'Cannot put'//colName//' into out cat.',&
              ZeroError,Status)
      
     endif
  enddo
  deallocate(tabi,tabr,tabr4)

  OutCatPtr = 0

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!=========================
END SUBROUTINE WriteOutCat
!=========================
 
!#########################################
! SUBROUTINE OF MIMOSA_BASE1_MODULE  
!#########################################

!==============================================================
SUBROUTINE ReadIndexInfo(ShdIdxPtr,EBins,Chans,Status)
!===============================================================
USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE2_MODULE
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER                       :: ShdIdxPtr
REAL(kind=4),dimension(:,:),pointer     :: EBins
INTEGER ,dimension(:,:),pointer :: Chans
INTEGER                       :: Status
!LOCAL VARIABLES
INTEGER  ::nrows,iok,tnum,enum,n,k
REAL(kind=4),dimension(:),pointer :: tabe1,tabe2,tabt1,tabt2,tabind,&
                                     tabc1,tabc2
CHARACTER(len=20) ::procname

Status = ISDC_OK
procName = ' ReadIndexInfo'
if(DebugMode.eq.3)&
  call Message(procName,'begin ',ZeroError,Status)

if(associated(EBins))deallocate(EBins)
if(associated(Chans))deallocate(Chans)

status = dal_table_get_num_rows(ShdIdxPtr,nrows,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' dal_table_get_num_rows problem ' ,&
        IsdcExitStatus,Status)
   return
endif



allocate(tabe1(nrows),tabe2(nrows),tabt1(nrows),&
         tabc1(nrows),tabc2(nrows),tabt2(nrows),tabind(nrows),stat=iok)
IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',&
        AllocError,Status)
   return
endif
tabe1(:) = -1.
tabe2(:) = -1.
tabt1(:) = -1.
tabt2(:) = -1.
tabc1(:) = -1.
tabc2(:) = -1.
tabind(:) = 0.

Status = DAL_TABLE_GET_COL(ShdIdxPtr,'E_MIN',&
        1,DAL_FLOAT,nrows,addrof(tabe1(1)),Status)

if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' DAL_TABLE_GET_COL problem for E_MIN ' ,&
        IsdcExitStatus,Status)
   return
endif
Status = DAL_TABLE_GET_COL(ShdIdxPtr,'E_MAX',&
        1,DAL_FLOAT,nrows,addrof(tabe2(1)),Status)

if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' DAL_TABLE_GET_COL problem for E_MAX ' ,&
        IsdcExitStatus,Status)
   return
endif

Status = DAL_TABLE_GET_COL(ShdIdxPtr,'CHANMIN',&
        1,DAL_FLOAT,nrows,addrof(tabc1(1)),Status)

if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' DAL_TABLE_GET_COL problem for CHANMIN ' ,&
        IsdcExitStatus,Status)
   return
endif
Status = DAL_TABLE_GET_COL(ShdIdxPtr,'CHANMAX',&
        1,DAL_FLOAT,nrows,addrof(tabc2(1)),Status)

if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' DAL_TABLE_GET_COL problem for CHANMAX ' ,&
        IsdcExitStatus,Status)
   return
endif


tabt1(:) = 0.

Status = DAL_TABLE_GET_COL(ShdIdxPtr,'TFIRST',&
        1,DAL_FLOAT,nrows,addrof(tabt1(1)),Status)


if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' DAL_TABLE_GET_COL problem for TFIRST ' ,&
        IsdcExitStatus,Status)
   return
endif

Status = DAL_TABLE_GET_COL(ShdIdxPtr,'TLAST',&
        1,DAL_FLOAT,nrows,addrof(tabt2(1)),Status)

if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' DAL_TABLE_GET_COL problem for TLAST ' ,&
        IsdcExitStatus,Status)
   return
endif

call VerifBins(nrows,tabe1,tabe2,tabind,Ebins,enum,'ENERGY',Status)
if(Status.ne.ISDC_OK)return


allocate(Chans(enum,2),stat=iok)  
 IF (iok /= 0) then
   call MESSAGE(procName,'Allocation problem',&
        AllocError,Status)
   return
endif



k=0
do n=1, nrows  
   if(tabind(n)==1)then
      k=k+1
      Chans(k,1) = tabc1(n)
      Chans(k,2) = tabc2(n)
   endif
enddo

    
deallocate(tabe1,tabe2,tabt1,tabt2,tabc1,tabc2,tabind)


if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!============================
END SUBROUTINE ReadIndexInfo
!============================






!........................................................
subroutine PutCol(Ptr,nsou,tab,colname,colnr,&
                  SourceList,Status)
!........................................................
USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_IMAGING_PARAM
IMPLICIT NONE
!INPUT VARIABLES
INTEGER                   :: Ptr,nsou,colnr
REAL(kind=4),dimension(:),pointer :: tab
REAL(kind=4),dimension(:,:),pointer :: SourceList
CHARACTER(len=*)                  :: colname
!OUTPUT VARIABLES
INTEGER             :: Status
!LOCAL VARIABLES
INTEGER   :: ns
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'PutCol'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)

do ns=1,nsou
    tab(ns) = SourceList(ns,colnr)
enddo 
if(colname=='FIN_RD_ERR')then
    if(FitOffset==0)then
       tab = tab*PixAngDim
    endif
endif
Status = dal_table_put_col(ptr,colname,1,DAL_FLOAT,nsou,&
         addrof(tab(1)),Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
             ' dal_table_put_col problem for '//colname,&
              IsdcExitStatus,Status)
  return
endif

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!.....................
end subroutine PutCol
!......................



!.........................................
subroutine UpdateOutCat(type,Scw,OutBand,&
           nsou,SourceList,alpha,delta,flux,fluxerr,snr,locerr,Status)
!........................................
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE
!INPUT VARIABLES
INTEGER             ::type,Scw,OutBand, nsou
REAL(kind=4),dimension(:,:),pointer :: SourceList
REAL(kind=4),dimension(:),pointer::  alpha,delta,flux,fluxerr,snr,locerr

!OUTPUT VARIABLES
INTEGER             :: Status
!LOCAL VARIABLES
REAL (kind=4), dimension(:,:), pointer :: temp
REAL (kind=4), dimension(:), pointer :: temp1,temp4
REAL (kind=4), dimension(:,:), pointer :: temp2,temp3
INTEGER   :: ns,k,nr,inr,outsou,iok,nn,i
REAL(kind=4)     :: dist,radidentif,mindist
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'UpdateOutCat'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)

radidentif = PixelAngDim*CleanParTab(6)
outsou = size(OutSourceCat,1)
!SourceNumber - number of sources in input catalogue

do ns=1,nsou

  ! k - number in input catalogue
  k = int(SourceList(ns,1))
   
   if(k.eq.0) then
      !nonidentified source
      nr = SourceNumber+1
      mindist = radidentif
      inr = 0
      !search for registered sources

      do while(nr .le.outsou)
         if(OutSourceCat(nr,1).gt.0)then
           ! source already seen in FOV
            dist =  sqrt((OutSourceCat(nr,2)/OutSourceCat(nr,1)-alpha(ns))**2+&
                    (OutSourceCat(nr,3)/OutSourceCat(nr,1)-delta(ns))**2)
           if(dist .lt.radidentif)then
              if(dist .lt.mindist) then
                 mindist = dist
                 inr = nr
              endif 
           endif  
         endif 
        nr = nr+1
     enddo

      if(inr.ne.0)then
         ! not catalogue source but already found
         nr = inr
      else
       ! new source
      
      if(associated(OutSourceCat))then
         allocate(temp(outsou,NumSouDes),temp1(outsou),&
              temp2(outsou,OutEnergyNumber),&
              temp3(outsou,OutEnergyNumber),temp4(outsou),stat=iok)
         if(iok.ne.0)then
            call MESSAGE (procName,' aloocation problem ',&
                 AllocError,Status)
            return
         endif
         temp = OutSourceCat
         temp2 = OutSourceCatFlux
         temp3 = OutSourceCatFluxErr
         temp1 = OutSourceCatSig 
         temp4 = OutSourceCatLocerr 
         deallocate(OutSourceCat,OutSourceCatFlux,OutSourceCatSig,&
              OutSourceCatLocerr)
      endif !associated(OutSourceCat)

      outsou = outsou+1
      allocate(OutSourceCat(outsou,NumSouDes),&
               OutSourceCatFlux(outsou,OutEnergyNumber),&
               OutSourceCatFluxErr(outsou,OutEnergyNumber),&
               OutSourceCatSig(outsou),&
               OutSourceCatLocerr(outsou),stat=iok)
      if(iok.ne.0)then
         call MESSAGE (procName,' aloocation problem ',&
                      AllocError,Status)
         return
      endif
      OutSourceCat(:,:) = 0.
      OutSourceCatFlux(:,:) = 0.
       OutSourceCatFluxErr(:,:) = 0.
      OutSourceCatSig(:) = 0.
      OutSourceCatLocerr(:) = 0.
      if(associated(OutSourceCat))then
         OutSourceCat(1:outsou-1,1:NumSouDes) = &
              temp(1:outsou-1,1:NumSouDes)
         OutSourceCatFlux(1:outsou-1,:) = &
              temp2(1:outsou-1,:)
         OutSourceCatFluxErr(1:outsou-1,:) = &
              temp3(1:outsou-1,:)
         OutSourceCatSig(1:outsou-1) = &
              temp1(1:outsou-1)
         OutSourceCatLocerr(1:outsou-1) = &
              temp4(1:outsou-1)
         deallocate(temp,temp1,temp2,temp3,temp4)
      endif
      nr = outsou
      endif !new source
  else
  !input catalogue source
   nr= InScwCat(10,k)
  endif 

 OutSourceCat(nr,1) = OutSourceCat(nr,1)+1.

 ! for the moment unweighted mean : sum calculation  here
 ! of alpha,delta,flux,fluxerr,locerr
 ! will be divised by number of measures in WriteOutCat

 OutSourceCat(nr,2) = OutSourceCat(nr,2)+alpha(ns)
 OutSourceCat(nr,3) =  OutSourceCat(nr,3)+delta(ns)
 nn = 4+type*OutEnergyNumber

 OutSourceCat(nr,nn+OutBand) = OutSourceCat(nr,nn+OutBand)+ 1.
 OutSourceCatFlux(nr,Outband) = OutSourceCatFlux(nr,Outband)+flux(ns)
 OutSourceCatFluxErr(nr,Outband) = OutSourceCatFluxErr(nr,Outband)+&
      fluxerr(ns)
 OutSourceCatLocerr(nr) = OutSourceCatLocerr(nr)+&
      locerr(ns)
 if(type==0)then ! scw source list
    OutSourceCatSig(nr) = OutSourceCatSig(nr)+snr(ns)**2
 else ! og source list
    if(OutSourceCatSig(nr).eq.0.)then
       OutSourceCatSig(nr) = OutSourceCatSig(nr)+snr(ns)**2
    endif
    !  include mosaicking snr
    !   only if mosa source 
 endif
! ns - input number of source
! nr - number in OutSourceCat
enddo


if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!......................
end subroutine UpdateOutCat
!......................
 
!........................................................
SUBROUTINE NewIndexMember&
          (IdxPtr,OutBand,StrucName,MemberPtr,Status)
!........................................................

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER             :: IdxPtr,OutBand
CHARACTER(len=*)            :: StrucName
!OUTPUT VARIABLES
INTEGER             :: MemberPtr
INTEGER             :: Status

!LOCAL VARIABLES
REAL(kind=8)   :: val
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'newIndexMember'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)


!creates new index member
Status = dal3gen_index_create_member(IdxPtr,&
          StrucName,'&',MemberPtr,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
   ' dal3gen_index_create_member problem for'//StrucName//' creating ',&
              IsdcExitStatus,Status)
  return
endif

Status = COMMON_STAMP_OBJECT(Memberptr,StampComment,Status)
if(Status.ne.ISDC_OK)then
    call WAR_MESSAGE(procName,'COMMON_STAMP_OBJECT problem  ',0,Status)
endif



!adds energy band info
val = OutEnergyBands(OutBand,1)
Status = dal_attribute_put_real(MemberPtr,"E_MIN",val,&
         "&","&",Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
             ' dal_attribute_put_real problem for E_MIN',&
              IsdcExitStatus,Status)
             return
endif
end0="&\0"
end1="&\0"

val = OutEnergyBands(OutBand,2)
Status = dal_attribute_put_real(MemberPtr,"E_MAX",val,&
         "&","&",Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
             ' dal_attribute_put_real problem for E_MAX',&
              IsdcExitStatus,Status)
             return
endif
end0="&\0"
end1="&\0"

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!............................
END SUBROUTINE NewIndexMember
!............................


!........................................................
SUBROUTINE AddIndexMemberArray&
          (ObjPtr,IdxPtr,scw,OutBand,ImaFlag,ArrName,ArrType,ArrTypeVal,Arr,&
           MemberPtr,Status)
!........................................................
USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE2_MODULE

IMPLICIT NONE

!INPUT VARIABLES
INTEGER                     ::ObjPtr, IdxPtr,scw,OutBand
!OutBand - number of output energy band
INTEGER                             :: ImaFlag
! if yes then write attitude info into the header
! and some additional keywords
CHARACTER(len=*)                    :: ArrName,ArrType,ArrTypeVal
REAL(kind=4),dimension(:,:),pointer :: Arr 
!OUTPUT VARIABLES
INTEGER                     :: MemberPtr
INTEGER                     :: Status

!LOCAL VARIABLES
INTEGER  :: idim,jdim,type,numAxes
INTEGER,dimension(1:DAL_MAX_ARRAY_DIMENSION)::axes 
INTEGER,dimension(2)::startVal,endVal
CHARACTER (len=20)  ::colname
CHARACTER(len=20) :: procName

Status = ISDC_OK
procName = 'AddIndexMemberArray'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)

!creates new index member
!adds energy band info
call  NewIndexMember&
          (IdxPtr,OutBand,ArrName,MemberPtr,Status)
if(Status.ne.ISDC_OK) return



!filles keyword in header
call WriteAttributes(ObjPtr,MemberPtr,scw,OutBand,&
     ImaFlag,ArrType,ArrTypeVal,&
           Status)
if(Status.ne.ISDC_OK) return

idim = SIZE(Arr,1)
jdim = SIZE(Arr,2)

if(ImaFlag.gt.0)then
  ! add attitude info
  call  WriteAttiHeader(imaflag,scw,idim,jdim,MemberPtr,Status)
  if(Status.ne.ISDC_OK) return
endif

!array structure verification

Status = dal_array_get_struct(MemberPtr,type,numAxes,axes,&
            Status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
       ' dal_array_get_struct problem ',&
       IsdcExitStatus,Status)
   return
endif
if((type.ne.DAL_FLOAT).or.(numAxes.ne.2))then
    call MESSAGE(procName,&
        ' type or numaxes of array invalid ',&
       InvalArrayError,Status)
     return
else
 if(axes(1).ne.0)then
   !predefined array size
    if((axes(1).ne.idim).or.(axes(2).ne.jdim))then
       call MESSAGE(procName,&
       ' array size invalid ',&
       InvalArrayError,Status)
      return
    endif
 else
   !free array size
     
     axes(1) = idim
     axes(2) = jdim
     Status = dal_array_mod_struct(MemberPtr,&
                      type,numAxes,axes,Status)
     if(Status.ne.ISDC_OK) then
        IsdcExitStatus=Status
        call MESSAGE(procName,&
       ' dal_array_mod_struct problem',IsdcExitStatus,Status)
       return
     endif

endif
endif
 
!writing
startVal = 1
endVal(1) = axes(1)
endVal(2) = axes(2)
! Status = dal_array_set_null_values(MemberPtr,numAxes,startVal,endVal,status) 
Status = dal_array_put_section&
         (MemberPtr,numAxes,startVal,endVal,type,&
          addrof(Arr(1,1)),status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
  call MESSAGE(procName,&
       ' dal_array_put_section problem ',&
       IsdcExitStatus,Status)
  return
endif

if(ProjType==0)then

   colname='LONGPOLE'
   Status = dal_attribute_delete (MemberPtr,colName,Status)
   status=0
   if(Status.ne.ISDC_OK) then
       status=0
       write(str250,*) 'cannot delete attribute '//colname
      call MESSAGE(procName,str250,ZeroError,Status)
      endif
   colname='LATPOLE'
   Status = dal_attribute_delete (MemberPtr,colName,Status)
   if(Status.ne.ISDC_OK) then
      status=0
      write(str250,*) 'cannot delete attribute '//colname
      call MESSAGE(procName,str250,ZeroError,Status)
   endif
endif

!UPDATING OF INDEX TABLE
Status = dal3gen_index_update(MemberPtr,IdxPtr,Status)

select case(Status)
case(ISDC_OK)
    Status = 0

case(DAL3GEN_INDEX_DOESNT_MATCH)
    call MESSAGE(procName,'DAL3GEN_INDEX_DOESNT_MATCH error '&
                // ArrName//ArrType//ArrTypeVal,&
         IndexError,Status)
    return

case(DAL3GEN_INDEX_KEY_NOT_FOUND)
   if(MissKeyStat==0)&
     call WAR_MESSAGE(procName,&
          'DAL3GEN_INDEX_KEY_NOT_FOUND error   ',0,Status)
   MissKeyStat=MissKeyStat+1
   Status = 0
case(DAL3GEN_MEMBER_NOT_ATTACHED)
   call MESSAGE(procName,'DAL3GEN_MEMBER_NOT_ATTACHED error for'&
                // ArrName//ArrType//ArrTypeVal,&
         IndexError,Status)
    return
case default
   call MESSAGE(procName,'unrecognized index update  error for'&
                // ArrName//ArrType//ArrTypeVal,&
         IndexError,Status)
    return
end select




if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!.................................
END SUBROUTINE AddIndexMemberArray
!.................................


!........................................................
SUBROUTINE CreateClear(ParentPtr,ElemType,ClearFlag,&
                       StrucNameExt,StrucNameInt,FitsName,ElemPtr,&
                       Status)
!........................................................

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE

IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
INTEGER  ,Intent(in)   :: ElemType
! if ElemType == 0 then element else index
INTEGER  ,Intent(in)   :: ClearFlag
! 0 - no clearing of existing element ( or index)
! 1 - clears rows(members) of element(index) if exist
character(len=*),Intent(in)   :: StrucNameExt
! for ordinary element its structure name
! for index            its index structure name
character(len=*),Intent(in)   :: StrucNameInt 
! for ordinary element not used
! for index            its member structure name
character(len=*),Intent(in)   :: FitsName  
! name of fits file of the element
!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
!LOCAL VARIABLES
character(len=300) :: name
character(len=100) :: str100
character(len=30)  :: tplName
character(len=20)  :: procName,str
INTEGER   :: len,num,n,ptr,nFound,num1,num2
LOGICAL            :: exist

Status = ISDC_OK
procName = 'Createclear'
if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)

!VERIF WHETHER ELEMENT ALREADY EXIST
Status = dal_object_find_element(ParentPtr,StrucNameExt,&
                                   ElemPtr,Status)
if(Status.ne.ISDC_OK)then
  exist = .false.
  Status = 0
else
  exist = .true.
endif

if(.not.exist) then
   !CREATION OF NEW ELEMENT
   !number of children before
   Status=dal_element_get_num_children(ParentPtr,num1,Status)
    if(Status.ne.ISDC_OK)then
          call MESSAGE(procName,&
               ' dal_element_get_num_children problem ',&
              FindElemError,Status)
          return
   endif
   ! parent element location
   Status=dal_element_get_location(ParentPtr,name,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
           ' dal_element_get_location problem ',&
           IsdcExitStatus,Status)
      return
   endif
   len=index(name,'/',back=.true.)
   name = name(1:len)
   name = name(1:len_trim(name))//FitsName
   
   ! element creation
   tplName = StrucNameExt(1:len_trim(StrucNameExt))//'.tpl'
   status = dal_object_create(name,&
             DAL_DISK,tplName,ElemPtr,Status)
   if(Status .ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,' dal_object_create problem ',&
             IsdcExitStatus,Status)
      return
   endif
   !attach new element
   Status = dal3gen_element_attach(ElemPtr,ParentPtr,&
                                      Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
              ' dal3gen_element_attach problem ',&
             IsdcExitStatus,Status)
      return
   endif
    !number of children after
   Status=dal_element_get_num_children(ParentPtr,num2,Status)
    if(Status.ne.ISDC_OK)then
       IsdcExitStatus=Status
       call MESSAGE(procName,&
            ' dal_element_get_num_children problem ',&
            IsdcExitStatus,Status)
       return
   endif
   if(num2.ne.(num1+1))then
      write(str250,*,err=200)&
      ' Number of elements after element attaching: ',&
      num2,' not correct. Before was : ',num1
      goto 201
200   str250 = 'Number of elements after element attaching not correct'&
                    //errorstr//' 200'
201      call MESSAGE(procName,str250,IsdcProcError,Status)
      return
   endif
else
!ELEMENT ALREADY EXIST
if(ClearFlag ==1 )then
!clearing of existing element

if(ElemType ==0)then
   !ordinary element
   Status = dal_table_get_num_rows(ElemPtr,num,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
       call MESSAGE(procName,&
            'dal_table_get_num_rows  problem',&
             IsdcExitStatus,Status)
       return
   endif
   if(num.ne.0)then
   !rows deleting in ordinary element
      Status = dal_table_del_rows(ElemPtr,1,num,Status)
      if(Status.ne.ISDC_OK)then
         IsdcExitStatus=Status
         call MESSAGE(procName,&
              'dal_table_del_rows problem ',&
              IsdcExitStatus,Status)
         return
      endif
      Status = dal_table_get_num_rows(ElemPtr,num,Status)
      if(Status.ne.ISDC_OK)then
         IsdcExitStatus=Status
          call MESSAGE(procName,&
               'dal_table_get_num_rows  problem',&
               IsdcExitStatus,Status)
          return
      endif    
     if(num.ne.0)then
         call MESSAGE(procName,&
             'element not cleared correctly ',&
              ClearElemError,Status)
     return
     endif
  endif !rows deleting in ordinary element
else  ! element exists and of index type
   Status = dal_table_get_num_rows(ElemPtr,num,Status)
   ! num - number of members in the index
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
          'dal_table_get_num_rows  problem',&
           IsdcExitStatus,Status)
     return
  endif
  if(num.ne.0)then
      !index clearing
      do n=1,num
         call OpenIndexMember&
             (ElemPtr,StrucNameInt,1,ptr,nFound,Status)
         if(Status.ne.ISDC_OK)return
         write(str,'(I7)',err=202)nfound
         goto 203
202      str = '*'&
                    //errorstr//' 202'
203      if(nFound.ne.1)then
            call MESSAGE(procName,&
                'number of members'//str//' .ne. 1 ',&
                 NumMembersError,Status)
            return
         endif
         Status = dal_element_detach(ptr,ElemPtr,Status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,&
                 ' dal_element_detach problem  ',&
                  IsdcExitStatus,Status)
            return
         endif
         Status = dal_object_close(ptr,DAL_DELETE,Status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,&
                ' dal_object_close problem  ',&
                  IsdcExitStatus,Status)
           return
        endif
      enddo ! end of index clearing
   Status = dal_table_get_num_rows(ElemPtr,num,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
          'dal_table_get_num_rows  problem',&
            IsdcExitStatus,Status)
     return
   endif
  if(num.ne.0)then
     call MESSAGE(procName,&
          'index not cleared correctly  problem',&
           ClearIndexError,Status)
     return
 endif
 endif !index clearing
endif ! element exists and of index type
endif ! clearing of existing element
endif !element exist

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!..........................
END SUBROUTINE CreateClear
!..........................



!........................................................
SUBROUTINE GetClearElement(ParentPtr,StrucName,ElemPtr,&
                       Status)
!........................................................

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
character(len=*),Intent(in)   :: StrucName 

!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
!LOCAL VARIABLES

character(len=20)  :: procName
INTEGER   :: num

Status = ISDC_OK
procName = ' GetClearElement'
if(DebugMode.eq.3)&
  call Message(procName,'begin',ZeroError,Status)

!VERIF WHETHER ELEMENT ALREADY EXIST
Status = dal_object_find_element(ParentPtr,StrucName,&
                                   ElemPtr,Status)
if(Status.ne.ISDC_OK)then
  call MESSAGE(procName,&
               ' Element '//StrucName//' not found',&
              FindElemError,Status)
  return
endif
Status = dal_table_get_num_rows(ElemPtr,num,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
            'dal_table_get_num_rows  problem for '//StrucName,&
             IsdcExitStatus,Status)
   return
endif
if(num.ne.0)then
!rows deleting 
   Status = dal_table_del_rows(ElemPtr,1,num,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
     call MESSAGE(procName,&
              'dal_table_del_rows problem for'//StrucName,&
              IsdcExitStatus,Status)
     return
   endif
   Status = dal_table_get_num_rows(ElemPtr,num,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
               'dal_table_get_num_rows  problem 2 for '//StrucName,&
               IsdcExitStatus,Status)
      return
   endif    
   if(num.ne.0)then
      call MESSAGE(procName,&
             'element '//StrucName//' not cleared correctly ',&
              ClearElemError,Status)
      return
   endif

endif 

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!..........................
END SUBROUTINE  GetClearElement
!..........................




!........................................................
SUBROUTINE ClearDetachElement(ParentPtr,StrucName,ElemPtr,&
                       Status)
!........................................................

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
character(len=*),Intent(in)   :: StrucName 

!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
!LOCAL VARIABLES
character(len=20)  :: procName
INTEGER   :: num

Status = ISDC_OK
procName = ' ClearDetachElement'
if(DebugMode.eq.3)&
  call Message(procName,'begin',ZeroError,Status)

!VERIF WHETHER ELEMENT ALREADY EXIST
Status = dal_object_find_element(ParentPtr,StrucName,&
                                   ElemPtr,Status)
if(Status.ne.ISDC_OK)then
  Status = ISDC_OK
  return
endif

Status = dal_table_get_num_rows(ElemPtr,num,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
            'dal_table_get_num_rows  problem for '//StrucName,&
             IsdcExitStatus,Status)
   return
endif
if(num.ne.0)then
!rows deleting 
   Status = dal_table_del_rows(ElemPtr,1,num,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
     call MESSAGE(procName,&
              'dal_table_del_rows problem for'//StrucName,&
              IsdcExitStatus,Status)
     return
   endif
   Status = dal_table_get_num_rows(ElemPtr,num,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
               'dal_table_get_num_rows  problem 2 for '//StrucName,&
               IsdcExitStatus,Status)
      return
   endif    
   if(num.ne.0)then
      call MESSAGE(procName,&
             'element '//StrucName//' not cleared correctly ',&
              ClearElemError,Status)
      return
   endif

endif 

Status = dal_element_detach(ElemPtr,Parentptr,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' dal_element_detach problem  ',&
        IsdcExitStatus,Status)
   return
endif
Status = dal_object_close(elemptr,DAL_DELETE,Status)



if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!..........................
END SUBROUTINE ClearDetachElement
!..........................


!........................................................
SUBROUTINE CreateClearIndex(ParentPtr,StrucNameExt,StrucNameInt,&
                            FitsName,ElemPtr,Status)
!........................................................

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
character(len=*),Intent(in)   :: StrucNameExt
! for ordinary element its structure name
! for index            its index structure name
character(len=*),Intent(in)   :: StrucNameInt 
! for ordinary element not used
! for index            its member structure name
character(len=*),Intent(in)   :: FitsName  
! name of fits file of the element
!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
!LOCAL VARIABLES
character(len=300) :: name
character(len=100) :: str100
character(len=30)  :: tplName
character(len=20)  :: procName,str
INTEGER   :: len,num,n,ptr,nFound,num1,num2
LOGICAL            :: exist

Status = ISDC_OK
procName = 'CreateclearIndex'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)

!VERIF WHETHER ELEMENT ALREADY EXIST
Status = dal_object_find_element(ParentPtr,StrucNameExt,&
                                   ElemPtr,Status)
if(Status.ne.ISDC_OK)then
  exist = .false.
  Status = 0
else
  exist = .true.
endif

if(.not.exist) then
   !CREATION OF NEW ELEMENT
   !number of children before
   Status=dal_element_get_num_children(ParentPtr,num1,Status)
    if(Status.ne.ISDC_OK)then
       IsdcExitStatus=Status
       call MESSAGE(procName,&
            ' dal_element_get_num_children problem ',&
            IsdcExitStatus,Status)
       return
   endif
   ! parent element location
   name = ' '
   Status=dal_element_get_location(ParentPtr,name,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
           ' dal_element_get_location problem ',&
           IsdcExitStatus,Status)
      return
   endif
   len=index(name,'/',back=.true.)
   name = name(1:len)
   name = name(1:len_trim(name))//ScwResultsPath(1:len_trim(ScwResultsPath))//FitsName
  
   ! element creation
   tplName = StrucNameExt(1:len_trim(StrucNameExt))//'.tpl'

   status = dal_object_create(name,&
             DAL_DISK,tplName,ElemPtr,Status)
   if(Status .ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,' dal_object_create problem '//name//tplname,&
             IsdcExitStatus,Status)
      return
   endif
   !attach new element
   Status = dal3gen_element_attach(ElemPtr,ParentPtr,&
                                      Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
              ' dal3gen_element_attach problem ',&
             IsdcExitStatus,Status)
      return
   endif
    !number of children after
   Status=dal_element_get_num_children(ParentPtr,num2,Status)
    if(Status.ne.ISDC_OK)then
       IsdcExitStatus=Status
          call MESSAGE(procName,&
               ' dal_element_get_num_children problem ',&
              IsdcExitStatus,Status)
          return
   endif
   if(num2.ne.(num1+1))then
      write(str250,*,err=204)&
      ' Number of elements after element attaching: ',&
      num2,' not correct. Before was : ',num1
      goto 205
204   str250 = 'Number of elements after element attaching  not correct'&
                    //errorstr//' 204'
205   call MESSAGE(procName,str250,IsdcProcError,Status)
      return
   endif
else
  !ELEMENT ALREADY EXIST clearing 
  Status = dal_table_get_num_rows(ElemPtr,num,Status)
   ! num - number of members in the index
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
          'dal_table_get_num_rows  problem',&
           IsdcExitStatus,Status)
     return
  endif
  if(num.ne.0)then
     write(str250,*)' Clearing of :',num,' elements'
     call message(procname,str250,Zeroerror,Status)
      !index clearing
      do n=1,num
         call OpenIndexMember&
             (ElemPtr,StrucNameInt,1,ptr,nFound,Status)
         if(Status.ne.ISDC_OK)return
         write(str,'(I6)',err=206)nfound
         goto 207
206       str = '*'&
                    //errorstr//' 207'
207         if(nFound.ne.1)then
            call MESSAGE(procName,&
                'number of members'//str//' .ne. 1 ',&
                 NumMembersError,Status)
            return
         endif
         Status = dal_element_detach(ptr,ElemPtr,Status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,&
                 ' dal_element_detach problem  ',&
                  IsdcExitStatus,Status)
            return
         endif
         Status = dal_object_close(ptr,DAL_DELETE,Status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,&
                ' dal_object_close problem  ',&
                  IsdcExitStatus,Status)
           return
        endif
      enddo ! end of index clearing


! commented 14.03.05 : after detaching from index
!  num still equal number of elements

!!$   Status = dal_table_get_num_rows(ElemPtr,num,Status)
!!$   if(Status.ne.ISDC_OK)then
!!$      call MESSAGE(procName,&
!!$          'dal_table_get_num_rows  problem',&
!!$            IsdcExitStatus,Status)
!!$     return
!!$   endif
!!$  if(num.ne.0)then
!!$     call MESSAGE(procName,&
!!$          'index not cleared correctly  problem',&
!!$           ClearIndexError,Status)
!!$     return
!!$ endif
 endif !index clearing
endif !index exist

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!..........................
END SUBROUTINE CreateClearIndex
!..........................

!........................................................
SUBROUTINE GetClearIndex(ParentPtr,StrucNameExt,StrucNameInt,ElemPtr,Status)
!........................................................

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  ,Intent(in)   :: ParentPtr
character(len=*),Intent(in)   :: StrucNameExt
! index structure name
character(len=*),Intent(in)   :: StrucNameInt
! index member structure name
!OUTPUT VARIABLES
INTEGER ,Intent(out)  :: ElemPtr
INTEGER               :: Status
!LOCAL VARIABLES
character(len=20)  :: procName,str
INTEGER   :: num,n,ptr,nFound

Status = ISDC_OK
procName = 'GetClearIndex'
if(DebugMode.eq.3)&
  call Message(procName,' begin',ZeroError,Status)


!VERIF WHETHER ELEMENT ALREADY EXIST
Status = dal_object_find_element(ParentPtr,StrucNameExt,&
                                   ElemPtr,Status)
if(Status.ne.ISDC_OK)then
   if(.not.FastOpenFlag)then
      ! element should be already initialised by CommonPreparePars
      call MESSAGE(procName,&
           ' dal_object_find_element problem for '//StrucNameExt,&
           FindElemError,Status)
   else
      ! element not found but it is OK with FastOpenFlag
      Status = FindElemError
   endif
   return
endif

Status = dal_table_get_num_rows(ElemPtr,num,Status)
! num - number of members in the index
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
      call MESSAGE(procName,&
          'dal_table_get_num_rows  problem',&
           IsdcExitStatus,Status)
     return
 endif

if(num.ne.0)then
      !index clearing
      do n=1,num
         call OpenIndexMember&
             (ElemPtr,StrucNameInt,1,ptr,nFound,Status)
         if(Status.ne.ISDC_OK)return
         write(str,'(I6)',err=208)nfound
         goto 209
208       str = '*'&
                    //errorstr//' 208'
209      if(nFound.ne.1)then
            call MESSAGE(procName,&
                'number of members'//str//' .ne. 1 ',&
                 NumMembersError,Status)
            return
         endif
         Status = dal_element_detach(ptr,ElemPtr,Status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,&
                 ' dal_element_detach problem  ',&
                  IsdcExitStatus,Status)
            return
         endif
         Status = dal_object_close(ptr,DAL_DELETE,Status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,&
                ' dal_object_close problem  ',&
                  IsdcExitStatus,Status)
           return
        endif
      enddo ! end of index clearing
   Status = dal_table_get_num_rows(ElemPtr,num,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
          'dal_table_get_num_rows  problem',&
            IsdcExitStatus,Status)
     return
   endif
  if(num.ne.0)then
     call MESSAGE(procName,&
          'index not cleared correctly  problem',&
           ClearIndexError,Status)
     return
 endif
 endif !index clearing

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!..........................
END SUBROUTINE GetClearIndex
!..........................

!...........................................
SUBROUTINE OpenIndexMember&
           (IdxPtr,MemberName,MemberNum,MemberPtr,&
             nFound,Status)
!...........................................

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   :: IdxPtr       ! pointer to index
CHARACTER(len=*)  :: MemberName   ! MEMBER NAME
INTEGER   :: MemberNum    ! member row number

!OUTPUT VARIABLES
INTEGER   :: MemberPtr !pointer to member
INTEGER   :: nFound
INTEGER   :: Status

!LOCAL VARIABLES
INTEGER   :: num
character(len=10) ::row,rowpom
character(len=20) :: procName

Status = ISDC_OK
procName = ' OpenIndexMember'
if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)

status = dal3gen_index_get_num_members(IdxPtr,MemberName,num,status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   str250 = ' dal3gen_index_get_num_members problem for'//MemberName
   call message(procname,str250,IsdcExitStatus,Status)
   return
endif

write(row,'(i8)',err=210)MemberNum
goto 211
210 write(str250,*) ' cannot write MemberNum ',MemberNum,' into selection line'
    call message(procname,str250,StringWriteError,Status)
    return

211 row = adjustl(row)
   
    rowpom = row(1:len_trim(row))
    row = '#row=='//rowpom
status = dal3gen_index_find_member(IdxPtr,&
            MemberName,row,nFound,MemberPtr,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
                ' dal3gen_index_find_member problem for'//MemberName//' selection :'//row,&
                  IsdcExitStatus,Status)
   return
endif


if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!..................................
end subroutine OpenIndexMember
!..................................




!..................................
SUBROUTINE ReadArray(ArrPtr,iDim,jDim,Arr,Status)
!..................................

USE ISDC
USE MIMOSA_CONTROL_MODULE
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

procName = 'ReadArray'
Status = ISDC_OK
if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)

 
! array size verification
Status = dal_array_get_struct(ArrPtr,type,numAxes,axes,&
            Status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,' dal_array_get_struct problem ',&
       IsdcExitStatus,Status)
   return
endif

if((type.ne.DAL_FLOAT).or.(numAxes.ne.2).or.&
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
Status = dal_array_get_section(ArrPtr,numAxes,startVal,endVal,&
         type,numValues,addrof(Arr(1,1)),status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,' dal_array_get_section problem ',&
         IsdcExitStatus,Status)
   return
endif


if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!.........................
END SUBROUTINE ReadArray
!.........................

!..................................
SUBROUTINE ReadUnknownArray(ArrPtr,iDim,jDim,Arr,Status)
!..................................

USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER :: ArrPtr,iDim,jDim
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: Arr
INTEGER            :: Status

!LOCAL VARIABLES
INTEGER,dimension(1:DAL_MAX_ARRAY_DIMENSION)::axes
INTEGER :: type,numAxes,numValues
INTEGER,dimension(2)::startVal,endVal
CHARACTER(len=20) :: procName

procName = 'ReadUnknownArray'
Status = ISDC_OK
if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)

 
! array size verification
Status = dal_array_get_struct(ArrPtr,type,numAxes,axes,&
            Status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,' dal_array_get_struct problem ',&
       IsdcExitStatus,Status)
   return
endif

if((type.ne.DAL_FLOAT).or.(numAxes.ne.2))then
       call MESSAGE(procName,' array invalid ',&
       InvalArrayError,Status)
   return
endif

idim = axes(1)
jdim = axes(2)

if(associated(arr))deallocate(Arr)
allocate(Arr(idim,jdim),stat=status)
if(Status.ne.ISDC_OK) then
   call MESSAGE(procName,&
        ' allocation problem ',AllocError,Status)
   return
endif
 
!reading
startVal = 1
endVal(1) = axes(1)
endVal(2) = axes(2)
numValues = axes(1)*axes(2)
Status = dal_array_get_section(ArrPtr,numAxes,startVal,endVal,&
         type,numValues,addrof(Arr(1,1)),status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,' dal_array_get_section problem ',&
         IsdcExitStatus,Status)
   return
endif


if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!.........................
END SUBROUTINE ReadUnknownArray
!.........................

!====================================================
SUBROUTINE ReadShdAttributes(MemberNum,MemberPtr,eMinVal,&
            eMaxVal,Status)
!====================================================

USE ISDC
USE DAL3GEN_F90_API
USE MIMOSA_USE_MODULE  
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   :: MemberNum,MemberPtr       ! pointer to index member
!OUTPUT VARIABLES
REAL(kind=4)      :: eMinVal,EMaxVal
INTEGER   :: Status

!LOCAL VARIABLES
REAL(kind=8)      ::val
INTEGER   :: col,timeproblem,errstat
character(len=DAL_MAX_ATTRIBUTE_SIZE) :: str30
character(len=300) :: errstr
character(len=20) :: procName,colName


Status = ISDC_OK
procName = 'ReadShdAttributes'


errstr = ' Reading problem for :'
errstat = 0
if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)

!ENERGY BAND READING
do col=1,2 
   select case(col)
   case(1) !E_MIN
      colName = 'E_MIN'
   case(2)  !E_MAX
      colName ='E_MAX'
   end select
   
   Status = dal_attribute_get_real&
        (MemberPtr,colName,val,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      IsdcExitStatus=Status
      write(str250,*,err=212) ' dal_attribute_get_real problem for ',colName
      goto 213
212   str250 = 'dal_attribute_get_real problem'&
                    //errorstr//' 212'
213   call MESSAGE(procName,str250,IsdcExitStatus,Status)
      return
   endif
 
   end0="&\0"
   end1="&\0"
 
   select case(col)
   case(1) !E_MIN
      eMinVal = val
      eMin = eMinVal
   case(2)  !E_MAX
      eMaxVal = val
      eMax = eMaxVal
   end select
enddo


!OTHER ATTRIBUTES READ FROM cexp
do col=6,8
   colname =  ShdTimeColNames(col)
   !ONTIME,EXPOSURE,DEADC
   
   Status = dal_attribute_get_real&
        (MemberPtr,colName,val,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      errstat = errstat+1
      errstr =errstr(1:len_trim(errstr))//' '//colname(1:len_trim(colname))//'( set to 1)'
      Status = ISDC_OK
      val = 1.
   endif
 
   end0="&\0"
   end1="&\0"
 
   select case(col)
   case(6) 
      shdontime = val
   case(7) 
      shdexposure = val
   case(8)
      shddeadc = val
   end select
enddo




!TIME READING
timeproblem=0
do col=1,5  
   colName =  ShdTimeColNames(col) 
   !TFIRST,TLAST,TELAPSE,TSTART,TSTOP
  
   Status = dal_attribute_get_real&
        (MemberPtr,colName,val,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      IsdcExitStatus=Status
      if(col.le.3)then !TFIRST,TLAST,TELAPSE reading error
         if(TimeInfo==0)then ! coherence with attitude data
            ! times can be set to scw times
            if(timeproblem == 0)then
               errstat = errstat+1
               errstr =errstr(1:len_trim(errstr))//' '//&
                    colname(1:len_trim(colname))//&
                    '( set to Scw value)'
               Status = ISDC_OK
               timeproblem=1
            endif!timeproblem == 0
         else ! cannot reset TFIRST,TLAST,TELAPSE
            write(str250,*,err=215) &
                 ' dal_attribute_get_real problem 2 for ',colName,MemberNum
            goto 216
215         str250 = 'dal_attribute_get_real problem 2'&
                    //errorstr//' 215'
216         call MESSAGE(procName,str250,IsdcExitStatus,Status)
            return
         endif
      else  !TSTOP,TSTART reading error
         if(timeproblem < 2)then
             errstat = errstat+1
             errstr =errstr(1:len_trim(errstr))//' '//&
                    colname(1:len_trim(colname))//&
                    '( set to Scw value)'
             Status = ISDC_OK
            
            timeproblem = 2 
         endif
        
      endif !TSTOP,TSTART reading error
   endif ! reading error
  end0="&\0"
  end1="&\0"

  select case(col)
  case(1)  !TFIRST
     if(timeproblem.ne.1)then
        shdtfirst = val
     else
        shdtfirst = scwtstart
     endif
  case(2)  !TLAST
     if(timeproblem.ne.1)then
        shdtlast = val
     else
        shdtlast = scwtstop
     endif
  case(3)  !TELAPSE
     if(timeproblem.ne.1)then
        shdtelapse = val
     else
        shdtelapse = scwtelapse
     endif
  case(4) !TSTART
     if(timeproblem <2)then
        shdtstart = val
     else
        shdtstart = scwtstart
     endif
  case(5) !TSTOP
     if(timeproblem <2)then
        shdtstop = val
     else
        shdtstop = scwtstop
     endif
   end select

enddo

colname = 'BKGPARAM'
Status = dal_attribute_get_char(MemberPtr,colName,BkgFile,end0,end1,Status)
if(Status.ne.ISDC_OK) then
   call war_MESSAGE(procName,&
        'Cannot read BKGPARAM - will be set to none ',0,&
        Status)
   BkgFile='NONE'
endif

end0="&\0"
end1="&\0"


colName = 'SHD_TYPE'
str30 = ' '
Status = dal_attribute_get_char(MemberPtr,colName,str30,end0,end1,Status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' dal_attribute_get_char problem for'//colName,IsdcExitStatus,&
        Status)
   return
endif

end0="&\0"
end1="&\0"

shd_type = str30(1 :len_trim(str30))

if(errstat >0)then
   call war_message(procname,errstr,0,status)
endif
if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!================================
end subroutine ReadShdAttributes
!================================


!===========================================================
SUBROUTINE ReadImaAttributes(MemberPtr,eMinVal,&
            eMaxVal,Status)
!===========================================================

USE ISDC
USE DAL3GEN_F90_API
USE MIMOSA_USE_MODULE  
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   ::MemberPtr       ! pointer to index member
!OUTPUT VARIABLES
REAL(kind=4)      :: eMinVal,EMaxVal
INTEGER   :: Status

!LOCAL VARIABLES
REAL(kind=8)      ::val
INTEGER   :: col,valint
character(len=20) :: procName,colName
character(len=DAL_MAX_ATTRIBUTE_SIZE) :: str30

Status = ISDC_OK
procName = ' ReadImaAttributes'

if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)

!ENERGY BAND READING
do col=1,2 
   select case(col)
   case(1) !E_MIN
      colName = 'E_MIN'
   case(2)  !E_MAX
      colName ='E_MAX'
   end select
   
   Status = dal_attribute_get_real&
        (MemberPtr,colName,val,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      IsdcExitStatus=Status
      write(str250,*,err=217) ' dal_attribute_get_real problem for ',colName
      goto 218
217   str250 = 'dal_attribute_get_real problem'&
                    //errorstr//' 217'
218   call MESSAGE(procName,str250,IsdcExitStatus,Status)
      return
   endif
 
   end0="&\0"
   end1="&\0"
 
   select case(col)
   case(1) !E_MIN
      eMinVal = val
      eMin = eMinVal
   case(2)  !E_MAX
      eMaxVal = val
      eMax = eMaxVal
   end select
enddo



do col=1,2   
   select case(col)
   case(1) !CHANMIN
      colName = 'CHANMIN'
   case(2) !CHANMAX
      colName = 'CHANMAX'
   end select
      Status = dal_attribute_get_int&
           (MemberPtr,colName,valint,end0,end1,Status)
      if(Status.ne.ISDC_OK) then
         call WAR_MESSAGE(procName,&
              ' dal_attribute_get_int problem for '//colName//' set to 0',0,Status)
         valint = 0.
      endif
 
      end0="&\0"
      end1="&\0"

      select case(col)
      case(1) !CHANMIN
         chanmin=valint 
      case(2) !CHANMAX
         chanmax=valint
      end select
   enddo



 colName = 'IMATYPE'

 str30 = ' '
 Status = dal_attribute_get_char&
      (MemberPtr,colName,str30,end0,end1,Status)
 if(Status.ne.ISDC_OK) then
    call WAR_MESSAGE(procName,&
         ' dal_attribute_get_char problem for'//colName//'set to 0 str',0,&
         Status)
    str30 = " 0"
 endif

 end0="&\0"
 end1="&\0"

 shd_type = str30(1 :len_trim(str30))

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!============================
end subroutine ReadImaAttributes
!============================


!====================================================
SUBROUTINE ReadScwAttributes(ScwPtr,Status)
!====================================================

USE ISDC
USE DAL3GEN_F90_API
USE MIMOSA_USE_MODULE  
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   :: ScwPtr       ! pointer to scw
INTEGER   :: Status

!LOCAL VARIABLES

REAL(kind=8)      ::val
INTEGER   :: col,valint
character(len=100):: str100
character(len=20) :: procName,colName
character(len=DAL_MAX_ATTRIBUTE_SIZE) :: str30

Status = ISDC_OK
procName = ' ReadScwAttributes'

if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)

do col=1,3   
   
   select case(col)
   case(1)  !TSTART
      colName ='TSTART'
   case(2)  !TSTOP
      colName ='TSTOP'
   case(3)  !TELAPSE
      colName ='TELAPSE'
   case(4)  !TEXPOSURE
      colName ='TEXPOSURE'
   end select
   
    Status = dal_attribute_get_real&
            (ScwPtr,colName,val,end0,end1,Status)
  if(Status.ne.ISDC_OK) then
     IsdcExitStatus=Status
      if(col.le.2)then !fatal
         write(str250,*,err=219) ' dal_attribute_get_real problem for ',colName
         goto 220
219      str250 = 'dal_attribute_get_real problem '&
                    //errorstr//' 219'
220      call MESSAGE(procName,str250,IsdcExitStatus,Status)
         return
      else           !non fatal
          write(str250,*,err=221) ' dal_attribute_get_real problem for ',colName,&
                          'set to 0.'
          goto 222
221       str250 = ' dal_attribute_get_real problem - set to 0.'&
                    //errorstr//' 221'
222       call WAR_MESSAGE(procName,str250,0,Status)
         val = 0.
      endif
  endif

  end0="&\0"
  end1="&\0"

   select case(col)
   case(1)  !TSTART
      scwtstart = val
   case(2)  !TSTOP
      scwtstop = val
   case(3)  !TELAPSE
      scwtelapse = val
   case(4)  !EXPOSURE
      scwexposure = val
   end select

enddo

 !REVOL
 colName = 'REVOL'
 
  Status = dal_attribute_get_int&
            (ScwPtr,colName,valint,end0,end1,Status)
  if(Status.ne.ISDC_OK) then
       call WAR_MESSAGE(procName,&
           ' dal_attribute_get_int problem for '//colName//' set to 0',0,Status)
       valint = 0.
  endif
  end0="&\0"
  end1="&\0"

 revol=valint
 

do col=1,3
   select case(col)
   case(1) !EXPID
     colName = 'EXPID'
   case(2) !SWID
      colName = 'SWID'
   case(3) !SW_TYPE
      colName = 'SW_TYPE'
  end select
  str30 = ' '
   Status = dal_attribute_get_char&
            (ScwPtr,colName,str30,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      call WAR_MESSAGE(procName,&
           ' dal_attribute_get_char problem for'//colName//'set to 0 str',0,&
            Status)
      str30 = " 0"
   endif

   end0="&\0"
   end1="&\0"

   select case(col)
   case(1) !EXPID
      expid = str30(1 :len_trim(str30))
   case(2) !SWID
      swid = str30(1 :len_trim(str30))
   case(3) !SW_TYPE
      sw_type = str30(1 :len_trim(str30))
  
  end select


enddo

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!============================
end subroutine ReadScwAttributes
!============================


!############################################
! SUBROUTINES FROM MIMOSA_BASE2_MODULE
!############################################

!..................................................
SUBROUTINE  WriteAttiHeader(imatype,num,nx,ny,ImaPtr,Status)
!..................................................
USE ISDC
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER  :: imatype ! 1 - scw, 2 - map
INTEGER  :: num    !pointing number
INTEGER  :: nx,ny  ! image dimensions
INTEGER  :: ImaPtr !image pointer
!OUTPUT VARIABLES
INTEGER  :: Status

!LOCAL VARIABLES
INTEGER  :: col,ix,iy
REAL(kind=4)      :: alpha,delta,cenx,ceny,l,b
!SPR 4357 declaration below changes to real*8 
REAL(KIND=8)      ::th,th_rad,st0,ct0 !ANTICLOCKWISE NORTH ANGLE
REAL(KIND=8)      ::incr1,incr2 ! PIXEL INCREMENT IN DEG
REAL(kind=8)      :: val
CHARACTER(len=20) :: ctyp, procName,colName


procName = 'WriteAttiHeader'
Status = ISDC_OK
if(DebugMode.eq.3)&
  call Message(procName,' begin ',ZeroError,Status)


!  RA DEC OF POINTING       
alpha = pointing_table(num,1)
delta = pointing_table(num,2)
l=reaL(LX(num))
b=real(BX(num))

th = pointing_table(num,3)
th_rad = th*sdeg_to_rad
ct0 = cos(th_rad)
st0 = sin(th_rad)

ix = (1+nx)/2
iy = (1+ny)/2



if(ImaType==1)then
   cenx = ix+Xdisi !- 0.729076
   ceny =iy+xDisj  ! +1.60576
elsE
   cenx = ix+Mapdisi !- 0.729076
   ceny =iy+MapDisj !+1.60576
endif

!PIXEL INCREMENT IN DEG
!incr1 = IOSphere*srad_to_deg*atan(pixel_ratio)
incr1 = srad_to_deg*atan(pixel_ratio)
incr2 = srad_to_deg*atan(pixel_ratio)


do col=1,8
   select case(col)
   case(1) !CRVAL1 - REFERENCE VALUE
       colName = 'CRVAL1'
       if(EquaGal==0)then
          val = alpha
       else
          val = l
       endif
       
   case(2)!CRPIX1 -  X COORDINATE OF IMAGE CENTER
       colName = 'CRPIX1'
       val = cenx
  
   case(3)! CRVAL2 : REFERENCE VALUE
       colName = 'CRVAL2'
       if(EquaGal==0)then
          val = delta
       else
          val = b
       endif
   case(4)!CRPIX2  - Y COORDINATE OF IMAGE CENTER
       colName = 'CRPIX2'
       val = ceny
 
   case(5) !CD1_1 CDELT1*cos(CROTA2)
       colName = 'CD1_1   '
      
       val=-incr1*ct0
      
   case(6) !CD1_2  -CDELT2*sin(CROTA2)
       colName = 'CD1_2   '
     
       val=-incr2*st0
      
   case(7) !CD2_1 CDELT1*sin(CROTA2)
       colName = 'CD2_1   '
       
       val=-incr1*st0
       
   case(8) !CD2_2 CDELT2*cos(CROTA2)
       colName = 'CD2_2   '
    
       val=incr2*ct0
    
   end select
if(val==0.0)val = 0.
Status = dal_attribute_put_real(ImaPtr,colName,val,"&","&",Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
             ' Cannot put attribute  '//colName,&
              IsdcExitStatus,Status)
   return
endif
end0="&\0"
end1="&\0"
enddo


! DESCRIPTION OF THE FIRST AXIS
! CTYPE1 : x AXIS TYPE
if(imatype==1)then
   ! scw image always in tan projection
   ctyp = 'RA---TAN'
else
   !final map can be in car projection
   select case (ProjType)
   case (0)
      if(EquaGal==0)then
         ctyp = 'RA---CAR' 
      else
         ctyp = 'GLON-CAR'
      endif
   case (1)
      ctyp = 'RA---TAN'          
   end select
endif

Status = dal_attribute_put_char(ImaPtr,"CTYPE1",ctyp,"&","&",Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,'Cannot put attribute CTYPE1 = '//ctyp,&
                IsdcExitStatus,Status)
   return
endif

end0="&\0"
end1="&\0"

! DESCRIPTION OF THE SECOND AXIS
! CTYPE 2 : y AXIS TYPE
if(imatype==1)then
   ! scw image always in tan projection
   ctyp = 'DEC--TAN'  
else
   !final map can be in car projection
   select case (ProjType)
   case (0)
      if(EquaGal==0)then
         ctyp = 'DEC--CAR' 
      else
         ctyp = 'GLAT-CAR' 
      endif
   case (1)
      ctyp = 'DEC--TAN'         
   end select
endif
Status = dal_attribute_put_char(ImaPtr,"CTYPE2",ctyp,"&","&",Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,'Cannot put attribute CTYPE2 = '//ctyp,&
                IsdcExitStatus,Status)
   return
endif

end0="&\0"
end1="&\0"


if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!.................................
END SUBROUTINE  WriteAttiHeader
!.................................

!======================================================
SUBROUTINE WriteAttributes(ObjPtr,MemberPtr,scw,OutBand,&
     ImaFlag,ImaType,ImaTypeVal,&
           Status)
!====================================================

USE ISDC
USE DAL3GEN_F90_API
USE MIMOSA_USE_MODULE  
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER   ::   ObjPtr,MemberPtr ,scw      ! pointer to index member
INTEGER   :: OutBand         !out energy band number
INTEGER           :: ImaFlag
! 0 - groupped shd
! 1 - image
! 2 - carte
character(len=*)  :: imaType,imaTypeVal
!OUTPUT VARIABLES
INTEGER   :: Status

!LOCAL VARIABLES

INTEGER   :: col,valint,k,k0
REAL(kind=8)      ::val
character(len=100) :: str100
character(len=15)  :: str15
character(len=20) :: colName,procName
character(len=30) :: str30
character(len=19) ::strutc0,strutc


Status = ISDC_OK
procName = ' WriteAttributes'
if(DebugMode.eq.3)&
  call Message(procName,'begin ',ZeroError,Status)



end0="&\0"
end1="&\0"


!!$Status = dal3gen_attribute_copy(ObjPtr,MemberPtr,'OGID',Status)
!!$if(Status.ne.ISDC_OK) then
!!$   call WAR_MESSAGE(procName,&
!!$        'cannot copy keyword OGID to the object ',0,Status)
!!$endif

strutc0 = '0000-00-00T00:00:00'

do col=1,2
   select case(col)
   case(1)  !CHANMIN
      valint = chanmin
      colName ='CHANMIN'
   case(2)  !CHANMAX
      valint = chanmax
      colName ='CHANMAX'
 
   end select
   Status = dal_attribute_put_int&
            (MemberPtr,colName,valint,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      call WAR_MESSAGE(procName,&
           ' dal_attribute_put_int problem for '//colname,0,Status)
   endif
   end0="&\0"
   end1="&\0"

enddo



select case(ImaTypeVal(1:2))
   case('IN')
     str15 = 'counts/sec'
    case('DE')
     str15 = 'counts/sec'
    case('VA')
     str15 = '(counts/sec)**2'
    case('RE')
     str15 = 'counts/sec'
    case('SI')
     str15 = 'no units'
     case('EF')
     str15 = 'no units'
      case('EX')
     str15 = 'sec'
end select


do col=1,7
   val = ogtimes(col)
   colName = ogTimeColNames(col)
   Status = dal_attribute_put_real&
        (MemberPtr,colName,val,"&","&",Status)
   if(Status.ne.ISDC_OK) then
      write(str250,*,err=231) &
           ' dal_attribute_put _real problem for a carte',colName
      goto 232
231   str250 = 'dal_attribute_put _real problem for a carte'&
           //errorstr//' 231'
232   call WAR_MESSAGE(procName,str250,0,Status)
   endif
   end0="&\0"
   end1="&\0"
enddo
do col=1,3
   select case(col)
  
   case(1) ! E_MIN
      val = OutEnergyBands(OutBand,1)
      colname = 'E_MIN'
   case(2) ! E_MAX
      val = OutEnergyBands(OutBand,2)
      colname = 'E_MAX'
   case(3)  !E_MEAN
      val = 0.5*(OutEnergyBands(OutBand,1)+OutEnergyBands(OutBand,2))
      colName ='E_MEAN'
   end select

  Status = dal_attribute_put_real&
            (MemberPtr,colName,val,"&","&",Status)
  if(Status.ne.ISDC_OK) then
      write(str250,*,err=233) &
      ' dal_attribute_put _real problem for ',colName
      goto 234
233   str250 = '  dal_attribute_put _real problem'&
                    //errorstr//' 233'
234   call WAR_MESSAGE(procName,str250,0,Status)
  endif
  end0="&\0"
  end1="&\0"
enddo

  colname = 'MOSASPR'
  Status = dal_attribute_put_int&
            (MemberPtr,colName,PixSpread,"&","&",Status)
  if(Status.ne.ISDC_OK) then
       write(str250,*,err=433) &
      ' dal_attribute_put int problem for ',colName
      goto 434
433   str250 = '  dal_attribute_put _int problem'&
                    //errorstr//' 433'
434   call WAR_MESSAGE(procName,str250,0,Status)
endif
end0="&\0"
end1="&\0"




do col=1,2
   select case(col)
   case(1)  !CHANMIN
      valint = chanmin
     colName ='CHANMIN'
   case(2)  !CHANMAX
      valint = chanmax
     colName ='CHANMAX'
 
   end select
   Status = dal_attribute_put_int&
            (MemberPtr,colName,valint,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      call WAR_MESSAGE(procName,&
           ' dal_attribute_put_int problem for '//colname,0,Status)
   endif
   end0="&\0"
   end1="&\0"

enddo



select case(ImaFlag)
case(0) !group shd
    k0 = 1
    k=1
case(1) !image
    k0 = 1
    k=2
case(2) !carte
    k0 = 1
    k=4
end select
end0="&\0"
end1="&\0"

do col=k0,k
select case(col)
 
   case(1) !IMATYPE
      str30 =ImaTypeVal
      colName = ImaType
    case(2) !BUNIT
      str30 =str15
      colName = 'BUNIT'
    ! for carte only  
    case(3) !DATE-OBS
      colName = 'DATE-OBS'
      val = ogtimes(4)
      Status = dal3gen_convert_ijd2utc(val,0,strUtc,Status)
      if(Status.ne.ISDC_OK) then
         str30 = strutc0
         TimeConvertStat = TimeConvertStat + 1
        ! if(TimeConvertStat==1)then
            write(str250,*,err=235)'dal3gen_convert_ijd2utc problem for ',&
            colName,' time : ',val,' set to ',strutc0
            goto 236
235         str250 = 'dal3gen_convert_ijd2utc problem -time' &
                    //errorstr//' 235'
236            call WAR_MESSAGE(procName,str250,0,Status)
        ! endif
         Status = ISDC_OK
!!$         write(str250,*)'test print :bad conversion  DATE-OBS : ',val
!!$          call message(procname,str250,0,Status)
      else
           str30 = strUtc         
      endif
     case(4) !DATE-END
      colName = 'DATE-END'
      val = ogtimes(5)
      Status = dal3gen_convert_ijd2utc(val,0,strUtc,Status)
      if(Status.ne.ISDC_OK) then
           str30 = strutc0
           TimeConvertStat = TimeConvertStat + 1
           if(TimeConvertStat==1)then
             write(str250,*,err=237)'dal3gen_convert_ijd2utc problem for ',&
                 colName,' time : ',val,' set to ',strutc0
             goto 238
237          str250 = 'dal3gen_convert_ijd2utc problem-time'&
                    //errorstr//' 237'
238          call WAR_MESSAGE(procName,str250,0,Status)
           endif
           Status = ISDC_OK
      else
           str30 = strUtc
      endif
   end select

Status = dal_attribute_put_char&
            (MemberPtr,colName,str30,"&","&",Status)

if(Status.ne.ISDC_OK) then
      PutCharStat = PutCharStat+1
      if(PutCharStat.le.1)then
         call WAR_MESSAGE(procName,&
           ' dal_attribute_put_char problem for '//colName,0,&
            Status)
      else
        Status = ISDC_OK
      endif
endif
end0="&\0"
end1="&\0"

enddo

if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!==============================
end subroutine WriteAttributes
!==============================

!--------------------------------------------------------
SUBROUTINE VerifBins(nrows,tab1,tab2,tabind,bins,bnum,str,Status)
!---------------------------------------------------------
USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER  :: nrows,bnum,iok
CHARACTER(len = *) :: str
REAL(kind=4),dimension(:),pointer :: tab1,tab2,tabind
REAL(kind=4),dimension(:,:),pointer ::bins
INTEGER  :: Status

!LOCAL VARIABLES
INTEGER  :: n
REAL(kind=4) :: val1,val2,tol
REAL(kind=4),dimension(:,:),pointer :: pomtab
LOGICAL :: jest
CHARACTER(len=20) :: procname

Status = ISDC_OK
procName = 'VerifBins'
if(DebugMode.eq.3)&
  call Message(procName,'begin ',ZeroError,Status)

tol = 0.0001
! input band verif
do n=1,nrows
   if((tab1(n).lt.-tol).or.(tab2(n).lt.-tol))then
      write(str250,*,err=239)'For ',str,' bins , row ',n,', lower or upper limit < 0'
      goto 240
239    str250 = 'For ** bins , row ** lower or upper limit < 0'&
                    //errorstr//' 239'
240   call MESSAGE(procname,str250,DataError,Status)
     
      return
   else
      if(tab1(n).gt.tab2(n)+tol)then
         write(str250,*,err= 242)'For ',str,' bins , row ',n,', lower limit  ge upper limit'
         goto 243
242      str250 = 'For ** bins , row ** lower limit  ge upper limit'&
                    //errorstr//' 242'
243         call MESSAGE(procname,str250,DataError,Status)
        
         return
      endif
   endif
enddo

bnum = 0

val1 = -1.
val2 = -1.
tabind = 0.
jest = .true.

do while(jest) ! if new band 
   jest = .false.
   n=1
   do while((n <= nrows).and.(.not.jest))
      if(tab1(n).ge.0.)then
         val1 = tab1(n)
         val2 = tab2(n)
         jest=.true.
         tabind(n) = 1
      endif
  n=n+1
   enddo
   if(jest)then ! new band adding 
     
      do n=1,nrows ! erase same band
         if((tab1(n).eq.val1).and.(tab2(n).eq.val2))then
         
           
            tab1(n) = -1.
            tab2(n) = -1.
         endif
      enddo ! erase same band
         if(bnum > 0)then ! bins already exists - copy
            
            allocate(pomtab(bnum,2),stat=iok)
            IF (iok /= 0) then
               call MESSAGE(procName,'Allocation problem',&
                    AllocError,Status)
               return
            endif
            pomtab = bins
            deallocate(bins)
         endif ! bins already exists - copy
         bnum = bnum+1
         Allocate(bins(bnum,2),stat=iok)
           
         IF (iok /= 0) then
            call MESSAGE(procName,'Allocation problem',&
                 AllocError,Status)
            return
         endif
         if(bnum > 1)then 
            bins(1:bnum-1,1:2) = pomtab
            deallocate(pomtab)
         endif
         bins(bnum,1) = val1
         bins(bnum,2) = val2
         
    endif ! new band adding 
enddo !if new band


if(DebugMode.eq.3)&
  call Message(procName,'end ',ZeroError,Status)
!-------------------------
END SUBROUTINE VerifBins
!-------------------------



!------------------------------------
SUBROUTINE PutCatNull(OutCatPtr,n1,n2,Status)
!-----------------------------------
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT VARIABLES
INTEGER :: OutCatPtr,n1,n2,Status
!LOCAL VARIABLES
INTEGER :: col
CHARACTER(len=20) :: colname,procname

procname = 'PutCatNull'
Status = ISDC_OK

do col=1,22
   select case(col)
   case(1)
      colname = 'RA_OBJ'
   
   case(2)
      colname = 'DEC_OBJ'
   case(3)
      colname = 'DAY_ID'
   case(4)
      colname = 'CLASS'
   case(5)
      colname = 'ERR_RAD'
   case(6)
      colname = 'RELDIST'
   case(7)
      colname = 'SPA_MODL'
   case(8)
      colname = 'SPA_NPAR'
   case(9)
      colname = 'SPE_MODL'
   case(10)
      colname = 'SPE_NPAR'
   case(11)
      colname = 'VAR_MODL'
   case(12)
      colname = 'VAR_NPAR'
   case(13)
      colname = 'COMMENTS'
  case(14)
      colname = 'SPI_FLUX_1'
  case(15)
      colname = 'SPI_FLUX_2'
   case(16)
      colname = 'ISGR_FLUX_1'
  case(17)
      colname = 'ISGR_FLUX_2'
   case(18)
      colname = 'PICS_FLUX_1'
  case(19)
      colname = 'PICS_FLUX_2'
   case(20)
      colname = 'JEMX_FLUX_1'
  case(21)
      colname = 'JEMX_FLUX_2'
   case(22)
      colname = 'SEL_FLAG'
   end select


   Status = dal_table_set_null_values(OutCatPtr,colname,0,n1,n2,Status)
   if(Status.ne.ISDC_OK)then
      call war_message(procname,'Cannot put NULL into'//colname,0,status)
   endif

enddo

!---------------------------
END SUBROUTINE PutCatNull
!---------------------------



!=============================================================
SUBROUTINE CreateElem(OutElemFile,StrucName,ElemPtr,Status)
!===============================================================
USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER        :: ElemPtr
CHARACTER(len=*)       :: StrucName,OutElemFile
INTEGER               :: Status
character(len=20) :: procname

procname = 'CreateElem'
Status = ISDC_OK

if(len_trim(OutElemFile) > 0)then
   ! not empty DOL of the file

   status = dal_object_create(OutElemFile,&
        DAL_DISK,trim(StrucName)//'.tpl',ElemPtr,Status)
   if(Status.ne.ISDC_OK)then
      IsdcExitStatus=Status
      call MESSAGE(procName,&
           ' cannot create output element'//StrucName,&
           IsdcExitStatus,Status)
      return
   endif
   Status = COMMON_STAMP_OBJECT(ElemPtr,StampComment,Status)
   if(Status.ne.ISDC_OK)then
      call WAR_MESSAGE(procName,'COMMON_STAMP_OBJECT problem  ',0,Status)
   endif
endif
!========================
END SUBROUTINE CreateElem
!========================

!=============================================================
SUBROUTINE ReplaceElement(GrpPtr,ElemPtr,StrucName,&
                       OutElemFile,Outname,clearflag,Status)
!===============================================================

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE2_MODULE
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER        :: GrpPtr,ElemPtr
CHARACTER(len=*)       :: StrucName,OutElemFile,Outname
LOGICAL :: clearFlag
INTEGER               :: Status
!LOCAL VARIABLES
INTEGER        :: ptr,num,i
character(len=20)  :: procName
procName = 'ReplaceElement '
if(DebugMode.eq.3)&
  call Message(procName,'begin',ZeroError,Status)


if(len_trim(OutElemFile) > 0)then
   ! not empty DOL of the file

   if(ChangeGroupFlag)then
      ! detach old and attach new element to the OG

      if(ClearFlag)then
         ! old element cleaned and detached
         call ClearDetachElement(GrpPtr,StrucName,ptr,&
              Status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,&
                 ' cannot detach old output element  ',&
                 IsdcExitStatus,Status)
            return
         endif
      else
         ! old element detached only
         Status = dal_object_find_element(GrpPtr,StrucName,&
                                   ptr,Status)
         if(Status.ne.ISDC_OK)then
            ! old element not found
            Status = ISDC_OK
         else
            Status = dal_element_detach(ptr,Grpptr,Status)
            if(Status.ne.ISDC_OK)then
               IsdcExitStatus=Status
               call MESSAGE(procName,&
                    ' dal_element_detach problem  ',&
                    IsdcExitStatus,Status)
               return
            endif
            Status = dal_object_close(ptr,DAL_SAVE,Status)

            if(Status.ne.ISDC_OK)then
               IsdcExitStatus=Status
               call MESSAGE(procName,&
                    ' cannot close old output element  ',&
                    IsdcExitStatus,Status)
               return
            endif
         endif ! old element found
         endif !.not.ClearFlag
      Status = dal3gen_element_attach(ElemPtr,GrpPtr,&
                                      Status)
      if(Status.ne.ISDC_OK)then
         IsdcExitStatus=Status
         call MESSAGE(procName,&
              ' dal3gen_element_attach problem ',&
              IsdcExitStatus,Status)
         return
      endif

   endif !ChangeGroupFlag

 else !len_trim(OutCatFile) = 0

    !Opening/Clearing OF linked  OUTPUT element
    if(clearFlag)then
       call GetClearElement(GrpPtr,StrucName,ElemPtr,&
            Status)
       if(Status.ne.ISDC_OK)then
          if(Status == FindElemError)then
             ! no element in OG nor defined by the user
             ! standart will be created
             Status = ISDC_OK
             call CreateElem(OutName,StrucName,ElemPtr,Status)
             if(Status.ne.ISDC_OK)return
             call message(procname,' Created std element:'//OutName,0,status)
             Status = COMMON_STAMP_OBJECT(ElemPtr,StampComment,Status)
             if(Status.ne.ISDC_OK)then
                call WAR_MESSAGE(procName,&
                     'COMMON_STAMP_OBJECT problem  ',0,Status)
             endif

             if(ChangeGroupFlag)then
                Status = dal3gen_element_attach(ElemPtr,GrpPtr,&
                                      Status)
                if(Status.ne.ISDC_OK)then
                   IsdcExitStatus=Status
                   call MESSAGE(procName,&
                        ' dal3gen_element_attach problem ',&
                        IsdcExitStatus,Status)
                   return
                endif
             endif !ChangeGroupFlag
          else 
             return
          endif

       endif!Opening/Clearing OF linked element

    else
      call message(procname,' Searching for '//Strucname//' in OG ...',0,status)
      Status = dal_object_find_element(GrpPtr,Strucname,ElemPtr,Status) 
       if(Status.ne.ISDC_OK)then
          IsdcExitStatus=Status
            call MESSAGE(procName,&
                 'Element not found  in the OG : '//Strucname,&
                 IsdcExitStatus,Status)
         return
      endif 
   endif !not clearflag
endif

!=============================================================
END SUBROUTINE ReplaceElement
!===============================================================


!=============================================================
SUBROUTINE ReplaceIdx(GrpPtr,ElemPtr,StrucName,StrucNameInt,&
                       OutElemFile,OutName,clearflag,Status)
!===============================================================

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE2_MODULE
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER        :: GrpPtr,ElemPtr
CHARACTER(len=*)       :: StrucName,StrucNameInt,OutElemFile,OutName
LOGICAL :: clearFlag
INTEGER               :: Status
!LOCAL VARIABLES
INTEGER        :: ptr,num,i
character(len=20)  :: procName
procName = 'ReplaceIdx '
if(DebugMode.eq.3)&
  call Message(procName,'begin',ZeroError,Status)


if(len_trim(OutElemFile) > 0)then
   ! not empty DOL of the file

 
   if(ChangeGroupFlag)then

      call GetClearIndex(GrpPtr,StrucName,StrucNameInt,Ptr,Status)
      if(Status.eq.FindElemError)then
         ! nothing found in the OG, new elem can be attached
         Status = ISDC_OK
         Status = dal3gen_element_attach(ElemPtr,GrpPtr,Status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,&
                 ' dal3gen_element_attach problem ',&
                 IsdcExitStatus,Status)
            return
         endif
      else
         if(Status.ne.ISDC_OK)return
         ! found in the OG
         if(ElemPtr.ne.Ptr)then
            ! not the same file
            Status = dal_element_detach(Ptr,GrpPtr,Status)
            if(Status.ne.ISDC_OK)then
               IsdcExitStatus=Status
               call MESSAGE(procName,' dal_element_detach problem  ',&
                    IsdcExitStatus,Status)
               return
            endif
            Status = dal_object_close(ptr,DAL_DELETE,Status)
 
            if(Status.ne.ISDC_OK)then
               IsdcExitStatus=Status
               call MESSAGE(procName,' cannot detach old output element  ',&
                    IsdcExitStatus,Status)
               return
            endif
        
            Status = dal3gen_element_attach(ElemPtr,GrpPtr,&
                                      Status)
            if(Status.ne.ISDC_OK)then
               IsdcExitStatus=Status
               call MESSAGE(procName,&
                    ' dal3gen_element_attach problem ',&
                    IsdcExitStatus,Status)
               return
            endif
         endif! not the same file
     endif ! found in the OG
   endif !ChangeGroupFlag
 else !len_trim(OutCatFile) = 0

    !Opening/Clearing OF linked  OUTPUT element
    if(clearFlag)then
       call GetClearElement(GrpPtr,StrucName,ElemPtr,&
            Status)
         if(Status.ne.ISDC_OK)then
          if(Status == FindElemError)then
             ! no element in OG nor defined by the user
             ! standart will be created
             Status = ISDC_OK
             call CreateElem(OutName,StrucName,ElemPtr,Status)
             if(Status.ne.ISDC_OK)return
             call message(procname,' Created std element:'//OutName,0,status)
    
             Status = COMMON_STAMP_OBJECT(ElemPtr,StampComment,Status)
             if(Status.ne.ISDC_OK)then
                call WAR_MESSAGE(procName,&
                     'COMMON_STAMP_OBJECT problem  ',0,Status)
             endif
             if(ChangeGroupFlag)then
                Status = dal3gen_element_attach(ElemPtr,GrpPtr,&
                     Status)
                if(Status.ne.ISDC_OK)then
                   IsdcExitStatus=Status
                   call MESSAGE(procName,&
                        ' dal3gen_element_attach problem ',&
                        IsdcExitStatus,Status)
                   return
                endif
             endif !ChangeGroupFlag
          else !Status .ne. FindElemError
             return
          endif
          endif!Opening/Clearing OF linked element
    else
      Status = dal_object_find_element(GrpPtr,Strucname,ElemPtr,Status) 
       if(Status.ne.ISDC_OK)then
          IsdcExitStatus=Status
            call MESSAGE(procName,&
                 'Element not found  in the OG : '//Strucname,&
                 IsdcExitStatus,Status)
         return
      endif 
   endif !not clearflag
endif

!=============================================================
END SUBROUTINE ReplaceIdx
!===============================================================

!====================================
subroutine deco(scwcleaimage)
!=============================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FITS_DECLARATIONS
IMPLICIT NONE

!INPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: scwcleaimage
integer(kind=4) entier,inan,i,j,im,jm
real(kind=4)reel

equivalence (reel,entier)

nanfilter(:,:)=1

im = size(scwcleaimage,1)
jm = size(scwcleaimage,2)

inan =  -1

do i = 1,im
do j = 1,jm
reel = scwcleaimage(i,j)
if (entier.eq.inan) nanfilter(i,j) = 0
end do
end do

!call FITS_FILE(1,nanfilter,1,nom_out='nanfilter.fits')
!====================================
END SUBROUTINE deco
!====================================



!====================================
subroutine RealSize(scwcleaimage,imin,imax,jmin,jmax)
!=============================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE FITS_DECLARATIONS
IMPLICIT NONE

!INPUT VARIABLES
REAL,dimension(:,:),pointer :: scwcleaimage
integer ::imin,imax,jmin,jmax
integer ::idim,jdim
real :: suma


idim = size(scwcleaimage,1)
jdim= size(scwcleaimage,2)



imin=1
suma = sum(nanfilter(imin,:))
do while((suma==0.).and.(imin .lt.idim))
   imin=imin+1
   suma = sum(nanfilter(imin,:))
enddo
imax = idim
suma = sum(nanfilter(imax,:))
do while((suma==0.).and.(imax .gt.1))
   imax=imax-1
   suma = sum(nanfilter(imax,:))
enddo

jmin=1
suma = sum(nanfilter(:,jmin))
do while((suma==0.).and.(jmin .lt.jdim))
   jmin=jmin+1
   suma = sum(nanfilter(:,jmin))
enddo
jmax = jdim
suma = sum(nanfilter(:,jmax))
do while((suma==0.).and.(jmax .gt.1))
   jmax=jmax-1
   suma = sum(nanfilter(:,jmax))
enddo

!!$print *,' default map size :',idim,jdim
!!$print *,' real map size i :',imin,imax,sum(nanfilter(imin-1,:)),sum(nanfilter(imax+1,:))
!!$print *,' real map size j :',jmin,jmax,sum(nanfilter(:,jmin-1)),sum(nanfilter(:,jmax+1))
!====================================
END SUBROUTINE RealSize
!====================================

!=============================================================
SUBROUTINE EnergyBandInfo(GrpPtr,EBins,CHans,Status)
!===============================================================

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE1_MODULE
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer :: EBins
INTEGER,dimension(:,:),pointer ::  Chans
INTEGER         :: GrpPtr,Status
!LOCAL VARIABLES
INTEGER         :: i,iok
character(len=20)  :: procName

procName = 'EnergybandInfo'
Status = ISDC_OK
if(OutEnergyNumber==0)then
   call Message(procname,&
        'No INTENSITY extension in'//trim(ogFile),DataError,status)
   return
endif

call ReadIndexInfo(GrpPtr,EBins,CHans,Status)
if(Status.ne.ISDC_OK)return

InEnergyNumber=OutEnergyNumber
if(associated(OutEnergyBands))deallocate(OutEnergyBands)
allocate(OutEnergyBands(OutEnergyNumber,2),stat=iok)
if(iok.ne.0)then
      call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
endif
if(associated(InEnergyBands))deallocate(InEnergyBands)
allocate(InEnergyBands(InEnergyNumber,2),stat=iok)
if(iok.ne.0)then
      call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
endif
if(associated(InToOutBandNums))deallocate(InToOutBandNums)
allocate(InToOutBandNums(InEnergyNumber),stat=iok)
if(iok.ne.0)then
      call MESSAGE(procName,'Allocation problem',&
                   AllocError,Status)
      return
endif


InEnergyBands = ebins

OutEnergyBands = ebins

InToOutBandNums=(/ (i, i=1,InEnergyNumber) /)

!====================================
END SUBROUTINE EnergyBandInfo
!====================================

!=====================================
SUBROUTINE FixSouPosInCat(Status)
!=====================================

USE ISDC
USE DAL3GEN_F90_API  
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE

IMPLICIT NONE

!INPUT/OUTPUT VARIABLES

INTEGER         :: Status
!LOCAL VARIABLES
INTEGER         ::i,j,infixsou
CHARACTER (len=10) :: str10
character(len=20)  :: procName

procName = ' FixSouPosInCat'
Status = ISDC_OK

InCatFixed(:) = .false.
        

! fixing source positions in the fit - only for Scw fitting
! if FitPosFixed then all source positions will be fixed
if(FitMode)then
   if(.not.FitPosFixed)then
      ! not all source pos fixed
      if(.not.associated(InCatFixedList))then
            if(InFixSouNumber==-1)then
               ! all cat sources fixed pos - modified in hidden file
               InCatFixed(:) = .true.
               InFixSouNumber= SourceNumber
            endif
         else
            ! list of pix pos sources ( FitPosFixed=.false.)
            ! can be read from hidden file of input catalog
            ! hidden file has higher priority
            InFixSouNumber = size(InCatFixedList,1)
            j=0
            do i=1,InFixSouNumber
               inFixSou = InCatFixedList(i)
               if(inFixSou.le.SourceNumber)then
                  InCatFixed(inFixSou) = .true.
                  str250 = adjustl(InCatName(inFixSou))
                  str10 = str250(1:10)
                  write(str250,&
                       '("For cat. source no.",I4,", ",A10," position will be fixed in Scw fit")',&
                       err=333)inFixSou,str10
                  goto 334
333               str250='Some source positions will be fixed in fitting'  
334               call message(procname,str250,ZeroError,Status)
                  j = j+1
               endif
            enddo
            InFixSouNumber = j 

            deallocate(InCatFixedList)
         endif! list of pix pos sources
      endif! not all source pos fixed
   endif!FitMode

!====================================
END SUBROUTINE FixSouPosInCat
!====================================

!=====================================
SUBROUTINE Read1Mosa(MosaNumber,EBins,Chans,Status)
!=====================================
USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API
USE DAL3CAT_F90_API
USE DAL3HK_F90_API  !NP
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE3_MODULE

IMPLICIT NONE 

!input/output variables
INTEGER :: MosaNumber,Status
REAL(kind=4),dimension(:,:),pointer :: EBins
INTEGER,dimension(:,:),pointer ::  Chans
!LOCAL  VARIABLES
INTEGER     :: i,GrpPtr,mosaIdxPtr
CHARACTER(len=100) :: selection
CHARACTER (len=30) :: strucname
CHARACTER (len=20)  :: procName

procName = 'Read1Mosa'
Status = ISDC_OK

open(20,file='mimosa.mosa')
read(20,'(i5)')MosaNumber

!-----------------------------------------------
do i=1,MosaNumber
!----------------------------------------------  

   read(20,'(a)')OgFile
   MosaNames(i) = ogFile
  
   ogfile=MosaNames(i)
   if ((i/100)*100==i)then
      write(str250,*)i,mosanumber,trim(ogfile)
      call message(procname,str250,ZeroError,Status)
   endif


   Status = dal_object_open(trim(OgFile)//'[1]',GrpPtr,Status)   

   if(Status.ne.ISDC_OK)then
      if (i==1)then
         call MESSAGE(procName,'dal_object_open problem for'//trim(ogfile),&
              DataError,Status)
         return
      else
         write(str250,*)'Rejected input file ',i,' case 2'
          call MESSAGE(procName,str250,ZeroError,Status)
          ActiveScwArray(i) = 0.
          Status=ISDC_OK
      endif
   endif

   if(ActiveScwArray(i)==1)then
      strucname=MosaImaStrucName
      selection= "(IMATYPE .eq.'INTENSITY')"
      status = dal3gen_index_find_member(GrpPtr,&
           StrucName,selection,OutEnergyNumber,mosaIdxPtr,Status)

      if(Status==ISDC_OK)then
         mosaType(i) = 1 !mosaics
      else
         Status=ISDC_OK
         strucname=ImaStrucName
         selection= "(IMATYPE .eq.'INTENSITY')"
         status = dal3gen_index_find_member(GrpPtr,&
              StrucName,selection,OutEnergyNumber,mosaIdxPtr,Status)
         if(Status==ISDC_OK)then
            mosaType(i) = 0 !sky image
         else
            if(i==1)then
               IsdcExitStatus=Status
               call MESSAGE(procName,&
                    'dal3gen_index_get_num_members problem for'//trim(ogFile),&
                    IsdcExitStatus,Status)
               return
            else
               write(str250,*)'Rejected input file ',i,' case 2'
               call MESSAGE(procName,str250,ZeroError,Status)
               ActiveScwArray(i) = 0.
               Status=ISDC_OK
            endif
         endif
      endif
      if(ActiveScwArray(i)==1)then
         if (i==1) then 
            call EnergyBandInfo(GrpPtr,EBins,CHans,Status)
            if(Status.ne.ISDC_OK)then
               IsdcExitStatus=Status
               call MESSAGE(procName,&
                    ' Energy info problem',IsdcExitStatus,Status)
               return
            endif
         endif

         if(OutEnergyNumber>1)then
            Status = dal_object_close(mosaIdxptr,DAL_DELETE,Status)
            if(Status.ne.ISDC_OK)then
               IsdcExitStatus=Status
               call MESSAGE(procName,&
                    ' cannot close mosaIdxptr '//trim(ogFile),IsdcExitStatus,Status)
               return
            endif
         endif

      endif!ActiveScwArray case 2
      Status = dal_object_close(grpptr,DAL_SAVE,Status)

      if(Status.ne.ISDC_OK)then
         IsdcExitStatus=Status
         call MESSAGE(procName,&
              ' cannot close Grp ptr  '//trim(ogFile),IsdcExitStatus,Status)
         return
      endif
endif!Activescwarray case 1
enddo
close(20)

      
!=====================================
END SUBROUTINE Read1Mosa
!=====================================


!=====================================
SUBROUTINE Read1ReCallMosa(MosaNumber,EBins,Chans,Status)
!=====================================
USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API
USE DAL3CAT_F90_API
USE DAL3HK_F90_API  !NP
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE3_MODULE

IMPLICIT NONE 

!input/output variables
INTEGER :: MosaNumber,Status
REAL(kind=4),dimension(:,:),pointer :: EBins
INTEGER,dimension(:,:),pointer ::  Chans
!LOCAL  VARIABLES
INTEGER     :: i,GrpPtr,mosaIdxPtr
CHARACTER(len=100) :: selection
CHARACTER (len=30) :: strucname
CHARACTER (len=20)  :: procName

procName = 'Read1ReCallMosa'
Status = ISDC_OK

open(20,file='mimosa.mosa')
read(20,'(i5)')MosaNumber
do i=1,MosaNumber
   read(20,'(a)')OgFile
   MosaNames(i) = ogFile
enddo
close(20)

ogfile = mosanames(1)
Status = dal_object_open(trim(OgFile)//'[1]',GrpPtr,Status)   

if(Status.ne.ISDC_OK)then
   call MESSAGE(procName,'dal_object_open problem for'//trim(ogfile),&
        DataError,Status)
   return
endif

call EnergyBandInfo(GrpPtr,EBins,CHans,Status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' Energy info problem',IsdcExitStatus,Status)
   return
endif
Status = dal_object_close(grpptr,DAL_SAVE,Status)

if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,&
        ' cannot close Grp ptr  '//trim(ogFile),IsdcExitStatus,Status)
   return
endif

!=====================================
END SUBROUTINE Read1ReCallMosa
!=====================================

!=====================================
SUBROUTINE Read2Mosa(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle,incr,Status)
!=====================================
USE ISDC
USE DAL3GEN_F90_API
USE DAL3AUX_F90_API
USE DAL3CAT_F90_API
USE DAL3HK_F90_API  !NP
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_BASE1_MODULE
USE MIMOSA_BASE3_MODULE

IMPLICIT NONE 

!input/output variables
INTEGER :: MosaNumber,Status
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle
REAL(kind=8)  ::incr

!LOCAL  VARIABLES
INTEGER     :: i,GrpPtr,iDim,jDim,imin,imax,jmin,jmax
! pointer to output source catalogue
REAL(kind=4),dimension(:,:),pointer :: MapArray
CHARACTER (len=20)  :: procName

procName = 'Read2Mosa'
Status = ISDC_OK

do i=1,MosaNumber
   if(ActiveScwArray(i)==1)then
   ogfile=MosaNames(i)
   if ((i/100)*100==i)then
      write(str250,*)i,trim(ogfile)
      call message(procname,str250,ZeroError,Status)
   endif
   Status = dal_object_open(trim(OgFile)//'[2]',GrpPtr,Status)
   if(Status.ne.ISDC_OK)then
      call MESSAGE(procName,'dal_object_open problem for'//trim(ogfile),&
           DataError,Status)
      return
   endif
   call  ReadAttributesAndAtti(i,GrpPtr,RaX,DecX,RaZ,DecZ,posAngle,&
     incr,Status)
   if(Status.ne.ISDC_OK) then
      write(str250,*)'Rejected input file ',i,' case 3'
      call MESSAGE(procName,str250,ZeroError,Status)
      ActiveScwArray(i) = 0.
      Status=ISDC_OK
   else
      call ReadUnknownArray(grpptr,iDim,jDim,MapArray,Status)
      if(Status.ne.ISDC_OK) then
         write(str250,*)'Rejected input file ',i,' case 4'
         call MESSAGE(procName,str250,ZeroError,Status)
         ActiveScwArray(i) = 0.
         Status=ISDC_OK
      else
         allocate(nanfilter(idim,jdim),stat=status)
         if(Status.ne.ISDC_OK) then
            call MESSAGE(procName,&
                 ' allocation problem ',AllocError,Status)
            return
         endif
   
         call deco(MapArray)
         DefMapSize(i,1) = idim
         DefMapSize(i,2) = jdim

         call RealSize(maparray,imin,imax,jmin,jmax)
         RealMapSize(i,1,1)=imin
         RealMapSize(i,1,2)=imax
         RealMapSize(i,2,1)=jmin
         RealMapSize(i,2,2)=jmax

         deallocate(MapArray,nanfilter)
 
         Status = dal_object_close(GrpPtr,DAL_SAVE,status)
         if(Status.ne.ISDC_OK)then
            IsdcExitStatus=Status
            call MESSAGE(procName,'dal_object_close problem for'//trim(ogFile),&
            IsdcExitStatus,Status)
            return
         endif
      endif
   endif
endif!ActiveScwArray=1
enddo
!=====================================
END SUBROUTINE Read2Mosa
!=====================================


!=====================================
SUBROUTINE ReadAttributesAndAtti(imos,GrpPtr,RaX,DecX,RaZ,DecZ,posAngle,&
     incr,Status)
!=====================================

USE ISDC
USE DAL3GEN_F90_API  
USE DAL3AUX_F90_API
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE

IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER         :: imos,GrpPtr
! pointer to output source catalogue
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle
REAL(kind=8)  ::incr
INTEGER         :: Status
!LOCAL VARIABLES
INTEGER             :: i,j
REAL(kind=8)        :: ptab(2)
CHARACTER (len=20)  :: colname 
REAL(kind=8)        :: val
CHARACTER (len=10) :: str10
character(len=20)   :: procName

procName = ' ReadAttributesAndAtti'
Status = ISDC_OK

ptab(:) = 0.

do i=1,7
    select case(i)
    case(1)
       colname = 'TFIRST'
    case(2)
       colname = 'TLAST'
    case(3)
       colname = 'TELAPSE'
    case(4)
       colname = 'TSTART'
    case(5)
       colname = 'TSTOP'
    case(6)
       colname = 'ONTIME'
    case(7)
       colname = 'DEADC'
    end select
    Status = dal_attribute_get_real&
         (GrpPtr,colName,val,end0,end1,Status)
    if(Status.ne.ISDC_OK) then
       IsdcExitStatus=Status
       write(str250,*) ' dal_attribute_get_real problem for ',&
            colName,' for ',trim(ogFile)
       call MESSAGE(procName,str250,IsdcExitStatus,Status)
       return
    endif
   
    shdtimes(imos,i) = val
enddo


colname='CTYPE1'
Status = dal_attribute_get_char&
     (GrpPtr,colName,str10,end0,end1,Status)
if(Status.ne.ISDC_OK) then
   IsdcExitStatus=Status
   write(str250,*) ' dal_attribute_get_real problem for ',&
        colName,' for ',trim(ogFile)
   call MESSAGE(procName,str250,IsdcExitStatus,Status)
   return
endif
j=verify('GLON',str10)
InputProjType(imos)=-1
if(j ==0)then
   InputProjType(imos) = 0
else
   j=verify('RA---CAR',str10)
   if (j == 0)then
      InputProjType(imos) = 1
   else
      InputProjType(imos) = 2
   endif
endif

if(InputProjType(imos) < 0)then
       call MESSAGE(procName,'Unknown projection type for '//trim(ogFile),DataError,Status)
       return
 endif


do i=1,5
    select case(i)
    case(1)
       colname='CRVAL1'
    case(2)
       colname='CRVAL2'
    case(3)
       colname='ONTIME'
    case(4)
       colname='CD1_1   '
    case(5)
       colname='CD1_2   '
   end select

   Status = dal_attribute_get_real&
        (GrpPtr,colName,val,end0,end1,Status)
   if(Status.ne.ISDC_OK) then
      IsdcExitStatus=Status
      write(str250,*) ' dal_attribute_get_real problem for ',&
           colName,' for ',trim(ogFile)
      call MESSAGE(procName,str250,IsdcExitStatus,Status)
      return
   endif

    select case(i)
    case(1)
       if(InputProjType(imos)==0)then
          LX(imos) = val
       else
          RAX(imos) = val
       endif
    case(2)
       if(InputProjType(imos)==0)then
          bX(imos) = val
       else
          decX(imos) = val
       endif
       if(InputProjType(imos)==0)then
          call D_CONV_LB_AD(LX(imos),BX(imos),rax(imos),decx(imos))
       endif
    case(3)
       Duration(imos) = val
    case(4)
       ptab(1) = val
    case(5)
       ptab(2) = val
    END SELECT
enddo

ptab=-ptab/incr
if(ptab(2)==0.)then
   posangle(imos) = 0.
else
   if(ptab(1)==0.)then
      posangle(imos) = 90.
   else
      posangle(imos)=atan2(ptab(2),ptab(1))*180.0d0/acos(-1.0d0)
   endif
endif

status=DAL3AUX_calc_Z_Axis(RAX(imos),DECX(imos),PosAngle(imos),RAZ(imos),DECZ(imos),status)
if(Status.ne.ISDC_OK)then
   IsdcExitStatus=Status
   call MESSAGE(procName,'DAL3AUXcalcZAxis problem for'//trim(ogFile),&
        IsdcExitStatus,Status)
   return
endif
!====================================
END SUBROUTINE ReadAttributesAndAtti
!====================================

!====================================
SUBROUTINE WriteRecallInfo(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle)
!=====================================

USE ISDC
USE DAL3GEN_F90_API  
USE DAL3AUX_F90_API
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE

IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER :: MosaNumber
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle

!LOCAL VARIABLES
INTEGER :: i
character(len=20)   :: procName

procName = ' WriteRecallInfo '


open(10,file='mimosa.info',form='formatted')
write(10,'(I6,1x,I2)')mosanumber,OutEnergyNumber
do i=1,MosaNumber
      write(10,'(i5,1x,i1,1x,2(I1,1x),"|",2(i4,1x),4(i4,1x),"|",2(2(F7.1,1x),f10.1,1x),f5.3,"|",8(g22.16,1x))')i,&
        activescwarray(i),mosatype(i),InputProjType(i),&
         DefMapSize(i,1:2),RealMapSize(i,1:2,1:2) ,shdtimes(i,1:7),&
       Duration(i),Posangle(i),Rax(i),decX(i),Raz(i),decZ(i),LX(i),BX(i)
enddo
close(10)


!====================================
END SUBROUTINE WriteRecallInfo
!====================================


!====================================
SUBROUTINE ReadRecallInfo(MosaNumber,RaX,DecX,RaZ,DecZ,posAngle)
!=====================================

USE ISDC
USE DAL3GEN_F90_API  
USE DAL3AUX_F90_API
USE ATTI_DEFS
USE ATTI_DECLARATIONS
USE ATTI_INTERNAL
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_USE_MODULE

IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER :: MosaNumber
REAL(kind=8), dimension(:), pointer :: RaX,DecX,RaZ,DecZ,posAngle

!LOCAL VARIABLES
INTEGER :: i,k
character(len=1) :: str1
character(len=20)   :: procName

procName = ' ReadRecallInfo '


open(10,file='mimosa.info',form='formatted')
read(10,'(I6,1x,I2)')mosanumber,OutEnergyNumber
do i=1,MosaNumber
   read(10,'(i5,1x,i1,1x,2(I1,1x),a1,2(i4,1x),4(i4,1x),a1,2(2(F7.1,1x),f10.1,1x),f5.3,a1,8(g22.16,1x))')k,&
        activescwarray(i),mosatype(i),InputProjType(i),&
         str1,DefMapSize(i,1:2),RealMapSize(i,1:2,1:2) ,str1,shdtimes(i,1:7),str1,&
       Duration(i),Posangle(i),Rax(i),decX(i),Raz(i),decZ(i),LX(i),BX(i)
enddo
close(10)

!====================================
END SUBROUTINE ReadRecallInfo
!====================================
