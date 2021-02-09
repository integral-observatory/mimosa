
!********************************
MODULE MIMOSA_DC_MODULE
!********************************

INTERFACE


SUBROUTINE IbisMapSearch(RefFrameNumber,OutBand,SKYCLE,&
           MapExpo,VAR,SkySnr,SourceQualityFlag,Mean,&
           SouNumber,Status)
!--------------------------------------------------------------
IMPLICIT NONE

!INPUT / OUTPUT  VARIABLES
INTEGER                       :: RefFrameNumber,OutBand
REAL(kind=4), dimension(:,:), pointer :: SKYCLE,MapExpo,VAR
REAL(kind=4), dimension(:,:), pointer :: SkySnr
INTEGER, dimension(:), pointer :: SourceQualityFlag
REAL(kind=4)                          :: Mean
INTEGER                       :: SouNumber,ns
INTEGER                       :: Status
END SUBROUTINE IbisMapSearch


END INTERFACE

!********************************
END MODULE MIMOSA_DC_MODULE
!********************************

!******************************
MODULE MAPLOOP
!******************************

INTERFACE

SUBROUTINE IbisMapSearchLoop(ifident,isou,excessnum,excessflag,&
           RefFrameNumber,OutBand,SKYCLE,MapExpo,VAR,Status)
!..............................................................
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER                       :: RefFrameNumber,OutBand
REAL(kind=4), dimension(:,:), pointer :: SKYCLE,MapExpo,VAR
Logical ::ifident,excessFlag
INTEGER                       :: isou,excessnum,Status
END SUBROUTINE IbisMapSearchLoop
END INTERFACE
!******************************
END MODULE MAPLOOP
!******************************

!============================================================
SUBROUTINE IbisMapSearch(RefFrameNumber,OutBand,SKYCLE,&
           MapExpo,VAR,SkySnr,SourceQualityFlag,Mean,&
           SouNumber,Status)
!============================================================

USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_DECON_LIB_MOD
USE IBIS_IMAGING_PARAM
USE FIT_PAR_MODULE
USE OVER_FIT_SUB_MODULE
USE FITS_DECLARATIONS
USE SOURCE_IDENT1_MOD
USE SEARCH_CLEAN_AUX_MODULE
USE IbisMapSearch_DEFS
USE MAP_SEARCH_AUX_MODULE
USE MAPLOOP
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER                       :: RefFrameNumber,OutBand
REAL(kind=4), dimension(:,:), pointer :: SKYCLE,MapExpo,VAR
REAL(kind=4), dimension(:,:), pointer :: SkySnr
INTEGER, dimension(:), pointer :: SourceQualityFlag
REAL(kind=4)                          :: Mean
INTEGER                       :: SouNumber,ns
INTEGER                       :: Status

! Local variables
INTEGER,parameter :: excessnum=10


REAL(kind=4) :: area
INTEGER :: iok,isou,i
INTEGER :: j,souno,nextsou,searchcatsou
! for identif
REAL(kind=4) :: fitflux,locerr,fluxerr,sigma,localmax
INTEGER :: isn
INTEGER :: maxn,identsou,identsou1,mincatnumber1,inCatNumber
REAL(kind=4) ::flux,coeff,theorFluxErr
REAL(kind=4) :: bett_mean
REAL(kind=8) :: xcdble,ycdble,xsdble,ysdble
REAL(kind=8) ::fitfluxdble,locerrdble,fluxerrdble
Logical ::excessFlag,peak,doubleIdent,filterval
Logical :: catSource,QQche
character(len=20)::procName
INTEGER :: localStatus,loopnumber

Status=ISDC_OK
procName = 'IbisMapSearch'

coeff = 1.
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)

LocalStatus = 0


! Initialize check-printing variable
iPsfTest=DebugMode
dbw = DBLE(MaElPixDim)

call  MapAllocInitPrints(SkyCle,SourceQualityFlag,Status)

WHERE (VAR==0)
    filter=.false.
endwhere

Filter0 = filter

soupar(1)=3
snr = 0.
isou = 0
SouNumber = 0

!image stat + source search
FILTER2 = FILTER
snr(:,:)=0.
soupar(1)=0 ! no sources to be searched for
CALL IBIS_IMASN(SKYCLE,VAR=VAR,MASK=FILTER2,&
     soupar=soupar,imastat=imastat,sousnr=snr,&
     status=Status)
if(Status.ne.0)then
   if(status==1)then ! allocation problem in IMASN
      call MESSAGE(procname,&
           ' Allocation problem in IBIS_IMASN',AllocError,Status)
      return
   else
      str250 = &
           ' 0 active pixels in the mosaicked image : end of source analysis'
      call WAR_MESSAGE(procName,str250,1,Status)
         localStatus = 3
         goto 103
      endif
   endif

bett_mean=0.
!!$call BETTER_MEAN(1,0,isky,jsky,skycle,filter,&
!!$                 fanMap,bett_mean,Status)
!!$if(Status.ne.ISDC_OK)return
rmean=bett_mean 
if(rmean == 0.)rmean= imastat(3)
if(rmean==0.)then
  call War_MESSAGE(procName,&
       'image mean ==0 no source analysis will be done ',1,&
        Status)
   localStatus = 2
   goto 103
endif
Mean = rmean

excessFlag = .true.
isou = 0

if( nSouToSearch >0)then
   ! excludes very small count pixels
   if(rmean.ne.0.)BORDERS = abs((rmean-skycle)/rmean)
   where(abs(1.-borders).lt.0.01)
      filter = .false.
      SKYCLE = 0.
   endwhere
   BORDERS = 0.
   WHERE ((filter))BORDERS=SKYCLE/ SQRT(VAR)
   where(BORDERS.lt.1.)filter=.false.

endif

 SKYMOD=0. !SNR
 WHERE ((filter)) SKYMOD=(SKYCLE-rmean)/ SQRT(VAR)

searchcatsou = 1
Do while(searchcatsou .le. MapSourceNumber)
   peak = .false.
   xc=inscwcat(1,searchcatsou)+REAL(icsky)
   yc=inscwcat(2,searchcatsou)+REAL(jcsky)
   xs=xc ; ys=yc ; is=NINT(xs) ; js=NINT(ys)
   localmax=maxval(SKYMOD(is-1:is+1,js-1:js+1))
   peak = Filter(is,js).and.(localmax > siglim1)
   if(peak)then
      MinCatNumber=searchcatsou
      call  IbisMapSearchLoop(.false.,isou,excessnum,excessflag,RefFrameNumber,&
           OutBand,SKYCLE,MapExpo,VAR,Status) 
   endif
   searchcatsou = searchcatsou+1
enddo

! SOURCE SEARCH for nSouToSearch sources
! ----------------------------------------------------
SOURCELOOP : DO while((excessFlag).and.(isou.lt.nSouToSearch))
! ----------------------------------------------------
  

   peak = .false. ! no good peak detected in the image
   QQche = .true. ! something peak-like in the image
  
   loopnumber = 0
   do while((.not.peak).and.(QQche).and.(loopnumber < 3))
      loopnumber = loopnumber+1
      ! searching until good peak found or no more excesses

      !hide found sources positions and forbidden positions
      FILTER2=FILTER
      snr(:,:)=0.
      soupar(1)=excessnum ! number of source to search for 
      CALL IBIS_IMASN(SKYMOD,VAR=VAR,MASK=FILTER2,&
           soupar=soupar,imastat=imastat,sousnr=snr,&
           status=Status)
      if(Status.ne.0)then
         if(status==1)then ! allocation problem in IMASN
            call MESSAGE(procname,&
                 ' Allocation problem in IBIS_IMASN',AllocError,Status)
            return
         else
            call WAR_MESSAGE(procName,&
              '0 active pixels in the deconvolved image :  no  source analysis ',1,&
              Status)
            localStatus = 3
            goto 103
         endif
      endif ! Status ne 0 after IBIS_IMASN

      !good peak detection
      call MapPeakIdent(OutBand,MapExpo,CatSource,Peak,IdentSou,&
           QQche,Status)
   enddo ! searching until good peak found or no more excesses


   IF (peak) THEN 
       call  IbisMapSearchLoop(.true.,isou,excessnum,excessflag,RefFrameNumber,&
           OutBand,SKYCLE,MapExpo,VAR,Status) 
       if(Status==3)goto 103
       if(Status.ne.ISDC_OK)return 

   ELSE
      ! No source found
      excessFlag = .false.
   ENDIF

ENDDO SOURCELOOP

SouNumber = isou ! isou - number of sources found
!quality assessment
do isn=1,SouNumber
  if(sourceList(isn)==0)then ! non-identified source
       SourceQualityFlag(isn) = 4
  else
     if(tempSnr(isou,2) > 1) then
        SourceQualityFlag(isn)  = 2 ! several sources at peak positions
     else
        if(sourceList(isn) > 0) then ! identified source
           SourceQualityFlag(isn) = 0
        else ! source identified at ghost position
           SourceQualityFlag(isn) = 1
        endif
    endif
  endif
enddo
do isn=1,SouNumber
 souno= abs(sourceList(isn))
 if(souno> 0) then
    doubleIdent = .false.
    do j = isn+1,SouNumber
       nextsou = abs(sourceList(j))
       if((isn.ne.j).and.(souno==nextsou))then
          doubleIdent = .true.
       endif
    enddo
    if(doubleIdent)SourceQualityFlag(isn) = 3
 endif
enddo

!----------------------------------
! POSSIBLE EXIT FROM SEARCHING LOOP
!-----------------------------------
103  SouNumber = isou 


!final mean calculation out of sources
snr(:,:)=0.
soupar(1)=0 ! 0 source to search for 
            ! all sources peak hidden in FILTER
CALL IBIS_IMASN(SKYCLE,VAR=VAR,MASK=FILTER,&
     soupar=soupar,imastat=imastat,sousnr=snr,&
     status=Status)

if(Status.ne.ISDC_OK)then
   imastat(3) = 0
   Status = 0
   call message(procname,&
        'Cannot calculate final mean with IBIS_IMASN',ZeroError,status)
else
   rmean = imastat(3)
endif

!!$call BETTER_MEAN(1,0,isky,jsky,skycle,filter,&
!!$     fanMap,bett_mean,Status)
!!$if(Status.ne.ISDC_OK)then
!!$   str250 = ' Cannot calculate final Mean'
!!$   call message(procname,str250,ZeroError,Status)
!!$   Status = 0
!!$else
!!$   rmean=bett_mean 
!!$endif
Mean = rmean



! normalisations to the highest source flux
!skycle = skycle*coeff
!var = var*coeff**2

!Final mean calculation

!Mean  = Mean*coeff



str250 = '  Flux upper limits for all FOV sources ...'
call message(procname,str250,Zeroerror,status)
! UPPER LIMITS of FLUX
if(SouNumber > 0)then
   do i=1,SouNumber
      if(sourcelist(i)>0)then ! catalog source
         AllSouMapFlux(sourcelist(i),OutBand) = FluxLowLimit
      endif
   enddo
endif
do isou = 1,ScwSourceNumber
   ns = Inscwcat(10,isou)
   if(AllSouMapFlux(ns,OutBand).eq.0.)then
      ! source not found in the MAP 
      xc=InScwCat(1,isou)+REAL(icsky)
      yc=InScwCat(2,isou)+REAL(jcsky)
      xs=xc ; ys=yc ; is=NINT(xs) ; js=NINT(ys)
      if((is>0).and.(is<isky+1).and.(js >0).and.(js < jsky+1))then
         if(filter0(is,js))then
            ! source in deconvolved zone
            AllSouMapFlux(ns,OutBand)=max(0.,SkyCle(is,js)/mapexpo(is,js))
            
           if(VAR(is,js) > 0.)then
               AllSouMapSnr(ns,OutBand)=&
                    max(0.,(skycle(is,js))/SQRT(VAR(is,js)))
            endif 
!!$print *,'isou,is,js,AllSouMapFlux',&
!!$                 isou,is,js,AllSouMapFlux(ns,:)
         endif ! source in deconvolved zone
      endif ! source in map FOV 
   endif ! source not found in the MAP 
enddo 
     


IF(ASSOCIATED(SkySnr))DEALLOCATE (SkySnr)
if(SouNumber > 0)then
   allocate(SkySnr(1:SouNumber,tempsnrsize),stat=iok)
   IF (iok /= 0) then
      call MESSAGE(procName,&
        'Allocation problem',AllocError,Status)
      return
   endif

   !OUTPUT LIST OF SOURCES
   SkySnr(1:SouNumber,:) = tempSnr(1:SouNumber,:) 
endif !SouNumber > 0

if((DebugMode > 0).and.(SouNumber > 0))then
   maxn = min(isou,20)
   write(str250, '(" Found    :",25(f5.0,2x))',err=311) Tempsnr(1:maxn,1)
   goto 312
311 str250 = errorstr//' 311'
312 call Message(procName,str250,ZeroError,Status)
   write(str250,'(" Count/pix :",25(f5.2,2x))',err=313)Tempsnr(1:maxn,6)
   goto 314
313 str250 = errorstr//' 313'
314 call Message(procName,str250,ZeroError,Status)
endif 


IF(ASSOCIATED(BORDERS))  DEALLOCATE(BORDERS)
IF(ASSOCIATED(catSouFlux))     DEALLOCATE(catSouFlux)
IF(ASSOCIATED(dista))   DEALLOCATE(dista)
IF(ASSOCIATED(Fanmap))  DEALLOCATE(fanmap)
IF(ASSOCIATED(FILTER0))  DEALLOCATE(FILTER0)
IF(ASSOCIATED(FILTER))  DEALLOCATE(FILTER)
IF(ASSOCIATED(FILTER2)) DEALLOCATE(FILTER2)
IF(ASSOCIATED(imastat))  DEALLOCATE(imastat)
IF(ASSOCIATED(index)) DEALLOCATE(index)
IF(ASSOCIATED(SKYMOD))  DEALLOCATE(SKYMOD)
IF(ASSOCIATED(snr))     DEALLOCATE(snr)
IF(ASSOCIATED(sourceList))  DEALLOCATE(sourceList)
IF(ASSOCIATED(soupar))   DEALLOCATE(soupar)
IF(ASSOCIATED(tempSnr)) DEALLOCATE(tempSnr)


!''''''''''''''''''''''''''''''
END SUBROUTINE IbisMapSearch
!'''''''''''''''''''''''''''''''




!********************************
MODULE MIMOSA_DECLEAN_MODULE
!********************************


INTERFACE
SUBROUTINE InMapSearching(RefFrameNumber,OutBand,Map,&
                         MapExpo,mapVar,&
                         MapSignif,MapSnr,&
                         SourceQualityFlag,SouNumber,Status)
!----------------------------------------------------------------------


IMPLICIT NONE

!INPUT /  OUTPUTVARIABLES
INTEGER                         :: RefFrameNumber,OutBand
REAL(kind=4)   , dimension(:,:) ,pointer :: Map,MapExpo
REAL(kind=4)   , dimension(:,:) ,pointer :: MapVar
REAL(kind=4)   , dimension(:,:) ,pointer :: MapSignif
REAL(kind=4), dimension(:,:), pointer    :: MapSnr
INTEGER, dimension(:), pointer    :: SourceQualityFlag
INTEGER                         :: SouNumber,Status

END SUBROUTINE InMapSearching

END INTERFACE

!*************************************
END MODULE MIMOSA_DECLEAN_MODULE
!*************************************


!====================================================================
SUBROUTINE InMapSearching(RefFrameNumber,OutBand,Map,&
                     MapExpo,mapVar,MapSignif,&
                     MapSnr,SourceQualityFlag,SouNumber,Status)
!====================================================================

USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_DC_MODULE
USE SEARCH_CLEAN_AUX_MODULE
USE FITS_DECLARATIONS
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER                         :: RefFrameNumber,OutBand
REAL(kind=4)   , dimension(:,:) ,pointer :: Map,MapExpo
REAL(kind=4)   , dimension(:,:) ,pointer :: MapVar
REAL(kind=4)   , dimension(:,:) ,pointer :: MapSignif
REAL(kind=4), dimension(:,:), pointer    :: MapSnr
INTEGER, dimension(:), pointer :: SourceQualityFlag
INTEGER                         :: SouNumber,Status

!LOCAL VARIABLES
REAL(kind=4)       :: mean
CHARACTER(len=20)  :: procName

Status = ISDC_OK
procName='InMapSearching'

if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)



!NO MAP  CLEANING - only source cearching
CleanParTab(2) = -abs(CleanParTab(2))
!SOURCE SEARCHING RADIUS 
CleanParTab(4)= PeakMinHeight+0.15  ! Minimum peak height
SouNumber = 0



CALL IbisMapSearch(RefFrameNumber,OutBand,Map,MapExpo,&
                   MapVar,MapSnr,SourceQualityFlag,&
                   mean,SouNumber,Status)
if(Status.ne.ISDC_OK)then
       call MESSAGE(procName,'Problem in IbisMapSearch',&
                 CleanError,Status)
       return
endif 
CleanParTab(2) = abs(CleanParTab(2))  

!significance array
MapSignif = 0.
!22.10.04 mean deleted here
where(MapVar .gt.0.)MapSignif = Map/sqrt(MapVar) 


IF (DebugMode.gt.0) THEN
  call fits_file(RefFrameNumber,mapsignif,1,nom_out='map_clean_signif_image.fits')
 
endif




!=========================
END SUBROUTINE InMapSearching
!=========================



!============================================================
SUBROUTINE IbisMapSearchLoop(ifident,isou,excessnum,excessflag,&
           RefFrameNumber,OutBand,SKYCLE,MapExpo,VAR,Status)
!============================================================

USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_DECON_LIB_MOD
USE IBIS_IMAGING_PARAM
USE FIT_PAR_MODULE
USE OVER_FIT_SUB_MODULE
USE FITS_DECLARATIONS
USE SOURCE_IDENT1_MOD
USE SEARCH_CLEAN_AUX_MODULE
USE IbisMapSearch_DEFS
USE MAP_SEARCH_AUX_MODULE
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER                       :: RefFrameNumber,OutBand
REAL(kind=4), dimension(:,:), pointer :: SKYCLE,MapExpo,VAR
Logical ::ifident,excessFlag
INTEGER                       :: isou,excessnum,Status

! Local variables
INTEGER :: i,j,identsou,iloop,identsou1,mincatnumber1,inCatNumber
REAL(kind=4) :: fitflux,locerr,fluxerr,sigma
REAL(kind=4) ::flux,theorFluxErr
REAL(kind=8) :: xcdble,ycdble,xsdble,ysdble
REAL(kind=8) ::fitfluxdble,locerrdble,fluxerrdble
Logical ::peak,doubleIdent
Logical :: catSource,QQche
character(len=20)::procName

Status=ISDC_OK
procName = 'IbisMapSearchLoop'

iloop = 0
! SOURCE SEARCH for nSouToSearch sources
! ----------------------------------------------------
SOURCELOOP : DO while((excessFlag).and.(iloop .le.excessnum).and.(isou.lt.nSouToSearch))
! ----------------------------------------------------
  
   if(ifident) then
      !good peak detection
      call MapPeakIdent(OutBand,MapExpo,CatSource,Peak,IdentSou,&
           QQche,Status)

   else
      peak=.true.
      Identsou=1
      iloop = excessnum
   endif

 write(str250,*)'CatSource,Peak,IdentSou : ',  CatSource,Peak,IdentSou
call message(procname,str250,0,status)
   IF (peak) THEN 

      ! EXCESS FOUND
      excessFlag = .true.
      fitflux = 0.
      locerr  =  0.
     
      theorFluxErr =0.
      if(var(is,js) >0.)theorFluxErr = sqrt(var(is,js))
      if(IdentSou ==1)then
         inCatNumber = Int(InScwCat(10,MinCatNumber))
      else
         inCatNumber = 0
      endif

      FitFlux = Skycle(is,js)

      ! FITTING to get fine source location
      if(FitMode)then  
         ! fixed pos. non permitted in mosaicked fit
         !  more than one source fit possible in case of fail
         xcdble = dble(xc)
         ycdble = dble(yc)
         xsdble = xcdble
         ysdble = ycdble

         CALL OVER_FIT(0,inCatNumber,1,1,SKYCLE,VAR,FILTER,&
              xcdble,ycdble,xsdble,ysdble,fitfluxdble,&
              locerrdble,fluxerrdble,MapFitRes)

         xs = real(xsdble)
         ys = real(ysdble)
         fitflux = real(fitfluxdble)
         locerr = real(locerrdble)
         fluxerr= real(fluxerrdble)

         if(MapFitRes.ne.0)then
            !bad fit result
            call WAR_MESSAGE(procName,' Failure of fit procedure',0,Status)
            xs = xc;ys = yc
            is=NINT(xs) ; js=NINT(ys)
            theorFluxErr =nanf
            
         else

            ! good fit result
            is=NINT(xs) ; js=NINT(ys)
            if(.not.FILTER(is,js))then 
               ! forbidden zone
               xs = xc
               ys = yc
               locerr = nanf
               fluxerr = nanf
               theorFluxErr =nanf
               is=NINT(xs) ; js=NINT(ys)
               fitflux = Skycle(is,js)
            else
                ! good fit result and in admissible zone
                ! Source Identification at the fitted position
               call Identification(xs-icsky-Mapdisi,&
                    ys-jcsky-Mapdisj,radidentif,dista,&
                    MinCatNumber1,IdentSou1)

               if(IdentSou1==0)then
                  ! no source identified at the fitted pos
                  if(IdentSou > 0)then ! source identified at raw pos.
                     call WAR_MESSAGE(procName,'Failure of source recognition at fitted pos ',0,Status)
                     locerr = nanf
                     fluxerr = nanf
                     theorFluxErr = 0.
                      xs = xc
                      ys = yc
                      is=NINT(xs) ; js=NINT(ys)
                      fitflux = Skycle(is,js)
                   endif! source identified at fitted pos.
               else ! source identified at fitted position
                  theorFluxErr =0.
                  if(var(is,js) >0.)theorFluxErr = sqrt(var(is,js))
                  MinCatNumber = MinCatNumber1
                  IdentSou = IdentSou1
               endif


            endif! good fit result and in admissible zone
         endif ! good fit result
      endif !fitting
   
     
       if(VAR(is,js) > 0.)then
         sigma=fitflux/SQRT(VAR(is,js)) 

!!$         write(str250,*)'is,js , raw and fit fitflu ,raw and fitsnr ',&
!!$         is,js,skycle(is,js),fitflux,skycle(is,js)/sqrt(var(is,js)),sigma
!!$         call MESSAGE(procName,str250,0,Status)
      else
          sigma = 0.
          write(str250,'(" cannot calculate sigma source on pos. is="&
               &,i5," js=",i5)',err=309)is,js
          goto 310
309       str250 = 'cannot calculate sigma source'&
                    //errorstr//' 309'
310       call WAR_MESSAGE(procName,str250,0,Status)
       endif

      ! Set filter to hide the source in next loop
      irad=INT(soupar(2))
      do i=is-irad,is+irad
         do j=js-irad,js+irad
             FILTER(i,j) = .false.
          enddo
       enddo
     
      
      isou = isou+1
     
      
      call  MapResults(RefFrameNumber,isou,IdentSou,flux,fitflux,locerr,&
           fluxerr,theorFluxErr,sigma,mapexpo,Status)


      ! tempsnr(1:tempsnrsize=11,isou) contents :
      !  1  number in INScwCat
      !  2 > 0 if catalogue source 
      !  3,4, distance in pixel units from the detector centre
      !       ( X telescope axis)
      !  5  counts/sec  ( will be modified by spectral fit)
      !  6  counts/pixel    "
      !  7  fitted flux
      !  8 location error
      !  9  distance from catalogue source if identified
      ! 10  significance
      ! 11  flux err  
   ELSE
      ! No source found
      
      excessFlag = .false.

   ENDIF
iloop = iloop+1
ENDDO SOURCELOOP

!''''''''''''''''''''''''''''''
END SUBROUTINE IbisMapSearchLoop
!'''''''''''''''''''''''''''''''
