

!***********************************
MODULE MIMOSA_WORK_MODULE
!***********************************

INTERFACE

SUBROUTINE MIMOSA_work(Status)
!-----------------------------------
 
    IMPLICIT NONE
    INTEGER::Status

END SUBROUTINE MIMOSA_work


END INTERFACE
!***********************************
END MODULE MIMOSA_WORK_MODULE
!***********************************


!***********************************
MODULE WORK_MODULE
!***********************************
INTERFACE

SUBROUTINE Dealloc_Scw_Arrays
!--------------------------------
IMPLICIT NONE
END SUBROUTINE Dealloc_Scw_Arrays


END INTERFACE


!***********************************
END MODULE WORK_MODULE
!***********************************

!##########################
!SUBROUTINES CODE        
!##########################




!====================================
SUBROUTINE MIMOSA_work(Status)
!====================================

USE ISDC
USE DAL3AUX_F90_API  
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE MIMOSA_AUX_MODULE
USE MIMOSA_INIT_MODULE
USE MIMOSA_BASE_MODULE
USE MIMOSA_DECLEAN_MODULE

USE IBIS_IMAGING_PARAM
USE ATTI_DECLARATIONS
USE ATTI_DEFS
USE FITS_DECLARATIONS
USE WORK_MODULE
IMPLICIT NONE

!OUTPUT VARIABLES
INTEGER :: Status 
 
!LOCAL VARIABLES
! OG / Scw POINTERS
!""""""""""""""""""""
INTEGER :: grpPtr         ! pointer to OG
INTEGER :: MosaImaIdxPtr  ! pointer to mosaicked image list 
INTEGER :: MosaSouIdxPtr  ! pointer to mosaicked image source list
INTEGER :: scwPtr         ! scw pointer

INTEGER :: imaIdxPtr      ! scw image index ptr
INTEGER :: souIdxPtr      ! scw source list index ptr
INTEGER :: outCatPtr      ! pointer to output source catalogue
! DETECTOR ARRAYS    
!""""""""""""""""""""

! Scw image arrays : raw,cleaned,significance,residual,variance
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
REAL(kind=4),dimension(:,:),pointer :: ScwImage,ScwVar,&
                          ScwCleaImage,ScwResid,ScwSignif


! OG MAP ARRAYS   : raw,cleaned,significance,residual,variance
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
REAL(kind=4),dimension(:,:),pointer :: map,mapExpo,&
                                       mapSignif,mapVar,pMapVar

! ATTITUDE VARIABLES
!""""""""""""""""""""
! arrays of telescope axes (X,z) position for each pointing
REAL(KIND=8),DIMENSION(:),POINTER::raX,decX,raZ,decZ,posAngle
INTEGER :: attiNumPoi !number of pointings+1
!last position is the map centre  




! Scw and OG source list from cleaning
REAL      (kind=4), dimension(:,:), pointer ::  skySnr,mapSnr 
INTEGER, dimension(:), pointer   ::  SourceQualityFlag

! Scw image sizes and centres
INTEGER :: iSkySize,jSkySize,im,jm!,iSkySizepoprz,jSkySizepoprz,
REAL      (kind=4) :: skyCentreI,skyCentreJ!,mem,mem_increase
!OG image sizes and centres
INTEGER :: iMapSize,jMapSize
REAL      (kind=4) :: mapCentreI,mapCentreJ



!CONTROL VARIABLES
!""""""""""""""""""""

INTEGER   :: refFrameNumber ! last position = map centre

INTEGER  :: SouNumber

! off-axis correction variables
REAL(kind=8),dimension(:,:),pointer,save :: OffCorrBins  
INTEGER  :: OffCorrIdxPtr,OffCorrBinNumber

!AUXILIARY VARIABLES
!""""""""""""""""""""
LOGICAL :: okScw,okExpoMap
real(kind=4),dimension(:,:),pointer :: xyscwlist,xymaplist
REAL      (kind=4) :: multifactor,moy,val,meandc,ta,tb,mosatime
INTEGER  ::  inBand,outBand,statusClose,status1
INTEGER  ::  MosaNumber,mosa
INTEGER  ::  iok,i,j,imod,nso,iso,is,js,k,meansum
character(len=80) :: str80
character(len=20) :: str20
character(len=15) :: str15,strp
character(len=2) :: strno
character(len=1):: str1
character(len=2) :: str2
character(len=200) :: filename
character(len=20)::  procName   


! Set imaging parameters
CALL IBIS_IMAGING_INIT(DetType,DebugMode)
! Coding Power (to be computed) for the IBIS mask

procName= 'MIMOSA_CONTROL '
statusClose = ISDC_OK
status1 = ISDC_OK
Status = ISDC_OK

!PAR FILE READING

call ReadParFile(Status)
if(Status.ne.ISDC_OK) return

if(JestOffCorrIndex)then
   call  ReadOffCorrIndex(OffCorrIdxPtr,OffCorrBinNumber,&
     OffCorrBins,Status)
   if(Status.ne.ISDC_OK) then
      call war_message(procname,&
           ' Off-axis correction will not be done',Zeroerror,Status)
     JestOffCorrIndex = .false.
   else
      JestOffCorrIndex = .true.
   endif
else
   call war_message(procname,&
        ' No off-axis correction will be done',Zeroerror,Status)
endif

if(JestOffCorrIndex)then
   call message(procname,' Off-axis correction included',ZeroError,Status)
endif

call InitWork(Status)
if(Status.ne.ISDC_OK) return

call INIT

call OgOpen(MosaNumber,MosaImaIdxPtr,MosaSouIdxPtr,&
            OutCatPtr,raX,decX,raZ,decZ,posAngle,Status)
if(Status.ne.ISDC_OK) goto 100
call message(procname,'Data information written into recall.info',0,status)
call message(procname,'Subsequent runs can be done with ReCall=1',0,status)

if(DoPart2==0) goto 100



!ATTITUDE INITIALISATION
call message(procname,'Attitude treatment...',ZeroError,Status)

call  AttiInitialise(MosaNumber,raX,decX,raZ,decZ,posAngle,&
                          AttiNumPoi,iMapSize,jMapSize,Status)
if(Status.ne.ISDC_OK)goto 100


! Scw image centers and OG map centers
!!$skyCentreI=(iSkySize+1)/2
!!$skyCentreJ=(jSkySize+1)/2

mapCentreI=(iMapSize+1)/2
mapCentreJ=(jMapSize+1)/2



! Read source list,param arrays init 
 call CleanStaInit(raX,decX,raZ,decZ,&
                          AttiNumPoi,nToClean,Status)
if(Status.ne.ISDC_OK)goto 100

! SOME COMMENTS
call Description(MosaNumber,Rax,Decx,imapsize,jmapsize)


! Reference frame is : the last attitude 
refFrameNumber = attiNumPoi


call MESSAGE(procName,'---------------------------------------------',ZeroError,Status) 
write(str250,*,err=264)OutEnergyNumber,' maps will be created'
goto 265
264 str250 = ' '
265 call MESSAGE(procName,str250,ZeroError,Status) 
write(str250,*,err=380)' Map size : ',imapsize,'x',jmapsize
goto 381
380 str250 = ' '
381 call MESSAGE(procName,str250,ZeroError,Status) 


ALLOCATE(map(iMapSize,jMapSize),mapExpo(iMapSize,jMapSize),&
         mapSignif(iMapSize,jMapSize),&
         mapVar(iMapSize,jMapSize),pMapVar(iMapSize,jMapSize),&
        offaxiscorr(400,400),stat=iok)
 if(iok.ne.0)then
   call MESSAGE(procName,&
        'Allocation problem ',AllocError,Status)
   goto 100
endif



TimeConvertStat = 100

Multifactor=1.


! possible one energy band
OneBandMode = EnergyBand

!one scw mode
!OneScwMode=5

outband_loop:do outBand = 1,OutEnergyNumber
 if((OneBandMode==0).or.&
           ( OneBandMode== outBand).or.&
           ((OneBandMode==-1).and.&
           (outBand.ge.BandStart).and.(outBand.le.BandStop))) then
      map(:,:)       =  0.0
      mapExpo(:,:)   = 0.
      pMapVar(:,:)   = 0.
      mapVar(:,:)    =  0.0
      mapSignif(:,:) = 0.0

      
     
      TotalTime = 0.

      if(JestOffCorrIndex)then
         call ReadOffCorrBand(OutBand, OffCorrIdxPtr,&
              OffCorrBinNumber, OffCorrBins, OffAxisCorr,Status)
         if(Status.ne.ISDC_OK)then
            Status = ISDC_OK
            call war_message(procname,&
                 ' Cannot create off-axis correction array',0,Status)
            OffAxisCorr(:,:) = 1.
         endif
              
      endif
    
!!$
!!$      !----------------------------
!!$      ! LOOP ON POINTINGS          
!!$      !----------------------------
      write(str250,*,err=266)' => Map no.',OutBand,' : ',MosaNumber,'mosaics images reading and mosaicking'
      goto 267
266   str250 = errorstr//' 266'
267   call MESSAGE(procName,str250,ZeroError,Status) 
!open(13,file='atti.txt')
     
      poi_loop :do mosa = 1,MosaNumber
      if(ActiveScwArray(mosa)==0) then
         write(str250,'(a,I5,a,a)')'REJECTED input image  no.:',mosa,' ',trim(mosanames(mosa))
         call message(procname,str250,ZeroError,Status)
      else
      if((OneScwMode==0).or.(OneScwMode==mosa))then
!!$         write(str250,'(a,I2,a,I5,a,a)')'E band :',outband,', Projection of input image  no.:',mosa,' ',trim(mosanames(mosa))
!!$         call message(procname,str250,ZeroError,Status)

         if(mosaType(mosa)==0)then
            Rebintype=17 ! individual Scw
         else
            Rebintype=19 !Mosaicks
         endif
         call ReadMosaBandImage(mosa,outBand,iSkySize,jSkySize,&
                            scwCleaImage,scwVar,scwExpo, Status)
              
         if(Status.ne.ISDC_OK)goto 100


        if(DebugMode >0)then
            call FITS_FILE(RefFrameNumber,scwcleaimage,1,nom_out='image.fits') 
            call FITS_FILE(RefFrameNumber,nanfilter,1,nom_out='nanfilter.fits')
         endif
         mosaTime = Duration(mosa)


         ! no time weighting
             
         TotalTime = TotalTime+mosaTime
    

         SkyCentreI=(DefMapSize(mosa,1)+1)/2
         SkyCentreJ=(DefMapSize(mosa,2)+1)/2
!!$
!!$
!           useful thing scw xy coordinates fot cat sources
!           call ScwCoordinates(iSkySize,jSkySize,&
!                     iEdge,jEdge,ScW,xyscwlist,Status)   
!!$
!!$           ! ADDS THE IMAGE TO THE FINAL MAP map FOR A GIVEN 
!!$           ! OUTPUT ENERGY BAND outBand 
!!$           ! mf,is,tc,rot_type
if((mosaType(mosa)==0).and.(detmean >0))then
   Rebintype=17 ! individual Scw
else
  Rebintype=19 !Mosaicks
endif 
! RebinType
!  5 - PixSpread=0
!  7 - PixSpread=1, no covar
!  8 - PixSpread=1, error sum , correl = 1
!  9 - PixSpread=1, covar 0.75
! 17 - Pixspread=1, covar from covartab
! 19 - PixSpread=1, covar 0.75, 0.75**2
! 10 - LS fit - covar from PSF shape zone +/-2
! 11 - Chi2 fit - covar from PSF shape zone +/-2
! 12 - LS fit - covar from PSF shape zone +/-1
! 13 - Chi2 fit - covar from PSF shape zone +/-1
! 14 - LS fit , psf centered on pix proj , zone +/- 2
! 15 - Chi2 fit , psf centered on pix proj , zone +/- 2
! 16 - PixSpread with norm a_i^2 , covar 0.75

select case(RebinType)
case(5)
   str250 = 'RebinType=5 : PixSpread=0 '
case(7)
   str250 = 'RebinType=7 : PixSpread=1, no covariance assumed'
case(8)
   str250 = 'RebinType=8 : PixSpread=1, error summation ,thus  correlation = 1'
case(9)
   str250 = 'RebinType=9 : PixSpread=1, covar = 0.75'
case(17)
   str250 = ' RebinType=17 - Pixspread=1, covar from covartab'
case(19)
   str250 = ' RebinType=19 : PixSpread=1, covar 0.75, 0.75**2'
case(10)
   str250 = 'RebinType=10 : LS fit - covar from PSF shape, zone +/-2'
case(11)
   str250 = 'RebinType=11 : Chi2 fit - covar from PSF shape, zone +/-2'
case(12)
   str250 = 'RebinType=12 : LS fit - covar from PSF shape, zone +/-1'
case(13)
   str250 = 'RebinType=13 : Chi2 fit - covar from PSF shape, zone +/-1'
case(14)
   str250 = 'RebinType=14 : LS fit , psf centered on pix proj , zone +/- 2'
case(15)
   str250 = 'RebinType=15 : Chi2 fit , psf centered on pix proj , zone +/- 2'
case(16)
   str250 = 'RebinType=16 :PixSpread with norm a_i^2 , covar 0.75'
end select

select case(MosaType(mosa))
case(0)
   str15='Input Scw'
case(1)
   str15 = 'Input Mosa'
end select
select case(InputProjType(mosa))
case(0)
   str20 = ', proj : GAL-CAR, '
case(1)
   str20 = ', proj : RA-CAR, '
case(2)
   str20 = ', proj : RA-TAN, '
end select
if(projType==0)then
   if(EquaGal==1)then
      strp='GAL-CAR, '
   else
      strp='RA-CAR, '
   endif
else
   strp='RA-TAN, '
endif

!!call message(procname,trim(str15)//trim(str20)//'Output proj :'//trim(strp)//trim(str250),ZeroError,Status)


         call ROTAT (mf=MagnifFactor,is=ImaCarte,& 
              tc=ProjType,rot_Type=RebinType,&
              idim1=iSkySize,jdim1=jSkySize,&
              idim2=iMapSize,jdim2=jMapSize,&
              ci1= skycentreI,cj1=skyCentreJ,&
              ci2=mapCentreI,cj2=mapCentreJ,&
              num1=mosa,num2=refFrameNumber,&
              image= scwCleaImage, scwExposition=scwExpo,&
              imagevar=scwVar,&
              carte=map,cartevar=mapVar,normvar = pMapvar,&
              exposure=mapExpo, &
              time=real(duration(mosa)),info=DebugMode,&
              detmeanval = detmean,status=Status)
        

         if(Status > 0)then
            write(str250,*)&
                 'Projection problems occured -'&
                 //'projection angle too big for mosa ',mosa
            call War_MESSAGE(procName,str250,0,Status)
         endif

         if(Status.ne.0)then
            call MESSAGE(procName,&
                 'Rotation problem ',RotatError,Status)
            goto 100
         endif
          
         deallocate(nanfilter)



!!$! status = 1  when  x< 0 
!!$!       The arrival direction is below the detector plane.
!!$!       In this case the direction cosines are conserved but 
!!$!       the resulting point (y,z) is placed outside the radius 
!!$!       yz_radius.To recover the original position do
!!$!             n= dsqrt(y**2+z**2)
!!$!             y = y*yz_radius/n
!!$!             z = z*yz_radius/n
!!$! status = 2 when x near 0
!!$!      This means a very great direction angle with respect to 
!!$!      the line-of-sight-axis X

     endif !oneScwMode
endif !ActiveScwArray
 ENDDO poi_loop

!close(13)    


 
      TotalTime2 =TotalTime**2 

      select case(RebinType)
      case(0) ! no flux spread no weighting
         where(mapExpo.gt.0.)
            map = map/mapExpo
            mapVar = mapVar/mapExpo/mapExpo
           
         endwhere
      case(1,2) !Bevington,modif
          where(mapvar >0.)
             map = map/mapvar
            mapvar = 1./mapvar
           
         endwhere
      case(4) !Willmore corr
         where(pMapVar >0.)
            map   = map/pMapVar
            mapvar = mapvar/pMapVar**2
            
         endwhere
      case(5)  ! no flux spread , variance weighting
         where(pmapvar >0.)
            map = map/pmapvar
            mapvar =  mapvar/pmapvar**2
          endwhere
!!$          if(MapFluxReNorm)then
!!$             where(MapPointVar >0.)
!!$                MapPoint = MapPoint/MapPointVar
!!$             endwhere
!!$          endif
       case(6) !Willmore 
         where(pMapVar >0.)
            map   = map/pMapVar
            mapvar = mapvar**2/pMapVar**2
            
         endwhere
       case(7) !Willmore with covar tab
         where(pMapVar >0.)
            map   = map/pMapVar
            mapvar = mapvar/pMapVar**2
         endwhere
!!$         if(MapFluxReNorm)then
!!$            where(MapPointVar >0.)
!!$               MapPoint = MapPoint/MapPointVar
!!$            endwhere
!!$         endif
      case(8,9,10,11,12,13,14,15,16,17,19) !
         where(pMapVar >0.)
            map   = map/pMapVar
            mapvar = mapvar/pMapVar**2
         endwhere
!!$         if(MapFluxReNorm)then
!!$            where(MapPointVar >0.)
!!$               MapPoint = MapPoint/MapPointVar
!!$            endwhere
!!$         endif
      end select

      if(DebugMode.gt.0)then
         call FITS_FILE(RefFrameNumber,map,1,nom_out='map.fits')
         call FITS_FILE(RefFrameNumber,mapexpo,1,nom_out='exposure.fits')
         call FITS_FILE(RefFrameNumber,mapvar,1,nom_out='mapvar.fits')
         where(mapvar >0.)
            mapsignif = map/sqrt(mapvar)
         endwhere
         call FITS_FILE(RefFrameNumber,mapsignif,1,nom_out='mapsignif.fits')
      endif


      ! FINDS SOURCES INSIDE THE MAP FOR A GIVEN ENERGY BAND
      write(str250,*)'ANALYSIS OF ENERGY BAND ',outband
      call MESSAGE(procName,str250,ZeroError,Status)

      call MapFovSources(iMapSize,jMapSize,outBand,&
           iDetDim,jDetDim,Status)

      if(Status.ne.0)goto 100

      ! CLEANES THE FINAL MAP 
           
    
      call MESSAGE(procName,&
           '  => Map analysis and searching for sources...'&
           ,ZeroError,Status) 
      call InMapSearching&
           (RefFrameNumber,outBand,Map,MapExpo,mapVar,MapSignif,&
           MapSnr,SourceQualityFlag,SouNumber,Status)
      if(Status.ne.ISDC_OK)goto 100
      If(SouNumber > 0)then
         write(str250,*,err=268)'      Found ',SouNumber,'sources'
         goto 269
268      str250 = errorstr//' 268'
269      call MESSAGE(procName,str250,ZeroError,Status) 
      else
         call MESSAGE(procName,&
              'Attention! : no sources found in the total FOV',&
              ZeroError,Status)
      endif
  
       

      ! WRITE TO OG THE FINAL CLEANED MAP FOR A GIVEN
      !ENERGY BAND

       !TFIRST,TLAST,TELAPSE,TSTART,TSTOP,ONTIME,DEADC
       ogtimes(1) =   minval(shdtimes(1:MosaNumber,1))  !TFIRST
       ogtimes(2) =   maxval(shdtimes(1:MosaNumber,2))  !TLAST
       ogtimes(3) =  sum(shdTimes(1:MosaNumber,3))     !TELAPSE
       ogtimes(4)  =   minval(shdtimes(1:MosaNumber,4)) !TSTART
       ogtimes(5)  =    maxval(shdtimes(1:MosaNumber,5))!TSTOP
       ogtimes(6) = totaltime                          !ONTIME
       if(totaltime >0.)&                              !DEADC
       ogtimes(7) = sum(shdtimes(1:MosaNumber,7)*shdtimes(1:MosaNumber,6))/&
                    totaltime

!!$          write(str250,*)'test print :TFIRST,TLAST,TELAPSE,TSTART,TSTOP,ONTIME,DEADC :  ',ogtimes 
!!$          call message(procname,str250,0,Status)
       


!!$       Map = Map/TotalTime
!!$       MapVar  = MapVar/TotalTime/TotalTime
       

       call WriteImages(1,refFrameNumber,outBand,&
                        grpPtr,mosaImaIdxPtr,Map,mapVar,&
                        MapExpo,MapSignif,Status)

       if(DebugMode.gt.0)then
         call FITS_FILE(RefFrameNumber,mapsignif,1,nom_out='mapsignif_po.fits')
      
      endif




      if(Status.ne.ISDC_OK) goto 100

      ! write OG source list
      !if(SouNumber > 0)then
         call WriteSourceList(1,refFrameNumber,SouNumber,&
              grpPtr,mosaSouIdxPtr,outBand,mapSnr,Status)

         if(Status.ne.ISDC_OK) goto 100
      !endif
!==========================================
!END OF  MAIN Loop on output energy bands 
!==========================================
endif
enddo outband_loop


call UpdateIdx(mosaSouIdxPtr,Status)
if(Status.ne.ISDC_OK)goto 100


!WRITING OUTPUT SOURCE CATALOGUE
call WriteOutCat(OutCatPtr,Status)
if(Status.ne.ISDC_OK)goto 100

 call MESSAGE(procName,'Og closing',ZeroError,statusClose)
100    statusClose = dal_object_close(grpPtr,DAL_SAVE,statusClose)

if(statusClose.ne.ISDC_OK)then
  call MESSAGE(procName,&
   ' Problem in closing of Og',IsdcProcError,statusClose)
endif

IF(Status.eq.0)then
  if(statusClose.ne.0)then
    Status = statusClose
  endif
endif

if(associated(ImaStatPar))deallocate(ImaStatPar)
if(associated(BkgParTab))deallocate(BkgParTab)

if(associated(TimeBins))deallocate(TimeBins)


if(associated(InEnergyBands))deallocate(InEnergyBands)
if(associated(InToOutBandNums))deallocate(InToOutBandNums)
if(associated(OutEnergyBands))deallocate(OutEnergyBands)
IF(ASSOCIATED(posAngle))DEALLOCATE(posAngle)
IF(ASSOCIATED(raX))DEALLOCATE(raX)
IF(ASSOCIATED(decX))DEALLOCATE(decX)
IF(ASSOCIATED(raZ))DEALLOCATE(raZ)
IF(ASSOCIATED(decZ))DEALLOCATE(decZ)
IF(ASSOCIATED(MapExpo))       DEALLOCATE (MapExpo)
IF(ASSOCIATED(Map))       DEALLOCATE (Map)
IF(ASSOCIATED(MapVar))    DEALLOCATE (MapVar)
IF(ASSOCIATED(pMapVar))    DEALLOCATE (pMapVar)
IF(ASSOCIATED(MapSignif)) DEALLOCATE (MapSignif)
IF(ASSOCIATED(InScwCat))  DEALLOCATE (InScwCat)
IF(ASSOCIATED(InSimCat))  DEALLOCATE (InSimCat)
IF(ASSOCIATED(InSourceList))DEALLOCATE (InSourceList)
IF(ASSOCIATED(FilterImage)) DEALLOCATE(FilterImage)
IF(ASSOCIATED(FilterSky)) DEALLOCATE(FilterSky)
IF(ASSOCIATED(skySnr))    DEALLOCATE (skySnr)
IF(ASSOCIATED(mapSnr))    DEALLOCATE (mapSnr)
IF(ASSOCIATED(cleanParTab))   DEALLOCATE (cleanParTab)
IF(ASSOCIATED(Duration))  DEALLOCATE (Duration)
IF(ASSOCIATED(ShdTimes))  DEALLOCATE (ShdTimes)
IF(ASSOCIATED(ScwImage))  DEALLOCATE(ScwImage) 
IF(ASSOCIATED(ScwCleaImage)) DEALLOCATE(ScwCleaImage)
IF(ASSOCIATED(ScwResid))  DEALLOCATE(ScwResid)
IF(ASSOCIATED(ScwNoSou))  DEALLOCATE(ScwNoSou)
IF(ASSOCIATED(ScwSignif)) DEALLOCATE(ScwSignif)
IF(ASSOCIATED(ScwVar))    DEALLOCATE(ScwVar)
IF(ASSOCIATED(ScwExpo))   DEALLOCATE(ScwExpo)
IF(ASSOCIATED(ThickCorr)) DEALLOCATE(ThickCorr) 
IF(ASSOCIATED(OffAxisCorr)) DEALLOCATE(OffAxisCorr) 
IF(ASSOCIATED(OutSourceCat))DEALLOCATE(OutSourceCat)
IF(ASSOCIATED(OutSourceCatFlux))DEALLOCATE(OutSourceCatFlux)
IF(ASSOCIATED(OutSourceCatLocErr))DEALLOCATE(OutSourceCatLocErr)
IF(ASSOCIATED(OutSourceCatFluxErr))DEALLOCATE(OutSourceCatFluxErr)
IF(ASSOCIATED(OutSourceCatSig))DEALLOCATE(OutSourceCatSig)
IF(ASSOCIATED(InCatAlpha))DEALLOCATE(InCatAlpha)
IF(ASSOCIATED(InCatDelta))DEALLOCATE(InCatDelta)
IF(ASSOCIATED(InCatId))DEALLOCATE(InCatId)
IF(ASSOCIATED(InCatName))DEALLOCATE(InCatName)
IF(ASSOCIATED(InCatFlux))DEALLOCATE(InCatFlux)
IF(ASSOCIATED(FitDecTab)) DEALLOCATE (FitDecTab)
IF(ASSOCIATED(FitPosX)) DEALLOCATE (FitPosX)
IF(ASSOCIATED(AllSouMapFlux))DEALLOCATE(AllSouMapFlux)
IF(ASSOCIATED(AllSouMapSnr))DEALLOCATE(AllSouMapSnr)
IF(ASSOCIATED(AllSouList))DEALLOCATE(AllSouList)
IF(ASSOCIATED(RawFluxestab))DEALLOCATE (RawFluxesTab)
IF(ASSOCIATED(point_rot_mat))DEALLOCATE (point_rot_mat)
IF(ASSOCIATED(pointing_table))DEALLOCATE (pointing_table)
IF(ASSOCIATED(InCatFixed))DEALLOCATE (InCatFixed)
IF(ASSOCIATED( OffCorrBins))DEALLOCATE ( OffCorrBins)
IF(ASSOCIATED( ActiveScwArray))DEALLOCATE (ActiveScwArray )

if (CopyCatStat.ne.0)then
    write(str250,*,err=270)'Input Catalogue Copying  Statistics: ',CopyCatStat
    goto 271
270 str250 = errorstr//' 270'
271 call MESSAGE(procName,str250,ZeroError,Status1)
endif
if (TimeInfo.ne.0)then
    write(str250,*,err=272)'Time problems in header : Statistics: ',TimeInfo
    goto 273
272 str250 = errorstr//' 272'
273 call MESSAGE(procName,str250,ZeroError,Status1)
endif
if (MissKeyStat.ne.0)then
    write(str250,*,err=274)' Missing Keywords Statistics: ',MissKeyStat
    goto 275
274 str250 = errorstr//' 274'
275 call MESSAGE(procName,str250,ZeroError,Status1)
endif
if (WriteSourceListStat.ne.0)then
    write(str250,*,err=276)' Statistics of Source list writing problems: ',WriteSourceListStat
    goto 277
276 str250 = errorstr//' 276'
277 call MESSAGE(procName,str250,ZeroError,Status1)
endif
if(WriteOutCatStat.ne.0)then
    write(str250,*,err=278)' Statistics of Output Catalogue writing problems: ',WriteOutCatStat
    goto 279
278 str250 = errorstr//' 278'
279 call MESSAGE(procName,str250,ZeroError,Status1)
endif 
if(EmptyEffiStat.ne.0)then
    write(str250,*,err=280)' Statistics of empty shd effi : ',EmptyEffiStat
    goto 281
280 str250 = errorstr//' 280'
281    call MESSAGE(procName,str250,ZeroError,Status1)
endif    
if(EmptyShdStat.ne.0)then
    write(str250,*,err=282)' Statistics of empty shd  : ',EmptyShdStat
    goto 283
282 str250 = errorstr//' 282'
283 call MESSAGE(procName,str250,ZeroError,Status1)
endif  
if(EmptyShdVarStat.ne.0)then
    write(str250,*,err=284)' Statistics of empty shd var : ',EmptyShdVarStat
    goto 285
284 str250 = errorstr//' 284'
285 call MESSAGE(procName,str250,ZeroError,Status1)
endif 


!===============================
END SUBROUTINE MIMOSA_work
!===============================




!====================================
SUBROUTINE Dealloc_Scw_Arrays
!====================================

USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
USE MIMOSA_USE_MODULE
USE IBIS_DECON_LIB_MOD

IMPLICIT NONE


if(associated(CatSourceModel))deallocate(CatSourceModel)
if(associated(CatSourceModelNumero))deallocate(CatSourceModelNumero)
if(associated(CatSourceModelSNR))deallocate(CatSourceModelSNR)
if(associated(ijd_stop))deallocate(ijd_stop)


!====================================
ENDSUBROUTINE Dealloc_Scw_Arrays
!====================================
