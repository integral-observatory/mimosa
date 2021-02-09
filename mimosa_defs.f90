

!*********************************
MODULE MIMOSA_CONTROL_MODULE
!*********************************
USE ISDC
IMPLICIT NONE 
!-------------------------------------------- 
!  MAIN CONTROL VARIABLES  DEFINED ONLY HERE   
!--------------------------------------------
INTEGER,parameter :: WorkMode=0!main control variable
!0 - default ISDC work mode
!1 - other work mode
! Attention :if 0 then control variables marked here by (WM) 
!  are set to their default values in ConfigMask
! control variables marked by (PF) are defined from .par file
INTEGER,parameter :: SaclayMode = 0
! 0 - ISDC mode
! 1 - Saclay tests - some parameters can be redefined
LOGICAL,parameter :: DEPISTAGE= .false.
! if .true. then some special printing added
INTEGER,parameter :: HiddenParFile = 1
! 0 - ISDC modea
! 1 - Saclay OSA mode - hidden file
Logical,save         :: FitMode = .true.
! .false. - no fitting 
! .true.  - fitting - adequate procedures should be given (see declean)
! reset to false by CalibMode if souDetDist > 0
INTEGER,save :: FitType = 1
! 0 - PSF FIT
! 1 - Gaussian fit
Logical,save         :: FitPosFixed = .false.
Integer,save         :: FitOffset = 1
! if 1 then fit offset in deg put to fin_rd_error rather than theor psle error

INTEGER ,save        :: VerifEnergyBandMode = 0
! 0 - input/output energy band verif - do not reset by WorkMode=0
! 1 - input energy band set to output ones

INTEGER,parameter :: CalibMode = 0
! 0 - definitive data ( infinite source distance)
! 1 - calibration mode ( finite source distance read from .par file)

INTEGER  :: IsdcExitStatus

!-------------------------------------------------------
!  CONTROL VARIABLES PREDEFINED HERE. THEY WILL BE RESET
!  TO DEFAULT VALUES in ConfigMask subroutine 
!  IF WorkMode = 0. THEY WILL BE REDEFINED IF 
!  ParaFile = 1  - read from .para file  in ConfigMask
!-------------------------------------------------------
INTEGER, save     :: ParaFile = 0     !(WM)
! 0 -  default value
!      control variables will not be redefined
! 1 -  control variables will be redefined - read from
!      .para file in ConfigMask subroutine
!----------------------
! MAIN PRINT/INFO LEVEL
!----------------------
INTEGER ,save        :: DebugMode=0    !(WM)  
!print control
! 0 - no prints except error messages and some controls
!     - default value
! 1 - auxiliary prints in logfile
! 2 - copy of the prints on the screen 
!     and some fits file creation
! 3 - tracing mode

!-----------------------
! ATTITUDE CONTROLS     
!-----------------------
INTEGER,save :: AttitudeTimeFlag = 2
! 0 - old attitude time procedure
! 1 - NP attitude time procedure
! 2 - BON attitude time procedure
INTEGER, save     :: IOSphere =1     !(WM)  
! 1 - images seen from the outside of celestial sphere : 
! alpha increases to the right
! -1 - images seen from theinside of celestial sphere :
!  alpha increases to the left - default value
REAL(kind=4) :: XDisi = 0.5                   !(WM)  
REAL(kind=4) :: XDisj = 0.5                   !(WM)  
! telescope X axis displacent with respect to the 
! central pixel of the image
REAL(kind=4) :: MapDisi = 0.0                 !(WM)  
REAL(kind=4) :: MapDisj = 0.0                  !(WM)  
! Map centre displacent with respect to the 
! central pixel of the image
! AUXILIARY VARIABLE                                           
INTEGER,parameter :: MaxScwNumber = 4
! used only in attitude/duration simulation - only in this unit
!..............................................................

INTEGER ,save        :: AttiMode = 0   !(WM)  
! 0  - attitude data read from OG
! -1 -  attitude date read from auxiliary file
! 1 - attitude date taken from UserRa(Dec)X(Z) arrays
REAL(kind=4),dimension(MaxScwNumber),save ::UserRAX =(/125.,135.,145.,0./)
REAL(kind=4),dimension(MaxScwNumber),save ::UserDecX=(/-30.,-40.,0.,0./)
REAL(kind=4),dimension(MaxScwNumber),save ::UserRAZ =(/125.,135.,145.,0./)
REAL(kind=4),dimension(MaxScwNumber),save ::UserDecZ=(/60.,50.,90.,0./)

!-----------------------
!REBINNING CONTROLS
!-----------------------
INTEGER  ,save       :: CleanMode=1   !(WM)  
!1 - cleaning in Scw
!-1   - no cleaning 
INTEGER  ,save       :: MapMode=0  
! 0 map from cleaned images
! 1 map from raw images
! 2 map from resid images

INTEGER    ,save     :: ProjType = 0   !(WM)  
! projection type in the final sky map
! 1 TAN-TAN
! 2 TAN-CAR
INTEGER    ,save     :: ImaCarte = 1   
! projection direction
! 0 image-to_sky
! 1 sky-to-image
INTEGER    ,save     :: RebinType = 0  !(WM)  
! 0 - no spread
! 1 - willmore
! 4 - pixel division  - default value
REAL(kind=4)     ,save    :: ImageMagnifFactor = 1.,MagnifFactor = 1.,&
                             PointMagnifFactor=1.       !(WM)  
!ImageMagnifFactor - image
!MagnifFactor,PointMagnifFactor - carte
INTEGER  ,save           :: MaxMapISize = 4320
INTEGER  ,save           :: MaxMapJSize = 2160
! ratio image/map/point map  pixel /detector pixel size
INTEGER  ,save           :: iSkydim=400,jSkydim=400 !(WM)
! deconvolved image size - should be even

!-----------------------------
! Science Window time CONTROLS
!-----------------------------
INTEGER  ,save       :: ScwDurationMode = 1 !(WM)  
! auxiliary variables for source simulations
! redefinition of Scw duration
! 0 - Scw duration from OG - default value
! 1 - Scw Duration takem from ScwDurationArray
!..............................................................

INTEGER :: scwi
REAL(kind=4),dimension(MaxScwNumber) ,save  :: ScwDurationArray=&
             (/(1.,scwi=1,MaxScwNumber)/)
REAL(kind=4)      ,save      ::PredefBKG=22. 
REAL(kind=4)      ,save      ::PredefFlux=64. !cts/pixel
REAL(kind=4)      ,save      ::TimeUnit=1000.
! PredefBKG and PredefFlux are use for SimulMode = 1 to simulate shd
! Attention : effective source flux is : PredefFlux*Duration(scw)/Timeunit
REAL(kind=4)      ,save      ::deltaLIMIT=90.
!--------------------------
! RANDOM GENERATOR CONTROLS
!--------------------------
INTEGER ,save        :: GRAIN_INIT=1    !(WM)  
! 0 - save random routines state 
! 1 - can be reset after first run with GRAIN_INIT=1 - default value
INTEGER ,save        :: RepetableGeneratorMode=1   !(WM)  
! 0 - random initialisation - default value
! 1 - previous sequence - default value

!-------------------------------
! GENERAL FUNCTIONALITY CONTROLS
!-------------------------------
LOGICAL,save                 :: EffiWeiFlag = .true.
! if tru then weighted deconvolution
LOGICAL,save                 :: IsGoodScwFound=.false.

LOGICAL,save                 :: FastOpenFlag = .false.
! if true my own opening of the OG
LOGICAL,save                 :: ChangeGroupFlag = .false.
LOGICAL,save                 :: FFTInitFlag = .true.
! changed to false after first FFT call
LOGICAL,save                 :: NormalFlag = .true.
! if true one normal array BAL & WEI created per Scw
LOGICAL,save                 :: OKNormArrays = .false.
REAL(kind=4),parameter       :: EffiEnergyThreshold = 63.
LOGICAL,save                 :: OneModelMode = .false.
!if true : in principle one model per source & ScW 
LOGICAL,save                 :: NegModelsFlag = .false.
! if false : no negative models creation
LOGICAL,save                 :: JestOffCorrIndex = .false.
LOGICAL,save                 :: RawFitFluxVerif = .false.
! if true then fit rejected if fit flux < rawflux
LOGICAL,save                 :: RawFluxEstim = .false.
! raw flux estimation added in the _res file if true.
INTEGER ,save        :: BiasCorrection = 1

LOGICAL,save                 :: NoSouMap = .false.
!Image without sources at 4th extention of image index if true
LOGICAL,save                 :: ResidMap = .true.
! Residual map at 4th extension of ima idx
LOGICAL,save                 :: ExpoMap = .false.
! true exposition map at 4 or 5 th extention of image index
! and 4 ext mosa ima index if true
LOGICAL,save                 :: OneExpoMap = .false.
! true exposition map at the end of ima idx
Logical ,save                 ::OneExpoMapNow = .false.
! auxiliary for OneExpoMap  mode
LOGICAL,save                 :: SnrCorr = .false.
! SNR image corrected to mean=0,var=1 if true
INTEGER ,save        :: ScwFitMode = 0 !(PF) for SaclayMode
! if SaclayMode ==0 then ScwFitMode is defined only here
! 0 - no Chi2 fit in Scw
! 1 - Chi2 fit in Scw
! 2 - LS fit in Scw
INTEGER ,save        :: EmptyMode = 0   !(WM)  
! 0 default value,no action
! 1 no FFT,deconv,cleaning
INTEGER ,save        :: OneScwMode=0   !(WM)  
! -1 - Scw in the range scwstart,scwstop treated
! 0 - all Scw in OG are treated - default value
! n - only one (n) scw treated
INTEGER ,save        :: ScwStart=0,ScwStop=0
INTEGER ,dimension(5) :: ScwExcluTab
INTEGER ,save        :: ScwExcluNumber
INTEGER ,save        :: OneBandMode = 0 !(WM)  

! 0 all band treated
! n - only nth band treated
! -1 band in the range BandStart,BandStop treated
INTEGER ,save        :: BandStart=0,BandStop=0
INTEGER  ,save       :: MapOnlyMode = 0 !(WM)  
! 0 - Scw and Map treatment - default value
! 1 - only Scw part treatment 
! 2 - only Map part treatment
! 3 - Scw( only shd reading) + Map
! 0 only Scw Part treatment
INTEGER  ,save       :: PixSpread = 1
! 0 - no spread
! 1 - flux spread in mosaicking
INTEGER  ,save       :: Recall = 0
INTEGER  ,save       :: DoPart2 = 1
INTEGER  ,save       :: EnergyBand = 0
! parameter read from .par file
INTEGER  ,save       :: UserProj = 0
!0 - GAL_CAR
!1 - RA_CAR
!2 - RA_TAN
INTEGER,dimension(4),save :: MapParDec=(/1,1,1,0/)
! alpha,delta,radius,pixel size of mosaicked map
! 0 calculated by the program
! 1 can be changed by the user in .par file
Real(kind=4),dimension(4),save :: MapParVal=(/0.,0.,0.,1./)
Real(kind=4),save :: MaxDiffRadius = 360.
Real(kind=4),save :: EffiTolPar = 0.05
! minimum efficiency pixel 
!-------------------------
! SOURCE SEARCH CONTROLS
!------------------------
INTEGER  ,parameter :: MaxSouInGenCat=2000
INTEGER  ,save       :: SearchMode =0    !(PF)
! read from .par file
!  source  searching/cleaning
! 0 - searching for all significant excesses (long !) 
!     default value
! 1 - searching for all catalog sources only in the individual images
! 2 - searching for a given number of sources 
!     = nToclean ( from .par file) 
! 3 - searching for all catalog sources in the individual images
!     and next for a given number of sources 
!     = nToclean ( from .par file) 

INTEGER  ,save       :: CatSouSearch = 0
! 0 - search for excesses
! 1 - search for catalog sources only in Scw images
! CatSouSearch =1 if SearchMode =1
! 2 - search for catalog sources n Scw images + excesses
INTEGER  ,save       :: CatSouSearchLim = 5
INTEGER    ,save     :: nToClean=1    !(PF)
! read from .par file 
! sign - => only searching,no cleaning
! number of significant excess to be searched for 
! meaningfull only for SearchMode==2
Logical ,save :: ReplaceSourcePeak = .true.
Logical ,save ::AddSources = .true.
INTEGER ,save        :: UserMapMode = 1
! 0 raw
! 1 cleaned
! 2 resid
! modified only in hidden file
! works only with CleanMode = 0,1
Logical ,save :: MapFluxReNorm = .false.
! if true then MapPoint used to normalise flux in the mosaic
! possibly modified  in hidden file
INTEGER ,save        :: FantomAnalysis = 1  !(WM)  
! 0 - fantom analysis off 
! 1 - fantom analysis on - default value
REAL(kind=4 ),parameter          :: DefPeakMinHeight = 6.0
! minimum peak height in sigma for catalogue source
REAL(kind=4 ),parameter          :: DefNonIdentPeakMinHeight = 7. 
! min source SNR :
REAL(kind=4 ),save          :: PeakMinHeight = DefPeakMinHeight 
! minimum peak height in sigma for catalogue source
REAL(kind=4 ),save          :: NonIdentPeakMinHeight = DefNonIdentPeakMinHeight
! read from .par file

! minimum peak height in sigma for new source
REAL(kind=4 ),save          :: SouSearchRadius    = 1.0  !(WM)  
! source identification radius 
REAL(kind=4 ),save          :: FineTol=1.  !(WM)  
! precision increase for fine position source increase
!--------------------------
! SHD  SIMULATION CONTROLS
!--------------------------
INTEGER   ,save      :: SimulMode=0     !(PF)
!   read from .par file 
!0  no simulation  - shd read from ScW
!1  image simulation 
!2  image read from ordinary file

! Controls used  if SimulMode=1
INTEGER   ,save      :: SourceCat=  3     !(WM)  
! 0 - input catalog from OG - default value
! 1 - aux. source list file 
! 2 - random source positions , flux = predefFlux, number = RanSouNumber
! 3 - random but repetable source position , random flux round predefFlux
! 4 - only one source, position given by souX,souY( set RanSouNumber=1)
! read from .para
REAL(kind=4)      ,save      :: souX = -36.9231,souY = 69.83173  
! meaningful only for SourceCat=3 - source position
INTEGER   ,save      :: DrawReadSource = 0
! meaningfull only for SourceCat = 2
! 0  - source positions drawn and written into MIMOSA.sou
! 1  - source positions read from MIMOSA.sou
INTEGER  ,save       :: FluxMode=1      !(WM)  
! auxiliary variables for source simulations
! redefinition of source flux
!0 -  flux calculated from OG data
!1 -  flux predefined - current default value
!2 -   "      "       and  only one source simulated
!     that indicated by nToSim 
INTEGER    ,save     :: nToSim=3
!number of sources to be simulated , meaningfull only for FluxMode==2
INTEGER   ,save      :: RanSouNumber=10
! number of sources to simulate in the SourceCat=2 mode
REAL(kind=4)      ,save      :: DrawZone=10.
! width of draw zone in deg in the SourceCat=2 mode
INTEGER   ,save      :: DeadZonesFilling =1    !(WM)  
!0 -  in simulated image dead zones are filling with 0
!1 -  in simulated image dead zones are filling with detector mean

!2 -                     no dead zones
INTEGER,save         :: SourcePattern = 0  !(WM)  
! 0 drawn in simulations (meaningful only for SimulMode = 1)
! 1 added "       "          "                   "          



!--------------------------
! SOURCE MODEL CONTROLS
!--------------------------
LOGICAL,save :: FineCleaning = .false.
INTEGER,save                 :: ModelType = 3
! 0  - IBIS_SIMPHO
! 1 -  Model0
! 2  - Simulator
! 3  - AnalModel   (WM)
INTEGER,save                 :: SimulModelType = 2
! 0  - IBIS_SIMPHO - make sense only for infinite source distance
! 1  - Model0
! 2  - Simulator (WM)
! 3 - analmodel
INTEGER   ,save      :: ModelPhotons=1000000 !(WM)
! number of photons drawn in source modelling


!--------------------------
! DECONVOLUTION CONTROLS
!--------------------------
INTEGER           , save     :: FFTArrays = 0 !(WM)
!0 - FFT arrays created in the exec  - default value
!1 - FFT arrays created  and written
!2 - FFT arrays read
INTEGER   ,save      :: AreaMode=1 !(WM)
! 0 - image/variance normalized by effective area old formula
! 1 - image/variance normalized by effective area new formula
!-1 - no normalisation
INTEGER ,save        ::DetType=1 !(WM)
! -n integral step=n 
! 0  129
! 1 ISGRI  real - default value
! 2 Picsit real
! replaces iDet 
INTEGER  ,save          :: DeconType =1 !(WM)
! 1 wei = 0 in dead zones ; noisy pixels set to effi corr. mean

!---------------------------
!THICKNESS EFFECT CORRECTION
!---------------------------
INTEGER,parameter :: ThickBandNumber = 5
REAL(kind=4),dimension(ThickBandNumber,2) :: ThickBand 
REAL(kind=4),dimension(ThickBandNumber,4)    ::  ThickX ,ThickY
INTEGER , save        :: ThickCorrection = 1
! 0 - no ethickness effect correction
! 1 - yes

!--------------------------
! OTHER EFFEXTS CORRECTION
!--------------------------
INTEGER,parameter :: OtherEffectBandNumber = 5
REAL(kind=4),dimension(OtherEffectBandNumber,2) :: OtherEffectBand 
REAL(kind=4),dimension(OtherEffectBandNumber,2)    :: OtherEffectXY 
INTEGER , save        :: OtherEffectCorrection = 1
! 0 - no ethickness effect correction
! 1 - yes


REAL(kind=4),save :: VertMaskAtten =0.
REAL      (kind=4),parameter :: OnAxisOpenFraction=8192.4
! Absorption approx. type
INTEGER,save :: AbsNormType = 0

!----------------------
!MISCELANOUS CONTROLS
!---------------------
Real(kind=4),save :: NormSou = 390.
Real(kind=4),save :: MaxFlux = 10000.
! standarized source peak at 1c/pixel
CHARACTER(len=Pil_LineSize),save       :: OutType
CHARACTER(len=10),Parameter   ::OutTypeName =  'OutType'
! SHD TYPE NAMES
CHARACTER(len=10),dimension(3),Parameter::ShdTypeName= (/'BIN_I','BIN_S','BIN_T' /) 
REAL(kind=4)     ,parameter          :: ActiveDetArea=(1.-0.2888)*3686.072
! Total active detector area in cm^2
!WARNING CONTROLS
!STATISTICS OF  EMPTY 0 SHD EFFI,SHD,SHD VAR
INTEGER , save        :: EmptyEffiStat = 0
INTEGER , save        :: EmptyShdStat = 0
INTEGER , save        :: EmptyShdVarStat = 0
INTEGER , save        :: MissKeyStat = 0
!STATISTICS OF MISSING KEYWORDS IN INDEX APPEND
INTEGER , save        :: TimeInfo = 0
!STATISTICS OF Duration info
INTEGER , save        :: TimeConvertStat = 0
!STATISTICS error in time converting from IJD to UTC
INTEGER , save        :: PutCharStat = 0
!STATISTICS error in char attributes putting
INTEGER , save        :: WriteSourceListStat = 0
!STATISTICS error in  source list writing
INTEGER , save        ::  WriteOutCatStat=0
!STATISTICS error in  output cat writing
INTEGER , save        ::  CopyCatStat=0
!STATISTICS error in  copying cat
INTEGER ,SAVE        ::iReb=1 
!subpixelisation of the detector(meaningful only for SimulType = 1)
INTEGER ,Parameter    :: SourceIdLength=16
REAL(kind=4)    ,Parameter    :: ObtUnitToSec = 61.0/1000000.0
REAL(kind=4)    ,Parameter    :: SecInDays=1./60./60./24.
REAL(kind=4)    ,Parameter    :: DaysToSec = 24.*60.*60.
REAL(kind=4)    ,Parameter    :: PixelAngDim=5.1/60.
INTEGER ,Parameter    :: ImaTypeNum = 4
INTEGER ,Parameter    :: ShdTypeNum = 3
CHARACTER(len=30),save :: errorstr='+ Write string error'
!NAMES OF PARAMS IN .par FILE
INTEGER ,Parameter    :: ParFileNum =7
CHARACTER(len=9),dimension(ParFileNum) :: ParFileNames=&
     (/'mask     ','deco     ','tungAtt  ','aluAtt   ','leadAtt  ','covrMod  ','corrDol  '/)
INTEGER ,Parameter    :: PointingTypeNum = 3
CHARACTER(len=8),dimension(0:PointingTypeNum) :: PointingTypeArray=&
 (/'ANY     ','POINTING','SLEW    ','OTHER   '/)
INTEGER ,Save         :: PointingType=0
CHARACTER(len=30),save :: StampComment = '  '
! STRUCTURE NAMES
CHARACTER(len=30),Parameter::MaskName            ='ISGR-MASK-MOD'
CHARACTER(len=30),Parameter::ProjDecoName        ='ISGR-DECO-MOD'
CHARACTER(len=30),Parameter::AttName             ='ISGR-ATTN-MOD'
! MOSAICKED IMAGE
CHARACTER(len=30),Parameter::MosaImaIdxStrucName ='ISGR-MOSA-IMA-IDX'
CHARACTER(len=30),Parameter::MosaImaStrucName    ='ISGR-MOSA-IMA'
CHARACTER(len=30),Parameter::MosaImaIdxFitsName  ='isgri_mosa_ima.fits'
!  SOURCE LIST
CHARACTER(len=30),Parameter::MosaSouIdxStrucName ='ISGR-MOSA-RES-IDX'
CHARACTER(len=30),Parameter::MosaSouStrucName    ='ISGR-MOSA-RES'
CHARACTER(len=30),Parameter::MosaSouIdxFitsName  ='isgri_mosa_res.fits'
!  SCIENCE WINDOW
!CHARACTER(len=30),Parameter::ScwStrucName    ='GNRL-SWG3-GRP'
CHARACTER(len=30),Parameter::ScwStrucName    ='GNRL-SCWG-GRP'
!  SHD
CHARACTER(len=30),Parameter::ShdIdxStrucName ='ISGR-CEXP-SHD-IDX'
CHARACTER(len=30),Parameter::ShdStrucName    ='ISGR-CEXP-SHD'
CHARACTER(len=30),Parameter::ShdIdxFitsName  ='isgri_cexp_shd.fits'
!  GROUPPED SHD
CHARACTER(len=30),Parameter::GrouShdIdxStrucName ='ISGR-GROU-SHD-IDX'
CHARACTER(len=30),Parameter::GrouShdStrucName    ='ISGR-GROU-SHD'
CHARACTER(len=30),Parameter::GrouShdIdxFitsName  ='isgri_grou_shd.fits'
!  IMAGE
CHARACTER(len=30),Parameter::ImaIdxStrucName ='ISGR-SKY.-IMA-IDX'
CHARACTER(len=30),Parameter::ImaStrucName    ='ISGR-SKY.-IMA'
CHARACTER(len=30),Parameter::ImaIdxFitsName  ='isgri_sky_ima.fits'
!  SOURCE LIST
CHARACTER(len=30),Parameter::SouIdxStrucName ='ISGR-SKY.-RES-IDX'
CHARACTER(len=30),Parameter::SouStrucName    ='ISGR-SKY.-RES'
CHARACTER(len=30),Parameter::SouIdxFitsName  ='isgri_sky_res.fits'
CHARACTER(len=30),Parameter::OutCatStrucName = 'ISGR-SRCL-RES'
CHARACTER(len=30),Parameter::OutCatFitsName  = 'isgri_srcl_res.fits'
CHARACTER(len=30),Parameter::InCatStrucName = 'ISGR-SRCL-CAT'
CHARACTER(len=30),Parameter::CovStrucName = 'ISGR-COVR-MOD'
!ATTENUATION length file
CHARACTER(len=30),Parameter:: AttStrucName = 'ISGR-ATTN-MOD'
CHARACTER(len=30),Parameter:: OffCorrStrucName = 'ISGR-EFFI-MOD'
CHARACTER(len=30),Parameter:: OffCorrIdxStrucName = 'ISGR-EFFI-MOD-IDX'


!ERRORS CODES

INTEGER,Parameter  :: ZeroError       = 0 
INTEGER,Parameter  :: ErrorBegin      = -122200
INTEGER,Parameter  :: AllocError      = ErrorBegin-1
INTEGER,Parameter  :: CenterError     = ErrorBegin-2
INTEGER,Parameter  :: DeconInitError  = ErrorBegin-3
INTEGER,Parameter  :: AttiError       = ErrorBegin-4
INTEGER,Parameter  :: InCatError      = ErrorBegin-5
INTEGER,Parameter  :: ParFileError    = ErrorBegin-6
INTEGER,Parameter  :: DeconError      = ErrorBegin-7
INTEGER,Parameter  :: CleanError      = ErrorBegin-8
INTEGER,Parameter  :: OutBandError    = ErrorBegin-9
INTEGER,Parameter  :: NumMembersError = ErrorBegin-10
INTEGER,Parameter  :: MemberSelError  = ErrorBegin-11
INTEGER,Parameter  :: InvalArrayError = ErrorBegin-12
INTEGER,Parameter  :: GroupShdError   = ErrorBegin-13
INTEGER,Parameter  :: OpenError       = ErrorBegin-14
INTEGER,Parameter  :: OutCatError     = ErrorBegin-15
INTEGER,Parameter  :: SouIndexError   = ErrorBegin-16
INTEGER,Parameter  :: ClearIndexError = ErrorBegin-17
INTEGER,Parameter  :: ClearElemError  = ErrorBegin-18
INTEGER,Parameter  :: FindElemError   = ErrorBegin-19
INTEGER,Parameter  :: IndexError      = ErrorBegin-20
INTEGER,Parameter  :: DeconLibError   = ErrorBegin-21
INTEGER,Parameter  :: ReflectionError = ErrorBegin-22
INTEGER,Parameter  :: FitError        = ErrorBegin-23
INTEGER,Parameter  :: AnalModelError  = ErrorBegin-26
INTEGER,Parameter  :: RotatError      = ErrorBegin-30
INTEGER,Parameter  :: DataError       = ErrorBegin-31
INTEGER,Parameter  :: AttenError      = ErrorBegin-32
INTEGER,Parameter  :: InternalError   = ErrorBegin-35
INTEGER,Parameter  :: StringWriteError = ErrorBegin-36
INTEGER,Parameter  :: OffCorrError    = ErrorBegin-37
INTEGER,Parameter  :: ISDCProcError   = ErrorBegin-40
INTEGER,Parameter  :: OtherError      = ErrorBegin-49

CHARACTER(len=Pil_lineSize) :: MaskFile,ProjDecoFile,&
                               TungFile,LeadFile,AluFile,CovrFile,OffCorrFile
CHARACTER(len=Pil_lineSize) :: OgFile,InCatFile,OutCatFile,OutMosImaFile,OutMosResFile
!!$INTEGER,parameter ::MaxMosaNum=10000
!!$CHARACTER(len=Pil_lineSize),dimension(MaxMosaNum) :: MosaNames
!!$INTEGER, dimension(MaxMosaNum,2) :: DefMapSize
!!$INTEGER, dimension(MaxMosaNum,2,2) :: RealMapSize
!!$INTEGER, dimension(MaxMosaNum) :: MosaType,InputProjType


CHARACTER(len=Pil_lineSize),dimension(:),pointer :: MosaNames
INTEGER, dimension(:,:),pointer :: DefMapSize
INTEGER, dimension(:,:,:),pointer :: RealMapSize
INTEGER, dimension(:),pointer :: MosaType,InputProjType


INTEGER,save :: EquaGal=1 ! projection type
REAL(kind=4),    save :: SouDetDist = 0.,SouTY = 0.,SouZ = 0.

!OneModel flag definitions
INTEGER,save  :: CatSourceModelMaxNumber=50
INTEGER,save :: CatSourceModelNumber
REAL(kind=4)::t1,t2
!*************************************
END MODULE MIMOSA_CONTROL_MODULE
!*************************************

!********************************************************************
MODULE MIMOSA_GLOBVAR_MODULE
!********************************************************************
USE ISDC
REAL :: detmean
REAL(kind=8),dimension(2,2,2,2) :: covpattern

Real(KIND=8),dimension(-5:5,-5:5),save::ucoeff
Real(KIND=8),parameter :: psfsigmacoeff = 3.19295!2.5005

REAL(KIND=8),DIMENSION(:),POINTER::LX,BX

!OneModel flag definitions
REAL(kind=4) :: aproj,dproj,lproj,bproj
REAL(kind=4),dimension(:,:,:),pointer :: CatSourceModel
INTEGER,dimension(:)  ,pointer :: CatSourceModelNumero
REAL(kind=4),dimension(:)  ,pointer :: CatSourceModelSNR
REAL(kind=4),dimension(:,:),pointer :: nanfilter

! no source image
REAL(kind=4),dimension(:,:),pointer :: ScwNoSou
REAL(kind=4),dimension(:,:),pointer :: ScwExpo
REAL(kind=4),dimension(:,:),pointer :: OrgEfficiency
REAL(kind=4),dimension(:)  ,pointer :: RawFluxesTab
REAL(kind=4),dimension(:,:),pointer :: ThickCorr,OffAxisCorr
INTEGER,dimension(8,4) :: ModulesLimits
REAL(kind=4),dimension(8) :: ModuleMean
REAL(kind=4)                           :: nanf
REAL      (kind=4) :: DetModelOpenFraction,CurEffAreaCorr

!parameters for cleaning
REAL      (kind=4), dimension(:),   pointer,save :: CleanParTab
!bkg parameters - for simulations
REAL(kind=4), dimension(:), pointer ,save:: BkgParTab
! array for statistic
REAL(kind=4), dimension(:), pointer,save :: ImaStatPar
INTEGER  ,save  :: scwNumber   ! number of Scw in OG
INTEGER  ,save  :: pointScwNumber   ! number of pointing Scw in OG
INTEGER  ,save  :: activescwNumber 
! number of active SCw ( not empty and pointing)
INTEGER  ,save  :: shdNumber   ! number of shd in Scw
LOGICAL,dimension(:,:),pointer :: DeTActivePixels
! VARIABLES FOR INPUT CATALOGUE SOURCES
INTEGER,save ::  SourceNumber
!number of sources in the input catalogue
INTEGER,save ::  ScwSourceNumber,MapSourceNumber
!number of sources in the current Scw FOV
INTEGER,save ::  EnergySourceNumber
! number of sources in the current Scw FOV 
! and in the current energy band
REAL (kind=4), dimension(:,:), pointer ,save:: InScwCat,InSimCat
!input sources in the current Scw FOV 
REAL (kind=4), dimension(:,:), pointer,save :: CurEffArea
! and in the current energy band
REAL (kind=4), dimension(:,:), pointer,save ::OutSourceCat
!final source list
REAL (kind=4), dimension(:,:), pointer,save ::OutSourceCatFlux
REAL (kind=4), dimension(:), pointer,save ::OutSourceCatSig
REAL (kind=4), dimension(:,:), pointer,save ::OutSourceCatFluxErr
REAL (kind=4), dimension(:), pointer,save ::OutSourceCatLocErr
!final source list flux
REAL (kind=4), dimension(:,:), pointer,save :: AllSouMapFlux,AllSouMapSnr
REAL (kind=4), save :: FluxLowLimit = 1000.
!upperlimits of source flux in the mosaicks per source end energy band
INTEGER             ,save           ::NumSouDes
!length of description of the final source list
INTEGER,save :: CoorType
! 0 TAN
! 1 Ra Dec
INTEGER,save :: GhostPeakWidth = 1
REAL(kind=4),save :: BkgMoyPar,BkgCts,FluxNorm
! parameters meaningfull for simulations only
REAL(kind=4),dimension(:),pointer,save ::InCatAlpha,InCatDelta
REAL(kind=4),dimension(:,:),pointer,save ::InCatFlux
CHARACTER(len=16),dimension(:),pointer ,save :: InCatId
CHARACTER(len=100),dimension(:),pointer ,save :: InCatName
LOGICAL,dimension(:),pointer ,save :: InCatFixed
INTEGER             ,save           ::InFixSouNumber = 0
INTEGER,dimension(:),pointer ,save ::InCatFixedList
REAL(kind=8), dimension(:), pointer :: Duration
! time in sec of pointing
REAL(kind=4), dimension(:,:), pointer,save  :: InSourceList !auxiliay array
REAL(kind=4),save :: CodPow
!VARIABLES FOR ENERGY BANDS
INTEGER  ,save    :: InEnergyNumber! number of input energy bands
REAL(kind=4),dimension(:,:),pointer,save :: InEnergyBands
!LOWER AND UPPER LIMITS OF INPUT ENERGY BANDS   
INTEGER  ,save    :: OutEnergyNumber! number of output energy bands
REAL(kind=4),dimension(:,:),pointer,save :: OutEnergyBands  
!LOWER AND UPPER LIMITS OF OUTPUT ENERGY BANDS
INTEGER,dimension(:),pointer,save :: InToOutBandNums
!InToOutNums(i)=k indicates that input energy band i is within output energy 
! band k
REAL(kind=4),dimension(:,:),pointer,save :: TimeBins  


!LOWER AND UPPER LIMITS OF INPUT ENERGY/Time  BANDS   
INTEGER  ,save, dimension(:,:),pointer ::ChanNums
INTEGER  ,save    :: ReGroup 
! 0 - no output energy band defined in .par file- no regrouppment
! 1 - elsewhere ( verified and set in OgOpen)
LOGICAL           , dimension(:,:), pointer,save :: FilterSky
CHARACTER(len=250),dimension(1),save::Dols
CHARACTER(len=250),save :: ScwResultsPath
! Dols of input/output fits
LOGICAL,save:: OgChange
! set in ReadParFile ; .true. if Og changes saved
REAL(kind=4),save :: TotalTime,TotalTime2
! HEADERS PARAMETERS
REAL(kind=4)     ,save  :: eMin,eMax
!og times
INTEGER ,Parameter    :: OGTimeNum = 7
REAL(kind=8)     ,dimension(OGTimeNum) :: ogtimes
CHARACTER(len=8),dimension(OGTimeNum) ::OgTimeColNames=&
 (/'TFIRST  ','TLAST   ','TELAPSE ','TSTART  ','TSTOP   ','ONTIME  ','DEADC   '/)
!shd times
INTEGER ,Parameter    :: ShdTimeNum = 8
REAL(kind=8)     ,save  :: ShdTFirst,ShdTLast,ShdTElapse,ShdTStart,ShdTStop
REAL(kind=8),dimension(:,:),pointer :: shdtimes
CHARACTER(len=8),dimension(ShdTimeNum) :: ShdTimeColNames=&
     (/'TFIRST  ','TLAST   ','TELAPSE ','TSTART  ','TSTOP   ',&
       'ONTIME  ','EXPOSURE','DEADC   '/)
REAL(kind=8)     ,save  :: shdontime,shdexposure,shddeadc
!sc times
REAL(kind=8)     ,save  :: ScwTElapse,ScwTStart,ScwTStop,ScwExposure
INTEGER  ,save  :: revol,chanMin,chanMax
character(len=30),save  :: expid,swid,sw_type,shd_type
CHARACTER (len=DAL_MAX_STRING),save       :: bkgfile,offcorrfilename
REAL(kind=8) , dimension(:), pointer  :: IJD_start,IJD_stop
!FITTING IN SHD SPACE DECLARATIONS
REAL   (kind=4), save,dimension(:,:,:),pointer :: SourceModel
REAL (kind=4), dimension(:), pointer ,save:: FitDecTab
REAL (kind=4), dimension(:), pointer ,save:: FitPosX,FitPosY
REAL(kind=4)    ,Parameter    :: MaxMaxLVar = 1.e+8 
CHARACTER(100),save::end0="&\0",end1="&\0"
CHARACTER(len=250) :: str250
LOGICAL, save,dimension(:,:),pointer :: DeadZonePattern
INTEGER  ,dimension(:), pointer  :: ActiveScwArray
! 0 -slew
! 1 - pointing
! -1 - empty
TYPE ScwSourcesType
   INTEGER :: souNumber
   real(kind=4),dimension(:),pointer :: alpha,delta
   real(kind=4),dimension(:),pointer :: xmap,ymap

end type ScwSourcesType

type (ScwSourcesType),dimension(:,:),pointer :: AllSouList
INTEGER :: AttenNum
real(kind=4),dimension(:),pointer :: SliceEnergy,SliceWidth,TungTab,AluTab,LeadTab

real(kind=8),dimension(2,2) :: adtab
real(kind=4),dimension(:,:),pointer :: CovarTab,Covar1tab,DetMeanTab
real(kind=4),save :: CodedMaxDist = 180.

! current pixel dead and life zones used in modelling

REAL    (kind=8), save :: LifeZoneX,LifeZoneY
REAL    (kind=8), save :: DeadZoneLeft,DeadZoneRight,&
                          DeadZoneUpper,DeadZoneLower,&
                          LifeZoneLeftShift,LifeZoneRightShift,&
                          LifeZoneUpperShift,LifeZoneLowerShift
REAL    (kind=8), save :: LifeZoneXshift,LifeZoneYshift,&
                          PixLeftShift,PixRightShift,&
                          PixLowerShift,PixUpperShift

LOGICAL :: DecodMode 


!************************************   
END MODULE MIMOSA_GLOBVAR_MODULE
!************************************


!***********************************
MODULE IBIS_IMAGING_PARAM
!***********************************
      IMPLICIT NONE

      ! Mask
INTEGER, save :: iMaBPDim, jMaBPDim, iMasDim, jMasDim
INTEGER, save :: iMasPixDim, jMasPixDim
REAL    (kind=4), save :: MasPixDimI, MasPixDimJ
REAL    (kind=4), parameter :: MasRotAng = 0.  ! Mask Rot. Ang (deg)
REAL    (kind=4), parameter :: MaElBriSiz= 2.0 ! Ma.El. bridge size (mm)
  

! ISGRI Parameters
CHARACTER (len=15) :: Detector
INTEGER, parameter :: IsgrNumPix = 16384
INTEGER, parameter :: IsgrYDim=128, IsgrZDim=128

REAL    (kind=4), parameter :: IsgrPixSiz=4.6 !mm
REAL    (kind=4), parameter :: IsgrPixDZmm =   0.3 !mm
REAL    (kind=4), parameter :: IsgrPixEffSiz=IsgrPixSiz - 2.*IsgrPixDZmm
REAL    (kind=4), parameter :: IsgrEffNorm =IsgrPixEffSiz**2/ IsgrPixSiz**2

REAL    (kind=8), parameter :: dIsgrPixDZmm =   0.3 
REAL    (kind=8), parameter :: dIsgrPixSiz=4.6 
REAL    (kind=8), parameter :: DetPixDeadZone   = dIsgrPixDZmm/dIsgrPixSiz

REAL    (kind=8), parameter :: DetPixLifeZone = (dIsgrPixSiz-2.0d0* dIsgrPixDZmm)/dIsgrPixSiz

REAL    (kind=4), parameter :: IsgrPixThi=2.0



REAL    (kind=4), parameter :: IsgrRotAng=0.
REAL    (kind=4), parameter :: IsgrXDetCentre = 598./2.
REAL    (kind=4), parameter :: IsgrYDetCentre = 616.4/2.
! error corrected 614.4/2.
     
! Mask and pixels
REAL    (kind=4), parameter :: MaElSiz = 11.2  ! Mask El. thickness (mm)
REAL    (kind=8), parameter :: DMaElSiz = 11.2


REAL    (kind=8), parameter :: MaElThi  = 16.0 !(mm)
REAL    (kind=8), parameter :: MaElThiCm =  MaElThi/10.


REAL    (kind=4), save :: MaElPixDim= MaElSiz/IsgrPixSiz  ! Mask El. to Pixel Size Ratio
REAL    (kind=4), save :: MasDetDisDim=3202. ! Mask Detector distance (pix)
REAL    (kind=4), save :: MaElThiDim=16.   ! Mask El. Thickness (pix)
REAL    (kind=4), save :: MasDetRotAng=0. ! Mask-Detector Rot angle (deg)
REAL    (kind=4), save :: PixAngDim    ! Sky-proj pix dim (deg)      

! Detector 
INTEGER, save :: iDetDim, jDetDim  ! Total Detector dim (pix)
INTEGER, save :: iActiveDetDim, jActiveDetDim 
!!$REAL    (kind=4), save :: PixThiDim, PixEffDim ! Pix thick/eff dim (pix)
!!$REAL    (kind=4), save :: PixDZDim ,PixDzHalfDim        ! Pixels Dead Zones dim (pix)
REAL    (kind=4), parameter :: DetMaskDistCorr = 8.9699999999998
REAL    (kind=4), parameter :: IsgrHalfMasDis=   3194. + DetMaskDistCorr
REAL    (kind=4), parameter :: IsgrBotMasDis= IsgrHalfMasDis   -  MaElThi/2.
                                                     ! Det-top to Mask-bottom
! Constants
REAL    (kind=4), parameter :: pi=3.141592654    


!***********************************
END MODULE IBIS_IMAGING_PARAM
!***********************************
