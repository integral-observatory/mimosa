!
!|************************************************************************
!| FILE:         ii_skyimage_model0.F90
!| VERSION:      2.8 Aleksandra Gros
!| COMPONENT:    ii_skyimage_extract
!| DeveloperS:      Aleksandra Gros   SAp/CEA-Saclay 
!|               ola@discovery.saclay.cea.fr 
!| Notes:        This file contains analytic source model creation subroutine 
!|***********************************************************************  
!| MODULES :1. MODEL0_VALS - global constants and variables
!|          2. MODEL0_AUX - auxiliary subroutines 
!|          3. MODEL_PROC- auxiliary subroutines 
!|          4. SOURCE_MODEL - main subroutine
!|************************************************************************
!**************************
!**************************
MODULE MODEL0_VALS
!**************************
! effective start and stop of mask pattern in the array
INTEGER,parameter :: iMaskStart = 33,iMaskStop = 127,&
                               jMaskStart = 33,jMaskStop = 127
REAL     (kind=8),parameter :: FigSurfMinTol= 1.E-3,FigTol = 1.E-4
LOGICAL ,save     :: MaskDetRotationFlag =.false.
REAL     (kind=8) :: xsoushift,ysoushift
REAL     (kind=8) :: relpi,theta,thetarad,tgtheta
REAL     (kind=8) :: cost,sent,alpha,alphax,alphay
REAL     (kind=8) :: xmaskcentre,ymaskcentre
REAL     (kind=8) :: xdetcentre,ydetcentre

REAL     (kind=8),dimension(4,2) ::  xyCorner
INTEGER,dimension(4,2) ::  ijCorner
INTEGER,dimension(4,2),save :: wzor 

INTEGER,parameter :: pktmax=20
TYPE line
   REAL  (kind=8)  :: a,b 
   REAL  (kind=8),dimension(pktmax,2) :: pkty
   INTEGER,dimension(pktmax) :: inlist 
   INTEGER :: ilepktow
   INTEGER :: fill_gap_dummy
END TYPE line

INTEGER,parameter :: pktmax1=20
TYPE cors
   REAL  (kind=8)  :: a,b 
   REAL  (kind=8),dimension(pktmax,2) :: pkty
   INTEGER,dimension(pktmax) :: inlist 
  INTEGER :: ilepktow
  INTEGER :: fill_gap_dummy
END TYPE cors

TYPE(line),dimension(4) :: PixLines
TYPE(cors)              :: PixCorners

INTEGER,parameter ::  mpsize = 40
! list of singular points of the projected mask pixel
REAL     (kind=8),dimension(mpsize,2) :: points
! list of int singular points of the projected mask pixel
INTEGER,dimension(mpsize,2) :: intpoints
INTEGER,dimension(mpsize) :: typepoints
! list of point belonging to the given det pixel
INTEGER,parameter :: MaxFigCorners=10
REAL     (kind=8),dimension(MaxFigCorners,2) :: Figure 
! figure point order
INTEGER,dimension(MaxFigCorners) :: vertex ,vertex1
! mask pixel surface fraction projected onto given det pixel
REAL     (kind=8),dimension(mpsize,mpsize) :: DetPixFraction
!**************************
END MODULE MODEL0_VALS
!**************************

!***************************
 MODULE MODEL0_AUX
!***************************

INTERFACE

SUBROUTINE AddPoint_MODEL0(x,y,point)
!-------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER   :: point
REAL     (kind=8) :: x,y
END SUBROUTINE AddPoint_MODEL0


SUBROUTINE FigVerif_MODEL0(ile,surface,Status)
!--------------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: ile,Status
REAL     (kind=8) :: surface
END SUBROUTINE FigVerif_MODEL0

SUBROUTINE Triangle_MODEL0(surface)
!----------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL     (kind=8) :: surface
END SUBROUTINE Triangle_MODEL0

SUBROUTINE Quad_MODEL0(surface)
!--------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL     (kind=8) :: surface
END SUBROUTINE Quad_MODEL0

SUBROUTINE Penta_MODEL0(surface)
!--------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL     (kind=8) :: surface
END SUBROUTINE Penta_MODEL0


SUBROUTINE Hexa_MODEL0(surface)
!--------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL     (kind=8) :: surface
END SUBROUTINE Hexa_MODEL0




SUBROUTINE MAKEVERTEX_MODEL0(ile)
!--------------------------
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
INTEGER :: ile  
END SUBROUTINE MAKEVERTEX_MODEL0
END INTERFACE

!***************************
END  MODULE MODEL0_AUX
!***************************



!**************************
MODULE MODEL0_PROC
!**************************

INTERFACE

SUBROUTINE DetMaskPos_MODEL0(idetpix,jdetpix,&
     xmaskpos,ymaskpos,imaskpix,jmaskpix)
!----------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: idetpix,jdetpix,imaskpix,jmaskpix
REAL(kind=8)     :: xmaskpos,ymaskpos
END SUBROUTINE DetMaskPos_MODEL0

SUBROUTINE MaskDetPos_MODEL0(imaskpix,jmaskpix,&
     xdetpos,ydetpos,idetpix,jdetpix)
!-----------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER :: imaskpix,jmaskpix,idetpix,jdetpix
REAL(kind=8)     :: xdetpos,ydetpos
END SUBROUTINE MaskDetPos_MODEL0

SUBROUTINE  MinMaxMaskPixel_MODEL0(i1,j1,i2,j2,&
            imin,imax,jmin,jmax)
!---------------------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: i1,j1,i2,j2,imin,imax,jmin,jmax
END SUBROUTINE  MinMaxMaskPixel_MODEL0

!!$SUBROUTINE  MaskBlok_MODEL0(MASKP,iMaskPix,jMaskPix,jMaskMax,ile)
!!$!--------------------------------------------------------
!!$IMPLICIT NONE
!!$
!!$!INPUT/OUTPUT VARIABLES
!!$REAL      (kind=8), dimension(:,:), pointer :: MASKP
!!$INTEGER ::iMaskPix,jMaskPix,jMaskMax,ile 
!!$END SUBROUTINE  MaskBlok_MODEL0

SUBROUTINE MaskDetDefLines_MODEL0(imin,imax,jmin,jmax)
!-------------------------------------------
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER ::imin,imax,jmin,jmax 
END SUBROUTINE MaskDetDefLines_MODEL0


SUBROUTINE DetPixChoice_MODEL0(i1,i2,j1,j2,Status)
!-----------------------------
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: i1,i2,j1,j2,Status
END SUBROUTINE DetPixChoice_MODEL0


SUBROUTINE SET_DET_MODEL0(idetpix,jdetpix,Det,MASKP,&
     imaskpix,jmaskpix,coeff)
!--------------------------------------
IMPLICIT NONE

! INPUT/OUTPUT VARIABLES
INTEGER :: idetpix,jdetpix,imaskpix,jmaskpix
REAL     (kind=8), dimension(:,:),pointer :: DET,MASKP
REAL     (kind=8) :: coeff
END SUBROUTINE SET_DET_MODEL0



END INTERFACE
!**************************
END MODULE MODEL0_PROC
!**************************




!**************************
MODULE MODEL0_DECLARATIONS
!**************************


INTERFACE



SUBROUTINE Model0(MASKP,DET,xsou,ysou,status)
!---------------------------------------------
IMPLICIT NONE

! INTERFACE variables
REAL     (kind=8), dimension(:,:), pointer :: DET,MASKP
REAL     (kind=8)                          :: xsou,ysou
INTEGER :: status
END SUBROUTINE Model0

END INTERFACE
!**************************
END MODULE MODEL0_DECLARATIONS
!**************************


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!           SUBROUTINES                 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!##################################
!  subroutines from MODEL0_AUX
!##################################



!===================================
SUBROUTINE AddPoint_MODEL0(x,y,point)
!===================================
USE MODEL0_VALS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER   :: point
REAL     (kind=8) :: x,y
!LOCAL VARIABLES
REAL     (kind=8) :: xi,yi

xi = aint(x)
yi = aint(y)
points(point,1) = x
points(point,2) = y
intpoints(point,1) = xi
intpoints(point,2) = yi
if(xi.ne.x)then
   !not on vertical line
   if(yi.eq.y)then
      ! only  on horizonthal line
      typepoints(point) = 2
   endif
else
   !on vertical line
   if(yi.eq.y)then
      !  on vertical and horizonthal line
      typepoints(point) = 4
   else
      ! only on vertical line
      typepoints(point) = 1
   endif
endif
!=======================
END SUBROUTINE AddPoint_MODEL0
!=======================

!============================================
SUBROUTINE FigVerif_MODEL0(ile,surface,Status)
!============================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
!!!!!!!!USE II_SKYIMAGE_GLOBVAR_MODULE
!!!!!!!!USE II_SKYIMAGE_USE_MODULE
USE MODEL0_VALS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: ile,Status
REAL     (kind=8) :: surface
!LOCAL VARIABLES
INTEGER :: ki,kj
REAL     (kind=8) :: x1,x2,y1,y2
INTEGER,dimension(10) :: tab
INTEGER,dimension(1) :: tt1

CHARACTER(len=100)  :: str100
CHARACTER(len=20)                 :: procName

Status = ISDC_OK
procName = 'FigVerif_MODEL0'



if(ile > 2)then
   !figure size verif

   x1 = minval(figure(1:ile,1)) 
   x2 = maxval(figure(1:ile,1)) 
   y1 = minval(figure(1:ile,2)) 
   y2 = maxval(figure(1:ile,2)) 
   surface = (x2-x1)*(y2-y1)

   if(surface.lt.FigSurfMinTol)then
      ile=0
   else ! not too small figure

      ! same points searching
      !numbers will be marked as 0 in the table tab

      tab(1:ile) = 1

      do ki=1,ile ! first point 
         if(tab(ki) > 0)then ! not excluded
            do kj=ki+1,ile ! second point
               if(tab(kj) > 0)then
                  if(abs(figure(ki,1)-figure(kj,1)) < figtol)then
                     if(abs(figure(ki,2)-figure(kj,2))< figtol)then
                        tab(kj) = 0
                     endif
                  endif
               endif
            enddo
         endif
      enddo

      if(sum(tab(1:ile)) < ile)then ! same points found
         do while(tab(ile)==0) ! ostatnie podwojne pkty usuniete
            ile = ile-1
         enddo
         if(ile==0)then

              write(str100,*)' figure verif problem  1'
              call MESSAGE(procName,str100, AnalModelError,Status)
              return
         endif
         do while(minval(tab(1:ile))==0) ! deleting duplicated 
            tt1=minloc(tab(1:ile))
            ki= tt1(1)
            ! point no ki duplicated - deleting - but not the last one 
            kj=ile ! point to be placed on ki place
            do while((tab(kj)==0).and.(kj > 0))
               kj = kj-1
            enddo
            if(kj==0)then
                write(str100,*)' figure verif problem  2'
              call MESSAGE(procName,str100, AnalModelError,Status)
               return
            endif
            ! kj - point to be put on the ki place
            figure(ki,1:2) = figure(kj,1:2)
            if(kj==ile)then
               ile=ile-1
            else
               ! point ile to be put at kj place
               figure(kj,1:2)=figure(ile,1:2)
               ile=ile-1
            endif
            tab(ki) = 1
            do while(tab(ile)==0) ! ostatnie podwojne pkty usuniete
               ile = ile-1
            enddo
         enddo
      endif ! same points found - deleting
  endif ! not too small figure
endif ! ile > 2

!==========================
END SUBROUTINE FigVerif_MODEL0
!==========================





!==================================
SUBROUTINE Triangle_MODEL0(surface)
!==================================
!calculates Triangle_MODEL0 surface
USE MODEL0_VALS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL     (kind=8) :: surface
!LOCAL VARIABLES
REAL     (kind=8),dimension(3) :: bok
REAL     (kind=8) :: p
surface = 0.

bok = 0.


bok(1) = sqrt((figure(vertex(1),1)-figure(vertex(2),1))**2+&
     (figure(vertex(1),2)-figure(vertex(2),2))**2 )   

bok(2) = sqrt((figure(vertex(2),1)-figure(vertex(3),1))**2+&
     (figure(vertex(2),2)-figure(vertex(3),2))**2 )  
 
bok(3) = sqrt((figure(vertex(3),1)-figure(vertex(1),1))**2+&
     (figure(vertex(3),2)-figure(vertex(1),2))**2 )  
 
p = 0.5*(bok(1)+bok(2)+bok(3))

surface = p*(p-bok(1))*(p-bok(2))*(p-bok(3))
surface = sqrt(surface)
!=======================
END SUBROUTINE Triangle_MODEL0
!========================

!==================================
SUBROUTINE Quad_MODEL0(surface)
!==================================
!calculates Quad_MODEL0rilateral surface
USE MODEL0_VALS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL     (kind=8) :: surface
!LOCAL VARIABLES
REAL     (kind=8),dimension(5) :: bok
REAL     (kind=8) :: p1,p2,s1,s2
surface = 0.

bok = 0.


bok(1) = sqrt((figure(vertex(1),1)-figure(vertex(2),1))**2+&
     (figure(vertex(1),2)-figure(vertex(2),2))**2 )  
 
bok(2) = sqrt((figure(vertex(2),1)-figure(vertex(3),1))**2+&
     (figure(vertex(2),2)-figure(vertex(3),2))**2 ) 
  
bok(3) = sqrt((figure(vertex(3),1)-figure(vertex(4),1))**2+&
     (figure(vertex(3),2)-figure(vertex(4),2))**2 ) 
 
bok(4) = sqrt((figure(vertex(4),1)-figure(vertex(1),1))**2+&
     (figure(vertex(4),2)-figure(vertex(1),2))**2 ) 
 
bok(5) = sqrt((figure(vertex(1),1)-figure(vertex(3),1))**2+&
     (figure(vertex(1),2)-figure(vertex(3),2))**2 )  
  
p1 = 0.5*(bok(1)+bok(2)+bok(5))
p2= 0.5*(bok(3)+bok(4)+bok(5))

s1 = p1*(p1-bok(1))*(p1-bok(2))*(p1-bok(5))
s1 = sqrt(s1)

s2 = p2*(p2-bok(3))*(p2-bok(4))*(p2-bok(5))
s2 = sqrt(s2)

surface = s1+s2
!=======================
END SUBROUTINE Quad_MODEL0
!========================

!==================================
SUBROUTINE Penta_MODEL0(surface)
!==================================
USE MODEL0_VALS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL     (kind=8) :: surface
!LOCAL VARIABLES
REAL     (kind=8),dimension(7) :: bok
REAL     (kind=8) :: p1,p2,p3,s1,s2,s3
surface = 0.

bok = 0.


bok(1) = sqrt((figure(vertex(1),1)-figure(vertex(2),1))**2+&
     (figure(vertex(1),2)-figure(vertex(2),2))**2 )   
bok(2) = sqrt((figure(vertex(2),1)-figure(vertex(3),1))**2+&
     (figure(vertex(2),2)-figure(vertex(3),2))**2 )   
bok(3) = sqrt((figure(vertex(3),1)-figure(vertex(4),1))**2+&
     (figure(vertex(3),2)-figure(vertex(4),2))**2 )  
bok(4) = sqrt((figure(vertex(4),1)-figure(vertex(5),1))**2+&
     (figure(vertex(4),2)-figure(vertex(5),2))**2 )  
bok(5) = sqrt((figure(vertex(5),1)-figure(vertex(1),1))**2+&
     (figure(vertex(5),2)-figure(vertex(1),2))**2 ) 

bok(6) = sqrt((figure(vertex(1),1)-figure(vertex(3),1))**2+&
     (figure(vertex(1),2)-figure(vertex(3),2))**2 )  
bok(7) = sqrt((figure(vertex(1),1)-figure(vertex(4),1))**2+&
     (figure(vertex(1),2)-figure(vertex(4),2))**2 ) 

  
p1 = 0.5*(bok(1)+bok(2)+bok(6))
p2= 0.5*(bok(6)+bok(3)+bok(7))
p3= 0.5*(bok(7)+bok(4)+bok(5))

s1 = p1*(p1-bok(1))*(p1-bok(2))*(p1-bok(6))
s1 = sqrt(s1)

s2 = p2*(p2-bok(6))*(p2-bok(3))*(p2-bok(7))
s2 = sqrt(s2)

s3 = p3*(p3-bok(7))*(p3-bok(4))*(p3-bok(5))
s3 = sqrt(s3)

surface = s1+s2+s3
!=======================
END SUBROUTINE Penta_MODEL0
!========================


!==================================
SUBROUTINE Hexa_MODEL0(surface)
!==================================
USE MODEL0_VALS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
REAL     (kind=8) :: surface
!LOCAL VARIABLES
REAL     (kind=8),dimension(9) :: bok
REAL     (kind=8) :: p1,p2,p3,p4,s1,s2,s3,s4

surface = 0.
bok = 0.


bok(1) = sqrt((figure(vertex(1),1)-figure(vertex(2),1))**2+&
     (figure(vertex(1),2)-figure(vertex(2),2))**2 )   
bok(2) = sqrt((figure(vertex(2),1)-figure(vertex(3),1))**2+&
     (figure(vertex(2),2)-figure(vertex(3),2))**2 )   
bok(3) = sqrt((figure(vertex(3),1)-figure(vertex(4),1))**2+&
     (figure(vertex(3),2)-figure(vertex(4),2))**2 )  
bok(4) = sqrt((figure(vertex(4),1)-figure(vertex(5),1))**2+&
     (figure(vertex(4),2)-figure(vertex(5),2))**2 )  
bok(5) = sqrt((figure(vertex(5),1)-figure(vertex(6),1))**2+&
     (figure(vertex(5),2)-figure(vertex(6),2))**2 ) 
bok(6) = sqrt((figure(vertex(6),1)-figure(vertex(1),1))**2+&
     (figure(vertex(6),2)-figure(vertex(1),2))**2 ) 

bok(7) = sqrt((figure(vertex(1),1)-figure(vertex(3),1))**2+&
     (figure(vertex(1),2)-figure(vertex(3),2))**2 )  
bok(8) = sqrt((figure(vertex(1),1)-figure(vertex(4),1))**2+&
     (figure(vertex(1),2)-figure(vertex(4),2))**2 ) 

bok(9) = sqrt((figure(vertex(1),1)-figure(vertex(5),1))**2+&
     (figure(vertex(1),2)-figure(vertex(5),2))**2 ) 
  
p1 = 0.5*(bok(1)+bok(2)+bok(7))
p2= 0.5*(bok(7)+bok(3)+bok(8))
p3= 0.5*(bok(8)+bok(4)+bok(9))
p4= 0.5*(bok(9)+bok(5)+bok(6))

s1 = p1*(p1-bok(1))*(p1-bok(2))*(p1-bok(7))
s1 = sqrt(s1)

s2 = p2*(p2-bok(7))*(p2-bok(3))*(p2-bok(8))
s2 = sqrt(s2)

s3 = p3*(p3-bok(8))*(p3-bok(4))*(p3-bok(9))
s3 = sqrt(s3)

s4 = p4*(p4-bok(9))*(p4-bok(5))*(p4-bok(6))
s4 = sqrt(s4)

surface = s1+s2+s3+s4
!=======================
END SUBROUTINE Hexa_MODEL0
!========================



!=======================
SUBROUTINE MAKEVERTEX_MODEL0(ile)
!=======================
USE MODEL0_VALS
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
INTEGER :: ile 
!LOCAL VARIABLES
REAL   (kind=8)     ::  x1,x2,a,b,y1,y2
INTEGER :: i,j,i1,i2,gora,dol,maxi,ival
REAL   (kind=8)     ::  val
INTEGER ,dimension(4) :: gorny,dolny
!REAL(kind=8) ,dimension(2) :: tab

gora = 0
dol = 0

x1 = minval(figure(1:ile,1)) 
x2 = maxval(figure(1:ile,1))

do i=1,ile
   if(figure(i,1) == x1)i1=i
   if(figure(i,1) == x2) i2=i
enddo
y1 = figure(i1,2)
y2 = figure(i2,2)
a = (y2-y1)/(x2-x1)
b = (y1+y2-a*(x1+x2))/2.

do i=1,ile
   if((i.ne.i1).and.(i.ne.i2))then
      if(a*figure(i,1)+b.lt.figure(i,2))then
         gora = gora+1
         gorny(gora) = i
      else
         dol = dol+1
         dolny(dol) = i
      endif
   endif
enddo
do i=1,gora ! sortowanie gory wedlug rosnacych x-ow
   maxi = 0
   val = 10000000.
   do j=i,gora
      if(figure(gorny(j),1).lt.val)then
         val =figure(gorny(j),1)
         maxi = j
      endif
   enddo
   if(maxi.ne.i)then ! zamiana
     ! tab = figure(gorny(maxi),1:2)
     ! figure(gorny(maxi),1:2) = figure(gorny(i),1:2)
     ! figure(gorny(i),1:2) = tab

      ival = gorny(maxi)
      gorny(maxi) = gorny(i)
      gorny(i) = ival
   endif
enddo

do i=1,dol ! sortowanie dolu wedlug malejacych x-ow
   maxi = 0
   val = 0.
   do j=i,dol
      if(figure(dolny(j),1).gt.val)then
         val =figure(dolny(j),1)
         maxi = j
      endif
   enddo
   if(maxi.ne.i)then ! zamiana
    !  tab = figure(dolny(maxi),1:2)
    !  figure(dolny(maxi),1:2) = figure(dolny(i),1:2)
    !  figure(dolny(i),1:2) = tab

      ival =dolny(maxi)
      dolny(maxi) = dolny(i)
      dolny(i) = ival
   endif
enddo      

vertex(1) = i1
do i=1,gora
   vertex(1+i)= gorny(i)
enddo
vertex(2+gora) = i2
do i=1,dol
   vertex(2+gora+i) = dolny(i)
enddo

!==========================
END SUBROUTINE MAKEVERTEX_MODEL0
!==========================




!##################################
!  subroutines from MODEL0_PROC     
!##################################


!=========================================
SUBROUTINE DetMaskPos_MODEL0(idetpix,jdetpix,&
     xmaskpos,ymaskpos,imaskpix,jmaskpix)
!=========================================
USE MODEL0_VALS
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER :: idetpix,jdetpix,imaskpix,jmaskpix
REAL(kind=8)     :: xmaskpos,ymaskpos
!LOCAL VARIABLES
REAL(kind=8)     :: xdp,ydp,xm,ym

xdp=idetpix+xsoushift 
ydp=jdetpix+ysoushift
!  Mask coordinates
if(theta .ne. 0.)then
   xm=xdp*cost-ydp*sent
   ym=xdp*sent+ydp*cost
   xm=xm/relpi
   ym=ym/relpi
else
   xm=xdp/relpi
   ym=ydp/relpi
endif

xmaskpos = xm+xmaskcentre  
ymaskpos = ym+ymaskcentre  


imaskpix=int(xmaskpos)+1
jmaskpix=int(ymaskpos)+1

!==========================
END SUBROUTINE DetMaskPos_MODEL0
!==========================



!=========================================
SUBROUTINE MaskDetPos_MODEL0(imaskpix,jmaskpix,&
     xdetpos,ydetpos,idetpix,jdetpix)
!=========================================
USE MODEL0_VALS
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER :: imaskpix,jmaskpix,idetpix,jdetpix
REAL(kind=8)     :: xdetpos,ydetpos
!LOCAL VARIABLES
REAL(kind=8)     :: xdp,ydp,xm,ym

xdp=imaskpix -xmaskcentre
ydp=jmaskpix -ymaskcentre
!  detector coordinates
if(theta .ne. 0.)then
   xm= xdp*cost + ydp*sent
   ym=-xdp*sent + ydp*cost
   xm=xm*relpi
   ym=ym*relpi
else
   xm=xdp*relpi
   ym=ydp*relpi
endif

xdetpos = xm-xsoushift 
ydetpos = ym-ysoushift


idetpix=int(xdetpos)+1
jdetpix=int(ydetpos)+1

!==========================
END SUBROUTINE MaskDetPos_MODEL0
!==========================

!=========================================
SUBROUTINE  MinMaxMaskPixel_MODEL0(i1,j1,i2,j2,&
            imin,imax,jmin,jmax)
!=========================================
USE MODEL0_VALS
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: i1,j1,i2,j2,imin,imax,jmin,jmax
!LOCAL VARIABLES
INTEGER,dimension(4) :: ipix,jpix
REAL(kind=8) :: xm,ym

call DetMaskPos_MODEL0(i1-1,j1-1,xm,ym,ipix(1),jpix(1))
call DetMaskPos_MODEL0(i2+1,j2+1,xm,ym,ipix(2),jpix(2))
call DetMaskPos_MODEL0(i2+1,j1-1,xm,ym,ipix(3),jpix(3))
call DetMaskPos_MODEL0(i1-1,j2+1,xm,ym,ipix(4),jpix(4))

imin = int(minval(ipix))
imax = int(maxval(ipix))
jmin = int(minval(jpix))
jmax = int(maxval(jpix))
imin = max(iMaskStart,imin)
imax = min(iMaskStop,imax)
jmin = max(jMaskStart,jmin)
jmax = min(jMaskStop,jmax)

!=========================================
END SUBROUTINE  MinMaxMaskPixel_MODEL0
!=========================================

!!$!===============================================
!!$SUBROUTINE  MaskBlok_MODEL0(MASKP,iMaskPix,jMaskPix,jMaskMax,ile)
!!$!===============================================
!!$ 
!!$USE MODEL0_VALS
!!$USE MODEL0_AUX
!!$USE II_SKYIMAGE_GLOBVAR_MODULE
!!$USE II_SKYIMAGE_CONTROL_MODULE
!!$USE II_SKYIMAGE_USE_MODULE
!!$IMPLICIT NONE
!!$
!!$!INPUT/OUTPUT VARIABLES
!!$REAL      (kind=8), dimension(:,:), pointer :: MASKP
!!$INTEGER ::iMaskPix,jMaskPix,jMaskMax,ile 
!!$!LOCAL VARIABLES
!!$INTEGER :: k,type
!!$
!!$type = MASKP(imaskpix,jmaskpix)
!!$ile=0
!!$k=jmaskpix+1
!!$do while((jmaskpix .le.jMaskMax).and.&
!!$     (type.eq.MASKP(imaskpix,k)))
!!$   k = k+1
!!$   ile = ile+1
!!$enddo
!!$ile = min(ile,4)
!!$!=========================
!!$END SUBROUTINE  MaskBlok_MODEL0    
!!$!=========================


!===============================================
SUBROUTINE MaskDetDefLines_MODEL0(imin,imax,jmin,jmax)
!================================================

USE MODEL0_VALS
USE MODEL0_AUX
IMPLICIT NONE

!INPUT/OUTPUT VARIABLES
INTEGER :: imin,imax,jmin,jmax 
!LOCAL VARIABLES
INTEGER :: lnum,pocz,kon,j1,j2,i1,i2,veri,horj,ipktow
REAL(kind=8) :: x1,x2,y1,y2,dx,dy,x,y,a,b

!Calculates cuts points


!-------------------------------
do lnum=1,4  !line loop
!-------------------------------

   PixLines(lnum)%inlist(:) = 1 ! different point list

   select case(lnum)
      case(2,3) 
         PixLines(lnum)%inlist(1) = 0
      case(4)
        PixLines(lnum)%inlist(1:2) = 0
   end select
   

   pocz = wzor(lnum,1) ! numer poczatku i konca boku
   kon = wzor(lnum,2)

   PixLines(lnum)%ilepktow = 2 ! min number of singular points on the line
   ipktow = 2

   x1 = xycorner(pocz,1)
   y1 = xycorner(pocz,2)
   PixLines(lnum)%pkty(1,1) = x1
   PixLines(lnum)%pkty(1,2) = y1

   x2 = xycorner(kon,1)
   y2 = xycorner(kon,2)
   PixLines(lnum)%pkty(2,1) = x2
   PixLines(lnum)%pkty(2,2) = y2
 
   i1 =int(min(x1,x2))  
   i2 = int(max(x1,x2))
  
   j1 = int(min(y1,y2))
   j2 = int(max(y1,y2))
  

   !---------------------------------------------
   if(.not.MaskDetRotationFlag)then !no rotation
   !---------------------------------------------

      ! line coeffs - not defined in this case
      if(i1.ne.i2)then 
         ! projected mask pixel side is a horizonthal line
         !vertical lines cut it
         do veri=i1+1,i2
            ipktow = ipktow+1
            PixLines(lnum)%pkty(ipktow,1) = real(veri)
            PixLines(lnum)%pkty(ipktow,2) =  xycorner(pocz,2)
         enddo
      else
         ! projected mask pixel side is a vertical line
         ! horizonthal lines cut it
         if(j1.ne.j2)then 
            do horj = j1+1,j2
               ipktow = ipktow +1
               PixLines(lnum)%pkty(ipktow,1) = xycorner(pocz,1)
               PixLines(lnum)%pkty(ipktow,2) =real(horj)
            enddo
      endif
      endif
   !------------------------------------------
   else ! rotation between mask/detector plane
   !------------------------------------------

      ! line coeffs 
      a = (y2-y1)/(x2-x1) 
      b = (y1+y2-a*(x1+x2))/2.
      PixLines(lnum)%a = a
      PixLines(lnum)%b = b
      if(i1.ne.i2)then 
         !  vertical lines cut this side 
         do veri=i1+1,i2
            !y1 = a*real(veri)+b
            y = y2-(x2-real(veri))*a
            ipktow =  ipktow+1
            PixLines(lnum)%pkty(ipktow,1) = real(veri)
            PixLines(lnum)%pkty(ipktow,2) =y
         enddo
      endif
      if(j1.ne.j2)then 
         ! horizonthal lines cut this side 
         do horj=j1+1,j2
           ! x1 = (real(horj)-b)/a
            x = x2-(y2-real(horj))/a
            if(aint(x).ne.x)then 
               ! not the vertical cut at the same time
               ipktow =  ipktow+1
               PixLines(lnum)%pkty(ipktow,1) = x
               PixLines(lnum)%pkty(ipktow,2) =real(horj)
            endif
         enddo
      endif

endif ! rotation        
 PixLines(lnum)%ilepktow  = ipktow
!---------------------
enddo  !line loop     
!---------------------


PixCorners%ilepktow = 0
PixCorners%inlist(:) = 1
!counting of internal det pixel vertex
ipktow = 0


if((imin < imax).and.( jmin < jmax)) then 
   !---------------------------------------------
   if(.not.MaskDetRotationFlag)then ! no rotation
   !---------------------------------------------
      do veri=imin+1,imax
         do horj=jmin+1,jmax
            ipktow = ipktow +1
            PixCorners%pkty(ipktow,1)= veri
            PixCorners%pkty(ipktow,2)= horj
         enddo
      enddo
   !---------------------------------------------
   else !rotation between mask and detector plane
   !---------------------------------------------
      if(theta.gt.0)then ! anticlockwise rotation
         do veri=imin+1,imax
         do horj=jmin+1,jmax
            dx = dble(veri)
            dy = dble(horj)
            if(PixLines(2)%a*dx + PixLines(2)%b .gt. dy)then 
               ! below 2 pixel line
               if(PixLines(4)%a*dx + PixLines(4)%b .lt. dy)then  
                  !above 4 line
                  if(PixLines(1)%a*dx + PixLines(1)%b .gt. dy)then 
                     ! above 1 line
                     if(PixLines(3)%a*dx + PixLines(3)%b .lt. dy)then 
                        !below 3 line
                       ipktow = ipktow + 1
                       PixCorners%pkty(ipktow,1)= dx
                       PixCorners%pkty(ipktow,2)= dy
                    endif
                 endif
              endif
           endif
        enddo
        enddo
     else !clockwise rotation
        do veri=imin+1,imax
         do horj=jmin+1,jmax
            dx = dble(veri)
            dy = dble(horj)
            if(PixLines(2)%a*dx + PixLines(2)%b .gt.dy)then 
               ! below 2  line
               if(PixLines(4)%a*dx + PixLines(4)%b .lt.dy)then  
                  !above 4 line
                  if(PixLines(1)%a*dx + PixLines(1)%b .lt. dy)then 
                     ! below 1 line
                     if(PixLines(3)%a*dx + PixLines(3)%b .gt. dy)then 
                        !above 3 line
                       ipktow = ipktow + 1
                       PixCorners%pkty(ipktow,1)= dx
                       PixCorners%pkty(ipktow,2)= dy
                    endif
                 endif
              endif
           endif
        enddo
        enddo
     endif !clockwise rotation

endif !rotation
endif !possible pkt inside det pixel      

PixCorners%ilepktow = ipktow
!=========================================
END SUBROUTINE MaskDetDefLines_MODEL0
!=========================================






!========================================
SUBROUTINE DetPixChoice_MODEL0(i1,i2,j1,j2,Status)
!=========================================
!!!zmianan procedury
USE ISDC
USE MIMOSA_CONTROL_MODULE
!!!!!!!!USE II_SKYIMAGE_GLOBVAR_MODULE
!!!!!!!!USE II_SKYIMAGE_USE_MODULE
USE MODEL0_VALS
USE MODEL0_AUX
IMPLICIT NONE
!INPUT/OUTPUT VARIABLES
INTEGER :: i1,i2,j1,j2,status
!LOCAL VARIABLES
INTEGER :: lnum,poi,point,ip,jp,num,ile,id,jd
REAL(kind=8) :: x,y,surf
INTEGER,dimension(10) :: tabnum
CHARACTER(len=20)                 :: procName

Status = ISDC_OK
procName = 'DetPixChoice_MODEL0'


points(:,:) = 0.
intpoints(:,:) = 0
typepoints(:) = 0
DetPixFraction(:,:) = 0.

point = 0

!different point list creation for pix lines
do lnum=1,4
do poi=1,PixLines(lnum)%ilepktow
   if(PixLines(lnum)%inlist(poi)==1)then 
      ! nowy pkt dorzucenie do listy
!!!!!zmiana
      if(point == mpsize)then
         call MESSAGE(procName,' too many points added ', &
              AnalModelError,Status)
         status = 10
         return
      endif
      point = point+1
      x = PixLines(lnum)%pkty(poi,1)
      y = PixLines(lnum)%pkty(poi,2)
      call AddPoint_MODEL0(x,y,point)
   endif ! on list point 
enddo
enddo

! adding internal corners to the point list
do poi=1,PixCorners%ilepktow
   if(PixCorners%inlist(poi)==1)then 
      ! nowy pkt dorzucenie do listy
!!!!!zmiana
      if(point == mpsize)then
         call MESSAGE(procName,' too many points added ', &
              AnalModelError,Status)
         status = 10
         return
      endif
      point = point+1
      x = PixCorners%pkty(poi,1)
      y = PixCorners%pkty(poi,2)
      call AddPoint_MODEL0(x,y,point)
   endif ! on list point 

enddo

!............................................................
! searching for all points belonging to the det pixel (ip,jp)
!............................................................
id = 0
do ip=i1,i2 ! det pixel loop 
id = id+1
jd = 0
do jp=j1,j2 
   jd = jd+1
   ile=0
   do num=1,point
      select case(typepoints(num))
         case(0) ! internal point
            if(intpoints(num,1)==ip)then
               if(intpoints(num,2)==jp)then
                  if(ile ==10)then
                     call MESSAGE(procName,' too many points in tabnum ', &
                          AnalModelError,Status)
                     status = 10
                     return
                  endif
                  ile=ile+1
                  Figure(ile,1:2) = points(num,1:2)
                  tabnum(ile) = num
               endif
            endif
         case(1) ! point on vertical line
            if(intpoints(num,2)==jp)then
               if((intpoints(num,1)-1==ip).or.(intpoints(num,1)==ip))then
                   if(ile ==10)then
                     call MESSAGE(procName,' too many points in tabnum ', &
                          AnalModelError,Status)
                     status = 10
                     return
                  endif
                 ile=ile+1
                  Figure(ile,1:2) = points(num,1:2)
                  tabnum(ile) = num
               endif
            endif
         case(2) ! point on horizonthal line
            if(intpoints(num,1)==ip)then
               if((intpoints(num,2)-1==jp).or.(intpoints(num,2)==jp))then
                  if(ile ==10)then
                     call MESSAGE(procName,' too many points in tabnum ', &
                          AnalModelError,Status)
                     status = 10
                     return
                  endif
                  ile=ile+1
                  Figure(ile,1:2) = points(num,1:2)
                  tabnum(ile) = num
               endif
            endif
         case(4) ! point on two lines
            if((intpoints(num,1)-1==ip).or.(intpoints(num,1)==ip))then
               if((intpoints(num,2)-1==jp).or.(intpoints(num,2)==jp))then
                   if(ile ==10)then
                     call MESSAGE(procName,' too many points in tabnum ', &
                          AnalModelError,Status)
                     status = 10
                     return
                  endif
                  ile=ile+1
                  Figure(ile,1:2) = points(num,1:2)
                  tabnum(ile) = num
               endif
            endif
      end select
   enddo
 
   call FigVerif_MODEL0(ile,surf,Status)
   if(Status.ne.0)return

   if(ile > 0) then
      ! calcul of surface fraction belonging to the det pixel
  
      if(.not.MaskDetRotationFlag)then
         ! no rotation - possible partitioning on squares only
         if(ile.ne.4)then
            call MESSAGE(procName,' bad fractioning', AnalModelError,Status)
            return
         endif
         ! square fraction
         DetPixFraction(id,jd) =surf

      else ! rotation
         surf = 0.
         select case(ile)
         case(1)
            if(typepoints(tabnum(1)).ne.4)then
               call MESSAGE(procName,'  problem : only one internal point ', &
                    AnalModelError,Status)
               return
            endif
         case(2)
            if((typepoints(tabnum(1)).ne.4).or.&
                 (points(tabnum(1),1).ne.points(tabnum(2),1)).or.&
                 (points(tabnum(1),2).ne.points(tabnum(2),2)))then
                call MESSAGE(procName,' problem : only two internal points', &
                    AnalModelError,Status)
               return
            endif
         case(3)
            !Triangle_MODEL0
            vertex(1:3) = (/ 1,2,3/)
            call Triangle_MODEL0(surf)
         case(4)
            !Quad_MODEL0rilateral
            call MakeVertex_MODEL0(4)
            call Quad_MODEL0(surf)
         case(5)
            call MakeVertex_MODEL0(5)
            call Penta_MODEL0(surf)
         case(6)
            call  MakeVertex_MODEL0(6)
            call Hexa_MODEL0(surf)
         end select
         DetPixFraction(id,jd) =surf
 endif !rotation
endif ! ile > 0
enddo
enddo

!==========================
END SUBROUTINE DetPixChoice_MODEL0
!==========================




!==================================================
SUBROUTINE SET_DET_MODEL0(idetpix,jdetpix,Det,MASKP,&
     imaskpix,jmaskpix,coeff)
!==================================================
USE MODEL0_VALS
!!!!!!!!USE II_SKYIMAGE_CONTROL_MODULE
!!!!!!!!USE II_SKYIMAGE_GLOBVAR_MODULE
!!!!!!!!USE II_SKYIMAGE_USE_MODULE
IMPLICIT NONE

! INPUT/OUTPUT VARIABLES
INTEGER :: idetpix,jdetpix,imaskpix,jmaskpix                
REAL     (kind=8), dimension(:,:),pointer :: DET,MASKP
REAL     (kind=8) :: coeff
!LOCAL VARIABLES


if(MASKP(imaskpix,jmaskpix)==1)then
   Det(idetpix,jdetpix) = Det(idetpix,jdetpix) +coeff
else
   Det(idetpix,jdetpix) = Det(idetpix,jdetpix) +0.*coeff
endif
  

!======================
END SUBROUTINE SET_DET_MODEL0
!======================



!#########################################
! SUBROUTINES FROM SOURCE_MODEL MODULE    
!#########################################

!=======================================================
SUBROUTINE Model0(MASKP,DET,xsou,ysou,status)
!========================================================
USE ISDC
USE MIMOSA_CONTROL_MODULE
!!!!!!!!USE II_SKYIMAGE_GLOBVAR_MODULE
!!!!!!!!USE II_SKYIMAGE_USE_MODULE
USE IBIS_IMAGING_PARAM
USE MODEL0_VALS
USE MODEL0_PROC
IMPLICIT NONE

! INTERFACE variables
REAL     (kind=8), dimension(:,:), pointer :: DET,MASKP
REAL     (kind=8)  :: xsou,ysou
INTEGER :: status
! Local Var
INTEGER,parameter :: DetDimPlus = 100
REAL     (kind=8), dimension(:,:), pointer :: DetPlus
INTEGER :: imaskdim,jmaskdim,iPlusDim,jPlusDim
INTEGER :: idetedim,jdetedim,iok
INTEGER :: imaskpix,jmaskpix,idetpix,jdetpix
INTEGER :: ipix,jpix,ord,ktory,iaxmin,iaxmax,jaxmin,jaxmax
INTEGER :: iMaskmin,iMaskmax,jMaskmin,jMaskmax
INTEGER :: iDetPixMin,jDetPixMin,ip,jp,id,jd
REAL     (kind=8) :: xdetpos,ydetpos
REAL     (kind=8) :: surface,val
INTEGER,dimension(4) :: order 
CHARACTER(len=20)                 :: procName


Status = ISDC_OK
procName = 'Model0'
if(DebugMode.eq.3)&
  call Message(procName,' ',ZeroError,Status)


DET = 0.


! CALL IBIS_IMAGING_INIT before 

relpi = MaElPixDim
theta = MasDetRotAng
if(theta.ne.0.)MaskDetRotationFlag = .true.

tgtheta = tan(theta*pi/180.0d0)
imaskdim = SIZE(MASKP,1)      
jmaskdim = SIZE(MASKP,2)      
  
xmaskcentre = real(imaskdim)/2.
ymaskcentre = real(jmaskdim)/2.   

!original detector dimension
IDeteDim = SIZE(DET,1)         
JDeteDim = SIZE(DET,2)         

iPlusDim = IDeteDim +2*DetDimPlus
jPlusDim = JDeteDim +2*DetDimPlus

xdetcentre = real(IPlusDim)/2. 
ydetcentre= real(JPlusDim)/2. 

allocate( DetPlus(iPlusDim,jPlusDim),stat = iok)
if(iok.ne.0)then
  call MESSAGE(procName,'allocation  problem',&
                   AllocError,Status)
  return

endif
!CodedPlus(:,:) = 0.
DetPlus(:,:) = 0.


wzor(1,1) = 1
wzor(1,2) = 2
wzor(2,1) = 2
wzor(2,2) = 3
wzor(3,1) = 3
wzor(3,2) = 4
wzor(4,1) = 4
wzor(4,2) = 1


!  Rotation parameters
thetarad=(theta/60.)*(PI/180.)
cost=COS(thetarad)
sent=SIN(thetarad)

!  Angle of source direction in degrees
alpha=SQRT(xsou*xsou + ysou*ysou)*PixAngDim
alphax=xsou*PixAngDim
alphay=ysou*PixAngDim
           
!  Shift to -xs,-ys          
! source distance from detector centre in det units

xsoushift=xsou-xdetcentre
ysoushift=ysou-ydetcentre

order(1) = 1
order(2) = 2
order(3) = 4
order(4) = 3


! det projection maximum size
!takes into account effective mask size
!given in i(j)MaskStart,i(j)MaskStop

call MinMaxMaskPixel_MODEL0(DetDimPlus,DetDimPlus,&
iPlusDim-DetDimPlus,jPlusDim-DetDimPlus,&
iMaskmin,iMaskmax,jMaskmin,jMaskmax)


!---------------------------------------------------
do imaskpix =  iMaskmin,iMaskmax !mask pixel i  loop
do jMaskPix = jMaskMin,jMaskMax  !mask pixel j  loop
!---------------------------------------------------
 

!      call MaskBlok_MODEL0(MASKP,iMaskPix,jMaskPix,jMaskMax,ileWTypie)
!      if(ilewtypie > ilemax)ilemax = ileWTypie
      
!!$if((iMaskPix==64).and.(jmaskpix==97))then
!!$print *,imaskpix,jmaskpix
!!$endif
      ord = 0 !vertex number

      do ipix = imaskpix-1,imaskpix  ! in-pix i loop
         do jpix = jmaskpix-1,jmaskpix  ! in-pix j loop  
            call  MaskDetPos_MODEL0(ipix,jpix,&
                 xdetpos,ydetpos,idetpix,jdetpix)
            ord = ord+1
            ktory = order(ord)
            xyCorner(ktory,1) = xdetpos
            xyCorner(ktory,2) = ydetpos
            ijCorner(ktory,1) = idetpix !det pixel number
            ijCorner(ktory,2) = jdetpix ! int(x)+1
         enddo ! in-pix i loop
      enddo  ! in-pix j loop 

      iAxMin = minval(ijcorner(:,1))-1 !min max number on axis
      iAxMax = maxval(ijcorner(:,1))-1
      jAxMin = minval(ijcorner(:,2))-1
      jAxMax = maxval(ijcorner(:,2))-1

      iDetPixMin = iAxMin+1 ! min i,j det pixel number
      jDetPixMin = jAxMin+1
    
      if((iAxMin==iAxMax).and.(jAxMin==jAxMax))then 
         ! mask pixel inside one det pixel
         val = 1.0
         call SET_DET_MODEL0(iDetpixMin,jDetpixMin,DetPlus,MASKP,&
              iMaskPix,jMaskPix,val)
      
      else ! mask pixel line cut(s) in the det pixel
         call MaskDetDefLines_MODEL0(iAxMin,iAxMax,jAxMin,jAxMax)
         call DetPixChoice_MODEL0(iAxMin,iAxMax,jAxMin,jAxMax,Status)
         if(Status.ne.0)return
         
         id = 0
         do ip = iAxMin,iAxMax
            id = id+1
            jd = 0
            do jp = jAxMin,jAxMax
               jd = jd+1
               surface = DetPixFraction(id,jd)
               call SET_DET_MODEL0(ip+1,jp+1,DetPlus,MASKP,&
                    iMaskPix,jMaskPix,Surface)
               !CodedPlus(ip+1,jp+1) = 1.

!!$               if((ip+1.eq.26+100).and.(jp+1.eq.110+100))then
!!$                  print *,imaskpix,jmaskpix
!!$               endif
         
            enddo
         enddo

 endif ! mask pixel line cut(s) 
ENDDO !detector pixel i  loop
ENDDO !detector pixel j  loop

Det(:,:) = DetPlus(DetDimPlus+1:iPlusDim-DetDimPlus,&
     DetDimPlus+1:jPlusDim-DetDimPlus)



!!$CodedZone(:,:) = CodedPlus(DetDimPlus+1:iPlusDim-DetDimPlus,&
!!$     DetDimPlus+1:jPlusDim-DetDimPlus)


deallocate(DetPlus)


!=================================== 
END SUBROUTINE Model0
!==================================
