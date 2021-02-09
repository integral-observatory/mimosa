
!|                   
!********************
Module ATTI_DEFS
!|*******************
!|MODULE ATTI_DEFS
!|*******************
!| CONTAINS GLOBAL DECLARATIONS FOR ATTITUDE CALCULATIONS
!| 
!| important variables:  
!|--------------------- 
!| point_max -  MAXIMAL POINTING NUMBER
!Integer,parameter::point_max=3000
Real (KIND=8),Parameter :: axes_precision = 1.0d-3
Real (KIND=8),Parameter :: zero_double    = 1.0d-15
!| err - ERROR NUMBER
!|   1            = INVALID NUMBER OF POINTINGS
!|   10+DATA_TYPE = OUTSIDE RANGE FOR THE GIVEN DATA_TYPE
!|   15           = INVALID DATA_TYPE
!|   20+DATA_TYPE = INVALID AXIS DIRECTION

Integer :: err


! IF IOView = 1 then image created will be as seen from outside of the sphere
! IF IOViev = -1 then image created will be as seen from inside of the sphere
Integer,Save:: IOView=1

! Image and Carte tangent point displacement with respect
! to the array centre 
Real(kind=4),save::ImaTanPointDisX=0.0,ImaTanPointDisY=0.0,&
                   MapTanPointDisX=0.0,MapTanPointDisY=0.0

!| point_num -  NUMBER OF POINTINGS
Integer,Save:: point_num

!| point_rot_mat -  ARRAY OF ROTATION MATRICES FOR EACH POINTING
!Real(KIND=8),Dimension(point_max,3,3)::point_rot_mat
Real(KIND=8),Dimension(:,:,:),pointer ::point_rot_mat
!| pointing_table -  AUXILIARY ARRAY FOR COMMUNICATION WITH ibis_atti_lib
!Real(KIND=4),Dimension(point_max,3)::pointing_table
!SPR 4357 - change below
Real(KIND=8),Dimension(:,:),pointer::pointing_table
!| array containing exclusion zones in the image
Logical,Dimension(:,:),pointer:: FilterImage

!|  Data_Type -  INPUT DATA TYPE 
!|   1 :  Xtel(ra dec) , theta - clockwise angle between
!|                       the Nord direction and Ztel
!|   2 :  Xtel(ra dec) , Ztel(ra dec) 
!|   3 :  Xtel(x,y,z)  , Ztel(x,y,z) in Equatorial System
!|   4 :  xXtel(l2,b2) , Ztel(l2,b2) in Galactic System
!| ...........................................
Integer,Save ::  Data_Type

!| pixel_ratio -  RATIO MASK_PIXEL_DIM/DISTANCE_MASK_TO_DETECTOR
!| ........................................................
Real (KIND=8),Save :: pixel_ratio,pixel_ang_size

  

! EQUATORIAL COORDINATES OF Xtel  AND  Ztel  AXES
!........................................................
Real(kind=8) , Dimension(3) ::Xtel_eq,Ztel_eq

! Galactic COORDINATES OF Xtel  AND  Ztel  AXES
!........................................................
Real(kind=8)  ::l2Xtel,b2Xtel,l2Ztel,b2Ztel

! GALACTIC SYSTEM TYPE
! 1 : B1950.0 FK4
! 2 : J2000.0 FK5
!........................
Integer,Save :: gal_sys_Type

! ALPHA , DELTA (IN DEG. AND RAD.) OF Xtel and Ztel AXES
!.........................................................

Real (KIND=8) :: alphaXtel_deg,deltaXtel_deg
Real (KIND=8) :: alphaXtel_rad,deltaXtel_rad

Real (KIND=8) :: alphaZtel_deg,deltaZtel_deg
Real (KIND=8) :: alphaZtel_rad,deltaZtel_rad


! ANTICLOCKWISE DIRECTION ANGLE (IN DEG. AND RAD.) BETWEEN
! NORD DIRECTION ( IN FIELD OF VIEW) AND +Ztel AXIS
!..........................................................
Real (KIND=8) :: theta_deg,theta_rad


! ROTATION MATRIX BETWEEN EQUATORIAL AND TELESCOPE SYSTEM
!.........................................................
Real (KIND=8) , Dimension(3,3) :: rot_mat

! THE NORD DIRECTION IN THE TANGENT PLANE
Real(kind=8) , Dimension(3) :: Z_nord_dir
Real(kind=8) :: alphaZnord,deltaZnord 



!| message_level -  MESSAGE LEVEL
!| 0 - NO MESSAGES FROM LIBRARY ROUTINES
!| 1 -  DIAGNOSTIC MASSAGES
!|.......................................
Integer ,Save:: message_level=1

!| Interface_Type - INTERFACE TYPE
!| 0 - LABO INPUT FILE
!| 1 - ISDC INTERFACE
!|......................................
Integer ,Save:: Interface_Type=0

! AUXILIARY VARIABLES
!.....................................................
Real(kind=8) , Dimension(3) :: vector,vector1,vector2

Real (KIND=8),Save :: dpi,d360_rad,d270_rad,d180_rad,d90_rad
Real (KIND=8),Save :: deg_to_rad,rad_to_deg,x_tol
Real(KIND=4) ,Save :: spi,s360_rad,s270_rad,s180_rad,s90_rad
! SPR 4357 - change below
Real (KIND=8),Save :: sdeg_to_rad,srad_to_deg

Real (KIND=8) :: yz_radius,d_one,norm
! SLA FUNCTIONS INTERFACE
!.....................................

INTERFACE
! SLA FUNCTIONS DECLARATIONS

      Function sla_drange(a)
         Real  (kind=8)  :: a,sla_drange
      End Function sla_drange

      Function sla_dranrm(a)
         Real  (kind=8)  :: a,sla_dranrm
      End Function sla_dranrm

       Function sla_dsep(a,b,c,d)
         Real  (kind=8)  :: a,b,c,d,sla_dsep
      End Function sla_dsep

      Function sla_dbear(a,b,c,d)
         Real  (kind=8)  :: a,b,c,d,sla_dbear
      End Function sla_dbear

END INTERFACE

!********************
End Module ATTI_DEFS
!********************



!********************
MODULE ATTI_INTERNAL
!|****************************************************
!|MODULE ATTI_INTERNAL - CONTAINS AUXILIARY PROCEDURES
!|                       You can use it with the statement 
!|                       USE MODULE ATTI_INTERNAL
!|*****************************************************

CONTAINS

!|  A. TAN projection :  (alpha,delta) in deg. ==> tangent plane coordinates 
!|     Versions:                          
!|       1.  DOUBLE PRECISION OF CALL LIST PARAMEMETRS  
!|       2.  SINGLE PRECISION "   "    "     "          
!|       3.  POINTING NUMBER VERIFICATION AND ROTATION MATRIX SEARCHING
!|          
!|      D_CONVN_ADYZ: 1
!|      S_CONV_ADYZ : 2+3
!|      S_CONVN_ADYZ :2
!|                      
 



!....................................................................
Subroutine D_CONVN_ADYZ(alpha_deg,delta_deg,y,z,rot_matrix,pixel,status)
!....................................................................


Use ATTI_DEFS
Implicit None

! INPUT VARIABLES

Real(kind=8)     :: alpha_deg,delta_deg    ! RA DEC COORD. IN deg
Real(kind=8),Dimension(3,3),Intent(in)   :: rot_matrix ! rotation matrix
REAL(kind=8),Intent(in)::pixel
! OUTPUT VARIABLES

Real(kind=8),Intent(out) :: y,z                   ! TANGENT PLANE COORD.
Integer     ,Intent(out) :: status                

! status = 1  when  x< 0 
!       The arrival direction is below the detector plane.
!       In this case the direction cosines are conserved but 
!       the resulting point (y,z) is placed outside the radius 
!       yz_radius.To recover the original position do
!             n= dsqrt(y**2+z**2)
!             y = y*yz_radius/n
!             z = z*yz_radius/n
! status = 2 when x near 0
!      This means a very great direction angle with respect to 
!      the line-of-sight-axis X
! status = 0  otherwise


! AUXILIARY VARIABLES

Real(kind=8) :: a_rad,d_rad,n1

! NO VERIFICATION OF POINTING NUMBER
! ROTATION MATRIX GIVEN IN  rot_matrix
 
STATUS = 0

If((alpha_deg  .Lt.0.0d0).or.( alpha_deg .Gt. 360.0d0)) Then
   If(message_level ==1) Print *,' ATTENTION RA < 0 !!'
   alpha_deg = mod(alpha_deg+360.0d0,360.0d0)

Endif

If(Abs(delta_deg).Gt.90.0d0) Then
   If(message_level ==1) Print *,' ATTENTION abs(DEC) > 90 !!'
   delta_deg = mod(delta_deg,90.0d0)
Endif

a_rad = alpha_deg*deg_to_rad
d_rad = delta_deg*deg_to_rad

! transformation (a,d) --> (x,y,z) in Equatorial System
Call sla_dcs2c(a_rad,d_rad,vector)


! rotation (x,y,z) Equat. --> (x,y,z) Telescope
Call sla_dmxv(rot_matrix,vector,vector1)
!normalisation
Call sla_dvn(vector1,vector,norm)

! NORMALISATION FACTOR
! transformation to pixel units

n1 = vector(1)

If(n1.Le.0.d0)Then
    ! ARRIVAL DIRECTION BELOW THE DETECTOR PLANE
    STATUS = 1
    n1 =  (1.d0+10.0d-10 -n1**2) / yz_radius
Else
   If(vector(1).Le.x_tol) Then
   ! TOO BIG ARRIVAL ANGLE
   STATUS = 2
   Endif
Endif

n1=n1*pixel

vector = vector/n1


!  reflection of Y axis

y = IOView*vector(2)
z = vector(3)


!...........................
End Subroutine D_CONVN_ADYZ
!...........................

!....................................................................
Subroutine S_CONVN_ADYZ(alpha_deg,delta_deg,y,z,rot_matrix,pixel,status)
!....................................................................



Use ATTI_DEFS
Implicit None

! INPUT VARIABLES

Real(kind=4),Intent(in) :: alpha_deg,delta_deg      ! RA DEC COORD. IN deg
Real(kind=8),Dimension(3,3),Intent(in)::rot_matrix  ! ROTATION MATRIX
REAL(kind=4),Intent(in)::pixel
! OUTPUT VARIABLES

Real(kind=4),Intent(out) :: y,z                   ! TANGENT PLANE COORD.
Integer     ,Intent(out) :: status                

! status = 1  when  x< 0 
!       The arrival direction is below the detector plane.
!       In this case the direction cosines are conserved but 
!       the resulting point (y,z) is placed outside the radius 
!       yz_radius.To recover the original position do
!             n= dsqrt(y**2+z**2)
!             y = y*yz_radius/n
!             z = z*yz_radius/n
! status = 2 when x near 0
!      This means a very great direction angle with respect to 
!      the line-of-sight-axis X
! status = 0  otherwise


! AUXILIARY VARIABLES

Real(kind=8) ::a,d, a_rad,d_rad,n1


 
STATUS = 0
a=alpha_deg
d=delta_deg

If((a.Lt.0.0d0).or.(alpha_deg .Gt. 360.0d0)) Then
   If(message_level ==1) Print *,' ATTENTION RA not in (0,360)' 
  a = mod(a,360.0d0)   
Endif

If(Abs(d).Gt.90) Then
   If(message_level ==1) Print *,' ATTENTION abs(DEC) > 90 !!'
   d = mod(d,90.0d0)
Endif

a_rad = a*deg_to_rad
d_rad = d*deg_to_rad

! transformation (a,d) --> (x,y,z) in Equatorial System
Call sla_dcs2c(a_rad,d_rad,vector)

! rotation (x,y,z) Equat. --> (x,y,z) Telescope
Call sla_dmxv(rot_matrix,vector,vector1)


! NORMALISATION FACTOR
! transformation to pixel units

n1 = vector1(1)

If(n1.Le.0.d0)Then
    ! ARRIVAL DIRECTION BELOW THE DETECTOR PLANE
    STATUS = 1
    n1 =  (1.d0+10.0d-10 -n1**2) / yz_radius
Else
   If(vector1(1).Le.x_tol) Then
   ! TOO BIG ARRIVAL ANGLE
   STATUS = 2
   Endif
Endif

n1=n1*pixel

vector1 = vector1/n1


!  reflection of Y axis
! RETOUR TO THE SINGLE PRECISION
y = IOView*vector1(2)
z = vector1(3)


!..........................
End Subroutine S_CONVN_ADYZ
!..........................

!........................................................
Subroutine D_CONV_LB_AD(l_deg,b_deg,alpha_deg,delta_deg)
!........................................................

Use ATTI_DEFS
Implicit None

!INPUT/OUTPUT VARIABLES
Real(kind=8)::l_deg,b_deg,alpha_deg,delta_deg
!LOCAL VARIABLES
Real(kind=8):: l_rad,b_rad,cosdltsinalf,sindlt,dlt,alf



l_rad = l_deg*deg_to_rad
b_rad = b_deg*deg_to_rad
call sla_galeq(l_rad,b_rad,alf,dlt)

alpha_deg=alf*rad_to_deg
delta_deg=dlt*rad_to_deg
!...........................
END Subroutine D_CONV_LB_AD
!...........................

!........................................................
Subroutine D_CONV_AD_LB(alpha_deg,delta_deg,l_deg,b_deg)
!........................................................

Use ATTI_DEFS
Implicit None

!INPUT/OUTPUT VARIABLES
Real(kind=8)::alpha_deg,delta_deg,l_deg,b_deg
!LOCAL VARIABLES
Real(kind=8):: dlt,alf,cosbsinl,sinb,b,l

dlt = delta_deg*deg_to_rad
alf = alpha_deg*deg_to_rad
call sla_eqgal(alf,dlt,l,b)

l_deg = l*rad_to_deg
b_deg = b*rad_to_deg
!...........................
END Subroutine D_CONV_AD_LB
!...........................


!........................................................
Subroutine S_CONV_LB_AD(l_deg,b_deg,alpha_deg,delta_deg)
!........................................................

Use ATTI_DEFS
Implicit None

!INPUT/OUTPUT VARIABLES
Real(kind=4)::l_deg,b_deg,alpha_deg,delta_deg
!LOCAL VARIABLES
Real(kind=8):: l_rad,b_rad,cosdltsinalf,sindlt,dlt,alf



l_rad = dble(l_deg)*deg_to_rad
b_rad = dble(b_deg)*deg_to_rad
call sla_galeq(l_rad,b_rad,alf,dlt)

alpha_deg=real(alf)*rad_to_deg
delta_deg=real(dlt)*rad_to_deg
!...........................
END Subroutine S_CONV_LB_AD
!...........................

!........................................................
Subroutine S_CONV_AD_LB(alpha_deg,delta_deg,l_deg,b_deg)
!........................................................

Use ATTI_DEFS
Implicit None

!INPUT/OUTPUT VARIABLES
Real(kind=4)::alpha_deg,delta_deg,l_deg,b_deg
!LOCAL VARIABLES
Real(kind=8):: dlt,alf,cosbsinl,sinb,b,l

dlt = dble(delta_deg)*deg_to_rad
alf = dble(alpha_deg)*deg_to_rad
call sla_eqgal(alf,dlt,l,b)

l_deg = real(l)*rad_to_deg
b_deg = real(b)*rad_to_deg
!...........................
END Subroutine S_CONV_AD_LB
!...........................

!|  B. TAN back-projection : tangent plane coordinates ==>  (alpha,delta) in deg.
!|     Versions:                          
!|        1.  DOUBLE PRECISION OF CALL LIST PARAMEMETRS  
!|        2.  SINGLE PRECISION "   "    "     "             
!|          
!|      D_CONVN_YZAD: 1
!|      S_CONVN_YZAD :2
!| 



!................................................................
Subroutine D_CONVN_YZAD(y,z,alpha_deg,delta_deg,rot_matrix,pixel)
!................................................................


Use ATTI_DEFS
Implicit None

! INPUT VARIABLES
Real(kind=8),Intent(in) :: y,z                     ! COORD. IN TANGENT PLANE
Real(kind=8),Dimension(3,3),Intent(in)::rot_matrix ! ROTATION MATRIX
REAL(kind=8),Intent(in)::pixel                     !pixel size
!OUTPUT VARIABLES
Real(kind=8),Intent(out) :: alpha_deg,delta_deg    ! RA DEC (deg) 


! AUXILIARY VARIABLES
Real(kind=8) :: a_rad,d_rad



! TRANSFORMATION TO A VECTOR OF THE DIRECTION IN THE SATELLITE SYSTEM
vector(1) = 1
vector(2) =  IOView*(y*pixel)
vector(3) = (z*pixel)

! BACK PROJECTION TO EQUATORIAL SYSTEM
Call sla_dimxv(rot_matrix,vector,vector1)
Call sla_dvn(vector1,vector2,norm)

! (a,d) --> (x,y,z)
Call sla_dcc2s(vector2,a_rad,d_rad)

! normalisation
a_rad = sla_dranrm(a_rad)
d_rad = sla_drange(d_rad)

alpha_deg = a_rad*rad_to_deg
delta_deg = d_rad*rad_to_deg

alpha_deg = MOD(alpha_deg+360.d0,360.0d0)



!..........................
END SUBROUTINE D_CONVN_YZAD
!...........................

!.............................................................
SUBROUTINE S_CONVN_YZAD(y,z,alpha_deg,delta_deg,rot_matrix,pixel)
!.............................................................

Use ATTI_DEFS
Implicit None

! INPUT VARIABLES
Real(kind=4),Intent(in) :: y,z                     ! COORD. IN TANGENT PLANE
Real(kind=8),Dimension(3,3),Intent(in)::rot_matrix ! ROTATION MATRIX
REAL(kind=4),Intent(in)::pixel
!OUTPUT VARIABLES
Real(kind=4),Intent(out) :: alpha_deg,delta_deg ! RA DEC (deg) 


! AUXILIARY VARIABLES
Real(kind=8) :: a_rad,d_rad


! TRANSFORMATION TO A VECTOR OF THE DIRECTION IN THE SATELLITE SYSTEM
vector(1) = 1
vector(2) = IOView*(y*pixel)
vector(3) = (z*pixel)

! BACK PROJECTION TO EQUATORIAL SYSTEM
Call sla_dimxv(rot_matrix,vector,vector1)


! (a,d) --> (x,y,z)
Call sla_dcc2s(vector1,a_rad,d_rad)

! normalisation
a_rad = sla_dranrm(a_rad)
d_rad = sla_drange(d_rad)

alpha_deg = a_rad*rad_to_deg
delta_deg = d_rad*rad_to_deg

alpha_deg = MOD(alpha_deg+360.d0,360.0d0)



!...........................
END SUBROUTINE S_CONVN_YZAD
!...........................

!|C. Cartesian projection/back-projection                                         
!|    Versions:                          
!|        1.  DOUBLE PRECISION OF CALL LIST PARAMEMETRS  
!|        2.  SINGLE PRECISION "   "    "     "          
!|          
!|     D_CAR_ADYZ   : 1
!|     S_CAR_ADYZ   : 2
!|     D_CAR_YZAD :   1
!|  
!...........................................................
SUBROUTINE D_CAR_ADYZ(alpha,delta,xn,yn,aref,dref,pixel_ang)
!...........................................................
Use ATTI_DEFS
Implicit None
!INPUT VARIABLES
Real(kind=8),Intent(in)::alpha,delta,aref,dref
REAL(kind=8),Intent(in)::pixel_ang !pixel ang size
!OUTPUT VARIABLES
Real(kind=8),Intent(out)::xn,yn

!Local variables
Real(kind=8) :: adiff,a,ar,zref

adiff = alpha-aref
xn =adiff
zref = mod(aref+180.0d0,360.0d0)
if(aref==180.0d0)zref=360.0d0
if((aref.ge.0).and.(aref.le.180.0d0))then
   !situation 1
  
   if(alpha .gt.zref)then
      xn = adiff-360.0d0
   endif
else
   !situation 2
   if(alpha.le.zref)then
      xn = 360.0d0+adiff
   endif
endif


xn = IOView*xn
xn = xn/pixel_ang

yn =( delta-dref)/(pixel_ang)

!........................
END SUBROUTINE D_CAR_ADYZ
!........................

!...................................................
SUBROUTINE D_CAR_YZAD(xn,yn,alpha,delta,aref,dref,pixel_ang)
!...................................................
Use ATTI_DEFS
Implicit None
!INPUT VARIABLES
Real(kind=8),Intent(in)::xn,yn,aref,dref
REAL(kind=8),Intent(in)::pixel_ang !pixelANGULAR size
!OUTPUT VARIABLES
Real(kind=8),Intent(out)::alpha,delta


alpha = Mod(IOView*xn*pixel_ang+aref+360.0d0,360.0d0)
!alpha =xn*pixel_ang+aref
delta = yn*pixel_ang+dref


!........................
END SUBROUTINE D_CAR_YZAD
!........................


!...................................................
SUBROUTINE S_CAR_ADYZ(alpha,delta,xn,yn,aref,dref,pixel_ang)
!...................................................
Use ATTI_DEFS
Implicit None
!INPUT VARIABLES
Real(kind=4),Intent(in)::alpha,delta,aref,dref
REAL(kind=4),Intent(in)::pixel_ang
!OUTPUT VARIABLES
Real(kind=4),Intent(out)::xn,yn
!Local variables
Real(kind=4) :: adiff,a,ar

!!$adiff = abs(alpha-aref)
!!$if(adiff.lt.180.0)then
!!$   xn = IOView*(alpha-aref)
!!$else
!!$   if(alpha.lt.180.0)then
!!$      a = alpha+180.0
!!$      ar= aref
!!$   else
!!$      a = alpha
!!$      ar = aref+180.0
!!$      endif
!!$xn = IOView*(a-ar)
!!$endif
xn = xn/pixel_ang
yn =( delta-dref)/(pixel_ang)
!........................
END SUBROUTINE S_CAR_ADYZ
!........................

!|  D. MORE GENERAL SUBROUTINES FOR PROJECTION/BACK PROJECTION 
!|     D_CONVN_YZYZ - TAN backprojection and projection
!|     PROJECTION TAN(CAR) - TAN/CAR backprojection and projection
!| 
!...................................................................
SUBROUTINE D_CONVN_YZYZ(y1,z1,y2,z2,rot_matrix,pixel1,pixel2,status)
!..................................................................
Use ATTI_DEFS
Implicit None
! INPUT VARIABLES
Real(kind=8) :: y1,z1            ! COORD.( IN THE IMAGE) TANGENT PLANE
Real(kind=8),Dimension(3,3)::rot_matrix ! ROTATION MATRIX
REAL(kind=8),Intent(in)::pixel1,pixel2 !pixel sizes
!OUTPUT VARIABLES
Real(kind=8),Intent(out) :: y2,z2    ! COORD.( IN THE CARTE) TANGENT PLANE

Integer    ,Intent(out)  :: status                

! status = 1  when  x< 0 
!       The arrival direction is below the detector plane.
!       In this case the direction cosines are conserved but 
!       the resulting point (y,z) is placed outside the radius 
!       yz_radius.To recover the original position do
!             n= dsqrt(y**2+z**2)
!             y = y*yz_radius/n
!             z = z*yz_radius/n
! status = 2 when x near 0
!      This means a very great direction angle with respect to 
!      the line-of-sight-axis X
! status = 0  otherwise

! AUXILIARY VARIABLES
Real(kind=8)::n1

vector(1) = 1.0d0
vector(2) =  IOView*(y1*pixel1)
vector(3) = (z1*pixel1)

Call sla_dmxv(rot_matrix,vector,vector1)
! NORMALISATION FACTOR
! transformation to pixel units

n1 = vector1(1)

If(n1.Le.0.d0)Then
    ! ARRIVAL DIRECTION BELOW THE DETECTOR PLANE
    STATUS = 1
    n1 =  (1.d0+10.0d-10 -n1**2) / yz_radius
Else
   If(vector1(1).Le.x_tol) Then
   ! TOO BIG ARRIVAL ANGLE
   STATUS = 2
   Endif
Endif

n1=n1*pixel2

vector1 = vector1/n1


! reflection of Y axis
! RETOUR TO THE SINGLE PRECISION
y2 = IOView*vector1(2)
z2 = vector1(3)
!..............................
END SUBROUTINE D_CONVN_YZYZ
!..............................

!.......................................................
SUBROUTINE PROJECTION(is,tc,pixel1,pixel2,&
                     pixel_ang2,aref,dref,lref,bref,&
                      rot_mat1,rot_mat3,rot_mat4,&
                      ri,rj,xn,yn,status)
!.......................................................

 Use ATTI_DEFS
! USE MIMOSA_GLOBVAR_MODULE
 USE MIMOSA_CONTROL_MODULE
 Implicit None
!INPUT VARIABLES
INTEGER,Intent(in)::is,tc
REAL(kind=8),Intent(in)::pixel1,pixel2,&
                         pixel_ang2,aref,dref,lref,bref
Real(kind=8),Dimension(3,3),Intent(in)::rot_mat1,rot_mat3,rot_mat4
REAL(kind=8),Intent(in)::ri,rj
!OUTPUT VARIABLES
REAL(kind=8),Intent(out)::xn,yn
INTEGER,Intent(out)::status

!AUXILIARY VARIABLES
REAL(kind=8)::ALPHA,DELTA,l,b

status = 0



 if(tc.eq.1)then
    !TAN-TAN
    if(is.eq.0) then

       Call D_CONVN_YZYZ(ri,rj,xn,yn,rot_mat3,&
                        pixel1,pixel2,status)
    else
       Call D_CONVN_YZYZ(ri,rj,xn,yn,rot_mat4,&
                        pixel2,pixel1,status)
    endif
 else
    !...-CAR
    if(is.eq.0) then !image-to_sky proj
       Call  D_CONVN_YZAD(ri,rj,alpha,delta,rot_mat1,pixel1)
       if(EquaGal==0)then
          !TAN--CAR
          Call D_CAR_ADYZ(alpha,delta,xn,yn,aref,dref,pixel_ang2)
       else
          !GAL--CAR
          call D_CONV_AD_LB(alpha,delta,l,b)
          Call D_CAR_ADYZ(l,b,xn,yn,lref,bref,pixel_ang2)
       endif
     else ! sky_to_image
         if(EquaGal==0)then
            !TAN--CAR
            Call D_CAR_YZAD(ri,rj,alpha,delta,&
                        aref,dref,pixel_ang2)
         else
            Call D_CAR_YZAD(ri,rj,l,b,&
                 lref,bref,pixel_ang2)
            call D_CONV_LB_AD(l,b,alpha,delta)
         endif
         Call D_CONVN_ADYZ(alpha,delta,xn,yn,&
                 rot_mat1,pixel1,status)
        
    endif
 endif
!aproj = alpha
!dproj=delta
!lproj=l
!bproj = b
!..........................
END SUBROUTINE PROJECTION
!..........................




!......................................
SUBROUTINE add_val0(is,im,i,j,in,jn,time,image,ScwExposition,imagevar,&
                     carte,cartevar,exposure)
!......................................
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,i,j,in,jn
REAL(kind=4) :: time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer,optional::imagevar
!INPUT/output VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer,optional::cartevar ,exposure
 if(is.eq.0) then
   carte(in,jn) = carte(in,jn)+image(i,j)
   if(im ==1)&
      cartevar(in,jn) = cartevar(in,jn)+imagevar(i,j)
   exposure(in,jn) = exposure(in,jn)+ScwExposition(i,j)
 else
   carte(i,j) = carte(i,j)+image(in,jn)
   if(im ==1)&
      cartevar(i,j) = cartevar(i,j)+imagevar(in,jn)
    exposure(i,j) = exposure(i,j)+ScwExposition(in,jn)
 endif
!......................
END SUBROUTINE add_val0
!......................



!.............................................
SUBROUTINE add_val1(is,im,i,j,in,jn,an,bn,time,image,ScwExposition,imagevar,&
                    carte,cartevar,exposure)
!.............................................
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) :: time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer,optional::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer,optional::cartevar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var

i2 = in+1
j2 = jn+1

if(is==0)then
!image-to-sky
   var = imagevar(i,j)
   if(var > 0.)then

      carte(i2,j2) = carte(i2,j2)+an*bn*image(i,j)/var
      carte(in,jn) = carte(in,jn)+(1-an)*(1-bn)*image(i,j)/var
      carte(in,j2) = carte(in,j2)+(1-an)*bn*image(i,j)/var
      carte(i2,jn) = carte(i2,jn)+an*(1-bn)*image(i,j)/var
  
      ! exposure
      exposure(i2,j2) = exposure(i2,j2)+an*bn*ScwExposition(i,j)
      exposure(in,jn) = exposure(in,jn)+(1-an)*(1-bn)*ScwExposition(i,j)
      exposure(in,j2) = exposure(in,j2)+(1-an)*bn*ScwExposition(i,j)
      exposure(i2,jn) = exposure(i2,jn)+an*(1-bn)*ScwExposition(i,j)
  

      ! variance
      a1 = an*bn
      a2 = (1-an)*(1-bn)
      a3 = (1-an)*bn
      a4 =  an*(1-bn)
      cartevar(i2,j2) = cartevar(i2,j2)+time**2*a1**2/var
      cartevar(in,jn) = cartevar(in,jn)+time**2*a2**2/var
      cartevar(in,j2) = cartevar(in,j2)+time**2*a3**2/var
      cartevar(i2,jn) = cartevar(i2,jn)+time**2*a4**2/var
  endif !var > 0.

else
!sky-to-image
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
      carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      cartevar(i,j) = cartevar(i,j)+ time**2*a1**2/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1=an*bn
      carte(i,j) =  carte(i,j)+a1*image(i2,j2)/var
      cartevar(i,j) = cartevar(i,j)+ time**2*a1**2/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      cartevar(i,j) = cartevar(i,j)+ time**2*a1**2/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      cartevar(i,j) = cartevar(i,j)+ time**2*a1**2/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif



endif

!......................
END SUBROUTINE add_val1
!......................



!.............................................
SUBROUTINE add_val01_point(is,i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,&
                           idim2,jdim2,i,j,in,jn,&
                           image,imagevar,carte,cartevar)
!.............................................
Implicit None
!INPUT/output VARIABLES 
INTEGER      ::     is,idim2,jdim2,i,j,in,jn
INTEGER      ::    i_dim_inp,j_dim_inp,i_dim_out,j_dim_out
REAL(kind=4),dimension(:,:),pointer::image,imagevar
Real(KIND=4)  ,dimension(:,:),pointer::carte,cartevar

!AUXILIARY VARIABLES
integer ::i1,j1,i2,j2,ia,ib,ja,jb,ii,jj,ipos,jpos
Real(KIND=4) :: maxvalue,var
if(is.eq.0)then
   i1=i;j1=j;i2=in;j2=jn
   ia = max(1,i1-1); ib = min(i1+1,i_dim_inp)
   ja = max(1,j1-1); jb = min(j1+1,j_dim_inp)
else
   i1=in;j1=jn;i2=i;j2=j
   ia = max(1,i1-1); ib = min(i1+1,i_dim_out)
   ja = max(1,j1-1); jb = min(j1+1,j_dim_out)
endif


maxvalue = image(ia,ja)
ipos = ia
jpos = ja
do ii=ia,ib
do jj=ja,jb
   if(image(ii,jj) > maxvalue)then
      maxvalue = image(ii,jj)
      ipos = ii
      jpos = jj
   endif
enddo
enddo
var = imagevar(ipos,jpos)
if(var > 0.)then
   carte(i2,j2) =carte(i2,j2)+maxvalue/var
   cartevar(i2,j2) =cartevar(i2,j2)+1./var
endif
!......................
END SUBROUTINE add_val01_point
!......................

!.............................................
SUBROUTINE add_val01(is,im,i,j,in,jn,time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im
INTEGER,Intent(in)::i,j,in,jn

REAL(kind=4) :: time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var

i2 = in+1
j2 = jn+1

if(is==0)then
!image-to-sky
   var = imagevar(i,j)
   if(var > 0.)then
      carte(in,jn) = carte(in,jn)+image(i,j)/var
      exposure(in,jn) = exposure(in,jn)+ScwExposition(i,j)
      cartevar(in,jn) = cartevar(in,jn)+1./var
      normvar(in,jn) = normvar(in,jn)+1./var
  endif !var > 0.

else
!sky-to-image
   var = imagevar(in,jn)
   if(var > 0.)then
      
      carte(i,j) =  carte(i,j)+image(in,jn)/var
      cartevar(i,j) = cartevar(i,j)+ 1./var
      normvar(i,j) = normvar(i,j)+ 1./var
      exposure(i,j) =  exposure(i,j)+ScwExposition(in,jn)
   endif



 


endif

!......................
END SUBROUTINE add_val01
!......................



!.............................................
SUBROUTINE add_val01_withtime(is,im,expo,i,j,in,jn,time,image,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,expo
INTEGER,Intent(in)::i,j,in,jn

REAL(kind=4) :: time   
REAL(kind=4),dimension(:,:),pointer::image
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var

i2 = in+1
j2 = jn+1

if(is==0)then
!image-to-sky
   var = imagevar(i,j)
   if(var > 0.)then
      carte(in,jn) = carte(in,jn)+image(i,j)/var
      exposure(in,jn) = exposure(in,jn)+time
      cartevar(in,jn) = cartevar(in,jn)+time**2/var
      normvar(in,jn) = normvar(in,jn)+time/var
  endif !var > 0.

else
!sky-to-image
   var = imagevar(in,jn)
   if(var > 0.)then
      
      carte(i,j) =  carte(i,j)+image(in,jn)/var
      cartevar(i,j) = cartevar(i,j)+ time**2/var
      normvar(i,j) = normvar(i,j)+ time/var


   endif

  

 ! exposure
  exposure(i,j) =  exposure(i,j)+time
 


endif

!......................
END SUBROUTINE add_val01_withtime
!......................

!.............................................
SUBROUTINE add_val2(is,im,i,j,in,jn,an,bn,time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) :: time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var

i2 = in+1
j2 = jn+1

if(is==0)then
!image-to-sky
   var = imagevar(i,j)
   if(var > 0.)then

      carte(i2,j2) = carte(i2,j2)+an*bn*image(i,j)/var
      carte(in,jn) = carte(in,jn)+(1-an)*(1-bn)*image(i,j)/var
      carte(in,j2) = carte(in,j2)+(1-an)*bn*image(i,j)/var
      carte(i2,jn) = carte(i2,jn)+an*(1-bn)*image(i,j)/var
  
      ! exposure
      exposure(i2,j2) = exposure(i2,j2)+an*bn*ScwExposition(i,j)
      exposure(in,jn) = exposure(in,jn)+(1-an)*(1-bn)*ScwExposition(i,j)
      exposure(in,j2) = exposure(in,j2)+(1-an)*bn*ScwExposition(i,j)
      exposure(i2,jn) = exposure(i2,jn)+an*(1-bn)*ScwExposition(i,j)
  

      ! variance
      a1 = an*bn
      a2 = (1-an)*(1-bn)
      a3 = (1-an)*bn
      a4 =  an*(1-bn)
      cartevar(i2,j2) = cartevar(i2,j2)+(time*a1)**2/var
      cartevar(in,jn) = cartevar(in,jn)+(time*a2)**2/var
      cartevar(in,j2) = cartevar(in,j2)+(time*a3)**2/var
      cartevar(i2,jn) = cartevar(i2,jn)+(time*a4)**2/var

      ! normalisation de variance 
      normvar(i2,j2) =  normvar(i2,j2)  +time*an*bn/var
      normvar(in,jn) = normvar(in,jn)  + time*(1-an)*(1-bn)/var
      normvar(in,j2) = normvar(in,j2)+time*(1-an)*bn/var
      normvar(i2,jn) = normvar(i2,jn)+time*an*(1-bn)/var
 
  endif !var > 0.

else
!sky-to-image
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)**2/var
      normvar(i,j) = normvar(i,j)+time*a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1=an*bn
      carte(i,j) =  carte(i,j)+a1*image(i2,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)**2/var
      normvar(i,j) = normvar(i,j)+ time*a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)**2/var
      normvar(i,j) = normvar(i,j)+ time* a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)**2/var
      normvar(i,j) = normvar(i,j)+  time*a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif



endif

!......................
END SUBROUTINE add_val2
!......................


!.............................................
SUBROUTINE add_val3(is,im,i,j,in,jn,an,bn,time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) :: time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer,optional::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer,optional::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var,dev

i2 = in+1
j2 = jn+1

if(is==0)then
!image-to-sky
   var = imagevar(i,j)
   if(var > 0.)then
      dev = sqrt(var)
      a1 = an*bn
      a2 = (1-an)*(1-bn)
      a3 = (1-an)*bn
      a4 =  an*(1-bn)
      carte(i2,j2) = carte(i2,j2)+a1*image(i,j)/var
      carte(in,jn) = carte(in,jn)+a2*image(i,j)/var
      carte(in,j2) = carte(in,j2)+a3*image(i,j)/var
      carte(i2,jn) = carte(i2,jn)+a4*image(i,j)/var
  
      ! exposure
      exposure(i2,j2) = exposure(i2,j2)+a1*ScwExposition(i,j)
      exposure(in,jn) = exposure(in,jn)+a2*ScwExposition(i,j)
      exposure(in,j2) = exposure(in,j2)+a3*ScwExposition(i,j)
      exposure(i2,jn) = exposure(i2,jn)+a4*ScwExposition(i,j)
  

      ! variance
     
      cartevar(i2,j2) = cartevar(i2,j2)+(time*a1)/dev
      cartevar(in,jn) = cartevar(in,jn)+(time*a2)/dev
      cartevar(in,j2) = cartevar(in,j2)+(time*a3)/dev
      cartevar(i2,jn) = cartevar(i2,jn)+(time*a4)/dev

      ! normalisation de variance 
      normvar(i2,j2) =  normvar(i2,j2)  + time*a1/var
      normvar(in,jn) = normvar(in,jn)   + time*a2/var
      normvar(in,j2) = normvar(in,j2)   + time*a3/var
      normvar(i2,jn) = normvar(i2,jn)   + time*a4/var
 
  endif !var > 0.

else
!sky-to-image
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)/sqrt(var)
      normvar(i,j) = normvar(i,j)+time*a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1=an*bn
      carte(i,j) =  carte(i,j)+a1*image(i2,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)/sqrt(var)
      normvar(i,j) = normvar(i,j)+ time*a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)/sqrt(var)
      normvar(i,j) = normvar(i,j)+ time* a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)/sqrt(var)
      normvar(i,j) = normvar(i,j)+  time*a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif

 
 


endif

!......................
END SUBROUTINE add_val3
!......................

!.............................................
SUBROUTINE add_val9(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k,ii1,jj1,ii2,jj2
REAL(kind=8)::a1,a2,a3,a4,var,dev,sumvar
REAL(kind=8),dimension(2,2) :: ab
i2 = in+1
j2 = jn+1

ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
  
 sumvar=0.0d0  
 
  
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       ab(1,1) = a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1= an*bn
      ab(2,2) =  a1/sqrt(var)
      carte(i,j) = carte(i,j)+ a1*image(i2,j2)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+ a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      ab(1,2) =  a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      sumvar = sumvar+  a1**2/var
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      ab(2,1) =a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif

  
cartevar(i,j) = cartevar(i,j)+sumvar
sumvar = 0.
do ii1=1,2
do jj1=1,2
   do ii2=1,2
   do jj2=1,2
      if((ii1.ne.ii2).or.(jj1.ne.jj2))then
         sumvar = sumvar+ ab(ii1,jj1)*ab(ii2,jj2)
      endif
   enddo
   enddo
enddo
enddo
cartevar(i,j) = cartevar(i,j)+sumvar*0.75
endif

!......................
END SUBROUTINE add_val9
!......................

!.............................................
SUBROUTINE add_val9a(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k,ii1,jj1,ii2,jj2
REAL(kind=8)::a1,a2,a3,a4,var,dev,sumvar
REAL(kind=8),dimension(2,2) :: ab

i2 = in+1
j2 = jn+1

ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
  
 sumvar=0.0d0  
 
  
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       ab(1,1) = a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1= an*bn
      ab(2,2) =  a1/sqrt(var)
      carte(i,j) = carte(i,j)+ a1*image(i2,j2)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+ a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      ab(1,2) =  a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      sumvar = sumvar+  a1**2/var
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      ab(2,1) =a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif

  
cartevar(i,j) = cartevar(i,j)+sumvar
sumvar = 0.
do ii1=1,2
do jj1=1,2
   do ii2=1,2
   do jj2=1,2
      if((ii1.ne.ii2).or.(jj1.ne.jj2))then
         sumvar = sumvar+ ab(ii1,jj1)*ab(ii2,jj2)*covpattern(ii1,jj1,ii2,jj2)
      endif
   enddo
   enddo
enddo
enddo
cartevar(i,j) = cartevar(i,j)+sumvar
endif

!......................
END SUBROUTINE add_val9a

!......................
!.............................................
SUBROUTINE add_val10(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::ii,jj,ii1,ii2,jj1,jj2,zak
REAL(kind=8)::m,sumdi,sumui,f,varcoeff,sumvar,sumexpo,var,sumerr,fac,p,q
REAL(kind=8),dimension(-2:2,-2:2) :: ab


ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
 
!LS fit of flux
sumdi = 0.
sumui = 0.
sumvar= 0.
sumexpo = 0.

zak=5
do ii=-zak,zak
do jj=-zak,zak
   ii1 = in+ii
   jj1 = jn+jj
   var=imagevar(ii1,jj1)
   if(var >0.)then
      sumdi = sumdi+image(ii1,jj1)*ucoeff(ii,jj)
      sumui = sumui+ucoeff(ii,jj)**2
      sumexpo=sumexpo+ScwExposition(ii1,jj1)
     
   endif
enddo
enddo

if(sumui > 0)then

   f = sumdi/sumui
   m=(an**2+bn**2)/ psfsigmacoeff
   m = exp(-m)
   f = f*m

   carte(i,j) =carte(i,j)+f*time
   normvar(i,j)= normvar(i,j)+time
   exposure(i,j) =exposure(i,j)+ sumexpo*sumui*m

   varcoeff = (time*m/sumui)**2
   sumvar = 0.

   do ii1=-zak,zak
   do jj1=-zak,zak
      do ii2=-zak,zak
      do jj2=-zak,zak
        ! if((ii1.ne.ii2).or.(jj1.ne.jj2))then
            ii=ii1-ii2
            jj=jj1-jj2
            if((abs(ii).le.zak).and.(abs(jj).le.zak))then
               fac =  sqrt(imagevar(in+ii1,jn+jj1)*imagevar(in+ii2,jn+jj2))
               if(fac > 0)then
                  q = ucoeff(ii1,jj1)*ucoeff(ii2,jj2)*ucoeff(ii,jj)*fac
                  sumvar = sumvar+ q
               endif
            endif
        ! endif
      enddo
      enddo
   enddo
   enddo
cartevar(i,j) = cartevar(i,j)+sumvar*varcoeff
endif
endif

!......................
END SUBROUTINE add_val10
!......................

!.............................................
SUBROUTINE add_val14(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::ii,jj,ii1,ii2,jj1,jj2,zak
REAL(kind=8)::m,sumdi,sumui,f,varcoeff,sumvar,sumexpo,var,sumerr,fac,p,q
REAL(kind=8)::x0,y0,u,uc
REAL(kind=8),dimension(-5:5,-5:5) :: ab


ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
 
!LS fit of flux
sumdi = 0.
sumui = 0.
sumvar= 0.
sumexpo = 0.
x0 = in+an
y0 = jn+bn
zak=2
do ii=-zak,zak
do jj=-zak,zak
   ii1 = in+ii
   jj1 = jn+jj
   var=imagevar(ii1,jj1)
   if(var >0.)then
      u=(ii1*1.-x0)**2+(jj1*1.-y0)**2
      uc = exp(-u/psfsigmacoeff)
      ab(ii,jj) = uc
      sumdi = sumdi+image(ii1,jj1)*uc
      sumui = sumui+uc**2
      sumexpo=sumexpo+ScwExposition(ii1,jj1)*uc
     
   endif
enddo
enddo

if(sumui > 0)then

   f = sumdi/sumui
!!$   m=(an**2+bn**2)/ psfsigmacoeff
!!$   m = exp(-m)
!!$   f = f*m

   carte(i,j) =carte(i,j)+f*time
   normvar(i,j)= normvar(i,j)+time
   exposure(i,j) =exposure(i,j)+ sumexpo

   varcoeff = (time/sumui)**2
   sumvar = 0.

   do ii1=-zak,zak
   do jj1=-zak,zak
      do ii2=-zak,zak
      do jj2=-zak,zak
        ! if((ii1.ne.ii2).or.(jj1.ne.jj2))then
            ii=ii1-ii2
            jj=jj1-jj2
            if((abs(ii).le.zak).and.(abs(jj).le.zak))then
               fac =  sqrt(imagevar(in+ii1,jn+jj1)*imagevar(in+ii2,jn+jj2))
               if(fac > 0)then
                  q = ab(ii1,jj1)*ab(ii2,jj2)*ucoeff(ii,jj)*fac
                  sumvar = sumvar+ q
               endif
            endif
        ! endif
      enddo
      enddo
   enddo
   enddo
cartevar(i,j) = cartevar(i,j)+sumvar*varcoeff
endif
endif

!......................
END SUBROUTINE add_val14
!......................

!.............................................
SUBROUTINE add_val11(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::ii,jj,ii1,ii2,jj1,jj2
REAL(kind=8)::m,sumdi,sumui,f,varcoeff,sumvar,sumexpo,var,sumerr,fac,p,q
REAL(kind=8),dimension(-2:2,-2:2) :: ab


ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
 
sumdi = 0.
sumui = 0.
sumvar= 0.
sumexpo = 0.
!CHI2 fit
do ii=-2,2
do jj=-2,2
   ii1 = in+ii
   jj1 = jn+jj
   var=imagevar(ii1,jj1)
   if(var >0.)then
      sumdi = sumdi+image(ii1,jj1)*ucoeff(ii,jj)/var
      sumui = sumui+ucoeff(ii,jj)**2/var
      sumexpo=sumexpo+ScwExposition(ii1,jj1)
     
   endif
enddo
enddo

if(sumui > 0)then

   f = sumdi/sumui
   m=(an**2+bn**2)/ psfsigmacoeff
   m = exp(-m)
   f = f*m

   carte(i,j) =carte(i,j)+f
   normvar(i,j)= normvar(i,j)+1
   exposure(i,j) =exposure(i,j)+ sumexpo*sumui*m

   varcoeff = m**2/sumui**2
  

   sumvar = 0.
   do ii1=-2,2
   do jj1=-2,2
      do ii2=-2,2
      do jj2=-2,2
        ! if((ii1.ne.ii2).or.(jj1.ne.jj2))then
            ii=ii1-ii2
            jj=jj1-jj2
            if((abs(ii).le.2).and.(abs(jj).le.2))then
               fac =  sqrt(imagevar(in+ii1,jn+jj1)*imagevar(in+ii2,jn+jj2))
               if(fac > 0)then
                  q = ucoeff(ii1,jj1)*ucoeff(ii2,jj2)*ucoeff(ii,jj)/fac
                  sumvar = sumvar+ q
               endif
            endif
        ! endif
      enddo
      enddo
   enddo
   enddo
cartevar(i,j) = cartevar(i,j)+sumvar*varcoeff
endif
endif

!......................
END SUBROUTINE add_val11
!......................

!.............................................
SUBROUTINE add_val12(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::ii,jj,ii1,ii2,jj1,jj2,zak
REAL(kind=8)::m,sumdi,sumui,f,varcoeff,sumvar,sumexpo,var,sumerr,fac,p,q
REAL(kind=8),dimension(-2:2,-2:2) :: ab


zak=1

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
 
!LS fit of flux
sumdi = 0.
sumui = 0.
sumvar= 0.
sumexpo = 0.

do ii=-zak,zak
do jj=-zak,zak
   ii1 = in+ii
   jj1 = jn+jj
   var=imagevar(ii1,jj1)
   if(var >0.)then
      sumdi = sumdi+image(ii1,jj1)*ucoeff(ii,jj)
      sumui = sumui+ucoeff(ii,jj)**2
      sumexpo=sumexpo+ScwExposition(ii1,jj1)
     
   endif
enddo
enddo

if(sumui > 0)then

   f = sumdi/sumui
   m=(an**2+bn**2)/ psfsigmacoeff
   m = exp(-m)
   f = f*m

   carte(i,j) =carte(i,j)+f*time
   normvar(i,j)= normvar(i,j)+time
   exposure(i,j) =exposure(i,j)+ sumexpo*sumui*m

   varcoeff = (time*m/sumui)**2
   sumvar = 0.

   do ii1=-zak,zak
   do jj1=-zak,zak
      do ii2=-zak,zak
      do jj2=-zak,zak
        ! if((ii1.ne.ii2).or.(jj1.ne.jj2))then
            ii=ii1-ii2
            jj=jj1-jj2
            if((abs(ii).le.zak).and.(abs(jj).le.zak))then
               fac =  sqrt(imagevar(in+ii1,jn+jj1)*imagevar(in+ii2,jn+jj2))
               if(fac > 0)then
                  q = ucoeff(ii1,jj1)*ucoeff(ii2,jj2)*ucoeff(ii,jj)*fac
                  sumvar = sumvar+ q
               endif
            endif
        ! endif
      enddo
      enddo
   enddo
   enddo
cartevar(i,j) = cartevar(i,j)+sumvar*varcoeff
endif
endif

!......................
END SUBROUTINE add_val12
!......................

!.............................................
SUBROUTINE add_val13(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::ii,jj,ii1,ii2,jj1,jj2,zak
REAL(kind=8)::m,sumdi,sumui,f,varcoeff,sumvar,sumexpo,var,sumerr,fac,p,q
REAL(kind=8),dimension(-2:2,-2:2) :: ab


ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
 
sumdi = 0.
sumui = 0.
sumvar= 0.
sumexpo = 0.
zak = 2
!CHI2 fit
do ii=-zak,zak
do jj=-zak,zak
   ii1 = in+ii
   jj1 = jn+jj
   var=imagevar(ii1,jj1)
   if(var >0.)then
      sumdi = sumdi+image(ii1,jj1)*ucoeff(ii,jj)/var
      sumui = sumui+ucoeff(ii,jj)**2/var
      sumexpo=sumexpo+ScwExposition(ii1,jj1)
     
   endif
enddo
enddo

if(sumui > 0)then

   f = sumdi/sumui
   m=(an**2+bn**2)/ psfsigmacoeff
   m = exp(-m)

   f = f*m

   carte(i,j) =carte(i,j)+f
   normvar(i,j)= normvar(i,j)+1
   exposure(i,j) =exposure(i,j)+ sumexpo*sumui*m

   varcoeff = m**2/sumui**2
  

   sumvar = 0.
   do ii1=-zak,zak
   do jj1=-zak,zak
      do ii2=-zak,zak
      do jj2=-zak,zak
        ! if((ii1.ne.ii2).or.(jj1.ne.jj2))then
            ii=ii1-ii2
            jj=jj1-jj2
            if((abs(ii).le.zak).and.(abs(jj).le.zak))then
               fac =  sqrt(imagevar(in+ii1,jn+jj1)*imagevar(in+ii2,jn+jj2))
               if(fac > 0)then
                  q = ucoeff(ii1,jj1)*ucoeff(ii2,jj2)*ucoeff(ii,jj)/fac
                  sumvar = sumvar+ q
               endif
            endif
        ! endif
      enddo
      enddo
   enddo
   enddo
cartevar(i,j) = cartevar(i,j)+sumvar*varcoeff
endif
endif

!......................
END SUBROUTINE add_val13
!......................


!.............................................
SUBROUTINE add_val15(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::ii,jj,ii1,ii2,jj1,jj2,zak
REAL(kind=8)::m,sumdi,sumui,f,varcoeff,sumvar,sumexpo,var,sumerr,fac,p,q
REAL(kind=8),dimension(-2:2,-2:2) :: ab
REAL(kind=8)::x0,y0,u,uc

ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
 
sumdi = 0.
sumui = 0.
sumvar= 0.
sumexpo = 0.
x0 = in+an
y0 = jn+bn

zak = 2
!CHI2 fit
do ii=-zak,zak
do jj=-zak,zak
   ii1 = in+ii
   jj1 = jn+jj
   var=imagevar(ii1,jj1)
   if(var >0.)then
      u=(ii1*1.-x0)**2+(jj1*1.-y0)**2
      uc = exp(-u/psfsigmacoeff)
      ab(ii,jj) = uc
      sumdi = sumdi+image(ii1,jj1)*uc/var
      sumui = sumui+uc**2/var
      sumexpo=sumexpo+ScwExposition(ii1,jj1)*uc
     
   endif
enddo
enddo

if(sumui > 0)then

   f = sumdi/sumui
!!$   m=(an**2+bn**2)/ psfsigmacoeff
!!$   m = exp(-m)
!!$
!!$   f = f*m

   carte(i,j) =carte(i,j)+f
   normvar(i,j)= normvar(i,j)+1
   exposure(i,j) =exposure(i,j)+ sumexpo*sumui

   varcoeff = 1./sumui**2
  

   sumvar = 0.
   do ii1=-zak,zak
   do jj1=-zak,zak
      do ii2=-zak,zak
      do jj2=-zak,zak
        ! if((ii1.ne.ii2).or.(jj1.ne.jj2))then
            ii=ii1-ii2
            jj=jj1-jj2
            if((abs(ii).le.zak).and.(abs(jj).le.zak))then
               fac =  sqrt(imagevar(in+ii1,jn+jj1)*imagevar(in+ii2,jn+jj2))
               if(fac > 0)then
                  q = ab(ii1,jj1)*ab(ii2,jj2)*ab(ii,jj)/fac
                  sumvar = sumvar+ q
               endif
            endif
        ! endif
      enddo
      enddo
   enddo
   enddo
cartevar(i,j) = cartevar(i,j)+sumvar*varcoeff
endif
endif

!......................
END SUBROUTINE add_val15
!......................

!.............................................
SUBROUTINE add_val8(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var,dev,sumvar
REAL(kind=8),dimension(2,2) :: ab
i2 = in+1
j2 = jn+1

ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
  
 sumvar=0.0d0  
 
  
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       ab(1,1) = a1/var
      carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      sumvar = sumvar+ a1/sqrt(var)
      normvar(i,j) = normvar(i,j)+a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1= an*bn
      ab(2,2) = a1/var
      carte(i,j) = carte(i,j)+ a1*image(i2,j2)/var
      sumvar = sumvar+ a1/sqrt(var)
      normvar(i,j) = normvar(i,j)+ a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      ab(1,2) = a1/var
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      sumvar = sumvar+ a1/sqrt(var)
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      ab(2,1) =a1/var
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      sumvar = sumvar+ a1/sqrt(var)
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif

  
cartevar(i,j) = cartevar(i,j)+sumvar**2

endif

!......................
END SUBROUTINE add_val8
!......................


!.............................................
SUBROUTINE add_val4(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var,dev
REAL(kind=8),dimension(2,2) :: ab
i2 = in+1
j2 = jn+1

ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
  
   
 
  
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       ab(1,1) = a1/var
      carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (a1)**2/var
      normvar(i,j) = normvar(i,j)+a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1= an*bn
      ab(2,2) = a1/var
      carte(i,j) =  carte(i,j)+a1*image(i2,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (a1)**2/var
      normvar(i,j) = normvar(i,j)+ a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      ab(1,2) = a1/var
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (a1)**2/var
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      ab(2,1) =a1/var
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (a1)**2/var
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif

  
!!$   cartevar(i,j) = cartevar(i,j) + 2.*detmeanval*&
!!$        (CovarTab(in,jn)*ab(1,1)*ab(1,2)/OffAxisCorr(in,jn) + &
!!$        CovarTab(in,jn)*ab(1,1)*ab(2,1)/OffAxisCorr(in,jn)  + &
!!$        CovarTab(i2,jn)*ab(2,1)*ab(2,2)/OffAxisCorr(i2,jn)  + &
!!$        CovarTab(in,j2)*ab(1,2)*ab(2,2)/OffAxisCorr(in,j2)  + &
!!$        Covar1Tab(in,jn)*ab(1,1)*ab(2,2)/OffAxisCorr(in,jn) + &
!!$        Covar1Tab(i2,jn)*ab(2,1)*ab(1,2)/OffAxisCorr(i2,jn))  &
!!$                           *1800./time**2
!!$  



endif

!!$if(carte(4321,576) .ne.0.)then
!!$   print *,'add_val4',i,j,carte(i,j),carte(4321,576)
!!$endif

!......................
END SUBROUTINE add_val4
!......................



!.............................................
SUBROUTINE add_val17(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var,dev
REAL(kind=8),dimension(2,2) :: ab
i2 = in+1
j2 = jn+1

ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
  
   

  
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       ab(1,1) = a1/var
      carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (a1)**2/var
      normvar(i,j) = normvar(i,j)+a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1= an*bn
      ab(2,2) = a1/var
      carte(i,j) =  carte(i,j)+a1*image(i2,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (a1)**2/var
      normvar(i,j) = normvar(i,j)+ a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      ab(1,2) = a1/var
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (a1)**2/var
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      ab(2,1) =a1/var
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (a1)**2/var
      normvar(i,j) = normvar(i,j)+  a1/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif

  
   cartevar(i,j) = cartevar(i,j) + 2.*detmeanval*&
        (CovarTab(in,jn)*ab(1,1)*ab(1,2)/OffAxisCorr(in,jn) + &
        CovarTab(in,jn)*ab(1,1)*ab(2,1)/OffAxisCorr(in,jn)  + &
        CovarTab(i2,jn)*ab(2,1)*ab(2,2)/OffAxisCorr(i2,jn)  + &
        CovarTab(in,j2)*ab(1,2)*ab(2,2)/OffAxisCorr(in,j2)  + &
        Covar1Tab(in,jn)*ab(1,1)*ab(2,2)/OffAxisCorr(in,jn) + &
        Covar1Tab(i2,jn)*ab(2,1)*ab(1,2)/OffAxisCorr(i2,jn))  &
                           *1800./time**2
  



endif



!......................
END SUBROUTINE add_val17
!......................

!.............................................
SUBROUTINE add_val16(num,detmeanval,is,im,i,j,in,jn,an,bn,&
                    time,image,ScwExposition,imagevar,&
                    carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image,ScwExposition
REAL(kind=4),dimension(:,:),pointer::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k,ii1,jj1,ii2,jj2
REAL(kind=8)::a1,a2,a3,a4,var,dev,sumvar
REAL(kind=8),dimension(2,2) :: ab
i2 = in+1
j2 = jn+1

ab = 0 ! SPR 4246 corrected in version 4.8

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
  
   
 
  
  
sumvar=0.0d0  
 
  
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       ab(1,1) = a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+a1**2/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,jn)
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1= an*bn
      ab(2,2) =  a1/sqrt(var)
      carte(i,j) = carte(i,j)+ a1*image(i2,j2)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+ a1**2/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,j2)
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      ab(1,2) =  a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      sumvar = sumvar+  a1**2/var
      normvar(i,j) = normvar(i,j)+  a1**2/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(in,j2)
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      ab(2,1) =a1/sqrt(var)
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      sumvar = sumvar+ a1**2/var
      normvar(i,j) = normvar(i,j)+  a1**2/var
      exposure(i,j) =  exposure(i,j)+a1*ScwExposition(i2,jn)
   endif

  
cartevar(i,j) = cartevar(i,j)+sumvar
sumvar = 0.
do ii1=1,2
do jj1=1,2
   do ii2=1,2
   do jj2=1,2
      if((ii1.ne.ii2).or.(jj1.ne.jj2))then
         sumvar = sumvar+ ab(ii1,jj1)*ab(ii2,jj2)
      endif
   enddo
   enddo
enddo
enddo
cartevar(i,j) = cartevar(i,j)+sumvar*0.75



endif



!......................
END SUBROUTINE add_val16
!......................

!.............................................
SUBROUTINE add_val4_withtime(num,detmeanval,is,im,expo,&
     i,j,in,jn,an,bn,time,image,imagevar,&
     carte,cartevar,normvar,exposure)
!.............................................
USE MIMOSA_GLOBVAR_MODULE
Implicit None
!INPUT VARIABLES 
INTEGER,Intent(in)::is,im,expo,num
INTEGER,Intent(in)::i,j,in,jn
REAL(kind=8),Intent(in)::an,bn
REAL(kind=4) ::detmeanval, time   
REAL(kind=4),dimension(:,:),pointer::image
REAL(kind=4),dimension(:,:),pointer,optional::imagevar,exposure
!OUTPUT VARIABLES
REAL(kind=4),dimension(:,:),pointer::carte
REAL(kind=4),dimension(:,:),pointer,optional::cartevar,normvar
!AUXILIARY VARIABLES
INTEGER::i2,j2,k
REAL(kind=8)::a1,a2,a3,a4,var,dev
REAL(kind=8),dimension(2,2) :: ab
i2 = in+1
j2 = jn+1

if(is==0)then
!image-to-sky
  print *,' error in covar method'
  return

else
!sky-to-image
  
   
 
  
   var = imagevar(in,jn)
   if(var > 0.)then
      a1 = (1-an)*(1-bn)
       ab(1,1) = a1*time/var
      carte(i,j) =  carte(i,j)+a1*image(in,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)**2/var
      normvar(i,j) = normvar(i,j)+time*a1/var
   endif

   var = imagevar(i2,j2)
   if(var > 0.)then
      a1= an*bn
      ab(2,2) = a1*time/var
      carte(i,j) =  carte(i,j)+a1*image(i2,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)**2/var
      normvar(i,j) = normvar(i,j)+ time*a1/var
   endif
   
   var = imagevar(in,j2)
   if(var > 0.)then
      a1=(1-an)*bn
      ab(1,2) = a1*time/var
      carte(i,j) =  carte(i,j)+a1*image(in,j2)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)**2/var
      normvar(i,j) = normvar(i,j)+ time* a1/var
   endif

   var = imagevar(i2,jn)
   if(var > 0.)then
      a1=an*(1-bn)
      ab(2,1) = a1*time/var
      carte(i,j) =  carte(i,j)+a1*image(i2,jn)/var
      cartevar(i,j) = cartevar(i,j)+ (time*a1)**2/var
      normvar(i,j) = normvar(i,j)+  time*a1/var
   endif

  
   cartevar(i,j) = cartevar(i,j) + 2.*detmeanval*&
        (CovarTab(in,jn)*ab(1,1)*ab(1,2) + CovarTab(in,jn)*ab(1,1)*ab(2,1)+&
        CovarTab(i2,jn)*ab(2,1)*ab(2,2)+CovarTab(in,j2)*ab(1,2)*ab(2,2)+&
        Covar1Tab(in,jn)*ab(1,1)*ab(2,2)+ Covar1Tab(i2,jn)*ab(2,1)*ab(1,2))/time
  
 ! exposure
  exposure(i,j) =  exposure(i,j)+ time
 


endif

!......................
END SUBROUTINE add_val4_withtime
!......................

!****************************
END MODULE ATTI_INTERNAL
!******************************



!***********************
MODULE ATTI_DECLARATIONS
!***********************
!SUBROUTINE DECLARATIONS
INTERFACE


SUBROUTINE ATTI_INIT(n,ax,dx,az,dz,posangle,status)
!------------------------------------


!SHOULD BE CALLED FIRST BEFORE USE OF ATTI_LIB SUBROUTINES
! FOR THE INTERFACE_TYPE==0
!  DATA FILE atti.data REQUIRED
! DATA FILE FORMAT IS THE FOLLOWING :
! 1st line : number of pointings
! next for each pointing :
! data_type
! next info about the telescope axis
! if data_type = 1 then  next_info is:
!                         alpha,delta,theta in deg
! if data_type = 2 then  next_info is:
!                         X(alpha,delta),Z(alpha,delta) in deg
! if data_type = 3 then  next_info is:
!                         X(x,y,z),Z(x,y,z) in deg
! if data_type = 3 then  next_info is:
!                         X(lII,bII),Z(lII,bII) 
!FOR THE INTERFACE_TYPE ==1 THE LIST OF ATTITUDES IN FORMAT 2
! IS GIVEN AS THE CALLING LIST
Use ATTI_DEFS
Implicit None
!INPUT VARIABLES
INTEGER ::n  ! NUMBER OF POINTINGS
REAL(KIND=8),DIMENSION(:),POINTER,OPTIONAL::ax,dx,az,dz,posangle
INTEGER :: Status
! RA,DEC OF TELESCOPE AXIS
End Subroutine ATTI_INIT


Subroutine ATTI_CONV_ADYZ(alpha_deg,delta_deg,y,z,num,proj,mf,status)
!-------------------------------------------------------------------
Use ATTI_DEFS
Implicit None
! INPUT VARIABLES :
Real(kind=4),Intent(in) :: alpha_deg,delta_deg   ! RA DEC COORD. IN deg
Integer     ,Intent(in) :: proj                  !projection type
Real(kind=4),Intent(in) :: mf
Integer     ,Intent(in) :: num                   ! POINTING NUMBER
!  OUTPUT VARIABLES
Real(kind=4),Intent(out) :: y,z                   ! TANGENT PLANE COORD.
Integer     ,Intent(out) :: status  
End Subroutine ATTI_CONV_ADYZ

Subroutine ATTI_CONV_YZAD(y,z,alpha_deg,delta_deg,num,proj,mf,status)
!-------------------------------------------------------------------
Use ATTI_DEFS
Implicit None
! INPUT VARIABLES
Real(kind=4),Intent(in) :: y,z      ! COORD. IN TANGENT PLANE
Integer     ,Intent(in) :: proj     ! projection type  
Real(kind=4),Intent(in) :: mf       !pixel_size/integral_pixel_size
Integer     ,Intent(in) :: num      ! POINTING NUMBER
!OUTPUT VARIABLES
Real(kind=4),Intent(out) :: alpha_deg,delta_deg ! RA DEC (deg) 
Integer     ,Intent(out) :: status
END SUBROUTINE ATTI_CONV_YZAD

SUBROUTINE ROTAT &
             (mf,is,tc,rot_Type,&
              idim1,jdim1,idim2,jdim2,ci1,cj1,ci2,cj2,&
             num1,num2, image,ScwExposition,imagevar,carte,cartevar,normvar,&
             exposure,surf,time,info,detmeanval,status)

!----------------------------------------------------------
USE ATTI_DEFS
USE ATTI_INTERNAL
 Implicit None

!INPUT VARIABLES

Real(kind=4),Intent(IN)::mf
! magnification factor carte_pixel_size/image_pixel_size
Integer   ::is
! 0 image-to-sky
! 1 inverse
Integer,Intent(IN)   :: rot_Type

! rot_type  :THE TYPE OF ROTATION
!
! rot_type       = 1 THEN :FLUX SPREAD (Willmore)
!                = 2 THEN :FLUX SPREAD (Approx)  
!                = 0 THEN :NO FLUX SPREAD
!                = 3 THEN :NO FLUX SPREAD : special peak conservation  

Integer,Intent(in) :: tc
! TYPE OF PROJECTION
! 1 - TAN-TAN
! 2 - TAN-CAR 





 Integer ,Intent(IN) :: num1,num2  ! IMAGE AND CARTE  POINTING NUMBER
 Integer ,Intent(IN) :: idim1,jdim1  ! IMAGE DIMENSIONS
 Real(KIND=4)  ,Intent(IN) :: ci1,cj1  ! IMAGE CENTRE (TANGENT POINT)
 Real(kind=4)  ,dimension(:,:),pointer :: image   ! IMAGE ARRAY
 Real(kind=4)  ,dimension(:,:),pointer :: ScwExposition ! scw exposure map
 Real(kind=4)  ,dimension(:,:),pointer,Optional :: imagevar  ! IMAGE VARIANCE ARRAY
  Real(kind=4)  ,dimension(:,:),pointer,Optional :: normvar
 Integer       ,Intent(IN) :: idim2,jdim2    ! SKY CARTE DIMENSIONS
 Real(KIND=4)  ,Intent(IN) :: ci2,cj2 ! SKY CARTE CENTRE (TANGENT POINT)  
Real(KIND=4) :: time,detmeanval
INTEGER                    :: info ! print level
! OUTPUT VARIABLES

Integer        :: STATUS
! status =   -1 if invalid pointing number or 
!                  invalid rotation type or 
!                  invalid weigt flag or
!                  invalid configuration of input/output files
! status =   1,2 if unusual situation occured in CONV_
!           ( in such a case a pixel contribution ignored)
!        =  0 othervise               

! INPUT/OUTPUT VARIABLES
Real(kind=4) ,dimension(:,:),pointer  :: carte   ! SKY ARRAY
Real(kind=4),dimension(:,:),pointer ,Optional :: cartevar,exposure ,surf  ! SKY VARIANCE ARRAY
END SUBROUTINE ROTAT


END INTERFACE
!***************************
END MODULE ATTI_DECLARATIONS
!***************************

!***************************
MODULE ATTI_ROTAT_MODULE
!***************************
INTERFACE
SUBROUTINE ROTAT_LOOP_ONE(mf,is,rot_Type,tc,&
            num1,num2,info,image_type,imin,imax,jmin,jmax,outside,&
            idim2,jdim2,i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,id,&
            image, ScwExposition,imagevar,normvar,&
            carte,cartevar ,exposure,surf,aref,dref,lref,bref,&
            ci1,cj1,ci2,cj2,detmeanval,time,&
            cb,ce,i_cen_inp,j_cen_inp,i_cen_out,j_cen_out,&
            pixel1,pixel2,pixel_ang1,pixel_ang2,&
            dispi_in,dispj_in, dispi_out,dispj_out,&
            rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4,status)
!------------------------------------------------------------------
 Implicit None
!INPUT/OUTPUT VARIABLES
Real(kind=4)::mf
Integer   ::is
Integer   :: rot_Type
Integer,Intent(in) :: tc
Integer  :: num1,num2,info,image_type,imin,imax,jmin,jmax,outside
INTEGER  ::idim2,jdim2,i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,id
Integer        :: STATUS
Real(kind=4)  ,dimension(:,:),pointer :: image, ScwExposition 
Real(kind=4)  ,dimension(:,:),pointer,Optional :: imagevar,normvar
Real(kind=4),dimension(:,:),pointer  :: carte  
Real(kind=4),dimension(:,:),pointer ,Optional :: cartevar ,exposure,surf
Real(kind=4) ::ci1,cj1,ci2,cj2,detmeanval,time
Real(KIND=8) ::aref,dref,lref,bref
Real(KIND=8) ::cb,ce,i_cen_inp,j_cen_inp,i_cen_out,j_cen_out
Real(KIND=8) :: pixel1,pixel2,pixel_ang1,pixel_ang2
Real(KIND=8) :: dispi_in,dispj_in, dispi_out,dispj_out
Real(kind=8),Dimension(3,3)::rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4 

!------------------
END SUBROUTINE ROTAT_LOOP_ONE
!------------------

END INTERFACE


!***************************
END MODULE ATTI_ROTAT_MODULE
!***************************

!################################
! CODE OF EXTERNAL SUBROUTINES 
!| ##############################
!| EXTERNAL SUBROUTINES 
!| ##############################


!===========================================================
Subroutine ATTI_INIT(n,ax,dx,az,dz,posangle,Status)
!===========================================================
!| ---------------------
!| Subroutine ATTI_INIT
!| ---------------------
!| SHOULD BE CALLED  BEFORE FIRST USE OF any ATTI_LIB SUBROUTINES
!|   FOR THE INTERFACE_TYPE==0
!|      DATA FILE atti.data REQUIRED
!|  DATA FILE FORMAT IS THE FOLLOWING :
!|  1st line : number of pointings
!|  next for each pointing :
!|  data_type
!|  next info about the telescope axis
!|  if data_type = 1 then  next_info is:
!|                          alpha,delta,theta in deg
!|  if data_type = 2 then  next_info is:
!|                          X(alpha,delta),Z(alpha,delta) in deg
!|  if data_type = 3 then  next_info is:
!|                          X(x,y,z),Z(x,y,z) in deg
!|  if data_type = 3 then  next_info is:
!|                          X(lII,bII),Z(lII,bII) 
!| FOR THE INTERFACE_TYPE ==1 THE LIST OF ATTITUDES IN FORMAT 2
!|  IS GIVEN AS THE CALLING LIST
Use ATTI_DEFS
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE 
USE MIMOSA_USE_MODULE
Implicit None
!INPUT VARIABLES
INTEGER ::n  ! NUMBER OF POINTINGS
REAL(KIND=8),DIMENSION(:),POINTER,OPTIONAL::ax,dx,az,dz,posangle
! RA,DEC OF TELESCOPE AXIS
INTEGER :: Status
!AUXILIARY VARIABLES
Integer poi,iok
character(len=20) :: procname

procname = 'ATTI_INIT'
Status = ISDC_OK

! SPR 4356 : err initialisation
err = 0

Call INIT

If(message_level==1) Print *,' SLA_INIT interface_type : ',Interface_Type

If (Interface_Type==0) Then

! LABO INPUT FILE
!-----------------
   If(message_level==1) Print *,' Labo input file required'

   Open(3,file='atti.data',status='old')

   Read (3,*) point_num

   N = POINT_NUM
   allocate(point_rot_mat(point_num,3,3),&
        pointing_table(point_num,3),stat=iok)
   if(iok.ne.0)then
      err=30
      return
   endif
   point_rot_mat(:,:,:) = 0.
   pointing_table(:,:) = 0.
!!$   If((point_num.Le.0).Or.(point_num.Gt.point_max)) Then
!!$      err = 1
!!$      If(message_level==1) Print *,' INVALID POINTING NUMBER ',point_num
!!$      Return 
!!$   Endif

   Do poi=1,point_num 
 
       Call Read_ATTI_Data
       If (err.Ne.0) Return
   
       Call LABO_TRANSFORM
       If (err.Ne.0) Return

       Call CREATE_ROTATION_MATRIX

       point_rot_mat(poi,1:3,1:3) = rot_mat
       pointing_table(poi,1) = alphaXtel_deg
       pointing_table(poi,2) = deltaXtel_deg
       pointing_table(poi,3) = theta_deg

   Enddo

   Close(3)
Else

! ISDC INTERFACE
!-----------------

   If(message_level==1) Print *,'  ISDC FORMAT '
   point_num=N

   allocate(point_rot_mat(point_num,3,3),&
           pointing_table(point_num,3),stat=iok)
  if(iok.ne.0)then
     err=30
     return
  endif
  point_rot_mat(:,:,:) = 0.
  pointing_table(:,:) = 0.


!!$   If((point_num.Le.0).Or.(point_num.Gt.point_max)) Then
!!$      err = 1
!!$      If(message_level==1) Print *,' INVALID POINTING NUMBER ',point_num
!!$      Return 
!!$   Endif



  data_type=2

   Do poi=1,point_num 
         If(message_level==1)Then
         Print *,' GIVEN X(RA,DEC), Z(RA,DEC)'
      Endif 

      alphaXtel_deg = ax(poi)
      deltaXtel_deg = dx(poi)
      alphaZtel_deg = az(poi)
      deltaZtel_deg = dz(poi)
      
      alphaXtel_deg = mod(alphaXtel_deg+360.0d0,+360.0d0)
      alphaZtel_deg = mod(alphaZtel_deg+360.0d0,+360.0d0)

     
      If ((deltaXtel_deg.Lt.-90.d0).or.(deltaXtel_deg.Gt.90.d0))&
          err = Data_Type
      If ((deltaZtel_deg.Lt.-90.d0).or.(deltaZtel_deg.Gt.90.d0))&
          err = Data_Type 
      If(err.Ne.0) Then
         write(str250,*)'  Attitude error 1 in Scw : ',poi
         call message(procname,str250,AttiError,Status)
           err = err+10
           If(message_level==1)&
               Print *,' ERROR IN DATA  ',err
           return
      Endif

      theta_deg = posangle(poi)
      !map theta angle should always be 0
      if (poi==point_num)data_type=1

      Call LABO_TRANSFORM
      If (err.Ne.0) then
          write(str250,*)'  Attitude error 2 in Scw : ',poi
         call message(procname,str250,AttiError,Status)
         Return
      endif

       Call CREATE_ROTATION_MATRIX

       point_rot_mat(poi,1:3,1:3) = rot_mat
       pointing_table(poi,1) = alphaXtel_deg
       pointing_table(poi,2) = deltaXtel_deg
       pointing_table(poi,3) = theta_deg
       
   Enddo
Endif


!==========================
End Subroutine ATTI_INIT
!==========================



!| --------------------------
!| Subroutine ATTI_CONV_ADYZ
!| --------------------------
!| TAN projection:
!| Conversion from (alpha,delta) in deg (sphere)to the tangent plane
!| coordinates  (y,z) in GNOMONIC (TAN) or CARTESIAN (CAR) projection 
!| INPUT VARIABLES :
!|   alpha_deg,delta_deg -  spherical coordinates in deg. 
!|   num - pointing number
!|   proj - projection type
!|   mf - magnification factor : ratio input_image_pixel_size/output_image_pixel_size
!| OUTPUT VARIABLES :  
!|   y,z -  TANGENT PLANE COORD.
!|   status  :             
!|   -1 when invalid pointing number
!|    1  when  x< 0 
!|        The arrival direction is below the detector plane.
!|        In this case the direction cosines are conserved but 
!|        the resulting point (y,z) is placed outside the radius 
!|        yz_radius.To recover the original position :
!|              n= dsqrt(y**2+z**2)
!|              y = y*yz_radius/n
!|              z = z*yz_radius/n
!|    2 when x near 0
!|       This means a very great direction angle with respect to 
!|       the line-of-sight-axis X
!|    0  otherwise
!| 
!============================================================
Subroutine ATTI_CONV_ADYZ(alpha_deg,delta_deg,y,z,num,proj,mf,status)
!============================================================
! CONVERSION FROM (RA,DEC) in deg IN THE EQUATORIAL SYSTEM
! TO THE TANGENT PLANE COORDINATES (y,z) IN GNOMONIC (TAN) 
! or CARTESIAN (CAR)PROJECTION     
Use ATTI_DEFS
Implicit None

! INPUT VARIABLES :

Real(kind=4),Intent(in) :: alpha_deg,delta_deg   ! RA DEC COORD. IN deg
Integer     ,Intent(in) :: proj                  !projection type
! 1 - TAN
! 2 - CAR
Real(kind=4),Intent(in) :: mf
Integer     ,Intent(in) :: num                   ! POINTING NUMBER

!  OUTPUT VARIABLES
Real(kind=4),Intent(out) :: y,z                   ! TANGENT PLANE COORD.

Integer     ,Intent(out) :: status  
! output status              
!   -1 when invalid pointing number or projection type
!   1  when  x< 0 
!        The arrival direction is below the detector plane.
!        In this case the direction cosines are conserved but 
!        the resulting point (y,z) is placed outside the radius 
!        yz_radius.To recover the original position :
!              n= dsqrt(y**2+z**2)
!              y = y*yz_radius/n
!              z = z*yz_radius/n
!   2 when x near 0
!       This means a very great direction angle with respect to 
!       the line-of-sight-axis X
!   0  otherwise


! AUXILIARY VARIABLES

Real(kind=8) :: a,d,a_rad,d_rad,n1,pixel_ang,ar
Real(kind=4) :: adiff


If((num.Lt.1).Or.(num.Gt.point_num)) Then
    STATUS = -1
    y = yz_radius
    z = yz_radius
    Return
Endif

If((proj.ne.0).and.(proj.ne.1))then
     STATUS = -1
     Return
Endif

STATUS = 0

a = alpha_deg
d = delta_deg

If((a.Lt.0.0d0).or.( a .Gt. 360.0)) Then
   If(message_level ==1) Print *,' ATTENTION RA not in (0,360) !!'
   a = mod(a,360.0d0)
Endif
If(Abs(d).Gt.90.0) Then
   If(message_level ==1) Print *,' ATTENTION abs(DEC) > 90 !!'
   d = mod(d,90.0d0)
Endif

if(proj==1)then ! TAN PROJECTION 
   rot_mat = point_rot_mat(num,1:3,1:3)  
   a_rad = a*deg_to_rad
   d_rad = d*deg_to_rad

   ! transformation (a,d) --> (x,y,z) in Equatorial System
   Call sla_dcs2c(a_rad,d_rad,vector)
   ! rotation (x,y,z) Equat. --> (x,y,z) Telescope
   Call sla_dmxv(rot_mat,vector,vector1)
   ! NORMALISATION FACTOR
   ! transformation to pixel units
   n1 = vector1(1)
   If(n1.Le.0.d0)Then
      ! ARRIVAL DIRECTION BELOW THE DETECTOR PLANE
      STATUS = 1
      n1 =  (1.d0+10.0d-10 -n1**2) / yz_radius
   Else
      If(vector1(1).Le.x_tol) Then
         ! TOO BIG ARRIVAL ANGLE
         STATUS = 2
      Endif
   Endif

   n1=n1*pixel_ratio*mf
   vector1 = vector1/n1
   !  reflection of Y axis
   ! RETOUR TO THE SINGLE PRECISION
   y =IOView* vector1(2)
   z = vector1(3)
else
   ! CAR PROJECTION
   adiff = abs(alpha_deg-pointing_table(num,1))
   if(adiff.lt.180.0)then
      y = IOView*(alpha_deg-pointing_table(num,1))
   else
      a = alpha_deg
      ar= pointing_table(num,1)
      if(alpha_deg.lt.180.0)then
         a =a+180.0d0
         else
         ar = ar+180.0d0
      endif
      y =IOView*( a-ar)
   endif
   pixel_ang=rad_to_deg*atan(pixel_ratio*mf)
   y = y/pixel_ang
   z =(delta_deg-pointing_table(num,2))/(pixel_ang)
endif
!=============================================
End Subroutine ATTI_CONV_ADYZ
!=============================================



!| --------------------------
!| Subroutine ATTI_CONV_YZAD
!| --------------------------
!|  TAN back-projection : 
!|  Conversion from the tangent plane corrdinates (y,z)
!|  IN GNOMONIC (TAN) PROJECTION  TO (alpha,delta) in deg
!|  INPUT VARIABLES : 
!|     y,z - tangent plane coordinates (referenced to its tangent point)
!|     proj - projection type
!|  OUTPUT VARIABLES
!|    alpha_deg,delta_deg -  spherical coordinates in deg. 
!|    num - pointing number
!|    mf - magnification factor : ratio input_image_pixel_size/output_image_pixel_size
!|    status
!|       -1 if invalid pointing number
!|        0  otherwise
!| 
!=========================================================
Subroutine ATTI_CONV_YZAD(y,z,alpha_deg,delta_deg,num,proj,mf,status)
!==========================================================
! CONVERSION FROM THE TANGENT PLANE COORDINATES (y,z)
! IN GNOMONIC (TAN) PROJECTION  TO (RA,DEC) in deg IN THE EQUATORIAL SYSTEM


Use ATTI_DEFS
Implicit None

! INPUT VARIABLES
Real(kind=4),Intent(in) :: y,z      ! COORD. IN TANGENT PLANE
Integer     ,Intent(in) :: proj     ! projection type  
Real(kind=4),Intent(in) :: mf       !pixel_size/integral_pixel_size
Integer     ,Intent(in) :: num      ! POINTING NUMBER

!OUTPUT VARIABLES
Real(kind=4),Intent(out) :: alpha_deg,delta_deg ! RA DEC (deg) 
Integer     ,Intent(out) :: status
! status = -1 if invalid pointing number or projection type
!        = 0  otherwise

! AUXILIARY VARIABLES
Real(kind=8) :: a_rad,d_rad,pixel_ang

If((num.Le.0).Or.(num.Gt.point_num)) Then
  STATUS = -1
  ALPHA_DEG = 0.0
  DELTA_DEG = 0.0
  Return
Endif
If((proj.ne.0).and.(proj.ne.1).and.(proj.ne.2))then
     STATUS = -1
     Return
Endif

STATUS = 0

if(proj==1) then  !TAN projection
   rot_mat = point_rot_mat(num,1:3,1:3)
   ! TRANSFORMATION TO A VECTOR OF THE DIRECTION IN THE SATELLITE SYSTEM
   vector(1) = 1
   vector(2) =  IOView*(y*pixel_ratio*mf)
   vector(3) = (z*pixel_ratio*mf)
   ! BACK PROJECTION TO EQUATORIAL SYSTEM
   Call sla_dimxv(rot_mat,vector,vector1)
   ! (a,d) --> (x,y,z)
   Call sla_dcc2s(vector1,a_rad,d_rad)
   ! normalisation
   a_rad = sla_dranrm(a_rad)
   d_rad = sla_drange(d_rad)
   alpha_deg = a_rad*rad_to_deg
   delta_deg = d_rad*rad_to_deg
   alpha_deg = MOD(alpha_deg+360.d0,360.0d0)
else
   pixel_ang=rad_to_deg*atan(pixel_ratio*mf)
   alpha_deg = Mod(IOView*y*pixel_ang+pointing_table(num,1)+360.0d0,360.0d0)
   delta_deg = z*pixel_ang+pointing_table(num,2)
endif

!=====================================================
End Subroutine ATTI_CONV_YZAD
!=====================================================
!| -----------------
!| SUBROUTINE ROTAT 
!| -----------------
!|  TAN_TAN/TAN_CAR PROJECTION OF AN input IMAGE/VARIANCE IMAGE
!|  onto an output SKY/SKY VARIANCE MAP 
!|  
!|  INPUT VARIABLES
!|     mf - magnification factor: carte_pixel_size/image_pixel_size
!|     is 
!|        0 image-to-sky projection
!|        1 inverse      projection 
!|     rot_Type - THE TYPE OF ROTATION
!|           1  :FLUX SPREAD (Willmore)
!|           2  :FLUX SPREAD (Approx)  
!|           0  :NO FLUX SPREAD
!|           3  :NO FLUX SPREAD : special peak conservation
!|     tc - TYPE OF PROJECTION
!|          1 - TAN-TAN projection
!|          2 - TAN-CAR projection
!|    num1,num2             - IMAGE AND CARTE  POINTING NUMBER
!|    idim1,jdim1           - IMAGE DIMENSIONS
!|    ci1,cj1               - IMAGE CENTRE (TANGENT POINT)
!|    image(idim1,jdim1)    - IMAGE ARRAY
!|    imagevar(idim1,jdim1) - IMAGE VARIANCE ARRAY
!|    iidim2,jdim2          -  SKY CARTE DIMENSIONS
!|    ci2,cj2               -  SKY CARTE CENTRE (TANGENT POINT)  
!|    info                  - prints level
!| 
!|  INPUT/OUTPUT VARIABLES
!|    carte(idim2,jdim2)    - SKY ARRAY
!|    cartevar(idim2,jdim2) - SKY VARIANCE ARRAY
!| 
!| OUTPUT VARIABLES
!|    STATUS :
!|  -1 if invalid pointing number or 
!|     invalid rotation type or 
!|      invalid weigt flag or
!|      invalid configuration of input/output files
!|  1,2 if unusual situation occured in CONV_
!|      ( in such a case a pixel contribution ignored)
!|  0 othervise               

!=========================================================
SUBROUTINE ROTAT &
             (mf,is,tc,rot_Type,&
              idim1,jdim1,idim2,jdim2,ci1,cj1,ci2,cj2,&
             num1,num2, image,ScwExposition,imagevar,carte,cartevar,normvar,&
             exposure,surf,time,info,detmeanval,status)
!=========================================================
USE ATTI_DEFS
USE ATTI_INTERNAL
USE ATTI_ROTAT_MODULE
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
 Implicit None

!INPUT VARIABLES

Real(kind=4),Intent(IN)::mf
! magnification factor carte_pixel_size/image_pixel_size
Integer   ::is
! 0 image-to-sky
! 1 inverse
Integer,Intent(IN)   :: rot_Type

! rot_type  :THE TYPE OF ROTATION
!
! rot_type       = 1 THEN :FLUX SPREAD (Willmore)
!                = 2 THEN :FLUX SPREAD (Approx)  
!                = 0 THEN :NO FLUX SPREAD
!                = 3 THEN :NO FLUX SPREAD - peak conservation


Integer,Intent(in) :: tc
! TYPE OF PROJECTION
! 1 - TAN-TAN
! 2 - TAN-CAR 

 Integer ,Intent(IN) :: num1,num2  ! IMAGE AND CARTE  POINTING NUMBER
 Integer ,Intent(IN) :: idim1,jdim1  ! IMAGE DIMENSIONS
 Real(KIND=4)  ,Intent(IN) :: ci1,cj1  ! IMAGE CENTRE (TANGENT POINT)
 Real(kind=4)  ,dimension(:,:),pointer :: image   ! IMAGE ARRAY
 Real(kind=4)  ,dimension(:,:),pointer :: ScwExposition ! scw exposure array
 Real(kind=4)  ,dimension(:,:),pointer,Optional :: imagevar   ! IMAGE VARIANCE ARRAY
  Real(kind=4)  ,dimension(:,:),pointer,Optional :: normvar
 Integer       ,Intent(IN) :: idim2,jdim2    ! SKY CARTE DIMENSIONS
 Real(KIND=4)  ,Intent(IN) :: ci2,cj2 ! SKY CARTE CENTRE (TANGENT POINT)  
Real(kind=4) :: detmeanval,time
INTEGER                   ::info ! print level
! OUTPUT VARIABLES

Integer        :: STATUS
! status =   -1 if invalid pointing number or 
!                  invalid rotation type or 
!                  invalid weigt flag or
!                  invalid configuration of input/output files
! status =   1,2 if unusual situation occured in CONV_
!           ( in such a case a pixel contribution ignored)
!        =  0 othervise               

! INPUT/OUTPUT VARIABLES
Real(kind=4),dimension(:,:),pointer  :: carte   ! SKY ARRAY
Real(kind=4),dimension(:,:),pointer ,Optional :: cartevar   ! SKY VARIANCE ARRAY
Real(kind=4),dimension(:,:),pointer ,Optional :: exposure,surf
! AUXILIARY VARIABLES

INTEGER::i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,OUTSIDE,i,j,k,itab(4),sitab(4),ntab(1)
INTEGER::id,expo,i1,j1,i2,j2,ia,ib,ja,jb,imin,imax,jmin,jmax
INTEGER::in,jn,stat,image_type,iloc(2),loci,locj
INTEGER::minin,maxin,minjn,maxjn,innermosa
LOGICAL :: surfstat,cond
Real(KIND=8) ::cb,ce,ri,rj,aref,dref,lref,bref,xn,yn,xn1,yn1,an,bn
Real(KIND=8) ::xnn,ynn,xnnp,ynnp,rip,rjp,surface
Real(KIND=8) ::kx,ky,i_cen_inp,j_cen_inp,i_cen_out,j_cen_out
Real(KIND=8) :: pixel1,pixel2,pixel_ang1,pixel_ang2
Real(KIND=8) :: dispi_in,dispj_in, dispi_out,dispj_out,mosaangsize
Real(kind=8),Dimension(3,3)::rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4 
Real(KIND=4)  ,dimension(idim2,jdim2) :: filltab
Real(KIND=4)  ,dimension(4,4) :: petitcarte,petitima,petitvar

!!$petitima = 0.
!!$petitcarte = 0
!!$petitvar = 0
!INPUT CONFIGURATION VERIFICATION

if(present(Surf))then
   surfStat = .true.
   surf = 0.
else
   surfStat = .false.
endif

i= size(image,1);j=size(image,2)
if((i.ne.idim1).or.(j.ne.jdim1))then
    status = -1
    If(message_level==1)Then
       Print *,' incorrect dimensions of input image :',&
       i,j,idim1,jdim1
 Endif
 Return
Endif

i= size(carte,1);j=size(carte,2)
if((i.ne.idim2).or.(j.ne.jdim2))then
    status = -1
    If(message_level==1)Then
       Print *,' incorrect dimensions of output carte : ',&
       i,j,idim2,jdim2    
 Endif
 Return
Endif

If((Present(imagevar)).And.(.Not.Present(cartevar)))Then
 status = -1
 If(message_level==1)Then
   Print *,' INP/OUT VARIANCE IMAGE MISSING'
 Endif
 Return
Endif

If(Present(imagevar))then
   if(.not.present(normvar))then
      if(rot_type.ne.3)then
         status = -1
         If(message_level==1)Then
            Print *,' INP/OUT Norm. VARIANCE IMAGE MISSING'
         Endif
         Return
      Endif
   endif

  i= size(imagevar,1);j=size(imagevar,2)
  if((i.ne.idim1).or.(j.ne.jdim1))then
    status = -1
    If(message_level==1)Then
       Print *,' incorrect dimensions of input image var :',&
        i,j,idim1,jdim1
    Endif
   Return
  Endif
   i= size(cartevar,1);j=size(cartevar,2)
  if((i.ne.idim2).or.(j.ne.jdim2))then
    status = -1
    If(message_level==1)Then
       Print *,' incorrect dimensions of output carte var ',&
         i,j,idim2,jdim2    
    Endif
   Return
  Endif
endif

!PROJECTION DIRECTION VERIFICATION
if((is.lt.0).or.(is.gt.1))then
 status = -1
 If(message_level==1)Then
   Print *,' BADPROJECTION  DIRECTION'
  Endif
 Return
Endif

!POINTING VERIFICATION
If((num1.Le.0).Or.(num1.Gt.point_num).Or.(num2.Le.0).Or.(num2.Gt.point_num))Then
 status = -1
 If(message_level==1)Then
   Print *,' BAD POINTING NUMBER'
 Endif
 Return
Endif

! ROTATION TYPE VERIFICATION
If((rot_Type.Lt.0).Or.(rot_Type.Gt.20))Then
  status = -1
   If(message_level==1)Then
   Print *,' BAD ROTATION TYPE'
 Endif
  Return
Endif

cond = (rot_type.ne.0).and.(rot_type.ne.3)
if(cond.and.&
     ((.not.Present(imagevar)).or.(.Not.Present(cartevar)).or.&
     (.Not.Present(exposure))  ))then
   status = -1
   If(message_level==1)Then
      Print *,' Willmore rebinning requires variances'
   Endif
   Return
Endif

! PROJECTION TYPE tc  VERIFICATION
If((tc.Lt.0).Or.(tc.Gt.2))Then
  status = -1
   If(message_level==1)Then
   Print *,' BAD PROJECTION TYPE'
 Endif
  Return
Endif
! MAGNIFICATION FACTOR VERIF
If(mf.Le.0.0)Then
 status = -1
   If(message_level==1)Then
   Print *,' CARTE/IMAGE RATIO INVALID'
 Endif
  Return
Endif

!If(mf.Lt.1.0)Then
! status = -1
!   If(message_level==1)Then
  
!   Print *,' CARTE/IMAGE RATIO uncorrect'
! Endif
!  Return
!Endif

if(info.gt.1)then
  Print *,' ROTATION/SUMMATION PROCEDURE'
  Print *,' CONFIGURATION:'
  Print *,'                projection   :'
  IF(IS.EQ.0) THEN
  PRINT *,'                             IMAGE-TO-SKY'
  ELSE
  PRINT *,'                             SKY-TO-IMAGE'
  ENDIF
  If(tc.Eq.1)Then
  Print *,'                             TAN - TAN'
  Else
  Print *,'                             TAN - CAR'
  Endif
 
  Print *,'                flux spread  :'
 Select Case(rot_Type)
  Case(0)
  Print *,'                             NO'
  Case(3)
  Print *,'                             SPECIAL - peak conservation'
  Case(1)
  Print *,'                             Willmore'
   Case(2)
  Print *,'                            Approx'
 

  End Select
  Print *,'                input        :'
  print *,'                              image'
  If(Present(imagevar)) Then
  Print *,'                              imagevar'
  Endif
  Print *,'                output        :'
  
  Print *,'                              carte'

  If(Present(cartevar)) Then
  Print *,'                              cartevar'
  Endif
  Print *,' CARTE_PIXEL_SIZE/IMAGE_PIXEL_SIZE : ',mf
endif

status = 0
outside = 0

if(.not.present(exposure))then
  expo = 0
else
  expo = 1
endif

if(.not.present(imagevar))then
  image_type = 0
else
  image_type = 1
endif

filltab = 0.

MapTanPointDisX =0.
MapTanPointDisY =0.

if(mosatype(num1)==1)then
   ImaTanPointDisX =0.
   ImaTanPointDisY=0.
else
   ImaTanPointDisX =0.5
   ImaTanPointDisY =0.5
endif


!!$! Input map shift
!!$if((idim1/2)*2==idim1)then
!!$   ImaTanPointDisX =0.5
!!$else
!!$   ImaTanPointDisX =0.
!!$endif
!!$if((jdim1/2)*2==jdim1)then
!!$   ImaTanPointDisY =0.5
!!$else
!!$   ImaTanPointDisY =0.
!!$endif
!!$!Output map shift
!!$if((idim2/2)*2==idim2)then
!!$   MapTanPointDisX =0.5
!!$else
!!$   MapTanPointDisX =0.
!!$endif
!!$if((jdim2/2)*2==jdim2)then
!!$   MapTanPointDisY =0.5
!!$else
!!$   MapTanPointDisY =0.
!!$endif




pixel1 = pixel_ratio
pixel2 = pixel_ratio*mf

pixel_ang1 = pixel_ang_size
pixel_ang2 = rad_to_deg*atan(pixel2)

rot_mat1 = point_rot_mat(num1,1:3,1:3)
rot_mat2 = point_rot_mat(num2,1:3,1:3)
rt = Transpose(rot_mat1)
rt2 = Transpose(rot_mat2)
Call sla_dmxm(rot_mat1,rt2,rot_mat4)
Call sla_dmxm(rot_mat2,rt,rot_mat3)
aref = POINTING_TABLE(num2,1)
dref = pointing_table(num2,2)
lref = LX(num2)
bref = BX(num2)

select case(rot_type)
CASE(0,3,5)
!NO FLUX SPREAD
  id= 1
  ce = 0.0
  cb = 0.0
 
case(1,4,6,7)
!FLUX SPREAD - SURFACE
  id = 0
  ce = 0.0
  cb = 0.0
 
case(2)
!FLUX SPREAD - APPROX
  id = 0
  ce = 0.0
  cb = 0.5
 

END SELECT

!IMAGE-TO-SKY-PROJECTION

! for transformation from pixel to telescope frame of the input image
dispi_in = ImaTanPointDisX 
dispj_in = ImaTanPointDisY
! for transformation from telescope to pixel frame of the output carte
dispi_out = MapTanPointDisX 
dispj_out = MapTanPointDisY

i_dim_inp = idim1
j_dim_inp = jdim1
i_dim_out = idim2
j_dim_out = jdim2
i_cen_inp = ci1
j_cen_inp = cj1
i_cen_out = ci2
j_cen_out = cj2
imin = 1
imax = i_dim_inp
jmin = 1
jmax = j_dim_inp

k=0
if(is.eq.1)then
   !SKY-TO-IMAGE PROJECTION
   imin = idim2
   imax = 1
   jmin = jdim2
   jmax = 1
   is = 0 ! temporary set to define projection limits
          ! will be reset just after
   ! calculate start and end of projected pixel
  
   do i1=1,4
!!$      select case (i1)
!!$      case(1)
!!$         i = RealMapSize(num1,1,1)  ; j=RealMapSize(num1,2,1)
!!$      case(2)
!!$         i = RealMapSize(num1,1,1)  ; j=RealMapSize(num1,2,2)
!!$      case(3)
!!$         i = RealMapSize(num1,1,2)  ; j=RealMapSize(num1,2,1)
!!$      case(4)
!!$         i = RealMapSize(num1,1,2)  ; j=RealMapSize(num1,2,2)
!!$      end select
       k=k+1
      select case (i1)
      case(1)
         i = 1  ; j=1
      case(2)
         i=1    ; j=jdim1
      case(3)
         i=idim1; j=1
      case(4)
         i=idim1; j=jdim1
      end select
    

      ri = Real(i)-cb-i_cen_inp - dispi_in 
      rj=Real(j)-cb-j_cen_inp   - dispj_in 

      call projection(is,tc,pixel1,pixel2,&
           pixel_ang2,aref,dref,lref,bref,&
           rot_mat1,rot_mat3,rot_mat4,&
           ri,rj,xn,yn,stat)
      if(stat==0)then
         xnn = xn-ce +i_cen_out- dispi_out
         ynn = yn-ce +j_cen_out- dispj_out
         in = Nint(xnn ) 
         jn = Nint(ynn )
         itab(k) = in
         if(imin > in)imin=in
         if(imax< in)imax = in
         if(jmin > jn)jmin=jn
         if(jmax < jn)jmax = jn

      endif
   enddo
 
   innermosa = 1
!!$   i = idim1/2
!!$   j = jdim1/2
!!$   ri = Real(i)-cb-i_cen_inp !+ dispi_in 
!!$   rj=Real(j)-cb-j_cen_inp   ! + dispj_in 
!!$   call projection(is,tc,pixel1,pixel2,&
!!$        pixel_ang2,aref,dref,lref,bref,&
!!$        rot_mat1,rot_mat3,rot_mat4,&
!!$        ri,rj,xn,yn,stat)
!!$   if(stat==0)then
!!$      xnn = xn-ce +i_cen_out- dispi_out
!!$      ynn = yn-ce +j_cen_out- dispj_out
!!$      in = Nint(xnn ) 
!!$      jn = Nint(ynn )
   mosaangsize = (imax-imin)*pixel_ang2
   if(mosaangsize > 180) then
         innermosa=0
         do i=1,4 
            ntab=minloc(itab)
            k=ntab(1)
            sitab(i)=itab(k)
            itab(k) = 20000
         enddo
      endif


   if ((imin==imax).or.(jmin==jmax))then
      imin = 1
      imax = idim2
      jmin = 1
      jmax = jdim2
   else
      if(imin <1)imin=1
      if(imax>idim2)imax = idim2
      if(jmin <1)jmin=1
      if(jmax>jdim2)jmax=jdim2

   endif
   ! for transformation from pixel to telescope frame of the input carte
   is = 1
   dispi_in = MapTanPointDisX 
   dispj_in = MapTanPointDisY
   ! for transformation from telescope to pixel frame of the output image
   dispi_out =  -ImaTanPointDisX 
   dispj_out =  - ImaTanPointDisY
   i_dim_inp = idim2
   j_dim_inp = jdim2
   i_dim_out = idim1
   j_dim_out = jdim1
   i_cen_inp = ci2
   j_cen_inp = cj2
   i_cen_out = ci1
   j_cen_out = cj1

endif

if(innermosa==1)then
   call ROTAT_LOOP_ONE (mf,is,rot_Type,tc,&
        num1,num2,info,image_type,imin,imax,jmin,jmax,outside,&
        idim2,jdim2,i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,id,&
        image, ScwExposition,imagevar,normvar,&
        carte,cartevar ,exposure,surf,aref,dref,lref,bref,&
        ci1,cj1,ci2,cj2,detmeanval,time,&
        cb,ce,i_cen_inp,j_cen_inp,i_cen_out,j_cen_out,&
        pixel1,pixel2,pixel_ang1,pixel_ang2,&
        dispi_in,dispj_in, dispi_out,dispj_out,&
        rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4,status)
else
!print *,'outer mosaicks'
   call ROTAT_LOOP_ONE (mf,is,rot_Type,tc,&
        num1,num2,info,image_type,1,sitab(2),jmin,jmax,outside,&
        idim2,jdim2,i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,id,&
        image, ScwExposition,imagevar,normvar,&
        carte,cartevar ,exposure,surf,aref,dref,lref,bref,&
        ci1,cj1,ci2,cj2,detmeanval,time,&
        cb,ce,i_cen_inp,j_cen_inp,i_cen_out,j_cen_out,&
        pixel1,pixel2,pixel_ang1,pixel_ang2,&
        dispi_in,dispj_in, dispi_out,dispj_out,&
        rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4,status)
    call ROTAT_LOOP_ONE (mf,is,rot_Type,tc,&
        num1,num2,info,image_type,sitab(3),idim2,jmin,jmax,outside,&
        idim2,jdim2,i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,id,&
        image, ScwExposition,imagevar,normvar,&
        carte,cartevar ,exposure,surf,aref,dref,lref,bref,&
        ci1,cj1,ci2,cj2,detmeanval,time,&
        cb,ce,i_cen_inp,j_cen_inp,i_cen_out,j_cen_out,&
        pixel1,pixel2,pixel_ang1,pixel_ang2,&
        dispi_in,dispj_in, dispi_out,dispj_out,&
        rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4,status)
endif

if(info.gt.1)then
   print *,'  rotat sum : ',sum(carte)
   if(outside .gt.0)then
      print *, 'Attention : ',outside ,' pixels outside a map'
   endif
endif

!print *,' minmax image zone treated :',minin,maxin,' ',minjn,maxjn
!============================
END SUBROUTINE ROTAT
!============================


!=========================================================
SUBROUTINE ROTAT_LOOP_ONE (mf,is,rot_Type,tc,&
            num1,num2,info,image_type,imin,imax,jmin,jmax,outside,&
            idim2,jdim2,i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,id,&
            image, ScwExposition,imagevar,normvar,&
            carte,cartevar ,exposure,surf,aref,dref,lref,bref,&
            ci1,cj1,ci2,cj2,detmeanval,time,&
            cb,ce,i_cen_inp,j_cen_inp,i_cen_out,j_cen_out,&
            pixel1,pixel2,pixel_ang1,pixel_ang2,&
            dispi_in,dispj_in, dispi_out,dispj_out,&
            rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4,status)
!=========================================================
USE ATTI_DEFS
USE ATTI_INTERNAL
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_GLOBVAR_MODULE
 Implicit None

!INPUT/OUTPUT VARIABLES

Real(kind=4)::mf
! magnification factor carte_pixel_size/image_pixel_size
Integer   ::is
! 0 image-to-sky
! 1 inverse
Integer   :: rot_Type
! rot_type       = 1 THEN :FLUX SPREAD (Willmore)
!                = 2 THEN :FLUX SPREAD (Approx)  
!                = 0 THEN :NO FLUX SPREAD
!                = 3 THEN :NO FLUX SPREAD - peak conservation
Integer,Intent(in) :: tc
! TYPE OF PROJECTION
! 1 - TAN-TAN
! 2 - TAN-CAR 

Integer  :: num1,num2,info,image_type,imin,imax,jmin,jmax,outside
INTEGER  ::idim2,jdim2,i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,id
Integer        :: STATUS
Real(kind=4)  ,dimension(:,:),pointer :: image, ScwExposition 
Real(kind=4)  ,dimension(:,:),pointer,Optional :: imagevar,normvar
Real(kind=4),dimension(:,:),pointer  :: carte  
Real(kind=4),dimension(:,:),pointer ,Optional :: cartevar ,exposure,surf
Real(kind=4) ::ci1,cj1,ci2,cj2,detmeanval,time
Real(KIND=8) ::aref,dref,lref,bref
Real(KIND=8) ::cb,ce,i_cen_inp,j_cen_inp,i_cen_out,j_cen_out
Real(KIND=8) :: pixel1,pixel2,pixel_ang1,pixel_ang2
Real(KIND=8) :: dispi_in,dispj_in, dispi_out,dispj_out
Real(kind=8),Dimension(3,3)::rot_mat1,rot_mat2,rot_mat3,rt,rt2,rot_mat4 

! AUXILIARY VARIABLES
INTEGER::i,j,in,jn,stat
Real(KIND=8) ::ri,rj,xn,yn,xn1,yn1,an,bn
Real(KIND=8) ::xnn,ynn,kx,ky





do j=jmin,jmax
   rj=Real(j)-cb-j_cen_inp! + dispj_in    
   do i=imin,imax
      ri = Real(i)-cb-i_cen_inp!+ dispi_in 

!!$      if((abs(i-519)<2).and.(abs(j-338)<2))then
!!$         print *,i,j
!!$      endif

!!$      if((i==324).and.(j==522))then
!!$         print *,i,j
!!$      endif

      call projection(is,tc,pixel1,pixel2,&
                      pixel_ang2,aref,dref,lref,bref,&
                      rot_mat1,rot_mat3,rot_mat4,&
                      ri,rj,xn,yn,stat)

      if(stat==0)then
          xnn = xn-ce +i_cen_out- dispi_out
          ynn = yn-ce +j_cen_out- dispj_out
         ! surface = (xnn-xnnp)*(ynn-ynnp)
          in = Nint(xnn ) 
          jn = Nint(ynn )
          If(MIN0(in,jn,i_dim_out-in+id,&
            j_dim_out-jn+id).Gt.0)Then
            If(((is==0).and.(nanfilter(i,j)==1)).or.&
              ((is==1).and.(nanfilter(in,jn)==1)))then

!!$if((in==203).and.(jn==155))then
!!$   print *, i,j,in,jn,image(in,jn)
!!$endif


            select case(rot_type)
            case(0) 
               ! NO FLUX SPREAD
               ! IMAGE/VARIANCE ROTATION
               call add_val0(is,image_type,&
                    i,j,in,jn,time,image,ScwExposition,imagevar,&
                     carte,cartevar,exposure)
           
             case(3) 
               ! NO FLUX SPREAD
               ! IMAGE/VARIANCE ROTATION
               ! source peak conservation
 
              call add_val01_point(is,&
                   i_dim_inp,j_dim_inp,i_dim_out,j_dim_out,idim2,&
                   jdim2,i,j,in,jn,&
                   image,imagevar,carte,cartevar)
            case(1)
               ! FLUX SPREAD- surface

               an = xnn-in
               if(an < 0.) then
                  in = in-1
                  an = xnn-in
               endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif
               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(FilterImage(i,j))).or.&
                       ((is==1).and.(FilterImage(in,jn))))then  
                     call add_val1(is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,exposure)
                  endif
               endif

            case(2)
               ! FLUX SPREAD- approx
               ! PIXEL BOTTOM LEFT CORNER in xn,yn
               ! PIXEL UPPER RIGHT CORNER
                 ky=Real(j)+cb-j_cen_inp
                 kx=Real(i)+cb-i_cen_inp
                 call projection(is,tc,pixel1,pixel2,&
                      pixel_ang2,aref,dref,lref,bref,&
                      rot_mat1,rot_mat3,rot_mat4,&
                      kx,ky,xn1,yn1,stat)
                  if(stat==0)then
                      kx = xn1+i_cen_out
                      ky = yn1+j_cen_out
                      an = (kx-in-0.5)*mf
                      If (an.Lt.0.0d0)an=0.0d0
                      bn = (ky-jn-0.5)*mf
                      If (bn.Lt.0.0d0)bn=0.0d0
                      call add_val1(is,image_type,&
                      i,j,in,jn,an,bn,time,&
                      image,ScwExposition,imagevar,&
                      carte,cartevar,exposure)
                   Else
                      status = stat
                   Endif 
          

           case(4) ! willmore corrected
              
              an = xnn-in
              if(an < 0.) then
                 in = in-1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif
               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(FilterImage(i,j))).or.&
                       ((is==1).and.(FilterImage(in,jn))))then  


                     call add_val2(is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,scwexposition,imagevar,&
                          carte,cartevar,normvar,exposure)
                  endif
               endif

            case(5) !no spread variance corrected
               
               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                   
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  

                     call add_val01(is,image_type,&
                          i,j,in,jn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)
                  endif
               endif

            case(6) ! Willmore var = sum(k_ierr_i)**2
              
              an = xnn-in
              if(an < 0.) then
                 in = in-1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif
               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(FilterImage(i,j))).or.&
                       ((is==1).and.(FilterImage(in,jn))))then  


                     call add_val3(is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)
                  endif
               endif


            case(7) ! Willmore correct var
              
              an = xnn-in
              if(an < 0.) then
                 in = in-1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  


                     call add_val4(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif


 case(8) ! sum errors
              
              an = xnn-in
              if(an < 0.) then
                 in = in-1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  

                     call add_val8(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif
 case(9) ! 0.75 covar
              
              an = xnn-in
              if(an < 0.) then
                 in = in-1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  
                     
                     call add_val9(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif

 case(17) ! Willmore correct var
              
              an = xnn-in
              if(an < 0.) then
                 in = in-1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  


                     call add_val17(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif


 case(19) ! 0.75 covar
              
              an = xnn-in
              if(an < 0.) then
                 in = in-1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  
                     
                     call add_val9a(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif
 case(10)!LS fit
              
              an = xnn-in
              if(an >0.5) then
                 in = in+1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn > 0.5) then
                  jn = jn+1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  
!!$ if((i==4321).and.(j==576))then
!!$    print*,in,jn
!!$endif
!!$ if((i==2194).and.(j==641))then
!!$    print*,in,jn
!!$endif
  if((abs(i-456) < 3).and.(j==400))then
   print*,i,j,in,jn
endif      
         
                     call add_val10(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif
 case(11)!chi2 fit
              
              an = xnn-in
              if(an >0.5) then
                 in = in+1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn > 0.5) then
                  jn = jn+1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  

                     call add_val11(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif

 case(12)!Ls fit +/1 1
              
              an = xnn-in
              if(an >0.5) then
                 in = in+1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn > 0.5) then
                  jn = jn+1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  

                     call add_val12(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif
 case(13)!Ls fit +/1 1
      
              an = xnn-in
              if(an >0.5) then
                 in = in+1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn > 0.5) then
                  jn = jn+1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  

                     call add_val13(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif

 case(14)!Ls fit +/1 1
      
              an = xnn-in
              if(an >0.5) then
                 in = in+1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn > 0.5) then
                  jn = jn+1
                  bn = ynn-jn
               endif
!!$ if((abs(i-456) < 3).and.(j==400))then
!!$   print*,i,j,in,jn
!!$endif    

               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  

                     call add_val14(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif

 case(15)!Ls fit +/1 1
      
              an = xnn-in
              if(an >0.5) then
                 in = in+1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn > 0.5) then
                  jn = jn+1
                  bn = ynn-jn
               endif
!!$ if((abs(i-456) < 3).and.(j==400))then
!!$   print*,i,j,in,jn
!!$endif    

               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  

                     call add_val15(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif

  case(16) 
              
              an = xnn-in
              if(an < 0.) then
                 in = in-1
                 an = xnn-in
              endif
               bn = ynn-jn
               if(bn < 0.) then
                  jn = jn-1
                  bn = ynn-jn
               endif


               If(MIN0(in,jn,i_dim_out-in+id,&
                    j_dim_out-jn+id).Gt.0)Then
                  If(((is==0).and.(nanfilter(i,j)==1)).or.&
                       ((is==1).and.(nanfilter(in,jn)==1)))then  


                     call add_val16(num1,detmeanval,is,image_type,&
                          i,j,in,jn,an,bn,time,&
                          image,ScwExposition,imagevar,&
                          carte,cartevar,normvar,exposure)


                  endif
               endif
    end select




           endif ! FilterImage
          else !MIN0
            if(is.eq.0)then
              if(info.gt.1)then
                  Print *,'ATTENTION  outside of the carte '
                  Print *,i,j,in,jn,idim2,jdim2
              endif
              outside = outside+1
           endif
        Endif  !MIN0
      Else
          status = 0
      Endif
  
 
   enddo ! I LOOP
enddo    ! J LOOP


if(info.gt.1)then
   print *,'  rotat sum : ',sum(carte)
   if(outside .gt.0)then
      print *, 'Attention : ',outside ,' pixels outside a map'
   endif
endif

!print *,' minmax image zone treated :',minin,maxin,' ',minjn,maxjn
!============================
END SUBROUTINE ROTAT_LOOP_ONE
!============================
!######################################
!INTERNAL SUBROUTINES WITHOUT VARIABLES
!######################################

!=============================================
Subroutine INIT
!=============================================

Use ATTI_DEFS
Use IBIS_IMAGING_PARAM
Implicit None

! DOUBLE PRECISION DEFINITIONS
dpi = 2.0D+0 * DASIN(1.0D+0)
rad_to_deg = 180.0d0/dpi
deg_to_rad=dpi/180.0d0
d360_rad=360.0d0*deg_to_rad
d180_rad=180.0d0*deg_to_rad
d270_rad=270.0d0*deg_to_rad
d90_rad=90.0d0*deg_to_rad

!SINGLE PRECISION DEFINITIONS
spi = dpi
sdeg_to_rad = deg_to_rad
srad_to_deg = rad_to_deg
s360_rad=d360_rad
s180_rad=d180_rad
s270_rad=d270_rad
s90_rad=d90_rad



! RADIUS FOR STOCKING OF UNUSUAL DIRECTIONS
yz_radius = 10d+4
!ARRIVAL ANGLE TOLERANCE
x_tol = 10d-4
d_one = 1.d0+10.0d-10
!pixel_ratio = 9.4/2500.0    ! SIGMA
pixel_ratio = IsgrPixSiz / IsgrHalfMasDis !- INTEGRAL
pixel_ang_size =  rad_to_deg*Atan(pixel_ratio)

If(message_level ==1) Then
   Print *,' INTEGRAL pixel definition  : ',pixel_ratio
Endif
!=============================================
End Subroutine INIT
!=============================================


!=============================================
Subroutine Read_ATTI_Data
!=============================================

Use ATTI_DEFS
Implicit None

err = 0


Read(3,*) Data_Type

Type:Select Case(Data_Type)

Case (1) Type

      If(message_level==1)Then
         Print *,' GIVEN X(RA,DEC), theta'
      Endif 
       
      Read(3,*) alphaXtel_deg,deltaXtel_deg,theta_deg

      If ((alphaXtel_deg.Lt.0.0d0).Or.(alphaXtel_deg.Gt.360.d0)&
     .Or.(deltaXtel_deg.Lt.-90.d0).Or.(deltaXtel_deg.Gt.90.d0)&
     .Or.(theta_deg.Lt.0.0d0).Or.(theta_deg.Gt.360.d0))&
          err = Data_Type

     
Case (2) Type

      If(message_level==1)Then
         Print *,' GIVEN X(RA,DEC), Z(RA,DEC)'
      Endif 

      Read(3,*) alphaXtel_deg,deltaXtel_deg,&
                alphaZtel_deg,deltaZtel_deg
      If ((alphaXtel_deg.Lt.0.0d0).Or.(alphaXtel_deg.Gt.360.d0)&
     .Or.(deltaXtel_deg.Lt.-90.d0).Or.(deltaXtel_deg.Gt.90.d0))&
          err = Data_Type
      If ((alphaZtel_deg.Lt.0.0d0).Or.(alphaZtel_deg.Gt.360.d0)&
     .Or.(deltaZtel_deg.Lt.-90.d0).Or.(deltaZtel_deg.Gt.90.d0))&
          err = Data_Type 

Case (3) Type
 
      If(message_level==1)Then
         Print *,' GIVEN X_eq, Z_eq'
      Endif 

      Read(3,*) vector

      Call sla_dvn(vector, Xtel_eq,norm)

       If(norm.Ne.0.d0) Then
        If(message_level==1)Then
          Print *,' Xtel_eq normalised '
        Endif
       Endif 

      Read(3,*)vector

      Call sla_dvn(vector,Ztel_eq,norm)

      If(norm.Ne.0.d0) Then
        If(message_level==1)Then
          Print *,'Ztel_eq normalised '
        Endif
       Endif 

Case (4) Type 

      If(message_level==1)Then
         Print *,' GIVEN X_gal, Z_gal'
      Endif 
  
      Read(3,*) gal_sys_Type
      
      Read(3,*) l2Xtel,b2Xtel
      Read(3,*) l2Ztel,b2Ztel
            
      If ((gal_sys_Type .Lt. 1)  .Or. (gal_sys_Type .Gt. 2) &
       .Or.(l2Xtel.Lt.0.0d0)  .Or. (l2Xtel.Gt.360.0d0)&
       .Or.(b2Xtel.Lt.-90.0d0).Or. (b2Xtel.Gt.90.0d0) &
       .Or.(l2Ztel.Lt.0.0d0)  .Or. (l2Ztel.Gt.360.0d0)&
       .Or.(b2Ztel.Lt.-90.0d0).Or. (b2Ztel.Gt.90.0d0))&
          err = Data_Type


Case default

      err = 5

End Select Type 

 If(err.Ne.0) Then
    err = err+10
    If(message_level==1)Then
           Print *,' ERROR IN DATA  ',err
    Endif
 Endif



!=============================================
End Subroutine Read_ATTI_Data
!=============================================
      
!=============================================
Subroutine CALCUL_THETA
!=============================================

! CALCULATES THE ANTICLOCWISE ANGLE THETA  FROM
! THE + NORD DIRECTION AXIS UNTIL THE +Z SATELLITE AXIS.
! THETA IS IN THE TANGENT PLANE (Y,Z) OF THE SATELLITE COORDINATE SYSTEM.
! THE + NORD DIRECTION AXIS IS DEFINED AS THE   ORTHOGONAL( X DIRECTION )
! PROJECTION OF THE POLAR AXIS OF THE EQUATORIAL SYSTEM ONTO THE TANGENT
!  PLANE (Y,Z).
! THE EQUATORIAL AND SATELLITE SYSTEM ARE TAKEN TO HAVE THE SAME ORIGIN
! AT THE CENTRE OF THE CELESTIAL SPHERE.

Use ATTI_DEFS
Implicit None


Real(kind=8) :: dir_angle

!CALCUL OF NORD DIRECTION AXIS Znord IN THE TANGENT PLANE

!!$if(deltaXtel_rad.ge.0.)then
!!$   alphaZnord = d180_rad+alphaXtel_rad
!!$   !NORMALISATION into [0 ,2*pi]
!!$   alphaZnord = sla_dranrm(alphaZnord)
!!$   deltaZnord =d90_rad-deltaXtel_rad
!!$   ! NORMALISATION INTO  [-pi, *pi]
!!$   deltaZnord = sla_drange(deltaZnord)
!!$else
!!$  alphaZnord = alphaXtel_rad
!!$  deltaZnord = d90_rad+deltaXtel_rad
!!$  ! NORMALISATION INTO  [-pi, *pi]
!!$  deltaZnord = sla_drange(deltaZnord)
!!$endif
!!$
!!$
!!$! THETA - ANGLE BETWEEN Znord AXIS AND Ztel axis
!!$theta_rad = sla_dsep(alphaZnord,deltaZnord,&
!!$                         alphaZtel_rad,deltaZtel_rad)
!!$dir_angle = sla_dbear(alphaZnord,deltaZnord,&
!!$                         alphaZtel_rad,deltaZtel_rad)
!!$!If(dir_angle.Gt.0) Then
!!$If(dir_angle.lt.0) Then
!!$    theta_rad = d360_rad-theta_rad
!!$Endif



!!$theta_rad = Mod(theta_rad,d360_rad)
!!$theta_deg = theta_rad*rad_to_deg

theta_rad = theta_deg*deg_to_rad

If(message_level==1)Then
         Print *,' theta = ',theta_deg
!         Print *,' position_angle = ',dir_angle
      Endif
!=============================================
End Subroutine CALCUL_THETA
!=============================================

!=============================================
Subroutine LABO_TRANSFORM
!=============================================
! AUXILIARY OPERATIONS
Use ATTI_DEFS
Implicit None


Real(kind=8) :: a,d

Type:Select Case(Data_Type)

Case (1) Type
   
     alphaXtel_rad = alphaXtel_deg*deg_to_rad
     deltaXtel_rad = deltaXtel_deg*deg_to_rad
     theta_rad = theta_deg*deg_to_rad
      
Case (2) Type
      
      alphaXtel_rad = alphaXtel_deg*deg_to_rad
      deltaXtel_rad = deltaXtel_deg*deg_to_rad

     
      alphaZtel_rad = alphaZtel_deg*deg_to_rad
      deltaZtel_rad = deltaZtel_deg*deg_to_rad

      a = sla_dsep(alphaXtel_rad,deltaXtel_rad,alphaZtel_rad,deltaZtel_rad)
      if(message_level==1)then
            print *,&
            ' axes precision : ',dabs(a-d90_rad),&
            ' required : ',axes_precision
      endif
       If(dabs(a-d90_rad).Gt.axes_precision) Then
          err = Data_Type
      Else
          Call CALCUL_THETA
      Endif

Case (3) Type

     ! (x,y,z) --> (a,d) for axis X
     Call sla_dcc2s(Xtel_eq,alphaXtel_rad,deltaXtel_rad)
     ! normalisation
     alphaXtel_rad = sla_dranrm(alphaXtel_rad)
     deltaXtel_rad = sla_drange(deltaXtel_rad)
     alphaXtel_deg = alphaXtel_rad*rad_to_deg
     deltaXtel_deg = deltaXtel_rad*rad_to_deg
     ! (x,y,z) --> (a,d) for axis Z
     Call sla_dcc2s(Ztel_eq,alphaZtel_rad,deltaZtel_rad)
     ! normalisation
     alphaZtel_rad = sla_dranrm(alphaZtel_rad)
     deltaZtel_rad = sla_drange(deltaZtel_rad)
     alphaZtel_deg = alphaZtel_rad*rad_to_deg
     deltaZtel_deg = deltaZtel_rad*rad_to_deg

     
      a = sla_dsep(alphaXtel_rad,deltaXtel_rad,alphaZtel_rad,deltaZtel_rad)
      if(message_level==1)then
            print *,&
            ' axes precision : ',dabs(a-d90_rad),&
            ' required : ',axes_precision
      endif
     If(dabs(a-d90_rad).Gt.axes_precision) Then
          err = Data_Type
      Else
          Call CALCUL_THETA
      Endif



Case (4) Type

     ! TRANSFORMATION :  GALACTIC TO EQUATORIAL SYSTEM

     sys: Select Case (gal_sys_Type)

     Case(1)sys


         Call sla_ge50(l2Xtel,b2Xtel,a,d)
         alphaXtel_rad = a
         deltaXtel_rad = d
         alphaXtel_deg = alphaXtel_rad*rad_to_deg
         deltaXtel_deg = deltaXtel_rad*rad_to_deg

         Call sla_ge50(l2Ztel,b2Ztel,a,d)
         alphaZtel_rad = a
         deltaZtel_rad = d
         alphaZtel_deg = alphaZtel_rad*rad_to_deg
         deltaZtel_deg = deltaZtel_rad*rad_to_deg

      Case(2)sys

         Call sla_galeq(l2Xtel,b2Xtel,a,d)
         alphaXtel_rad = a
         deltaXtel_rad = d
         alphaXtel_deg = alphaXtel_rad*rad_to_deg
         deltaXtel_deg = deltaXtel_rad*rad_to_deg

         Call sla_GALEQ(l2Ztel,b2Ztel,a,d)
         alphaZtel_rad = a
         deltaZtel_rad = d
         alphaZtel_deg = alphaZtel_rad*rad_to_deg
         deltaZtel_deg = deltaZtel_rad*rad_to_deg
      End Select sys

      a = sla_dsep(alphaXtel_rad,deltaXtel_rad,alphaZtel_rad,deltaZtel_rad)
       if(message_level==1)then
            print *,&
            ' axes precision : ',dabs(a-d90_rad),&
            ' required : ',axes_precision
      endif
      If(dabs(a-d90_rad).Gt.axes_precision) Then
          err = Data_Type
      Else
          Call CALCUL_THETA
      Endif



End Select Type

 If(err.Ne.0) Then
    err = err+20
    If(message_level==1)Then
           Print *,' ERROR IN DATA  ',err
    Endif
Else
   If(message_level==1)Then
   Print *,' Xtel(ra,dec) = ',alphaXtel_deg,'  ',deltaXtel_deg
   Print *,' theta = ',theta_deg
   Endif
Endif

!=============================================
End Subroutine LABO_TRANSFORM
!=============================================


!================================================================
 Subroutine CREATE_ROTATION_MATRIX
!================================================================
! ROTATION  MATRIX  G *X( IN EQUATORIAL S.) = Y(IN TELESCOdPE S.)
!...............................................................

Use ATTI_DEFS
Implicit None
REAL(kind=8)::dd
 dd = mod(d360_rad-deltaXtel_rad,d360_rad)

Call sla_deuler('321',&
                 alphaxtel_rad,dd,theta_rad,rot_mat)



!==========================================
 End Subroutine CREATE_ROTATION_MATRIX
!==========================================



