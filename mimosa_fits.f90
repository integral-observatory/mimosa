

!===========================
MODULE FITS_DEFS
!==========================

INTEGER,save :: BIT_PAR_PIXEL=-32 ! PIXEL LENGTH IN BITS
INTEGER,save :: naxis=2          ! NUMBER OF DIMENSIONS
INTEGER,save :: ieq = 2000       ! COORDINATE EPOCH
INTEGER,save :: lenh  = 36       ! FITS HEADER LENGTH
CHARACTER(len=20)::radesys='FK5              '!MEAN PLACE ETC
!==============================
END MODULE FITS_DEFS
!==============================

!==============================
MODULE FITS_DECLARATIONS
!==============================

INTERFACE


 
SUBROUTINE FITS_FILE(num,image,proj,nom_out,iprint,info)
!-------------------------------------------------
IMPLICIT NONE

!INPUT VARIABLES

INTEGER::num                        ! num - POINTING NUMBER
REAL(kind=4),dimension(:,:),pointer :: image !image array
INTEGER                     :: proj
! PROJECTIONtype  BETWEEN SPHERICAL(NATIVE) AND PLANE(WORLD) COORDINATES
CHARACTER(len=*)                    :: nom_out
! FITS FILE NAME
INTEGER,optional            :: iprint
CHARACTER(len=*),Optional           :: info
END SUBROUTINE FITS_FILE


SUBROUTINE FITS_FILE8(num,image,proj,nom_out,iprint,info)
!-------------------------------------------------
USE FITS_DEFS
IMPLICIT NONE

!INPUT VARIABLES

INTEGER::num                        ! num - POINTING NUMBER
REAL(kind=8),dimension(:,:):: image !image array
INTEGER :: proj
! PROJECTIONtype  BETWEEN SPHERICAL(NATIVE) AND PLANE(WORLD) COORDINATES
CHARACTER(len=*):: nom_out
! FITS FILE NAME
INTEGER,optional::iprint
CHARACTER(len=*),Optional::info
END SUBROUTINE FITS_FILE8

END INTERFACE 

!==============================
END MODULE FITS_DECLARATIONS
!==============================


!=============================================
SUBROUTINE FITS_FILE(num,image,proj,nom_out,iprint,info)
!=============================================
USE FITS_DEFS
IMPLICIT NONE

!INPUT VARIABLES
INTEGER::num  ! num - POINTING NUMBER
REAL(kind=4),dimension(:,:),pointer:: image
INTEGER :: proj
! PROJECTION BETWEEN SPHERICAL(NATIVE) AND PLANE(WORLD) COORDINATES
CHARACTER(len=*):: nom_out
! FITS FILE NAME
INTEGER,optional::iprint
CHARACTER(len=*),optional::info

!AUXILIARY VARIABLES
CHARACTER(len=80):: hline(36)
INTEGER::nx,ny


nx = size(image,1)
ny = size(image,2)

! 36*80 = 2880 - 36 LINES OF HEADER

 OPEN(90,FILE=nom_out,ACCESS='direct',RECL=2880)
CALL create_fits_header(num,nx,ny,proj,hline)
if(present(iprint))PRINT 101, hline
 101  FORMAT(A80)
WRITE(90,REC=1) hline
CALL wf_data(image,nx*ny,720)
CLOSE(90)
if(present(info))then
write(info,*)' file (',nx,',',ny,') : ',TRIM(nom_out),' created '
endif
!=======================
END SUBROUTINE FITS_FILE
!=======================

!=============================================
SUBROUTINE FITS_FILE8(num,image,proj,nom_out,iprint,info)
!=============================================
USE FITS_DEFS
IMPLICIT NONE

!INPUT VARIABLES
INTEGER::num  ! num - POINTING NUMBER
REAL(kind=8),dimension(:,:):: image
INTEGER :: proj
! PROJECTION BETWEEN SPHERICAL(NATIVE) AND PLANE(WORLD) COORDINATES
CHARACTER(len=*):: nom_out
! FITS FILE NAME
INTEGER,optional::iprint
CHARACTER(len=*),optional::info

!AUXILIARY VARIABLES
CHARACTER(len=80):: hline(36)
INTEGER::nx,ny


nx = size(image,1)
ny = size(image,2)
! 36*80 = 2880 - 36 LINES OF HEADER

OPEN(90,FILE=nom_out,ACCESS='direct',RECL=2880)
CALL create_fits_header(num,nx,ny,proj,hline)
if(present(iprint))PRINT 101, hline
 101  FORMAT(A80)
WRITE(90,REC=1) hline
CALL wf_data8(image,nx*ny,720)
CLOSE(90)
if(present(info))then
write(info,*)' file (',nx,',',ny,') : ',TRIM(nom_out),' created '
endif
!=======================
END SUBROUTINE FITS_FILE8
!=======================


!=====================================================
SUBROUTINE create_fits_header(num,nx,ny,proj,hline)
!=====================================================

USE FITS_DEFS
USE ATTI_DEFS
IMPLICIT NONE

!INPUT VARIABLES
INTEGER :: num       !POINTING NUMBER
INTEGER :: nx,ny       !IMAGE DIMENSIONS IN PIXEL
INTEGER :: proj      ! PROJECTION TYPE 
                     ! 0 - FOR CAR
                     ! 1 - FOR TAN 

!INPUT/OUTPUT VARIABLES
CHARACTER(len=80),dimension(lenh):: hline !HEADER
!AUXILIARY VARIABLES
REAL(kind=4)::cenx,ceny !IMAGE CENTRE
REAL(KIND=4)::alpha,delta ! POINTING RA DEC
REAL(KIND=4)::th,th_rad,st0,ct0 !ANTICLOCKWISE NORTH ANGLE
REAL(KIND=4)::incr ! PIXEL INCREMENT IN DEG
INTEGER ::i, nline  ! LINE NUMBER
CHARACTER(len=80)::dum,cur
CHARACTER(len=20):: ctyp
LOGICAL ,save::simple=.true.



 WRITE(dum,111)
 111  FORMAT(80x)


!  RA DEC OF POINTING 
      
alpha = pointing_table(num,1)
delta = pointing_table(num,2)

th = pointing_table(num,3)

!print *,' alpha,delta,th : ', alpha,delta,th
th_rad = th*sdeg_to_rad
ct0 = cos(th_rad)
st0 = sin(th_rad)
!print *,' cos,sin : ',ct0,st0
!IMAGE CENTRE
! ATTENTION : PIXEL I RUNS FROM I-0.5 TO I+0.5
! 0 IS AT THE CENTRE OF 0 PIXEL

cenx = 0.5 *(1+nx)
ceny = 0.5 *(1+ny)

!PIXEL INCREMENT IN DEG

incr = srad_to_deg*atan(pixel_ratio)



!      IF (ABS(ct0(1)) .GT. ABS(st0(1))) THEN
!         theta = - ASIND(st0(1))
!         IF (ct0(1) .LT. 0.)  theta = theta + 180.
!      ELSE
!         theta = ACOSD(ct0(1))
!         IF (st0(1) .GT. 0.)  theta = theta + 180.
!      ENDIF


!DATE
!      annee = an(1) - 1900

!      CALL idate(idtab)
!      idtab(3) = idtab(3) - 1900

nline = 1


 WRITE(hline(nline),101) 'SIMPLE  ',simple,dum

nline = nline+1     
cur = 'REAL*4' // dum
WRITE(hline(nline),102)'BITPIX  ',BIT_PAR_PIXEL,cur

 101  FORMAT(A8,'= ',L20,' / ',A47)
 102  FORMAT(A8,'= ',I20,' / ',A47)
 103  FORMAT(A8,'= ',A20,' / ',A47)
 104  FORMAT(A8,'= ',A20,' / ',A47)
 105  FORMAT(A8,'= ',F20.4,' / ',A47)
 106  FORMAT(A8,'= ',1f20.4,' / ',A47)
 107  FORMAT(A8,'= ',1PE20.4,' / ',A47)

!NAXIS
nline = nline+1
cur = 'NUMBER OF AXES' // dum
WRITE(hline(nline),102) 'NAXIS   ',naxis,cur

 
!NAXIS 1
nline = nline+1 
cur = 'IMAGE SIZE ALONG x AXIS ' // dum
WRITE(hline(nline),102) 'NAXIS1  ',nx,cur


!NAXIS 2
nline = nline+1
cur = 'IMAGE SIZE ALONG y AXIS ' // dum
WRITE(hline(nline),102) 'NAXIS2  ',ny,cur

!EQUINOX    
nline = nline+1
cur = 'COORDINATES in 2000.0' // dum
WRITE(hline(nline),102) 'EQUINOX ',ieq,cur

!RADESYS    
nline = nline+1
cur = 'MEAN PLACE, IAU 1984 system' // dum
WRITE(hline(nline),103) 'RADESYS ',radesys,cur

! DESCRIPTION OF THE FIRST AXIS


! CTYPE1 : x AXIS TYPE
select case (proj)
case (0)
nline = nline+1
ctyp = 'RA---CAR'          
cur = 'CARTESIAN PROJECTION' // dum
WRITE(hline(nline),104) 'CTYPE1  ',ctyp,cur

case (1)
nline = nline+1
ctyp = 'RA---TAN'          
cur = 'GNOMONIC PROJECTION' // dum
WRITE(hline(nline),104) 'CTYPE1  ',ctyp,cur
end select



! CRVAL1 : REFERENCE VALUE
nline = nline+1
cur = 'RIGHT ASCENSION OF POINTING DIRECTION' // dum
WRITE(hline(nline),105) 'CRVAL1  ',alpha,cur


!CRPIX1
nline = nline+1
cur = 'X COORDINATE OF IMAGE CENTER' // dum
WRITE(hline(nline),105) 'CRPIX1  ',cenx,cur


!CDELT1
nline = nline+1
cur = 'DEGREES PER PIXEL ALONG x AXIS ' // dum
WRITE(hline(nline),106) 'CDELT1  ',incr,cur

! DESCRIPTION OF THE SECOND AXIS


! CTYPE 2 :  AXIS TYPE

select case(proj)
case(2)
nline = nline+1
ctyp = 'DEC--CAR'          
cur = 'CARTESIAN PROJECTION' // dum
WRITE(hline(nline),104) 'CTYPE2  ',ctyp,cur

case(1)
nline = nline+1
ctyp = 'DEC--TAN'          
cur = 'GNOMONIC PROJECTION' // dum
WRITE(hline(nline),104) 'CTYPE2  ',ctyp,cur
end select


! CRVAL2 : REFERENCE VALUE   
nline = nline+1
cur = 'DECLINATION OF POINTING DIRECTION' // dum
WRITE(hline(nline),105) 'CRVAL2  ',delta,cur
 


    
!CRPIX2
nline = nline+1
cur = 'Y COORDINATE OF IMAGE CENTER' // dum
WRITE(hline(nline),105) 'CRPIX2  ',ceny,cur


!CDELT2
nline = nline+1
cur = 'DEGREES PER PIXEL ALONG y AXIS ' // dum
 WRITE(hline(nline),106) 'CDELT2  ',incr,cur
      

! ROTATION MATRIX

nline = nline+1
cur = ' ROTATION MATRIX : CDELT1*cos(CROTA2)' // dum
WRITE(hline(nline),107) 'CD1_1   ',-incr*ct0,cur

nline = nline+1
cur = '                 : -CDELT2*sin(CROTA2)' // dum
WRITE(hline(nline),107) 'CD1_2   ',-incr*st0,cur

nline = nline+1
cur = '                 : CDELT1*sin(CROTA2)' // dum
WRITE(hline(nline),107) 'CD2_1   ',-incr*st0,cur

nline = nline+1
cur = '                 : CDELT2*cos(CROTA2)' // dum
WRITE(hline(nline),107) 'CD2_2   ',incr*ct0,cur

!LONGPOLE
nline = nline+1
cur = ' NORTH POLE LATITUDE IN THE NATIVE SYSTEM' // dum

WRITE(hline(nline),104) 'LONGPOLE','180',cur 

 DO i = nline+1,lenh-1
     hline(i) = dum
 ENDDO

hline(lenh)='END'//dum
!do i=1,36
!  hline(i) = hline(i)(1:78)//'A'//' '
!enddo
END SUBROUTINE create_fits_header

!==================================
SUBROUTINE wf_data8(image,dim,len)
!===================================

USE FITS_DEFS
IMPLICIT NONE

INTEGER::dim
REAL(kind=8),dimension(dim)::image
INTEGER ::nrec,i,n0,n1,n      
INTEGER::len
! 720 = 2880/4
REAL(kind=4),dimension(len)::fin
REAL(kind=4),dimension(dim)::image4

image4 = image

fin = 0.0
 nrec = (dim-1)/len + 1
  
 DO i = 1,nrec-1
     n0 = (i-1) *len
     WRITE(90,REC=i+1) (image4(n),n=n0+1,n0+len)
 ENDDO

! Ecriture du dernier enregistrement complete par des 0.

      n0 = (nrec-1) *len
      n1 = len - dim + n0
      WRITE(90,REC=nrec+1) (image4(n),n=n0+1,dim),(fin(n),n=1,n1)

!==================================
END SUBROUTINE wf_data8
!==================================



!==================================
SUBROUTINE wf_data(image,dim,len)
!===================================

USE FITS_DEFS
IMPLICIT NONE

INTEGER::dim
REAL(kind=4),dimension(dim)::image
INTEGER ::nrec,i,n0,n1,n      
INTEGER::len
! 720 = 2880/4
REAL(kind=4),dimension(len)::fin

fin = 0.0
 nrec = (dim-1)/len + 1
  
 DO i = 1,nrec-1
     n0 = (i-1) *len
     WRITE(90,REC=i+1) (image(n),n=n0+1,n0+len)
 ENDDO

! Ecriture du dernier enregistrement complete par des 0.

      n0 = (nrec-1) *len
      n1 = len - dim + n0
      WRITE(90,REC=nrec+1) (image(n),n=n0+1,dim),(fin(n),n=1,n1)


END SUBROUTINE wf_data


