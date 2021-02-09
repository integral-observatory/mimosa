
!*************************************
MODULE MIMOSA_USE_MODULE
!*************************************

INTERFACE
SUBROUTINE MESSAGE(ProcName,Comment,Val,Status)
!------------------------------------------------
IMPLICIT NONE 

!INPUT VARIABLES
CHARACTER(len=*),Intent(in)::ProcName,Comment
INTEGER,Intent(in)::Val
!OUTPUT VARIABLES
INTEGER,Intent(out)::Status
END SUBROUTINE MESSAGE

SUBROUTINE WAR_MESSAGE(ProcName,Comment,warlevel,Status)
!-------------------------------------------------
IMPLICIT NONE 

!INPUT VARIABLES
CHARACTER(len=*),Intent(in)::ProcName,Comment
INTEGER,Intent(in)::warlevel
!OUTPUT VARIABLES
INTEGER,Intent(out)::Status
END SUBROUTINE WAR_MESSAGE

SUBROUTINE WAR_CONT_MESSAGE(ProcName,Comment,warlevel,Status)
!-------------------------------------------------
IMPLICIT NONE 

!INPUT VARIABLES
CHARACTER(len=*),Intent(in)::ProcName,Comment
INTEGER,Intent(in)::warlevel
!OUTPUT VARIABLES
INTEGER,Intent(out)::Status
END SUBROUTINE WAR_CONT_MESSAGE

FUNCTION GET_SEC()
!---------------
IMPLICIT NONE
REAL(kind=4) :: get_sec
END FUNCTION GET_SEC


END INTERFACE

!*************************************
END MODULE MIMOSA_USE_MODULE
!*************************************

!==============================================
SUBROUTINE MESSAGE(ProcName,Comment,Val,Status)
!==============================================

USE ISDC
USE MIMOSA_CONTROL_MODULE
IMPLICIT NONE 

!INPUT VARIABLES
CHARACTER(len=*),Intent(in)::ProcName,Comment
INTEGER,Intent(in)::Val
!OUTPUT VARIABLES
INTEGER,Intent(out)::Status

!LOCAL  VARIABLES
INTEGER::n,rilStat,k
CHARACTER(len=20) :: str20
CHARACTER(len=45) ::str45
CHARACTER(len=300) :: str300

rilStat =ISDC_OK

n=len_trim(ProcName)
if(n.le.20)then
   str20 = ProcName(1:n)//':'
else
   str20 = ProcName(1:19)//':'
endif



if(Val.eq.ZeroError)then
   str300 = str20(1:len_trim(str20))//Comment
   n = len_trim(str300)
   do while(n.gt.0)
     k=index(str300,' ',.true.)
     rilStat =  RIL_LOG_MESSAGE(null,log_0,str300(1:k))
     if(rilStat.ne.ISDC_OK)goto 101
     str300 = str300(k+1:n)
     n = len_trim(str300)
  enddo
  
else
  Status = val
  write(str300,*)str20(1:len_trim(str20)),' Status : ',Status 
  str300 = '!! : '//str45(1:len_trim(str45))//Comment
  n = len_trim(str300)
  do while(n.gt.0)
     k=index(str300,' ',.true.)
     if(Val.eq.IsdcProcError)then
        rilStat =  RIL_LOG_MESSAGE(null,error_1,str300(1:k))
         if(rilStat.ne.ISDC_OK)goto 101
     else
        rilStat =   RIL_LOG_MESSAGE(null,error_2,str300(1:k))
         if(rilStat.ne.ISDC_OK)goto 101
     endif
    
     str300 = str300(k+1:n)
     n = len_trim(str300)
  enddo


 
endif

101 Status=Val
!======================================
END SUBROUTINE MESSAGE     
!=====================================


!==============================================
SUBROUTINE WAR_MESSAGE(ProcName,Comment,warlevel,Status)
!==============================================

USE ISDC
USE MIMOSA_CONTROL_MODULE
IMPLICIT NONE 

!INPUT VARIABLES
CHARACTER(len=*),Intent(in)::ProcName,Comment
!warning level : between 0 and 3
INTEGER,Intent(in)::warlevel
!OUTPUT VARIABLES
INTEGER,Intent(out)::Status

!LOCAL  VARIABLES
INTEGER::rilStat,n,k,war,last
CHARACTER(len=20)::statStr
CHARACTER(len=300)::str300



if(Status.ne.ISDC_OK)then
   write(statStr,'(I10)')Status    
   Status=ISDC_OK
   str300 =ProcName(1:len_trim(ProcName))//': ATTENTION ! Status :  '//statStr//Comment 
else
   str300 =ProcName(1:len_trim(ProcName))//': ATTENTION ! '//Comment
endif

rilStat =ISDC_OK

war = warlevel
if(war > 3)war = 3
if(war < 0 )war = 0
war = WARNING_0+war

n = len_trim(str300)

if(n .le. 80)then
   rilStat =   RIL_LOG_MESSAGE(null,war,str300)
   if(rilStat.ne.ISDC_OK)goto 102
else
   do while(n.gt.0)
      
      if (n < 300)then
         last = min(80,n+1)
      else
         last = 80
      endif
      k = index(str300(1:last),' ',.true.)
      if( k == 0)k = last 
      rilStat =   RIL_LOG_MESSAGE(null,war,str300(1:k))
      if(rilStat.ne.ISDC_OK)goto 102
      str300 = str300(k+1:n)
      n = len_trim(str300)
   enddo

endif
102  continue
!======================================
END SUBROUTINE  WAR_MESSAGE     
!=====================================

!==============================================================
SUBROUTINE WAR_CONT_MESSAGE(ProcName,Comment,warlevel,Status)
!==============================================================

USE ISDC
USE MIMOSA_CONTROL_MODULE
IMPLICIT NONE 

!INPUT VARIABLES
CHARACTER(len=*),Intent(in)::ProcName,Comment
!warning level : between 0 and 3
INTEGER,Intent(in)::warlevel
!OUTPUT VARIABLES
INTEGER,Intent(out)::Status

!LOCAL  VARIABLES
INTEGER::rilStat,n,k,war,last
CHARACTER(len=300)::str300

Status=ISDC_OK
str300 =ProcName(1:len_trim(ProcName))//Comment 
n = len_trim(str300)

war = warlevel
if(war > 3)war = 3
if(war < 0 )war = 0
war =WARNING_0+ war

rilStat=ISDC_OK

if(n .le. 80)then
   rilStat =   RIL_LOG_MESSAGE(null,war,str300)
   if(rilStat.ne.ISDC_OK)goto 103
else
   do while(n.gt.0)
      
      if (n < 300)then
         last = min(80,n+1)
      else
         last = 80
      endif
      k = index(str300(1:last),' ',.true.)
      if( k == 0)k = last 
      rilStat =   RIL_LOG_MESSAGE(null,war,str300(1:k))
      if(rilStat.ne.ISDC_OK)goto 103
      str300 = str300(k+1:n)
      n = len_trim(str300)
   enddo

endif


103 continue

!======================================
END SUBROUTINE  WAR_CONT_MESSAGE     
!=====================================


!=================
FUNCTION GET_SEC()
!=================
IMPLICIT NONE
REAL(kind=4) :: get_sec

CHARACTER(len=8) :: date
CHARACTER(len=10) :: time
CHARACTER(len=10) :: zone
INTEGER,dimension(8) :: values

call date_and_time(date,time,zone,values)

get_sec= (values(5)*60. + values(6))*60. + values(7)+values(8)/1000.
!====================
END FUNCTION GET_SEC
!===================


