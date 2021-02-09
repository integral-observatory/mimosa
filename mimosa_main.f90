
!============================   
PROGRAM MIMOSA_main
!============================   
USE ISDC
USE MIMOSA_CONTROL_MODULE
USE MIMOSA_WORK_MODULE
Use MIMOSA_USE_MODULE

IMPLICIT NONE

 
   
INTEGER                      :: mode 
INTEGER                      :: status 
CHARACTER (len=30)                   :: title
INTEGER,dimension(8)         :: values1,values2
CHARACTER (len=10)                   :: version='1.3'
CHARACTER (len=5) :: zone
CHARACTER (len=10)                   :: str,time1,time2
CHARACTER (len=8)                    :: date
INTEGER                      :: rilStat
CHARACTER (len=DAL_MAX_STRING)       :: procName





status = ISDC_OK
procName = 'mimosa'
title = 'EXEC mimosa '

DO

      ! Common initialization

      mode = COMMON_INIT(procName, version)
      IF (mode /= ISDC_SINGLE_MODE) THEN
         status = mode      
         rilStat = RIL_LOG_MESSAGE (null, error_0, &
                   'Program not in the SINGLE mode')
         EXIT
      ENDIF
    
      rilStat = RIL_LOG_MESSAGE (null, log_0,title//' start ')
      if(WorkMode > 0) then
         call date_and_time(date,time1,zone,values1)
         rilStat = RIL_LOG_MESSAGE (null, log_0,&
              procName//' START TIME '//time1(1:2)//':'//time1(3:4))
      endif
      !******************************
      CALL MIMOSA_work(status)
      !******************************
      if(WorkMode > 0) then
         call date_and_time(date,time2,zone,values2)
         rilStat = RIL_LOG_MESSAGE (null, log_0,&
              procName//' START TIME '//time1(1:2)//':'//time1(3:4)&
              //' END TIME '//time2(1:2)//':'//time2(3:4))
      endif
      IF (status /= ISDC_OK) THEN
         rilStat = RIL_LOG_MESSAGE (null, error_0, &
                   'Program exit with error')
      ENDIF

      EXIT
ENDDO

write(str,'(I10)')status
rilStat = RIL_LOG_MESSAGE(null,log_0,title//'output status : '//str//&
                           'end of execution')
CALL COMMON_EXIT (status) 


!============================   
END PROGRAM MIMOSA_main
!============================   

