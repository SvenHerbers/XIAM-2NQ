C
C     This file contains system dependend subroutines and functions.
C     Modify it for your system/compiler!
C     the "myand" and "myor" functions are demanded by xiam,
C     "mysignal" and "mydate" are optional. 
C

      integer function myand(i1,i2)
C     some compiler don't know the generic "and",
C     use "iand" instead.
C     if neiher "and" nor "iand" is available you're in trouble here. 
      implicit none
      integer i1,i2
      myand=and(i1,i2)
      return
      end

      integer function myor(i1,i2)
C     some compiler don't know the generic "or",
C     use "ior" instead.
      implicit none
      integer i1,i2
      myor=or(i1,i2)
      return
      end

      subroutine mysignal()
      implicit none
      integer sig_stat
      common/sig_com/sig_stat
      external sig_func

C     for most unix systems signal(2) means SIGINT=control C (see e.g. /usr/include/signal.h) 
C     if signal is not known to your (unix) f77 compiler,
C     try _signal or signal_
C     on other machines (DOS, VMS), comment the next line out.
C
C     signal works with: 
C     AIX xlf
C
C     signal does not work with:
C     linux g77

C     call c_signal()
      return
      end

      subroutine mydate()
C     this subroutine writes the date and time to stdout 
C     other compiler may use other functions.
C     VMS use SYS$... functions.
      implicit none

C     If you find no corresponding commands for your compiler simple
C     comment it out. 
C     call dateput()
      return
      end
