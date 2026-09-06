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

