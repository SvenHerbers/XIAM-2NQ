C     -------------------------------------------------------
C     getbuf(gu,unit)
C     nx_ss(gu)
C     pr_ss(gu)
C     rd_ss(gu,ss)
C     len_ss(gu)
C     end_ss(gu)
C
C     -------------------------------------------------------------------- 
      integer function getbuf(gu,ui) 
C     read from unit ui, eleminates comments, and copies into buf
C     sets bi (buf index)  to 1
C     sets pointer bp(bi) (buf poiter) to 1 (beginnig of buf) 
C     sets pointer bl (buf length) to end of buf
C     readln returns the iostat of the last read statement
      implicit none 
      include 'mgetx.fi' 
      integer gu,ierr,i,ui,endp,maxl
      logical next,nocom,noquo
      character*500 instr 
c      character*1 CMT1a,CMT1b,CMT2,CMT3,QUOT
      integer  nx_ss
      external nx_ss

      maxl=len(buf)
      next =.true.              ! next line is to be read 
      nocom=.true.              ! no comment mode
      noquo=.true.              ! no quote mode
      bl=0 
      call fillsp(buf)
      do while (next.or.(.not.nocom)) 
        call fillsp(instr)
        next=.false. 
        read(unit=ui,fmt='(A)',iostat=ierr,err=10) instr 

C     get the pointer endp to the end of instr
        i=len(instr) 
        do while ((i.ge.1).and.(instr(i:i).le.' ')) 
          i=i-1 
        end do
        endp=i 
        i=1 
        do while (i .le. min(len(instr),maxl,endp)) 
          if (nocom) then
            if (instr(i:i) .eq. QUOT) then
              if (noquo) then
                noquo=.false.
              else
                noquo=.true.
              end if
            end if
            if (noquo) then
              next=.false.
            else
              next=.true.
            end if
          end if
          
C     if in quot modus, skip all comments
          if (.not.noquo) goto 30

C     eliminates TAB
          if (instr(i:i).eq.char(9)) instr(i:i)=' '
C     ignore all other nonprintable characters 
          if (instr(i:i).lt.' ') goto 20 

          if (nocom) then 
            if (instr(i:i) .eq. CMT1a) then 
              nocom=.false. 
            end if 
          end if 
          
          if (nocom) then 
            if (instr(i:i) .eq. CMT2) then 
              i=len(instr)      ! abbruch
              goto 20
            end if
          end if
          
          if (nocom) then 
            if (instr(i:i) .eq. CMT3) then 
              if ((i.lt.min(len(instr),maxl,endp)).and.
     $             ((instr(i+1:i+1).eq.CMT1a)
     $             .or.(instr(i+1:i+1).eq.CMT2)
     $             .or.(instr(i+1:i+1).eq.CMT3))) then
                i=i+1
              else
                next = .true. 
                i=len(instr)    ! abbruch
                goto 20
              end if
            end if
          end if

 30       continue
          if (nocom) then 
            bl=bl+1 
            buf(bl:bl)=instr(i:i) 
          end if 
          
          if (.not.nocom) then 
            if (instr(i:i) .eq. CMT1b) then 
              nocom=.true. 
            end if 
          end if
 20       continue
          i=i+1 
        end do                  ! i 
      end do                    ! read instr  
      bp=1
C     get the first non blank character
      do while ((buf(bp:bp).eq.' ').and.(bp.le.bl)) 
        bp=bp+1
      end do 
      bpf=bp 

      bpl=0
      bpsp=0
      getbuf = bl-bp+1
      if (rdldbg.gt.0) call writebuf(0)
      return 
 10   getbuf = ierr 
      call writebuf(0)
      return
      end 
      
C     -------------------------------------------------------------------- 
      integer function nx_ss(gu) 
C     delete leading spaces, commas and semicolons
C     set bpf (buf pointer first) to first character
C     set bpl (buf pointer last) to last character 
C     returns -1 if end of line reached
      implicit none 
      include 'mgetx.fi' 
      integer gu,p
      logical noquo
      
      p=bp
C     get the first non blank character
      do while ((buf(p:p).eq.' ').and.(p.le.bl)) 
        p=p+1
      end do 
      bpf=p 
      bp=p

C     scan over the string and stop when a blank of delimiter is reached
C     if a quote character is found, scan to the next quote ignoring all
C     delimiters and spaces
      
      noquo=.false.
      do while  (((buf(p:p).ne.' ')
     &     .and.(buf(p:p).ne.',') 
     &     .and.(buf(p:p).ne.';') 
     &     .and.(p.le.bl))
     $     .or.noquo)
        if (buf(p:p).eq.QUOT) then
          if (noquo) then
            noquo=.false.
          else
            noquo=.true.
          end if
        end if
        p=p+1     
      end do 
      bpl=p-1 
      
C     eliminate the following spaces 
      do while ((buf(p:p).eq.' ').and.(p.le.bl)) 
        p=p+1 
      end do 
      del_ss=0
C     and a following delimiter 
      if ((buf(p:p).eq.',').or.(buf(p:p).eq.';')) then 
        p=p+1
        del_ss=1
      end if
      bp=p
      nx_ss=bpl-bpf+1
      if (bpf.gt.bl) nx_ss=-1
      if (rdldbg.gt.2) then
        if (bpf.le.bpl) 
     $       write(0,'(3I4,3A)') p,bpf,bpl,'>>',buf(bpf:bpl),'<<'
        if ((bpf.gt.bpl).and.(bpf.le.bl)) 
     $       write(0,'(3I4,1A)') p,bpf,bpl,'>><<'
        if (bpf.gt.bl) 
     $       write(0,'(3I4,1A)') p,bpf,bpl,'>>EOF<<'
      end if
      return 
      end 

C     -------------------------------------------------------------------- 
      subroutine setdbg(i) 
      implicit none 
      integer i
      include 'mgetx.fi' 
      rdldbg=i
      return
      end

C     -------------------------------------------------------------------- 
      integer function rd_ss(gu,ss)
C     ReaD one SubString of the input buffer taking quotes into account 
C     and interprets a double quote '' as a single ' (fortran like)
C     returns the length of ss
      implicit none 
      integer gu,i,j,p
      character*(*) ss
      include 'mgetx.fi'

      call fillsp(ss)
      j=0
      p=bpf-1
      do i=bpf, bpl
        p=p+1
        if (buf(p:p).ne.QUOT) then
          j=j+1
          if (len(ss).ge.j) ss(j:j)=buf(p:p)
        else
          if (p.lt.bpl) then
            if (buf(p+1:p+1).eq.QUOT) then
              j=j+1
              p=p+1
              if (len(ss).ge.j) ss(j:j)=buf(p:p)
            end if
          end if
        end if
      end do
      rd_ss=j
      return
      end

C     -------------------------------------------------------------------- 
      integer function len_ss(gu)
C     returns the length of the next substring in the input buffer  
C     analog rd_ss
      implicit none 
      integer gu,i,j,p
      include 'mgetx.fi'
      j=0
      p=bpf-1
      do i=bpf, bpl
        p=p+1
        if (buf(p:p).ne.QUOT) then
          j=j+1
        else
          if (p.lt.bpl) then
            if (buf(p+1:p+1).eq.QUOT) then
              j=j+1
              p=p+1
            end if
          end if
        end if
      end do
      len_ss=j
      return
      end

C     -------------------------------------------------------------------- 
      logical function end_ss(gu)
      implicit none 
      integer gu
      include 'mgetx.fi'
      end_ss=.false.
      if (del_ss.ge.1) end_ss=.true.
      return
      end

C     -------------------------------------------------------------------- 
      logical function nx_end(gu)
      implicit none 
      integer gu
      include 'mgetx.fi'
      if (bp.lt.bl) then
        nx_end=.false.
      else
        nx_end=.true.
      end if
      return
      end

C     -------------------------------------------------------------------- 
      logical function is_end(gu)
      implicit none 
      integer gu
      include 'mgetx.fi'
      if (bpf.le.bl) then
        is_end=.false.
      else
        is_end=.true.
      end if
      return
      end

C     -------------------------------------------------------------------- 
      logical function getend(gu)
      implicit none 
      integer gu
      include 'mgetx.fi'
      if (bpf.le.bl) then
        getend=.false.
      else
        getend=.true.
      end if
      return
      end

C     --------------------------------------------------------------- 
      subroutine fillsp(var) 
      implicit none 
      character*(*) var
      integer  i
      do i=1, len(var)
        var(i:i)=' ' 
      end do 
      return
      end

C     ---------------------------------------------------------------------
      subroutine writebuf(UI)
      implicit none 
      include 'mgetx.fi' 
      integer ui
      if (bp.le.bl) write(UI,'(A)') buf(bp:bl)
      return
      end

C     -------------------------------------------------------------------- 
      block data getxdb   
      implicit none 
      include 'mgetx.fi' 
      data rdldbg /0/
      end 


C ---------------------------------------------------------------------
      logical function isdecno(c)
      implicit none
      character*(*) c
      if (((c(1:1).ge.'0').and.(c(1:1).le.'9')).or.
     $     ((c(2:2).ge.'0').and.(c(2:2).le.'9').and.
     $     ((c(1:1).eq.'-').or.(c(1:1).eq.'+').or.(c(1:1).eq.'.'))).or.
     $     ((c(3:3).ge.'0').and.(c(3:3).le.'9').and.(c(2:2).eq.'.').and.
     $     ((c(1:1).eq.'-').or.(c(1:1).eq.'+')))) then
        isdecno=.true.
      else
        isdecno=.false.
      end if
      return
      end

C ---------------------------------------------------------------------
      logical function isdualno(c)
      implicit none
      character*(*) c
      if (len(c).lt.2) then
        isdualno=.false.
        return
      end if
      if (((c(1:1).eq.'%').and.(c(2:2).eq.'0')).or.
     $    ((c(1:1).eq.'%').and.(c(2:2).eq.'1'))) then 
        isdualno=.true.
      else
        isdualno=.false.
      end if
      return
      end

C ---------------------------------------------------------------------
      logical function ishexno(c)
      implicit none
      character*(*) c
      if (len(c).lt.2) then
        ishexno=.false.
        return
      end if
      if (((c(1:1).eq.'0').and.(c(2:2).eq.'X')).or.
     $    ((c(1:1).eq.'0').and.(c(2:2).eq.'x'))) then 
        ishexno=.true.
      else
        ishexno=.false.
      end if
      return
      end

C ---------------------------------------------------------------------
      integer function isnumber(c)
      implicit none
      character*(*) c
      logical  isdualno, isdecno, ishexno, isvar
      external isdualno, isdecno, ishexno, isvar
      isnumber=0
      if (isdecno(c)) isnumber=1
      if (isdualno(c)) isnumber=2
      if (ishexno(c)) isnumber=3
      if (isvar(c)) isnumber=4
      return
      end

C ---------------------------------------------------------------------
      integer function isoper(c)
      implicit none
      character*(*) c
      if (len(c).lt.1) then
        isoper=0
        return
      end if
      isoper=0
      if ((c(1:1).eq.'+').and.(c(2:2).eq.' ')) isoper=1
      if ((c(1:1).eq.'-').and.(c(2:2).eq.' ')) isoper=2
      return
      end

C ---------------------------------------------------------------------
      logical function isvar(c)
      implicit none
      character*(*) c
      isvar=.false.
      if (c(1:1).eq.'$') isvar=.true.
      return
      end

C ---------------------------------------------------------------------
      integer function dualread(var)
      implicit none
      character*(*) var
      integer i,z,j 
      i=0
      z=1
      do j=len(var),2,-1
        if (var(j:j).eq.'1') i=i+z
        if ((var(j:j).eq.'1').or.(var(j:j).eq.'0')) then
          z=z*2
        else if (var(j:j).ne.' ') then
          write(0,'(A,A)') 'supposed binary: ',var
          stop 'ERROR: binary constant must be of 0 and 1'
        end if
      end do 
      dualread=i
      return
      end

C ---------------------------------------------------------------------
      integer function hexread(var)
      implicit none
      character*(*) var
      integer i,z,j
      logical hexfnd
      integer asc
      external asc
      i=0
      z=1
C     do loop down to 3 to ignore leading '0x'
      do j=len(var),3,-1
        if (var(j:j).ne.' ') then
          hexfnd=.false.
          if ((var(j:j).ge.'0').and.(var(j:j).le.'9')) then
             i=i+z*(asc(var(j:j))-asc('0'))
             hexfnd=.true.
          end if
          if ((var(j:j).ge.'a').and.(var(j:j).le.'f')) then
             i=i+z*(asc(var(j:j))-asc('a')+10)
             hexfnd=.true.
          end if
          if ((var(j:j).ge.'A').and.(var(j:j).le.'F')) then
             i=i+z*(asc(var(j:j))-asc('A')+10)
             hexfnd=.true.
          end if
          z=z*16
          if (.not.hexfnd) then
            write(0,'(A,A)') 'supposed hexadecimal: ',var
            stop 'ERROR: hex digit out of range !'
          end if
        end if
      end do 
      hexread=i
      return
      end

C ---------------------------------------------------------------------
      integer function asc(var)
      implicit none
      character*(*) var
      character*1 c
      integer i
      c=var(1:1)
      asc=ichar(c)
      return
      end

C ---------------------------------------------------------------------
      integer function ss2iv(c,i,vard,varc,ivar)
      implicit none
      character*(*) c
      integer       ivar
      character*(*) varc(ivar)
      real*8        vard(ivar)
      integer isno,i
      integer isnumber,dualread,hexread
      real*8  varread
      external isnumber,dualread,hexread
      external varread

      i=0
      isno=isnumber(c)
      if (isno.eq.1) read(c,'(I40)') i 
      if (isno.eq.2) i=dualread(c)
      if (isno.eq.3) i=hexread(c)
      if (isno.eq.4) i=varread(c,vard,varc,ivar)
      ss2iv=isno
      return
      end

C ---------------------------------------------------------------------
      integer function ss2dv(c,d,vard,varc,ivar)
      implicit none
      character*(*) c
      real*8        d
      integer       ivar
      character*(*) varc(ivar)
      real*8        vard(ivar)
      integer isno
      integer isnumber,dualread,hexread
      real*8  varread
      external isnumber,dualread,hexread
      external varread

      d=0.0
      isno=isnumber(c)
      if (isno.eq.1) read(c,'(F40.0)') d
      if (isno.eq.2) d=dualread(c)
      if (isno.eq.3) d=hexread(c)
      if (isno.eq.4) d=varread(c,vard,varc,ivar)
      ss2dv=isno
      return
      end

C ---------------------------------------------------------------------
      integer function ss2i(c,i)
      implicit none
      character*(*) c
      integer  isno,i
      integer  isnumber,dualread,hexread
      external isnumber,dualread,hexread

      i=0
      isno=isnumber(c)
      if (isno.eq.1) read(c,'(I40)') i 
      if (isno.eq.2) i=dualread(c)
      if (isno.eq.3) i=hexread(c)
      ss2i=isno
      return
      end

C ---------------------------------------------------------------------
      integer function ss2d(c,d)
C     converts a substring into a real number
C     returns the typ of number (isnumber)
C     returns zero if no valid number was found
      implicit none
      character*(*) c
      real*8        d
      integer isno
      integer isnumber,dualread,hexread
      external isnumber,dualread,hexread

      d=0.0
      isno=isnumber(c)
      if (isno.eq.1) read(c,'(F40.0)') d
      if (isno.eq.2) d=dualread(c)
      if (isno.eq.3) d=hexread(c)
      ss2d=isno
      return
      end

C ---------------------------------------------------------------------
      real*8 function varread(c,vard,varc,ivar)
      implicit none
      character*(*) c
      integer       ivar
      character*(*) varc(ivar)
      real*8        vard(ivar)
      integer  i,l
      integer  len_c
      external len_c

      l=len_c(c)
      do i=1, ivar
        if (varc(i)(1:l-1).eq.c(2:l)) then
          varread=vard(i)
          goto 10
        end if
      end do
      write(0,'(3A)') ' ERROR: Variable ',c(1:l),' not defined !'
      stop 
 10   continue
      return
      end

C- -------------------------------------------------------------------- 
      integer function len_c(c)
      implicit none 
      character*(*) c
      integer i
      do i=len(c),1,-1
        if (c(i:i).gt.' ') goto 10
      end do
 10   continue
      len_c=i
      return
      end

C ---------------------------------------------------------------------
      integer function getc(gu,c)
      implicit none
      integer  gu
      character*(*) c
      integer  l
      integer  nx_ss,rd_ss
      external nx_ss,rd_ss

      l=nx_ss(gu)
      if (l.lt.0) goto 10
      l=rd_ss(gu,c)
 10   continue
      getc=l
      return
      end

C ---------------------------------------------------------------------
      integer function getcd(gu,c,d)
      implicit none
      integer  gu
      character*(*) c
      real*8   d,dd
      integer  no,l
      logical  end_ss
      external end_ss
      integer  nx_ss,rd_ss,getd
      external nx_ss,rd_ss,getd

      l=0
      dd=d
      no=getd(gu,dd)
C      write(0,'(A,I3,$)') '  no-a',no
      if (no.le.0) then
        l=nx_ss(gu)
        l=rd_ss(gu,c)
        no=getd(gu,d)
C        write(0,'(A,I3,$)') '  no-b',no
      else
        d=dd
      end if
      getcd=no+10*l
C      write(0,'(A,I3)') '  getcd',no+10*l,no
      return
      end

C ---------------------------------------------------------------------
      integer function getd(gu,d)
C     gets a number (with calculations) out of buffer
C     if no number was found  : returns a zero 
C                               does not advance buffer poiter
C                               sets d to zero
C     returns a value .gt. zero if a number d was read. 
      implicit none
      integer gu
      real*8  d,dd
      include 'mgetx.fi' 
      integer  no,op
      integer  xgetop,xgetd
      logical  end_ss
      external xgetop,xgetd,end_ss

      no=0
      op=xgetop(gu)
      if (end_ss(gu).and.(op.gt.0)) goto 30

 10   continue
      no=xgetd(gu,dd)
      if ((no.le.0).and.(op.gt.0)) goto 30

      if (op.eq.0) d=dd 
      if (op.eq.1) d=d+dd 
      if (op.eq.2) d=d-dd 
      if (end_ss(gu)) goto 20
      
      op=xgetop(gu)
      if (op.gt.0) goto 10
      
 20   continue
      getd=no
      goto 40
 30   continue
      getd=5
      if (op.eq.2) d=-d
 40   continue
      return
      end

C ---------------------------------------------------------------------
      integer function geti(gu,i)
C     gets a number (with calculations) out of buffer
C     if no number was found  : returns a zero 
C                               does not advance buffer poiter
C                               sets d to zero
C     returns a value .gt. zero if a number d was read. 
      implicit none
      integer gu
      integer i,ii
      include 'mgetx.fi' 
      integer  no,op
      integer  xgetop,xgeti
      logical  end_ss
      external xgetop,xgeti,end_ss

      no=0
      op=xgetop(gu)
      if (end_ss(gu).and.(op.gt.0)) goto 30

 10   continue
      no=xgeti(gu,ii)
      if ((no.le.0).and.(op.gt.0)) goto 30

      if (op.eq.0) i=ii 
      if (op.eq.1) i=i+ii 
      if (op.eq.2) i=i-ii 
      if (end_ss(gu)) goto 20
      
      op=xgetop(gu)
      if (op.gt.0) goto 10
      
 20   continue
      geti=no
      goto 40
 30   continue
      geti=5
      if (op.eq.2) i=-i
 40   continue
      return
      end

C ---------------------------------------------------------------------
      integer function xgetop(gu)
      implicit none
      integer  gu,l,op
      character*(40) ss
      integer  m
      integer  nx_ss,s_mark,rd_ss,isoper
      external nx_ss,s_mark,rd_ss,isoper

      op=0
      m=s_mark(gu)
      l=nx_ss(gu)
      l=rd_ss(gu,ss)
c      write(0,'(2A)') 'xgetop:',ss(1:5)
      if (l.le.0) then
        call g_mark(gu,m)
      else
        op=isoper(ss)
        if (op.le.0) call g_mark(gu,m)
      end if
      xgetop=op
      return
      end

C ---------------------------------------------------------------------
      integer function xgetd(gu,d)
      implicit none
      integer  gu,l,no,m
      real*8   d
      character*(40) ss
      integer  nx_ss,s_mark,rd_ss,ss2d
      external nx_ss,s_mark,rd_ss,ss2d
      
      d=0.0
      no=0
      m=s_mark(gu)
      l=nx_ss(gu)
      l=rd_ss(gu,ss)
      if (l.le.0) then
        call g_mark(gu,m)
      else
        no=ss2d(ss,d)
c     write(0,'(2A)') 'xgetd: ',ss(1:5)
        if (no.le.0) call g_mark(gu,m)
      end if
      xgetd=no
      return
      end

C ---------------------------------------------------------------------
      integer function xgeti(gu,i)
      implicit none
      integer  gu,l,no,m,i
      character*(40) ss
      integer  nx_ss,s_mark,rd_ss,ss2i
      external nx_ss,s_mark,rd_ss,ss2i
      
      i=0
      no=0
      m=s_mark(gu)
      l=nx_ss(gu)
      l=rd_ss(gu,ss)
      if (l.le.0) then
        call g_mark(gu,m)
      else
        no=ss2i(ss,i)
        if (no.le.0) call g_mark(gu,m)
      end if
      xgeti=no
      return
      end

C     --------------------------------------------------------------- 
      subroutine pop_mark(gu) 
      implicit none 
      include 'mgetx.fi' 
      integer gu
      if (bpsp.gt.1) then
        bpsp=bpsp-1
        bp=bp_stack(bpsp)
      end if
      return 
      end 
      
C     --------------------------------------------------------------- 
      subroutine push_mark(gu) 
      implicit none 
      include 'mgetx.fi' 
      integer gu
      if (bpsp.ge.dimbpst) stop 'Buffer Pointer Stack exceeded !'
      bpsp=bpsp+1
      bp_stack(bpsp)=bp
      return 
      end 
      
C     --------------------------------------------------------------- 
      integer function s_mark(gu) 
      implicit none 
      include 'mgetx.fi' 
      integer gu
      s_mark=bp
      return 
      end 
      
C     --------------------------------------------------------------- 
      subroutine g_mark(gu,m) 
      implicit none 
      include 'mgetx.fi' 
      integer gu,m
      bp=m
      return 
      end 
      
C     ---------------------------------------------------------------------
      subroutine getln(gu,varc,vard,varx,DIMVAR)
      implicit none
      integer gu,DIMVAR
      character*(*) varc(DIMVAR)
      real*8        vard(DIMVAR)
      integer       varx(DIMVAR)
C     ..local
      integer lasthit,ilen,iq,ihit,qdef,no,l,strt,ovrwr
      integer skip! Herbers2024 - I had to add some conditions to add F1 and spin2 properly.
      real*8  dqn,oldv
      character*40 chrqn,ochrqn
C     ..externals
      logical  is_end
      integer  getd,nx_ss,pr_ss,rd_ss,s_mark
      external is_end
      external getd,nx_ss,pr_ss,rd_ss,s_mark

      call fillsp(ochrqn)
      call fillsp(chrqn)
      qdef=0
C      write(0,*) 'start getln'
      do while (.not.is_end(gu))
C        write(0,*) 'start getln loop'
        strt = s_mark(gu)
        no=getd(gu,dqn)
        if (no.le.0) then
          l=nx_ss(gu)
          l=rd_ss(gu,chrqn)
          ochrqn=chrqn
          strt = s_mark(gu)
        else
          chrqn=ochrqn
        end if
        if (chrqn(1:1).eq.' ') then
          qdef=qdef+1
          if (qdef.gt.DIMVAR) return
          chrqn=varc(qdef)
        end if

        ilen=1
        do while (chrqn(ilen:ilen).ne.' ')
          ilen=ilen+1
        end do
        ilen=ilen-1   ! up until here the start and end of the quantum number string is read. eg. ' F1up ' is read as 'F1up'
C        write(0,*) 'getln:',chrqn(1:ilen),no,ilen
        
        ihit=0 ! now it will be checked if it hits with a string in varc.
        ovrwr=0
 33     continue

        do iq=1,DIMVAR !all variables
         skip=0
         if ((ilen.eq.1).and.(varc(iq)(1:1).eq.'F')        ! condition to be able to handle F1 vs F (so that F does not trigger F1 to be read)
     $                   .and.(varc(iq)(2:2).eq.'1')) then ! condition to be able to handle F1 vs F (so that F does not trigger F1 to be read)
          skip=1
         end if
         if ((ilen.eq.4).and.(varc(iq)(1:4).eq.'spin')     ! condition to be able to handle spin2 vs spin (so that spin doe not trigger spin2 to be read)
     $                   .and.(varc(iq)(5:5).eq.'2')) then ! condition to be able to handle spin2 vs spin (so that spin doe not trigger spin2 to be read)
          skip=1
         end if
        
          if ((chrqn(1:ilen).eq.varc(iq)(1:ilen)).and.(skip.eq.0)) then !check if strings are identical to substrings of varq
C            write(0,*) 'getln:>>',varc(iq),'<<'
            if (varx(iq).eq.ovrwr) then
              ihit=ihit+1
              lasthit=iq
              call g_mark(gu,strt)
              oldv=vard(iq)
              no=getd(gu,vard(iq))
C              write(0,*) no,ihit,vard(iq)
              if (ovrwr.eq.1) write(0,*) 
     $          ' GETLN: overwriting ',chrqn(1:ilen)
     $             ,' from ',oldv,' to ',vard(iq)
              varx(iq)=ihit
            else
            if (varx(iq).gt.1) then
                ihit=ihit+1
                lasthit=iq
                call g_mark(gu,strt)
                no=getd(gu,vard(iq))
C                write(0,*) no,ihit,vard(iq)
                varx(iq)=varx(iq)-1
              end if
            end if
          end if
        end do !for reading F as F1 or F as fup and flo, it will not be identical though.
        if (ihit.eq.1) varx(lasthit)=1
        if ((ovrwr.eq.0).and.(ihit.eq.0)) then
          ovrwr=1
          goto 33
        end if
        if ((ovrwr.eq.1).and.(ihit.eq.0)) then
          write(0,'(2A)')
     $         'Input Error: undef character ',chrqn(1:ilen)
        end if
        if (ihit.eq.1) then
          call fillsp(ochrqn)
        end if
      end do
      return
      end

C     -------------------------------------------------------
      subroutine ind_ss(c,maxi,i1,i2,j2)
      implicit none
C     i1,i2 : index of array
C     j2 : len of substring of c without array information
C     maxi  : maximal arrayelement number

      integer maxi,i1,i2,j2
      character*(*) c
      character*(40) cc

      integer i,l,h1,h2,x1,x2
      l=len(c)
      h1=0
      h2=0
      x1=0
      x2=0
      do i=1, l
        if (c(i:i).eq.'(') then
          x1=i
          h1=h1+1
        end if
        if (c(i:i).eq.')') then
          x2=i
          h2=h2+1
        endif
        if (h2.gt.h1) stop ' index error 2'
        if (c(i:i).eq.' ') goto 10
      end do
 10   continue
      i1=1
      i2=maxi
      j2=i-1
      if (h1.ne.h2) stop ' index error 1'
      if ((h1.eq.0).and.(h2.eq.0)) return
      j2=x1-1
      if (x1.ge.(x2-1)) return

      do i=x1+1,x2-1
        if ((c(i:i).lt.'0').or.(c(i:i).gt.'9')) then
          if ((c(i:i).ne.'-').and.(c(i:i).ne.':'))
     $         stop 'index range - or :'
          goto 20
        end if
      end do
 20   continue
      
      if (i.gt.(x1+1)) then
        call fillsp(cc)
        cc=c(x1+1:i-1)
        read(cc,'(I40)') i1
      end if
      i2=i1
      if (i.eq.x2) goto 30
      i2=maxi
      if ((i+1).eq.x2) goto 30
      call fillsp(cc)
      cc=c(i+1:x2-1)
      read(cc,'(I40)') i2

 30   continue
      if ((i1.lt.1).or.(i2.lt.1)) stop ' index .lt. 1' 
      if ((i1.gt.maxi).or.(i2.gt.maxi)) stop ' index gt maxi' 
      return
      end

C     ---------------------------------------------------------------------
      subroutine getxln(gu,varc,vard,varx,DIMVAR,DIMP)
      implicit none
      integer gu,DIMVAR,DIMP
      character*(*) varc(DIMVAR)
      real*8        vard(DIMVAR,DIMP)
      integer       varx(DIMVAR,DIMP)
C     ..local
      integer lasthit,ilen,iq,ihit,qdef,no,l,strt
      integer id,i1,i2,ovrwr
      real*8  dqn,oldv
      character*40 chrqn,ochrqn,cq1,cq2!adding to characters for checking cq1 ending with space, cq2 ending with underscore _
C     ..externals
      logical  is_end
      integer  getd,nx_ss,pr_ss,rd_ss,s_mark
      external is_end
      external getd,nx_ss,pr_ss,rd_ss,s_mark

      call fillsp(ochrqn)
      call fillsp(chrqn)
      qdef=0
C      write(0,*) 'start getln'
      do while (.not.is_end(gu))
C        write(0,*) 'start getln loop'
        strt = s_mark(gu)
        no=getd(gu,dqn)
        if (no.le.0) then
          l=nx_ss(gu)
          l=rd_ss(gu,chrqn)
          ochrqn=chrqn
          strt = s_mark(gu)
        else
          chrqn=ochrqn
        end if
        if (chrqn(1:1).eq.' ') then
          qdef=qdef+1
          if (qdef.gt.DIMVAR) return
          chrqn=varc(qdef) !the defined parameters in iam.fi
        end if

        call ind_ss(chrqn,DIMP,i1,i2,ilen)

C        write(0,*) 'getln:',chrqn(1:ilen),no,ilen,i1,i2
        ilen=ilen+1 ! Herbers2024, Together with changing character*7 parstr(MAXPAR) 
C                   ! to character*8 parstr(MAXPAR) in iamdata.di, and adding an extra space to all parameter strings
C                   ! This should prevent reading parameters that are not suppose to be read (e.g. Dc3 would reald all parameters starting with Dc3, 
C                   ! which is not supposed to happen
        cq1=chrqn(1:ilen)!The underscore exception to the internal rotation parameters can still be read for all tops at once   
        cq2=chrqn(1:ilen)!Herbers 2024
        cq2(ilen:ilen)='_'    !
        if (cq1(ilen:ilen).eq.'(') then !Herbers 2024 ! the parentheses exception, to allow for different rotational constants.
        ilen=ilen-1
        cq1=chrqn(1:ilen)
        cq2=chrqn(1:ilen)
        end if        
        do id=i1,i2
          ihit=0
          ovrwr=0
 33       continue
          do iq=1,DIMVAR
            
            if ((cq1.eq.varc(iq)(1:ilen)).or.
     $          (cq2.eq.varc(iq)(1:ilen))) then ! Herbers 2024
C              write(0,*) 'getln:>>',varc(iq),'<<'
              if (varx(iq,id).eq.ovrwr) then
                ihit=ihit+1
                lasthit=iq
                call g_mark(gu,strt)
                oldv=vard(iq,id)
                no=getd(gu,vard(iq,id))
C                write(0,*) no,ihit,vard(iq,id)
                varx(iq,id)=ihit
                if (ovrwr.eq.1) write(0,*)
     $               'GETXLN: overwriting ',chrqn(1:ilen)
     $               ,' from ',oldv,' to ',vard(iq,id)
              else
                if (varx(iq,id).gt.1) then
                  ihit=ihit+1
                  lasthit=iq
                  call g_mark(gu,strt)
                  no=getd(gu,vard(iq,id))
C                  write(0,*) no,ihit,vard(iq,id)
                  varx(iq,id)=varx(iq,id)-1
                end if
              end if
            end if
          end do
          if (ihit.eq.1) varx(lasthit,id)=1
          if ((ovrwr.eq.0).and.(ihit.eq.0)) then
            ovrwr=1
            goto 33
          end if
          if ((ovrwr.eq.1).and.(ihit.eq.0)) then
            write(0,'(2A)')
     $           'Input Error: undef character ',chrqn(1:ilen)
          end if
          if (ihit.eq.1) then
            call fillsp(ochrqn)
          end if
        end do
      end do
      return
      end

