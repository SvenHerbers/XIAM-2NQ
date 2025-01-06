C ----------------------------------------------------------
      subroutine parinp(a,palc,pali,ifit,dfit,npar,nfit)

      implicit none
      include 'iam.fi'
      real*8  a(DIMPAR,DIMVB)
      real*8  palc(DIMFIT,-1:DIMPLC)
      integer pali(DIMFIT, 0:DIMPLC,2)
      integer ifit(DIMPAR,DIMVB), dfit(DIMFIT),npar, nfit

      integer i,j,iq,iv,ib,ndata,nsplit,ift,apri,xpri,sizeb,ffree
      real*8  indeg,inkj,pi
      real*8  msplit,efdata,m2split

C     local variables: header
      character*10 helpstr

C     local variables: reading list of control variables 
      real*8  dctl(DIMCPAR+DIMCINT)   
      integer cdone(DIMCPAR+DIMCINT)
      integer adj

C     local variables: reading list of parameters 
      integer padone(DIMPAR,DIMVB)

C     local variables: reading definition of gammas 
      real*8 dgam(2*DIMTOP)
      integer gdone(2*DIMTOP),is,it,itop
      character*10 symstr

C     local variables: reading linear comb. of  parameters to fit
      real*8 afit(DIMPAR,DIMVB)
      real*8 sum
      integer ix,gl,gm,oldix,df,no
      character*10 fitstr,spaz

C     local variables: reading list of transitions 
      real*8  dqno(DIMQC)
      integer qdone(DIMQC),ilr,ila,ilx,iqq,il
      character*40 fmtstr
      logical stpflg

      integer  getd,geti,getc,getbuf,myand,myor,s_mark,len_c
      logical getend
      external geti,getd,getbuf,getc,myand,myor,s_mark,len_c
      external getend
C
      include 'iamdata.fi'

      pi=dacos(-1.0d0)
      inkj=3.9903132D-04
      indeg=180.0d0/pi

      call fillsp(spaz)

      do i=1, DIMFIT
        do j=1, DIMPLC
          palc(i,j)=0.0d0
          pali(i,j,1)=0
          pali(i,j,2)=0
        end do
        dfit(i)=0
      end do

      do i=1, DIMPAR
        do j=1, DIMVB
          ifit(i,j)=0
        end do
      end do

C     Header of Input File
      do while (.true.)
        if (getbuf(gu,ui).le.0) goto 4
        call writebuf(6)
        call fillsp(helpstr)
        i=getc(gu,helpstr)
        if ((helpstr.eq."help").or.
     $      (helpstr.eq."Help").or.
     $      (helpstr.eq."HELP")) then
          write(*,'(A)') ' possible parameters: '
          do i=1,DIMPAR
            write(*,'(X,A)') parstr(i)
          end do
          write(*,'(/,A,I3)') ' DIMTOP ',DIMTOP
          write(*,'(A,I3)') ' DIMJ   ',DIMJ
          write(*,'(/,A)')
     $         ' See help file "xiam_v25.txt" for more information'
      stop 
        end if
      end do
 4    continue

      do i=1, DIMCINT+DIMCPAR
        cdone(i)=0
        dctl(i)=0.0d0
      end do
      do while (.true.)
        if (getbuf(gu,ui).le.0) goto 6
        call getln(gu,ctlstr,dctl,cdone,DIMCINT+DIMCPAR)               
      end do
 6    continue
      do i=1, DIMCINT
        ctlint(i)=0
      end do
      do i=DIMCINT+1, DIMCINT+DIMCPAR
        ctlpar(i)=0.0d0
      end do
      ctlint(C_MAXM)=8
      ctlint(C_NZYK )=1
      ctlint(C_PRI)=0 
      ctlint(C_XPR)=0
      ctlint(C_PRINT)=3
      ctlint(C_NFOLD)=3 
      ctlint(C_FITSC)=0
      ctlint(C_WOODS)=33
      ctlint(C_ADJF)=0
      ctlpar(C_EPS)=1.0d-12
      ctlint(C_DW)=3 !Herbers2024 default use exact NQC in fitting
C     ctlint(C_NQI)=0 !Herbers2024 default, printing of approximate hyperfine intensities is off.
      ctlint(C_SORT)=3
      ctlpar(C_DEFER)=1.0d-5
      ctlpar(C_CNVG)=0.999d0
      ctlpar(C_LMBDA)=0.00001d0
      ctlpar(C_FRQLO)=6.0d0
      ctlpar(C_QSUM)=1000. !Herbers2023 Default value for qsum is 1000 if not specified.  
      ctlpar(C_FRQUP)=40.0d0
      ctlpar(C_INTLM)=0.1d0
      ctlpar(C_TEMP)=273.0d0
      

      do i=1, DIMCINT
        if (cdone(i).ne.0) then
          ctlint(i)=int(dctl(i))
        end if  
      end do
      do i= DIMCINT+1, DIMCINT+DIMCPAR
        if (cdone(i).ne.0) then
          ctlpar(i)=dctl(i)
        end if  
      end do
      apri=0
      xpri=0
      if (ctlint(C_PRINT).eq.1) apri=1 
      if (ctlint(C_PRINT).eq.2) apri=1 + AP_TL
      if (ctlint(C_PRINT).ge.2) xpri=XP_CC
      if (ctlint(C_PRINT).eq.3) apri=1 + AP_TF
      if (ctlint(C_PRINT).eq.4) apri=2 + AP_TF
      if (ctlint(C_PRINT).eq.5) apri=2 + AP_TF + AP_TL
      if (ctlint(C_PRINT).eq.6) apri=2 + AP_TL + AP_LT
      ctlint(C_PRI)=myor(ctlint(C_PRI),apri)
      ctlint(C_XPR)=myor(ctlint(C_XPR),xpri)
      ctlint(C_NZYK)=max(ctlint(C_NZYK),ctlint(C_NCYCL))
      do i=1, DIMCINT
        if (mod(i,4).eq.1) write(*,*)
        if (ctlstr(i)(1:1).ne.'_')
     $       write(*,'(3X,A,3X,I5,$)') ctlstr(i),ctlint(i)
      end do
C      write(*,*)
      j=0
      do i= DIMCINT+1, DIMCINT+DIMCPAR
        j=j+1
        if (mod(j,3).eq.1) write(*,*)
        if (len(ctlstr(i)).ne.0) write(*,'(3X,A,3X,D13.6,$)')
     $       ctlstr(i),ctlpar(i)
      end do
      write(*,*)
      npar=DIMPRR+ctlint(C_NTOP)*DIMPIR

      if (ctlint(C_NTOP).eq.3) then
        if (myand(ctlint(C_ADJF),4).ne.0) write(*,*)
     $       ' I recommend to set adj to 4 for a three top molecule!'
      end if
      if (ctlint(C_NTOP).gt.DIMTOP) stop 'ERROR: ntop > DIMTOP'
      do itop=1, DIMTOP
        if ((2*ctlint(C_MAXM)+1).gt.DIMM) stop 'ERROR: 2maxm+1 > DIMM'
      end do
      if (ctlint(C_EVAL).ne.0) then
        write(*,'(A)') '\\ writing eigenvalues in file eval.out ' 
        open(20,file='eval.out',status='unknown')
      end if
      if (ctlint(C_DFRQ).ne.0) then
        write(*,'(A,A)') '\\ writing deviation of frequencies in',
     $       ' file dfreq.out'
        open(21,file='dfreq.out',status='unknown')
      end if

C     set the parameter name acoording to the reduction A or S
      if (ctlint(C_RED).eq.1) then
c        parstr(P_DKD)='R6      '
        write(*,*) 'Using Watson S Reduction '
      end if 
      if (ctlint(C_RED).eq.2) then
        parstr(P_DKD)='R6      '
        write(*,*) 'Using van Eijck-Typke Reduction '
      end if 
      if (ctlint(C_RED).eq.0) then
        write(*,*) 'Using Watson A Reduction '
      end if

      do itop=1, ctlint(C_NTOP)
        size(S_MAXM+itop)=ctlint(C_MAXM)
        size(S_MAXV+itop)=ctlint(C_MAXV)
      end do
        
C --------------
      write(*,*)
      do j=1, DIMVB
        do i=1, DIMPAR
          padone(i,j)=0
          a(i,j)=0.0d0
        end do
      end do
      do i=1,DIMPAR+1
        if (getbuf(gu,ui).lt.0) stop 'reading a'
        if (getend(gu)) goto 5
        call getxln(gu,parstr,a,padone,DIMPAR,DIMVB)
      end do
 5    continue

C     fill parameters of not used tops with zeros
      do ib=1, DIMVB
        do i=DIMPRR+DIMPIR*ctlint(C_NTOP)+1,DIMPAR
          a(i,ib)=0.0d0
          ifit(i,ib)=0
        end do
        if (ctlint(C_NTOP).le.1) then 
          a(P_FF,ib)=0.0d0
          a(P_VSS,ib)=0.0d0
          a(P_VCC,ib)=0.0d0
        end if
      end do

C     get a preliminary value of sizeb
      sizeb=1
      do i=1, DIMPRR+DIMPIR*ctlint(C_NTOP)
        do ib= DIMVB, 2, -1
          if ((a(i,ib).ne.a(i,(ib-1))).and.(a(i,ib).ne.0.0d0))
     $         sizeb=max(sizeb,ib)
        end do
      end do
      size(S_NB)=sizeb
      write(*,*) 'assumed sizeb',sizeb

      do ib=1, size(S_NB)
        adj=ctlint(C_ADJF)
        if ((ctlint(C_INTS).ne.0) .and. (a(P_MUX,ib).eq.0.0d0)
     $       .and. (a(P_MUY,ib).eq.0.0d0).and. (a(P_MUZ,ib).eq.0.0d0))
     $       stop 'ERROR: need mu_x mu_y or mu_z for intensities'
        if (a(P_BJ,ib).eq.0.0d0) stop 'ERROR: BJ can not be zero'
        if (ctlint(C_NTOP).ge.1) then
C          if (a(P_FF,j).eq.a(P1_F,j)) a(P_FF,j)=0.0d0
C          if (a(P1_F0,j).eq.a(P1_F,j)) a(P1_F0,j)=0.0d0
          if (padone(P_FF,ib).ge.DIMTOP)  a(P_FF,ib)=0.0d0
          if (padone(P1_F0,ib).ge.2) a(P1_F0,ib)=0.0d0
          if (padone(P2_F0,ib).ge.3) a(P2_F0,ib)=0.0d0
          if (((a(P1_ANGZ,ib).ne.0.0d0).or.(a(P1_ANGX,ib).ne.0.0d0))
     $         .and.
     $         (a(P1_BETA,ib).eq.0.0d0).and.(a(P1_GAMA,ib).eq.0.0d0)
     $         .and.(myand(adj,16).eq.0)) then 
            adj=myor(adj,16)
            write(*,'(A)') ' \\ set (adj  or 16)'
          end if
          if ((a(P1_F0,ib).ne.0.0d0).and.(a(P1_RHO,ib).eq.0.0d0)
     $         .and.(myand(adj,8).eq.0)) then 
            adj=myor(adj,8)
            write(*,'(A)') ' \\ set (adj  or  8)'
          end if
          if ((a(P1_F,ib).eq.0.0d0).and.(myand(adj,1).eq.0)) then 
            adj=myor(adj,1)
            write(*,'(A)') ' \\ set (adj  or  1)'
          end if
        end if
        if (ctlint(C_NTOP).ge.2) then
          if ((a(P_FF,ib).eq.0.0d0).and.(myand(adj,2).eq.0)) then 
            adj=myor(adj,2)
            write(*,'(A)') ' \\ set (adj  or  2)'
          end if
        end if

c        if ((a(P1_F0).eq.0.0d0).and.(a(P1_F).eq.0.0d0)
c     $       .and.(myand(ctlint(C_ADJF),1).eq.0)) then 
c          ctlint(C_ADJF)=myor(ctlint(C_ADJF),1)
c          write(*,'(A)') ' \\ set (adj  or  1)'
c        end if
        ctlnb(CB_ADJ,ib)=adj
        ctlnb(CB_WDS,ib)=ctlint(C_WOODS)
        call adjusta(a(1,ib),npar,adj)
        if (myand(adj,1).ne.0) write(*,*)
     $     '\\ adj  1: adjust F according to rho, beta and gamma'
        if (myand(adj,2).ne.0) write(*,*)
     $     '\\ adj  2: adjust F12 according to rho, beta and gamma'
        if (myand(adj,4).ne.0) write(*,*)
     $     '\\ adj  4: adjust F (one top case) and ignore F'' '
        if (myand(adj,8).ne.0) write(*,*)
     $     '\\ adj  8: adjust rho according to F0 = 1/(2 I_alpha)'
        if (myand(adj,16).ne.0) write(*,*)
     $     '\\ adj 16: adjust beta and gamma  according delta + epsil'
        write(*,'(A,I4)') ' new adj :',adj
      end do
      do ib=sizeb+1, DIMVB
        ctlnb(CB_ADJ,ib)=adj
        ctlnb(CB_WDS,ib)=ctlint(C_WOODS)
      end do

      call pra(a,a,ifit,3,0,0)
      write(*,*)
C ------
      do i=1,DIMFIT+1
        do is=1, DIMPAR
          do ib=1, DIMVB
            padone(is,ib)=0
            afit(is,ib)=0.0d0
          end do
        end do
        if (getbuf(gu,ui).lt.0) stop ' Error reading fit parameters '
        if (getend(gu)) goto 8
        dfit(i)=0
        df=0
        gl=getc(gu,fitstr)
        if (fitstr.eq.'dqx') df=-2 
        if (fitstr.eq.'dqu') df=-1
        if (fitstr.eq.'fit') df=1
        if ((i.gt.DIMFIT).and.(df.ne.0))
     $       stop ' maximum number (DIMFIT) of fit variables exceeded !'
        dfit(i)=df
        if (dfit(i).eq.0) goto 8

C     palc(i,0): stepwidth in % for differential quotient calc
        no=getd(gu,palc(i,0))
        if (no.le.0) palc(i,0)=0.1d0

C     palc(i,-1): factor to scale the new parameters in fit 
        no=getd(gu,palc(i,-1))
        if (no.le.0) palc(i,-1)=1.0d0 

        if ((dfit(i).gt.0))
     $       write(*,'(A,2D8.1,3X,$)') ' fit ',palc(i,0),palc(i,-1)
        if ((dfit(i).eq.-1))
     $       write(*,'(A,2D8.1,3X,$)') ' dqu ',palc(i,0),palc(i,-1) 
        if ((dfit(i).eq.-2))
     $       write(*,'(A,2D8.1,3X,$)') ' dqx ',palc(i,0),palc(i,-1) 

        call getxln(gu,parstr,afit,padone,DIMPAR,DIMVB)
        sum=0.0d0
        oldix=0
        ix=0
        do is=1, DIMPAR
          oldix=ix
          do ib=1, DIMVB
            if (padone(is,ib).eq.1) then 
              ifit(is,ib)=1
              if (afit(is,ib).eq.0.0d0) afit(is,ib)=1.0d0
              ix=ix+1
              if (ix.gt.DIMPLC) stop ' Dimension Error: DIMPLC'
              pali(i,ix,1)=is
              pali(i,ix,2)=ib
              palc(i,ix)=afit(is,ib)          
              sum=sum+dabs(a(is,ib))
            end if
          end do
          if ((ix-oldix).eq.DIMVB) then
            write(*,'(3X,1A,3X,1F6.2,$)')
     $           parstr(is),palc(i,ix)
          else
            do j=oldix+1, ix
              write(*,'(3X,2A,I1,2A,1F5.2,$)')
     $             parstr(is)(1:len_c(parstr(is)))
     $             ,'(',pali(i,j,2),')'
     $             ,spaz(1:(len(parstr(is))-len_c(parstr(is))+1))
     $             ,palc(i,j)
            end do
          end if
        end do
        write(*,*)
        pali(i,0,1)=ix
        do j=1, pali(i,0,1)
          if ((parfit(pali(i,j,1)).ne.0).and.(dfit(i).eq.1)) then
c           write(*,'(3A)')' Warning: fit changed to dqu!',
c    $           ' No analytic derivatives for ',parstr(j)
            dfit(i)=-1
          end if
        end do
        
       if (((ctlint(C_DW).eq.3).and.          ! Herbers2024 ! no analytical gradients for double well or quadrupole treatment $     
     $     (ctlint(C_SPIN).ne.0)).or.(ctlint(C_DW).eq.1)) then                            ! Herbers2024
C      write(0,*) dfit                                              ! Herbers2024
C      do i=1, DIMFIT                                               ! Herbers2024
         if (dfit(i).eq.1) then                                     ! Herbers2024
          write(0,*) 'Disabling analytical gradients for exact ',      ! Herbers2024
     $             'quadrupole treatment or DW (not ',                    ! Herbers2024
     $             'implemented). "fit" has same meaning as "dqu".' ! Herbers2024
          write(*,*) 'Disabling analytical gradients for exact ',      ! Herbers2024
     $             'quadrupole treatment (not ',                    ! Herbers2024
     $             'implemented). "fit" has same meaning as "dqu".' ! Herbers2024
          dfit(i) = -1       ! 1 means analytical gradient       ! Herbers2024
         end if                 ! -1 means dqu                      ! Herbers2024
C     end do                    ! -2 means dqx                      ! Herbers2024
      endif                                                         ! Herbers2024
        
        
        if ((dfit(i).eq.-1).and.(sum.eq.0.0d0)) then
             write(*,*) 'ERROR: dqu/x parameter can not be zero !' 
             stop 'ERROR: dqu/x parameter can not be zero !' 
        end if
      end do
 8    continue
      nfit=i-1
      write(*,*)
      do ib=1, DIMVB
        adj=ctlnb(CB_ADJ,ib)
        if (myand(adj,2).gt.0) then
          if (ifit(P_FF,ib).ne.0)
     $         stop ' Fit ERROR: can not fit F12: adjf = 2'
        end if
        do itop=1, ctlint(C_NTOP)
          ift=(itop-1)*DIMPIR
          if (myand(adj,1).gt.0) then
            if (ifit(P1_F +ift,ib).ne.0)
     $           stop ' Fit ERROR: can not fit F: adjf = 1'
          else
            if (ifit(P1_F0+ift,ib).ne.0)
     $           stop ' Fit ERROR: can not fit F0: adjf <> 1'
          end if
          if (myand(adj,8).gt.0) then
            if (ifit(P1_RHO+ift,ib).ne.0)
     $           stop ' Fit ERROR: can not fit rho: adjf = 8'
          else
            if (ifit(P1_F0+ift,ib).ne.0)
     $         stop ' Fit ERROR: can not fit F0: adjf <> 8'
          end if
          if (myand(adj,16).gt.0) then
            if ((ifit(P1_GAMA+ift,ib).ne.0)
     $           .or.(ifit(P1_BETA+ift,ib).ne.0))
     $           stop ' Fit ERROR: can not fit beta/gamma: adjf = 16'
          else
            if ((ifit(P1_ANGZ+ift,ib).ne.0)
     $           .or.(ifit(P1_ANGX+ift,ib).ne.0))
     $           stop ' Fit ERROR: can not fit delta/epsil: adjf <> 16'
          end if
        end do                  ! itop
      end do                    ! ib
      
C ------
      size(S_G)=0
      do i=1, DIMGAM 
        do itop=1, DIMTOP
          gamma(i,itop)=NaQN
        end do
        gamma(i,0)=0
      end do
      do i=1, 2*DIMTOP
        dgam(i)=0
      end do
      if (ctlint(C_NTOP).ne.0) then
        do i=1, DIMGAM+1 
          do itop=1, 2*DIMTOP
            gdone(itop)=0
          end do
          if (getbuf(gu,ui).lt.0) stop ' Error: reading gamma' 
          if (getend(gu)) goto 7
          if (i.gt.DIMGAM) stop ' to many lines reading S! '

C     symmetry-species-name beginning with '/' 
          gm=s_mark(gu)
          gl=getc(gu,symstr)
          call fillsp(qnostr(MAXQC+i))
          if (symstr(1:1).eq.'/') then
            qnostr(MAXQC+i)=symstr(1:len(qnostr(MAXQC+i))+1)
          else
            call g_mark(gu,gm)
          end if

          call getln(gu,gamstr,dgam,gdone,2*DIMTOP)
          do itop=1, ctlint(C_NTOP)
            if ((int(dgam(2*itop-1)).eq.0).and.(int(dgam(2*itop)).eq.0))
     $           gamma(i,itop)=0
            if ((int(dgam(2*itop-1)).ne.0).and.(int(dgam(2*itop)).eq.0))
     $           gamma(i,itop)=int(dgam(2*itop-1))
            if ((int(dgam(2*itop-1)).eq.0).and.(int(dgam(2*itop)).ne.0))
     $           gamma(i,itop)=int(dgam(2*itop))
            if ((int(dgam(2*itop-1)).ne.0).and.(int(dgam(2*itop)).ne.0))
     $           stop 'Error: use G or S as Keyword' 
          end do    
          size(S_G)=i           
        end do
 7      continue
        do is=1, size(S_G)
          if (qnostr(MAXQC+is)(1:1).ne.' ') 
     $         write(*,'(3X,A,$)') qnostr(MAXQC+is) 
          write(*,'(2X,A,4I4)') 'S  ',(gamma(is,it),it=1,ctlint(C_NTOP))
          do it=1,ctlint(C_NTOP)
            if (gamma(is,it).gt.DIMSIG) stop 'ERROR: sigma > DIMSIG'
          end do
        end do
      else
        gamma(1,1)=-999
      end if

C ----
      do i=1, DIMVV
        do ib=1, DIMVB
          do itop=1, DIMTOP
            qvv(i,itop,ib)=-1
          end do
        end do
      end do
      size(S_NB)=1
      if (ctlint(C_NTOP).ne.0) then
        write(*,*)
        do ib=1, DIMVB+1
          if (getbuf(gu,ui).lt.0) stop ' Error: reading qvv' 
          if (getend(gu)) goto 20
          do iv=1, DIMVV
            call fillsp(fitstr)
            gl=getc(gu,fitstr)
            if (fitstr(1:1).ne.'V') goto 19
            do itop=1, ctlint(C_NTOP)
              no=geti(gu,j)
              if (no.le.0) stop ' Error: V no. for each top necessary!'  
              qvv(iv,itop,ib)=j+1
            end do
          end do
 19       continue
          size(S_NB)=ib
        end do
 20     continue
        do ib=size(S_NB)+1, DIMVB
          do i=1, DIMVV
            do itop=1, DIMTOP
              qvv(i,itop,ib)=qvv(i,itop,size(S_NB))
            end do
          end do
        end do
          
        do ib=1, size(S_NB)
          do iv=1, DIMVV
            if (qvv(iv,1,ib).eq.-1) goto 21
            write(*,'(A,$)') '  V'
            do itop=1,ctlint(C_NTOP) 
              write(*,'(I3,$)') qvv(iv,itop,ib)-1  
            end do
          end do
 21       continue
          write(*,*)
        end do
      end if

C ----
      write(*,*)
      if (ctlint(C_NDATA).le.0) ctlint(C_NDATA)=DIMLIN
      do iq=1, DIMQC
        dqno(iq)=0.0d0
      end do 
C     initial values
      dqno( 8)=NOFIT ! Err
      dqno(15)=1     ! Bup
      dqno(16)=1     ! Blo
      dqno(17)=-1    ! Fup
      dqno(18)=-1    ! Flo
      dqno(31)=-1    ! F1up Herbers2024
      dqno(32)=-1    ! F1lo Herbers2024
      dqno(19)=0     ! Tup
      dqno(20)=0     ! Tlo
      dqno(11)=1     ! V1up
      dqno(12)=1     ! V1lo
      ndata   =0
      msplit  =0.0d0
      m2split =0.0d0
      nsplit  =0
      efdata  =0.0d0
      stpflg  =.false.
      do il=1,ctlint(C_NDATA)
        if (getbuf(gu,ui).lt.0) goto 10
        if (getend(gu)) goto 10
        dqno(23)=0  ! # 
        dqno(24)=0  ! &
        dqno(27)=0  ! diff
        do iq=1, DIMQC
          qdone(iq)=0
        end do
        call getln(gu,qnostr,dqno,qdone,DIMQC) ! READING THE LINE FREQUENCY
        
          if ((ctlint(C_SPIN2).eq.0)) then !This block is to read F as F1 in the single nucleus case, if it is input instead of F1.
           if (dqno(17).ne.-1) then
             dqno(31)=dqno(17)
             dqno(17)=-1
            end if
           if (dqno(18).ne.-1) then
             dqno(32)=dqno(18)
             dqno(18)=-1
           end if
          end if
        
        
        
        
        if ((qdone(25).ne.0).and.(qdone(21).ne.0)) !simplest solution find a different name... but which one? I will use nF1
     $       write(*,*) 'UP Tau overrides K',il
        if ((qdone(26).ne.0).and.(qdone(22).ne.0))
     $       write(*,*) 'LO Tau overrides K',il
        if (((qdone(2).ne.0).or.(qdone(3).ne.0)).and.(qdone(21).ne.0))
     $       write(*,*) 'UP Tau (t) overrides K- K+',il
        if (((qdone(5).ne.0).or.(qdone(6).ne.0)).and.(qdone(22).ne.0))
     $       write(*,*) 'LO Tau (t) overrides K- K+',il
        if ((int(dqno(1)).gt.DIMJ).or.(int(dqno(4)).gt.DIMJ)) then
          write(*,'(A,I4,A)')
     $         ' Warning Max. J exceeded. Line',il,' error'
          stpflg=.true.
        end if
        if (qdone(19).ne.0) then  ! Tup
          dqno(25)=0.0d0          ! Kup
          dqno(21)=0.0d0          ! tup
          qdone(25)=0           
          dqno( 2)=0.0d0          ! K- up
          dqno( 3)=0.0d0          ! K+ up
        end if
        if (qdone(20).ne.0) then  ! Tlo 
          dqno(26)=0.0d0          ! Klo
          dqno(22)=0.0d0          ! tlo
          qdone(26)=0
          dqno( 5)=0.0d0          ! K- lo
          dqno( 6)=0.0d0          ! K+ lo
        end if
        if (qdone(21).ne.0) then  ! tup
          dqno(25)=0.0d0          ! Kup
          dqno(19)=0.0d0          ! Tup
          qdone(25)=0           
          dqno( 2)=0.0d0          ! K- up
          dqno( 3)=0.0d0          ! K+ up
        end if
        if (qdone(22).ne.0) then  ! tlo 
          dqno(20)=0.0d0          ! Tlo
          dqno(26)=0.0d0          ! Klo
          qdone(26)=0
          dqno( 5)=0.0d0          ! K- lo
          dqno( 6)=0.0d0          ! K+ lo
        end if
        if (qdone(25).ne.0) then
          dqno( 2)=-10.0d0          ! K- up
          dqno( 3)=-10.0d0          ! K+ up
          dqno(21)=0.0d0            ! tup,Q_TJ
          qdone(21)=0
        end if
        if (qdone(26).ne.0) then 
          dqno( 5)=-20.0d0          ! K- lo
          dqno( 6)=-20.0d0          ! K+ lo
          dqno(22)=0.0d0            ! tlo,Q_TJ
          qdone(22)=0
        end if
C     copy the vector dqno to the transition list qlin
        do iq=1, DIMQC  
C          write(0,*) int(dqno(iq))
          if ((q_q(iq).gt.0).and.(qul(iq).gt.0)) then
            qlin(il,q_q(iq),qul(iq))=int(dqno(iq))
          end if
        end do
        if (qdone(7).ne.0) then   ! = (frequency)
          dln(il,LN_FREQ)=dqno(7)
          dln(il,LN_ERR)=ctlpar(C_DEFER)
          if (qdone(8).ne.0) then ! Err
            dln(il,LN_ERR)=dqno(8)
          end if
          ndata=ndata+1
        else   
          dln(il,LN_FREQ)=0.0d0
          dln(il,LN_ERR)=NOFIT
        end if
        if ((qdone(2).ne.0).or.(qdone(3).ne.0)) then 
          qlin(il,Q_TJ,Q_UP)=0
        end if
        if ((qdone(5).ne.0).or.(qdone(6).ne.0)) then 
          qlin(il,Q_TJ,Q_LO)=0
        end if

C    Check for symmetry-labels
        do i=1,DIMGAM
          if (qdone(MAXQC+i).ne.0) then
            if ((qdone(9).eq.0).and.(qdone(10).eq.0)) then
              qlin(il,Q_S,Q_UP)=i
              qlin(il,Q_S,Q_LO)=i
            else
              write(*,*)' Warning: S and Symmetry Label given (using S)'
            end if
          end if
        end do
        if ((qlin(il,Q_S,Q_UP).gt.size(S_G)).or.
     $       (qlin(il,Q_S,Q_LO).gt.size(S_G))) then
          write(*,*)'ERROR: Symmetry spezies not defined in line',il
          stpflg=.true.
        end if

C    Check for units (GHz default, MHZ and cm)
        if ((qdone(29).ne.0).and.(qdone(30).eq.0)) then
          dln(il,LN_FREQ)=dln(il,LN_FREQ)/1000.0d0
          if (qdone(8).ne.0) 
     $         dln(il,LN_ERR)=dln(il,LN_ERR)/1000.0d0
        end if
        if ((qdone(29).eq.0).and.(qdone(30).ne.0)) then
          dln(il,LN_FREQ)=dln(il,LN_FREQ)*29.9792458d0
          if (qdone(8).ne.0) 
     $         dln(il,LN_ERR)=dln(il,LN_ERR)*29.9792458d0
        end if
 
C     # references (23): ref lines
        ilr=int(dqno(23))
        ilr=ilr+il
        qlin(il,Q_REF,Q_UP)=ilr
        qlin(il,Q_REF,Q_LO)=ilr

C     splitting frequency as input (27)
        if (qdone(27).ne.0) then ! diff
          ilx=int(dqno(27))
          if (ilx.eq.0) then
            ilx=ilr
          else
            ilx=ilx+il
          end if
          if ((dln(il,LN_FREQ).ne.0.0d0).and.
     $         (dln(ilx,LN_FREQ).ne.0.0d0)) then
            dln(il,LN_FREQ)=dln(ilx,LN_FREQ)+dln(il,LN_FREQ)
          else
            write(0,'(A,2I4)') 'WARNING: diff frequency zero at',il,ilx
          end if
        end if

        if ((ilr.ne.il).and.(dln(il,LN_ERR).ne.NOFIT)) then
          msplit=msplit+abs(dln(il,LN_FREQ)-dln(ilr,LN_FREQ))
          m2split=m2split+(dln(il,LN_FREQ)-dln(ilr,LN_FREQ))**2
          nsplit=nsplit+1
        endif

        efdata=efdata+ctlpar(C_DEFER)**2/dln(il,LN_ERR)**2

C     Check for an average line (dqno(24))
        ila=int(dqno(24))
        ila=ila+il
        qlin(il,Q_AVG,Q_UP)=ila
        qlin(il,Q_AVG,Q_LO)=ila

        if (ila.ne.il) then
          if (dln(ila,LN_ERR).ne.NOFIT) then 
            write(*,'(A,I4,A,I4,A)')
     $           'Warning: average freq. in',il,': line',ila,' not used'
            ilr=qlin(ila,Q_REF,Q_UP)
            if (ilr.ne.ila) then
              msplit=msplit-abs(dln(ila,LN_FREQ)-dln(ilr,LN_FREQ))
              m2split=m2split-(dln(ila,LN_FREQ)-dln(ilr,LN_FREQ))**2
              nsplit=nsplit-1
            endif
            efdata=efdata-ctlpar(C_DEFER)**2/dln(ila,LN_ERR)**2
            dln(ila,LN_ERR)=NOFIT
            efdata=efdata+ctlpar(C_DEFER)**2/dln(ila,LN_ERR)**2
          end if
        end if

C     calc tau (Q_TJ)(dqno(21/22) from K- dqno(2/5) and K+ dqno(3/6)    
        if (((dqno(21).eq.0).and.(dqno(19).eq.0)).and.
     $       ((dqno(2).ge.0.0d0).or.(dqno(3).ge.0.0d0)))
     $       qlin(il,Q_TJ,Q_UP)=int(dqno(2))-int(dqno(3))+1
     $       +int(dqno(1))
        if (((dqno(22).eq.0).and.(dqno(20).eq.0)).and.
     $       ((dqno(5).ge.0.0d0).or.(dqno(6).ge.0.0d0))) 
     $       qlin(il,Q_TJ,Q_LO)=int(dqno(5))-int(dqno(6))+1
     $       +int(dqno(4))

C Q_TJ numbering of eigenvalues of one V
        if (qlin(il,Q_TJ,Q_UP).ne.0) then 
          qlin(il,Q_STAT,Q_UP)=myor(qlin(il,Q_STAT,Q_UP),2)
        end if
        if (qlin(il,Q_TJ,Q_LO).ne.0) then
          qlin(il,Q_STAT,Q_LO)=myor(qlin(il,Q_STAT,Q_LO),2)
        end if

C Q_T numbering of eigenvalues of the whole matrix (with several V''s)
        if (qlin(il,Q_T,Q_UP).ne.0) then 
          qlin(il,Q_STAT,Q_UP)=myor(qlin(il,Q_STAT,Q_UP),4)
        end if
        if (qlin(il,Q_T,Q_LO).ne.0) then
          qlin(il,Q_STAT,Q_LO)=myor(qlin(il,Q_STAT,Q_LO),4)
        end if

        if (dln(il,LN_ERR).ne.NOFIT) then
          qlin(il,Q_STAT,Q_UP)=myor(qlin(il,Q_STAT,Q_UP),16)
          qlin(il,Q_STAT,Q_LO)=myor(qlin(il,Q_STAT,Q_LO),16)
        end if
        write(fmtstr,'(A,I2,A,I2,A)')
     $       '(',DIMQLP,'I3,3X,',DIMQLP,'I3,F18.7,D12.2)'
        if (myand(ctlint(C_PRI),AP_IO).ne.0) write(*,fmtstr)
     $        (qlin(il,iqq,Q_UP),iqq=1,DIMQLP)
     $       ,(qlin(il,iqq,Q_LO),iqq=1,DIMQLP)
     $       ,dln(il,LN_FREQ),dln(il,LN_ERR)

        if (myand(ctlint(C_PRI),AP_IO).ne.0) then
          if (qlin(il,Q_S,Q_LO).eq.1) then
          write(0,'(2(3I3,2X))')
     $         qlin(il,Q_J,Q_UP)
     $         ,qlin(il,Q_TJ,Q_UP)/2
     $         ,(2*qlin(il,Q_J,Q_UP)-qlin(il,Q_TJ,Q_UP)+2)/2
     $         ,qlin(il,Q_J,Q_LO)
     $         ,qlin(il,Q_TJ,Q_LO)/2
     $         ,(2*qlin(il,Q_J,Q_LO)-qlin(il,Q_TJ,Q_LO)+2)/2
          write(0,'(I4,F15.4)')
     $         qlin(il,Q_S,Q_LO),dln(il,LN_FREQ)*1000.0d0
          end if
          if (qlin(il,Q_S,Q_LO).gt.1) then
          write(0,'(I4,F15.4)')
     $         qlin(il,Q_S,Q_LO),dln(il,LN_FREQ)*1000.0d0 !Frequencz is written here
          end if
        end if
        if ((qlin(il,Q_TJ,Q_UP).le.0).and.(qlin(il,Q_T,Q_UP).le.0)
     $       .and.(qdone(25).eq.0)) then
          write(0,'(A,i4)') 'No tau or K (up) in line',il
          stpflg=.true.
        end if
        if ((qlin(il,Q_TJ,Q_LO).le.0).and.(qlin(il,Q_T,Q_LO).le.0)
     $       .and.(qdone(26).eq.0)) then
          write(0,'(A,i4)') 'No tau or K (lo) in line',il
          stpflg=.true.
        end if
        if ((qlin(il,Q_TJ,Q_LO).gt.(2*qlin(il,Q_J,Q_LO)+1)).or.
     $       (qlin(il,Q_TJ,Q_UP).gt.(2*qlin(il,Q_J,Q_UP)+1))) then
          write(0,'(A,i4)') 'tau > 2J+1 at line',il
          stpflg=.true.
        end if
        if ((qlin(il,Q_V1,Q_LO).le.0).or.(qlin(il,Q_V1,Q_UP).le.0)) then
          write(0,'(A,i4)') 'V < 1 at line',il
          stpflg=.true.
        end if
        if ((qlin(il,Q_B,Q_LO).le.0).or.(qlin(il,Q_B,Q_UP).le.0)) then
          write(0,'(A,i4)') 'Error: B < 1 at line',il
          stpflg=.true.
        end if
        
        if ((qlin(il,Q_F,Q_LO).lt.-1).or.(qlin(il,Q_F,Q_UP).lt.-1))then
          write(0,'(A,i4)') 'Error: F < -1 at transition no.',il
          stpflg=.true.
        end if
       
        if (qlin(il,Q_F,Q_LO).ge.0) then
          if (abs(2*qlin(il,Q_J,Q_LO)-qlin(il,Q_F,Q_LO))
     $         .gt.(ctlint(C_SPIN)+ctlint(C_SPIN2))) then
            write(0,'(A,i4)') 'Error: |J-F| > I ! at transition no.',il
            stpflg=.true.
          end if
        end if
        if (qlin(il,Q_F,Q_UP).ge.0) then
          if (abs(2*qlin(il,Q_J,Q_UP)-qlin(il,Q_F,Q_UP))
     $         .gt.(ctlint(C_SPIN)+ctlint(C_SPIN2))) then
            write(0,'(A,i4)') 'Error: |J-F| > I at transition no.',il
            stpflg=.true.
          end if
        end if
        if ((qlin(il,Q_B,Q_LO).gt.DIMVB)
     $      .or.(qlin(il,Q_B,Q_UP).gt.DIMVB)) then
          write(0,'(A,i4)') 'B > DIM B' 
          stpflg=.true.
        end if
        if ((qlin(il,Q_B,Q_LO).gt.size(S_NB))
     $      .or.(qlin(il,Q_B,Q_UP).gt.size(S_NB))) then
          size(S_NB)=max(qlin(il,Q_B,Q_LO),qlin(il,Q_B,Q_UP)) ! Here it finds the maximum B in dataset.
        end if
      end do
 10   continue
        if (ctlint(C_DW).eq.1) then ! Sven 2024, for DW treatment.
         if (MODULO(size(S_NB),2).eq.1) then
            size(S_NB)=size(S_NB)+1
            write(0,*) 'For DW treatment, the maximum B quanum number', 
     $          ' is considered to be ', size(S_NB)
         end if
        end if
      if (stpflg) then
        write(*,*) 'INPUT ERROR(S)' 
        write(0,*) 'INPUT ERROR(S)'
        stop
      end if

      ctlint(C_NDATA)=il-1
      write(*,'(3(3X,A,I4),/,3X,A,F6.1)')
     $     ctlstr(C_NDATA),ctlint(C_NDATA)
     $    ,'Data Points',ndata,'Splittings',nsplit
     $    ,'Effective Data Points',efdata
      if (nsplit.ne.0) then
        msplit=msplit/nsplit
        m2split=sqrt(m2split/nsplit)
        write(*,'(3X,2(A,F12.6))') 'Mean Experimental Splitting:',msplit
     $       ,'   squared:',m2split
      end if

C     clean up pali and palc for not used ib values
      do i=1, nfit
        do j=1, pali(i,0,1)
          if (pali(i,j,2).gt.size(S_NB)) then
            pali(i,j,2)=0
            pali(i,j,1)=0
            palc(i,j)=0.0d0
          end if
        end do  
      end do
      do i=1, nfit
        do j=1, pali(i,0,1)
          if (pali(i,j,2).eq.0) goto 100
        end do
 100    continue
        ffree=j
        do j=ffree, pali(i,0,1)
          if (pali(i,j,2).gt.0) then
            pali(i,ffree,1)=pali(i,j,1)
            pali(i,ffree,2)=pali(i,j,2)
            palc(i,ffree)=palc(i,j)
            pali(i,j,1)=0
            pali(i,j,2)=0
            palc(i,j)=0.0d0
            ffree=ffree+1
          end if
        end do
        pali(i,0,1)=ffree-1
      end do

C- delete here
c        do ib=1, size(S_NB)
cc          do iv=1, DIMVV
c            if (qvv(iv,1,ib).eq.-1) goto 22
c            write(*,'(A,$)') '  V'
cc            do itop=1,ctlint(C_NTOP) 
c              write(*,'(I3,$)') qvv(iv,itop,ib)-1  
c            end do
c          end do
c 22       continue
c          write(*,*)
c        end do

      return
      end

C----------------------------------------------------------------------
      subroutine prrrp(a,da,covar,palc,pali,nfit,ib)
C     print rotational constants and errors
      implicit none
      include 'iam.fi'
      real*8  a(DIMPAR),da(DIMPAR),covar(DIMFIT,DIMFIT)
      real*8  palc(DIMFIT,-1:DIMPLC)
      integer pali(DIMFIT, 0:DIMPLC,2)
      integer nfit,ib

      integer ibj,ibk,ibd,i,j
      integer ixj,ixk,ixd
      real*8  cjk,cjd,ckd,kappa
      real*8  b(3),d(3),bb,dd
      character*3 c(3),cc
      logical noc

C     find correlation coeff between BJ BK and B- 
C     works only if no linear comb. of BJ BK B- is fitted
      ibj=0
      ibk=0
      ibd=0
      ixj=0
      ixk=0
      ixd=0
      do i=1, nfit
        do j=1, pali(i,0,1)
          if ((pali(i,j,1).eq.P_BJ).and.(pali(i,j,2).eq.ib)) then 
            ibj=i
            ixj=ixj+1
          end if
        end do
        do j=1, pali(i,0,1)
          if ((pali(i,j,1).eq.P_BK).and.(pali(i,j,2).eq.ib)) then
            ibk=i
            ixk=ixk+1
          end if
        end do
        do j=1, pali(i,0,1)
          if ((pali(i,j,1).eq.P_BD).and.(pali(i,j,2).eq.ib)) then
            ibd=i
            ixd=ixd+1
          end if
        end do
      end do
C     covar(j,i) with j<i only
      if ((ixj.eq.1).and.(ixk.eq.1)) then
        cjk=covar(min(ibj,ibk),max(ibj,ibk))
      else
        cjk=0.0d0
      end if
      if ((ixj.eq.1).and.(ixd.eq.1)) then
        cjd=covar(min(ibj,ibd),max(ibj,ibd))
      else
        cjd=0.0d0
      end if
      if ((ixk.eq.1).and.(ixd.eq.1)) then
        ckd=covar(min(ibk,ibd),max(ibk,ibd))
      else
        ckd=0.0d0
      end if
      noc=.true.
      if ((ixj.ne.1).or.(ixk.ne.1).or.(ixd.ne.1)) noc=.false.
      if (noc) then
        write(*,'(A)') ' Rotational Constants and Errors (in GHz)'
      else
        write(*,'(A)')
     $       ' Rotational Constants and approx. Errors (in GHz)'
      end if

      c(3)='B_z'
      b(3)=a(1)+a(2)
      d(3)=dsqrt(da(1)**2+da(2)**2+2.0d0*da(1)*da(2)*cjk) !changed for correct A constant error
      c(2)='B_y'
      b(2)=a(1)-a(3)
      d(2)=dsqrt(da(1)**2+da(3)**2-2.0d0*da(1)*da(3)*cjd)
      c(1)='B_x'
      b(1)=a(1)+a(3)
      d(1)=dsqrt(da(1)**2+da(3)**2+2.0d0*da(1)*da(3)*cjd)

      do i=1,3
        do j=1,3
          if (b(i).gt.b(j)) then 
            bb=b(i)
            b(i)=b(j)
            b(j)=bb
            dd=d(i)
            d(i)=d(j)
            d(j)=dd
            cc=c(i)
            c(i)=c(j)
            c(j)=cc
          end if
        end do
      end do
      
      write(*,'(2X,A7,1F18.9,1F18.9)') c(1),b(1),d(1)
      write(*,'(2X,A7,1F18.9,1F18.9)') c(2),b(2),d(2)
      write(*,'(2X,A7,1F18.9,1F18.9)') c(3),b(3),d(3)

      kappa=(2.d0*b(2)-b(1)-b(3))/(b(1)-b(3))
      write(*,'(A,F10.5)') ' Ray''s kappa ',kappa
      return
      end

C----------------------------------------------------------------------
      subroutine prpot(a,da,inkj,inkc,incm)
C     print potential paramaters
      implicit none
      include 'iam.fi'
      real*8  a(DIMPAR),da(DIMPAR),inkj,inkc,incm
      integer itop,ift
      include 'iamdata.fi'

      do itop=1, ctlint(C_NTOP)
        ift=DIMPIR*(itop-1)
        write(*,'(1X,2(A,F12.6,A,F8.6,A))')
     $       parstr(P1_VN1+ift)
     $       ,a(P1_VN1+ift)*inkj,' kj +/- ',da(P1_VN1+ift)*inkj
     $       ,' kj ','    '
     $       ,a(P1_VN1+ift)*inkc,' kcal +/- ',da(P1_VN1+ift)*inkc
     $       ,' kcal'
        write(*,'((A,F12.6,A,F8.4,A),A,F11.6)')
     $       '        '
     $       ,a(P1_VN1+ift)*incm,' cm +/- ',da(P1_VN1+ift)*incm
     $       ,' cm ','   s='
     $       ,4.0d0*a(P1_VN1+ift)/(9.0d0*a(P1_F))
      end do
      if (a(P1_VN2+ift).ne.0.0d0) then
      do itop=1, ctlint(C_NTOP)
        ift=DIMPIR*(itop-1)
        write(*,'(1X,2(A,F12.6,A,F8.6,A))')
     $       parstr(P1_VN2+ift)
     $       ,a(P1_VN2+ift)*inkj,' kj +/- ',da(P1_VN2+ift)*inkj
     $       ,' kj ','    '
     $       ,a(P1_VN2+ift)*inkc,' kcal +/- ',da(P1_VN2+ift)*inkc
     $       ,' kcal'
        write(*,'((A,F12.6,A,F8.4,A))')
     $       '        '
     $       ,a(P1_VN2+ift)*incm,' cm +/- ',da(P1_VN2+ift)*incm
     $       ,' cm '
      end do
      end if
      if (a(P_VSS).ne.0.0d0) then
      write(*,'(/,1X,2(A,F12.6,A,F8.6,A))')
     $     parstr(P_VSS)
     $     ,a(P_VSS)*inkj,' kj +/- ',da(P_VSS)*inkj
     $     ,' kj ','    '
     $     ,a(P_VSS)*inkc,' kcal +/- ',da(P_VSS)*inkc,' kcal'
      write(*,'((A,F12.6,A,F8.4,A))')
     $     '        '
     $     ,a(P_VSS)*incm,' cm +/- ',da(P_VSS)*incm
     $     ,' cm '
      end if
      if (a(P_VCC).ne.0.0d0) then
      write(*,'(1X,2(A,F12.6,A,F8.6,A))')
     $     parstr(P_VCC)
     $     ,a(P_VCC)*inkj,' kj +/- ',da(P_VCC)*inkj
     $     ,' kj ','    '
     $     ,a(P_VCC)*inkc,' kcal +/- ',da(P_VCC)*inkc,' kcal'
      write(*,'((A,F12.6,A,F8.4,A))')
     $     '        '
     $     ,a(P_VCC)*incm,' cm +/- ',da(P_VCC)*incm
     $     ,' cm '
      end if
      write(*,*)
      return
      end
C----------------------------------------------------------------------
      subroutine prirp(a,da,indeg)
      implicit none
      include 'iam.fi'
      real*8  a(DIMPAR),da(DIMPAR),indeg
      integer itop,ift,i
      real*8  lz,lx,ly,aclz,aclx,acly,f0calc,pi
      real*8  ax,ay,az,axerr,ayerr,azerr,f0err,iaerr,sq,sqd
      real*8  adjdev(2:5,3)
      integer myand
      external myand

      pi=dacos(-1.0d0)
      do itop=1, ctlint(C_NTOP)
        ift=DIMPIR*(itop-1)
        call recalf(a(P_BJ),a(P_BK),a(P_BD)
     $       ,a(P1_RHO+ift)
     $       ,a(P1_BETA+ift)
     $       ,a(P1_GAMA+ift)
     $       ,f0calc,lx,ly,lz)
        aclx=acos(lx)
        acly=acos(ly)
        aclz=acos(lz)
        i=1
        call recalf(a(P_BJ),a(P_BK),a(P_BD)
     $       ,a(P1_RHO+ift)+da(P1_RHO+ift)
     $       ,a(P1_BETA+ift)
     $       ,a(P1_GAMA+ift)
     $       ,adjdev(2,i),adjdev(3,i)
     $       ,adjdev(4,i),adjdev(5,i))
        adjdev(2,i)=abs(adjdev(2,i)-f0calc)
        adjdev(3,i)=abs(dacos(adjdev(3,i))-aclx) 
        adjdev(4,i)=abs(dacos(adjdev(4,i))-acly) 
        adjdev(5,i)=abs(dacos(adjdev(5,i))-aclz) 
        i=2
        call recalf(a(P_BJ),a(P_BK),a(P_BD)
     $       ,a(P1_RHO+ift)
     $       ,a(P1_BETA+ift)+da(P1_BETA+ift)
     $       ,a(P1_GAMA+ift)
     $       ,adjdev(2,i),adjdev(3,i)
     $       ,adjdev(4,i),adjdev(5,i))
        adjdev(2,i)=abs(adjdev(2,i)-f0calc)
        adjdev(3,i)=abs(dacos(adjdev(3,i))-aclx) 
        adjdev(4,i)=abs(dacos(adjdev(4,i))-acly) 
        adjdev(5,i)=abs(dacos(adjdev(5,i))-aclz) 
        i=3
        call recalf(a(P_BJ),a(P_BK),a(P_BD)
     $       ,a(P1_RHO+ift)
     $       ,a(P1_BETA+ift)
     $       ,a(P1_GAMA+ift)+da(P1_GAMA+ift)
     $       ,adjdev(2,i),adjdev(3,i)
     $       ,adjdev(4,i),adjdev(5,i))
        adjdev(2,i)=abs(adjdev(2,i)-f0calc)
        adjdev(3,i)=abs(dacos(adjdev(3,i))-aclx) 
        adjdev(4,i)=abs(dacos(adjdev(4,i))-acly) 
        adjdev(5,i)=abs(dacos(adjdev(5,i))-aclz) 
        
        if (myand(ctlint(C_ADJF),8).eq.0) then
          f0err=adjdev(2,1)+adjdev(2,2)+adjdev(2,3)
          iaerr=(505.37907d0/f0calc**2)
     $         *(adjdev(2,1)+adjdev(2,2)+adjdev(2,3))
          write(*,*)
     $         'Following values are derived from rho, beta, gamma'
          if (myand(ctlint(C_ADJF),1).eq.0) then
            f0err=0.0d0
            iaerr=0.0d0
          end if
        else
          f0err=da(P1_F0+ift)
          iaerr=(505.37907d0/f0calc**2)*da(P1_F0+ift)
        end if
        write(*,'(A,F18.9,F18.9)') ' F0(calc)'
     $       ,f0calc,f0err
        write(*,'(A,F18.9,F18.9)') ' I_alpha '
     $       ,505.37907d0/f0calc,iaerr

        if (myand(ctlint(C_ADJF),16).eq.0) then
          ax=aclx
          ay=acly
          az=aclz
          axerr=(adjdev(3,1)+adjdev(3,2)+adjdev(3,3))
          ayerr=(adjdev(4,1)+adjdev(4,2)+adjdev(4,3))
          azerr=(adjdev(5,1)+adjdev(5,2)+adjdev(5,3))
        else
          az=a(P1_ANGZ+ift)
          azerr=abs(da(P1_ANGZ+ift))
          sq=sqrt(1.0d0-cos(az)**2)
          sqd=sqrt(1.0d0-cos(az+azerr)**2)
          ax=acos(cos(a(P1_ANGX+ift))*sq)
          ay=acos(sin(a(P1_ANGX+ift))*sq)
          axerr=abs(ax-acos(cos(a(P1_ANGX+ift)+da(P1_ANGX+ift))*sq))
     $         +abs(ax-acos(cos(a(P1_ANGX+ift))*sqd))
          ayerr=abs(ay-acos(sin(a(P1_ANGX+ift)+da(P1_ANGX+ift))*sq))
     $         +abs(ay-acos(sin(a(P1_ANGX+ift))*sqd))
        end if

        write(*,'(A,3F14.4,A)') ' <(i,x)  <(i,y)  <(i,z)'
     $       ,ax*indeg,ay*indeg,az*indeg
        write(*,'(A,3F14.4)') 'd<(i,x) d<(i,y) d<(i,z)'
     $       ,axerr*indeg,ayerr*indeg,azerr*indeg
        write(*,*)
      end do
      return
      end

C     ---------------------------------------------------------------------
      subroutine funpr(rofit,ndat,chi2,wght)
C     print transition list
      implicit none
      include 'iam.fi'
      real*8  rofit
      real*8  chi2(DIMVB,DIMVV,2),wght(DIMVB,DIMVV,2)
      real*8  g2c
      real*8  maxdf ! stores frequency of maximum deviation line
      integer ndat(DIMVB,DIMVV,2)
      integer ii,ref,avg
      integer i,ib,iv
      integer fp1len,fp2len,fp3len,fp4len,fp5len,fp6len,fp7len,fp8len
      integer fp2blen!Herbers 2024, for F1
      integer fp4x
      real*8  reffreq,refcalc,fcalc,ffreq,maxd,dmaxd
      integer maxi,gam
      logical lfit
      character*1  errstr
      character*10 fp1astr
      character*10 fp1bstr
      character*10 fp1cstr
      character*10 fp1dstr
      character*10 fp1estr
      character*40 fp1str,lb1str
      character*40 fp2str,lb2str
      character*40 fp2bstr,lb2bstr !herbers2024
      character*40 fp3str,lb3str
      character*40 fp4str,lb4str
      character*60 fp5str,lb5str !Herbers2023 - needed some extra
      character*40 fp6str,lb6str
      character*40 fp7str,lb7str
      character*40 fp8str,lb8str
      character*40 fp9str,lb9str
      character*10 sym_lo,sym_up
      character*12 refstr
      integer myand
      external myand
      include 'iamdata.fi'
      g2c=1/29.9792458d0! To convert from GHz to cm-1

      do ib=1, DIMVB
        do iv=1, DIMVV
          do ii=1,2
            chi2(ib,iv,ii)=0.0d0
            wght(ib,iv,ii)=0.0d0
            ndat(ib,iv,ii)=0
          end do
        end do
      end do
      maxd=0.0d0

      lfit=.false.
      fp4x=0
      do i=1, size(S_G)
        do ii=1, len(qnostr(1)) 
          if (qnostr(MAXQC+i)(ii:ii).ne.' ') then
            if (fp4x.lt.ii) fp4x=ii
          end if
        end do
      end do

      do i=1, ctlint(C_NDATA)
        ib=qlin(i,Q_B,Q_UP)
        iv=qlin(i,Q_V1,Q_UP)
        errstr=' '
        if (((dnv(i,NV_ENG,Q_UP)-dnv(i,NV_ENG,Q_LO)).eq.0.0d0)
     $       .and.(dln(i,LN_FREQ).eq.0.0d0))  then
          write(*,'(A)') ' - - - - -'
          return
        end if
        do ii=1,10
          fp1astr(ii:ii)=' '
          fp1bstr(ii:ii)=' '
          fp1cstr(ii:ii)=' '
          fp1dstr(ii:ii)=' '
          fp1estr(ii:ii)=' '
          sym_lo(ii:ii)=' '
          sym_up(ii:ii)=' '
        end do

        do ii=1, 40
          fp1str(ii:ii)=' '
          lb1str(ii:ii)=' '
          fp2bstr(ii:ii)=' '
          fp2str(ii:ii)=' '
          lb2bstr(ii:ii)=' '
          lb2str(ii:ii)=' '
          fp3str(ii:ii)=' '
          lb3str(ii:ii)=' '
          fp4str(ii:ii)=' '
          lb4str(ii:ii)=' '
          fp5str(ii:ii)=' '
          lb5str(ii:ii)=' '
          fp6str(ii:ii)=' '
          lb6str(ii:ii)=' '
          fp7str(ii:ii)=' '
          lb7str(ii:ii)=' '
          fp8str(ii:ii)=' '
          lb8str(ii:ii)=' '
          fp9str(ii:ii)=' '
          lb9str(ii:ii)=' '
        end do
        
        
c     write(fp1str,'(A,2I3,A,2I4,A,2I3,A,I2)') 
c     $      'J',qlin(i,Q_J,Q_UP),qlin(i,Q_J,Q_LO)
c     $     ,' K',qlin(i,Q_K,Q_UP),qlin(i,Q_K,Q_LO)
c     $     ,' t',qlin(i,Q_TJ,Q_UP),qlin(i,Q_TJ,Q_LO)
c     $     ,' B',qlin(i,Q_B,Q_UP)

C     J K- K+ notation
c     if (size(S_NB).gt.1) then
c     write(lb1str,'(3A3,3A3,A,A2)')
c     $         'J','K-','K+','J','K-','K+','  ',' B'
c     write(fp1str,'(3I3,3I3,A,I2)') 
c     $         qlin(i,Q_J,Q_UP)
c     $         ,qlin(i,Q_TJ,Q_UP)/2
c     $         ,(2*qlin(i,Q_J,Q_UP)+2-qlin(i,Q_TJ,Q_UP))/2
c     $         ,qlin(i,Q_J,Q_LO)
c     $         ,qlin(i,Q_TJ,Q_LO)/2
c     $         ,(2*qlin(i,Q_J,Q_LO)+2-qlin(i,Q_TJ,Q_LO))/2
c     $         ,' B',qlin(i,Q_B,Q_UP)
c     fp1len=22
c     else
        
        write(fp1astr,'(I3)') 
     $       qlin(i,Q_J,Q_UP)
        write(fp1cstr,'(I3)') 
     $       qlin(i,Q_J,Q_LO)
        write(lb1str,'(3A3,3A3,A,A2)')
     $       'J','K-','K+','J','K-','K+','  ',' B'

        gam=qlin(i,Q_S,Q_UP)
        if (gam.ne.0) then
          if (.false.) then !(gamma(gam,0).eq.0) then !deactived the K printing, because I wouldnt know what to use it for. - Sven 2024
            write(fp1bstr,'(A3,I3)') ' K ',qlin(i,Q_K,Q_UP)
          else
            write(fp1bstr,'(2I3)') 
     $           qlin(i,Q_TJ,Q_UP)/2
     $           ,(2*qlin(i,Q_J,Q_UP)+2-qlin(i,Q_TJ,Q_UP))/2
          end if
        else
          write(fp1bstr,'(2I3)') 
     $         qlin(i,Q_TJ,Q_UP)/2
     $         ,(2*qlin(i,Q_J,Q_UP)+2-qlin(i,Q_TJ,Q_UP))/2
        end if

        gam=qlin(i,Q_S,Q_LO)
        if (gam.ne.0) then
          if (.false.) then !(gamma(gam,0).eq.0) then !deactived the K printing, because I wouldnt know what to use it for. - Sven 2024
            write(fp1dstr,'(A3,I3)') ' K ',qlin(i,Q_K,Q_LO)
          else
            write(fp1dstr,'(2I3)') 
     $           qlin(i,Q_TJ,Q_LO)/2
     $           ,(2*qlin(i,Q_J,Q_LO)+2-qlin(i,Q_TJ,Q_LO))/2
          end if
        else
          write(fp1dstr,'(2I3)') 
     $         qlin(i,Q_TJ,Q_LO)/2
     $         ,(2*qlin(i,Q_J,Q_LO)+2-qlin(i,Q_TJ,Q_LO))/2
        end if

        fp1len=19
        if (size(S_NB).gt.1) then
          write(fp1estr,'(A,I2)') 
     $          ' B',qlin(i,Q_B,Q_UP)
          fp1len=fp1len+4
        end if
        fp1str=':'//fp1astr(1:3)//fp1bstr(1:6)//fp1cstr(1:3)
     $       //fp1dstr(1:6)//fp1estr(1:4)

        fp2blen=1
        if (ctlint(C_SPIN2).ne.0) then
          write(lb2bstr,'(A,2A3)') '       ','F1 ','F1 '
          write(fp2bstr,'(A,2I3)')
     $         ' F1',qlin(i,Q_F1,Q_UP),qlin(i,Q_F1,Q_LO) ! this makes F quantum numbers printed.
          fp2blen=fp2blen+8 ! adding to F printing
        end if
        fp2len=1
C        write(0,*) ,fp1len,'$',fp1str
        if (ctlint(C_SPIN).ne.0) then
         if (ctlint(C_SPIN2).ne.0) then
          write(lb2str,'(A,2A3)') '  ','F','F'
          write(fp2str,'(A,2I3)')
     $         ' F',qlin(i,Q_F,Q_UP),qlin(i,Q_F,Q_LO) ! this makes F quantum numbers printed.
         else 
          write(lb2str,'(A,2A3)') '  ','F','F'
          write(fp2str,'(A,2I3)')
     $         ' F',qlin(i,Q_F1,Q_UP),qlin(i,Q_F1,Q_LO) ! use F1 as F in single nucleus case.
         end if
          fp2len=fp2len+7
        end if





        fp3len=1
        if ((ctlint(C_NTOP).ne.0).and.(size(S_VV).ne.1)) then
          if (qlin(i,Q_V1,Q_LO).eq.qlin(i,Q_V1,Q_up)) then
            write(fp3str,'(A,I2)')
     $           ' V',qlin(i,Q_V1,Q_UP)
            fp3len=2+2
          else
            write(fp3str,'(A,2I2)')
     $           ' V',qlin(i,Q_V1,Q_UP),qlin(i,Q_V1,Q_LO)
            fp3len=2+4
          end if
        end if 

        if (qlin(i,Q_S,Q_UP).eq.0) then
          write(sym_up,'(A)') 'rigid'
          write(sym_lo,'(A)') 'rigid'
        end if
        if (qlin(i,Q_S,Q_UP).lt.0) then
          write(sym_up,'(A)') 'rig-1'
          write(sym_lo,'(A)') 'rig-1'
        end if
        if (qlin(i,Q_S,Q_UP).gt.0) then
          write(sym_up,'(A)') 
     $         qnostr(MAXQC+qlin(i,Q_S,Q_UP))
          write(sym_lo,'(A)') 
     $         qnostr(MAXQC+qlin(i,Q_S,Q_LO))
        end if

        if ((ctlint(C_NTOP).ne.0)) then
          If (fp4x.eq.0) then
            if (qlin(i,Q_S,Q_LO).eq.qlin(i,Q_S,Q_UP)) then
              write(fp4str,'(A,I2)')
     $             ' S',qlin(i,Q_S,Q_UP)
              fp4len=4
              lb4str=' Sym'
            else
              write(fp4str,'(A,2I2)')
     $             ' S',qlin(i,Q_S,Q_UP),qlin(i,Q_S,Q_LO)
              fp4len=6
              lb4str='SymSym'
            end if
          else
            if (qlin(i,Q_S,Q_LO).eq.qlin(i,Q_S,Q_UP)) then
              fp4len=fp4x+1
              write(fp4str,'(A,A)') ' ',
     $             sym_up(1:fp4x)
              lb4str='Sym'
            else
              fp4len=2*fp4x+2
              write(fp4str,'(A,3A)') ' ',
     $             sym_up(1:fp4x),' ',sym_lo(1:fp4x)
              lb4str='Sym Sym'
            end if
          end if
        else
          fp4len=4
        end if

        avg=qlin(i,Q_AVG,Q_UP)
        fcalc=(dnv(avg,NV_ENG,Q_UP)-dnv(avg,NV_ENG,Q_LO)
     $       +dnv(i,NV_ENG,Q_UP)-dnv(i,NV_ENG,Q_LO))/2.0d0
        ffreq=dln(i,LN_FREQ)

        ref=qlin(i,Q_REF,Q_UP)
        refcalc=0.0d0
        reffreq=0.0d0
        if ((ref.ne.0).and.(ref.ne.i)) then
          refcalc=dnv(ref,NV_ENG,Q_UP)-dnv(ref,NV_ENG,Q_LO)
          reffreq=dln(ref,LN_FREQ)
        end if
        if (ctlint(C_PRINT).ne.6) then
        write(lb5str,'(A,A13,A16,A14)')
     $       ' ','   calc/GHz ','      diff/kHz  ','   obs/GHz '
        else
        write(lb5str,'(A,A18,A21,A17)')
     $       ' ','c. MW/GHz IR/cm-1 ','o.-c. MW/kHz IR/cm-1 '
     $       ,'o. MW/GHz IR/cm-1' 
        end if     
        if (dln(i,LN_ERR).ne.NOFIT) then
          lfit=.true.        
          fp5len=44 !Herbers2023 - some extra also chardef changed
          if (i.ne.avg) then ! Herbers2024 - changed diff to kHz scale
            if (abs(fcalc).gt.999.99) then
             if (ctlint(C_PRINT).ne.6) then
              write(fp5str,'(F14.6,F13.1,A1,F14.6)')  !Herbers2023
     $              fcalc-refcalc
     $             ,(ffreq-reffreq
     $             -fcalc+refcalc)*1000.0d3,'&'
     $             ,dln(i,LN_FREQ)-reffreq
             else 
               write(fp5str,'(F14.5,F13.5,A1,F14.5)')  !Herbers2023
     $              (fcalc-refcalc)*g2c
     $             ,(ffreq-reffreq
     $             -fcalc+refcalc)*g2c,'&'
     $             ,dln(i,LN_FREQ)*g2c-reffreq*g2c    
             end if     
            else
              write(fp5str,'(F14.9,F13.3,A1,F14.9)')  !Herbers2023
     $              fcalc-refcalc
     $             ,(ffreq-reffreq
     $             -fcalc+refcalc)*1000.0d3,'&'
     $             ,dln(i,LN_FREQ)-reffreq
            end if
          else
            if (abs(fcalc).gt.999.99) then
             if (ctlint(C_PRINT).ne.6) then
              write(fp5str,'(F14.6,F13.1,A1,F14.6)') !Herbers2023
     $             fcalc-refcalc
     $             ,(ffreq-reffreq
     $             -fcalc+refcalc)*1000.0d3,' '
     $             ,dln(i,LN_FREQ)-reffreq
              else 
               write(fp5str,'(F14.5,F13.5,A1,F14.5)') !Herbers2023
     $             (fcalc-refcalc)*g2c
     $             ,(ffreq-reffreq
     $             -fcalc+refcalc)*g2c,' '
     $             ,dln(i,LN_FREQ)*g2c-reffreq*g2c
             endif
            else
              write(fp5str,'(F14.9,F13.3,A1,F14.9)')  !Herbers2023
     $             fcalc-refcalc
     $             ,(ffreq-reffreq
     $             -fcalc+refcalc)*1000.0d3,' '
     $             ,dln(i,LN_FREQ)-reffreq
            end if
          end if
          if ((ref.ne.0).and.(ref.ne.i)) then
            ii=1
          else
            ii=2
          end if
          dmaxd=((ffreq-reffreq-fcalc+refcalc)/dln(i,LN_ERR))**2
          if (dmaxd.gt.maxd) then
            maxdf = dln(i,LN_FREQ)-reffreq
            maxd=dmaxd
            maxi=i
          end if
            
          chi2(ib,iv,ii)=chi2(ib,iv,ii)+ dmaxd
          wght(ib,iv,ii)=wght(ib,iv,ii)+ 1.0d0/dln(i,LN_ERR)**2
          ndat(ib,iv,ii)=ndat(ib,iv,ii)+ 1           

          if (rofit.ne.0.0d0) then
            fp6len=18
            write(fp6str,'(A,1D8.1,F5.2)') 
     $           ' Err'
     $           ,dln(i,LN_ERR)
     $           ,dln(i,LN_PSI)
          else
            fp6len=13
            write(fp6str,'(A,1D8.1)') 
     $           ' Err'
     $           ,dln(i,LN_ERR)
          end if
        else
          refstr='        --- '
          if (avg.ne.i) refstr='        --& '
          fp5len=44
          if ((ref.ne.0).and.(ref.ne.i)) then        
            if (abs(fcalc).gt.999.99) then
              write(fp5str,'(F12.4,F9.1,A1,F12.4)') 
     $             fcalc
     $             ,(-fcalc+refcalc)*1000.0d3,'#'
     $             ,reffreq-refcalc+fcalc
            else
              write(fp5str,'(F12.7,F9.3,A1,F12.7)') 
     $             fcalc
     $             ,(-fcalc+refcalc)*1000.0d3,'#'
     $             ,reffreq-refcalc+fcalc
            end if
          else
            if (abs(fcalc).gt.999.99) then
             if (ctlint(C_PRINT).ne.6) then
              write(fp5str,'(F12.4,A10,A12)') 
     $             fcalc
     $             ,'      --- '
     $             ,refstr
             else 
              write(fp5str,'(F12.5,A10,A12)') 
     $             fcalc*g2c
     $             ,'      --- '
     $             ,refstr             
             endif
            else
              write(fp5str,'(F12.7,A10,A12)') 
     $             fcalc
     $             ,'      --- '
     $             ,refstr
            end if
          end if
          if (rofit.ne.0.0d0) then
            fp6len=18
          else
            fp6len=13
          end if
          write(fp6str,'(A)') 
     $         '      ---  '
        end if

        fp7len=1
        if (ctlint(C_INTS).ne.0) then
          write(fp7str,'(A,1F7.3)') 
     $         ' log10_str**2',log10(dln(i,LN_INT)) !Herbers2023
          fp7len=20!Herbers2023
        end if

        write(lb8str,'(2A11)') ' 1.K 2.K % ',' 1.K 2.K % '
        write(fp8str,'(2(1X,2I3,1X,I3))')
     $       qlin(i,Q_K,Q_UP),qlin(i,Q_K2,Q_UP),qlin(i,Q_GK,Q_UP),
     $       qlin(i,Q_K,Q_LO),qlin(i,Q_K2,Q_LO),qlin(i,Q_GK,Q_LO)
        fp8len=22
        
        if ((ref.ne.0).and.(ref.ne.i)) then        
          write(fp9str,'(A1,I4)') 
     $         '#',ref-i
        end if

        if ((myand(qlin(i,Q_STAT,Q_UP),1).ne.1)
     $       .or.(myand(qlin(i,Q_STAT,Q_LO),1).ne.1))  errstr='*'
        if (myand(qlin(i,Q_STAT,Q_UP),8).ne.0)  errstr='x'
        if (i.eq.1) then
          if (myand(ctlint(C_PRI),AP_TE).ne.0) then
            write(*,'(A1,3X,10A)')
     $           errstr
     $           ,lb1str(1:fp1len)
     $           ,lb2bstr(1:fp2blen+3) !herbers 2024, added for F1
     $           ,lb2str(1:fp2len)
     $           ,lb3str(1:fp3len)
     $           ,lb4str(1:fp4len)
     $           ,lb5str(1:fp5len),lb6str(1:fp6len),lb9str(1:5)
     $           ,lb8str(1:fp8len),lb7str(1:fp7len)
          else
            write(*,'(A1,3X,9A)')
     $           errstr
     $           ,lb1str(1:fp1len)
     $           ,lb2bstr(1:fp2blen+3) !herbers 2024, added for F1
     $           ,lb2str(1:fp2len)
     $           ,lb3str(1:fp3len)
     $           ,lb4str(1:fp4len)
     $           ,lb5str(1:fp5len),lb6str(1:fp6len),lb9str(1:5)
     $           ,lb7str(1:fp7len)
          end if
        end if
        
        if (myand(ctlint(C_PRI),AP_TE).ne.0) then
          write(*,'(A1,I6,14A)')!Sven2023/2024
     $         errstr,i
     $         ,fp1str(1:fp1len)
     $         ,fp2bstr(1:fp2blen) !Herbers2024
     $         ,fp2str(1:fp2len)
     $         ,fp3str(1:fp3len)
     $         ,fp4str(1:fp4len)
     $         ,fp5str(1:fp5len),fp6str(1:fp6len),fp9str(1:5)
     $         ,fp8str(1:fp8len),fp7str(1:fp7len)
     $         ,' ',sym_lo(1:5),'-',sym_up(1:5)
        else
          write(*,'(A1,I6,A,13A)')!Sven2023/2924
     $         errstr,i
     $         ,fp1str(1:fp1len)
     $         ,fp2bstr(1:fp2blen) !Herbers2024
     $         ,fp2str(1:fp2len)
     $         ,fp3str(1:fp3len)
     $         ,fp4str(1:fp4len)
     $         ,fp5str(1:fp5len),fp6str(1:fp6len),fp9str(1:5)
     $         ,fp7str(1:fp7len)
     $         ,' ',sym_lo(1:5),'-',sym_up(1:5)
        end if
      end do                    ! ndata
      if (lfit) then
       if (ctlint(C_PRINT).ne.6) then
        write(*,'(A,I4,A,F17.3)') ' Maximum (obs-calc)/err in line'
     $       ,maxi,' '
     $       ,sqrt(maxd)*1000.d3*dln(maxi,LN_ERR)
       else
        if (maxdf.le.999) then
          write(*,'(A,I4,A,F17.3,A)') ' Maximum (obs-calc)/err in line'
     $       ,maxi,' '
     $       ,sqrt(maxd)*1000.d3*dln(maxi,LN_ERR), ' kHz'
        else
          write(*,'(A,I4,A,F17.5,A)') ' Maximum (obs-calc)/err in line'
     $       ,maxi,' '
     $       ,sqrt(maxd)*g2c*dln(maxi,LN_ERR), ' cm-1'
        end if
       end if
             
      end if
      return
      end 

C----------------------------------------------------------------------
      subroutine prpar(a)
      implicit none
      include 'iam.fi'
      real*8 a(DIMPAR,DIMVB)
      real*8 f
      real*8 bx,by,bz
      real*8 agam,abet,arho
      real*8 rhox1,rhoy1,rhoz1,rhox2,rhoy2,rhoz2
      real*8 bxeff,byeff,bzeff
      real*8 qxeff,qyeff,qzeff
      integer if1,it1,if2,it2,b

      do b=1, size(S_NB)
        write(*,'(A,I2)') ' B:',b
        bx=a(P_BJ,b)+a(P_BD,b)
        by=a(P_BJ,b)-a(P_BD,b)
        bz=a(P_BJ,b)+a(P_BK,b)
        bxeff=bx
        byeff=by
        bzeff=bz
        write(*,'(A,3F18.9)') ' Bx By Bz          ',bx,by,bz
        do it1=1, ctlint(C_NTOP)
          qxeff=0.0d0
          qyeff=0.0d0
          qzeff=0.0d0
          if1=(it1-1)*DIMPIR  
          abet=a(if1+P1_BETA,b)
          agam=a(if1+P1_GAMA,b)
          arho=a(if1+P1_RHO,b)
          rhoz1=dcos(abet)*arho
          rhox1=dsin(abet)*dcos(agam)*arho
          rhoy1=dsin(abet)*dsin(agam)*arho
          write(*,'(A,3F18.9)') ' rho_x rho_y rho_z ',rhox1,rhoy1,rhoz1
          do it2=1, ctlint(C_NTOP)
            if2=(it2-1)*DIMPIR
            abet=a(if2+P1_BETA,b)
            agam=a(if2+P1_GAMA,b)
            arho=a(if2+P1_RHO,b)
            rhoz2=dcos(abet)*arho
            rhox2=dsin(abet)*dcos(agam)*arho
            rhoy2=dsin(abet)*dsin(agam)*arho
            if (it1.ne.it2) then
              f=a(P_FF,b)
            else
              f=a(if1+P1_F,b)
            end if
            bxeff=bxeff+rhox1*rhox2*f
            byeff=byeff+rhoy1*rhoy2*f
            bzeff=bzeff+rhoz1*rhoz2*f
            qxeff=qxeff-rhox2*f
            qyeff=qyeff-rhoy2*f
            qzeff=qzeff-rhoz2*f
          end do
          write(*,'(A,3F18.9)') ' Qx Qy Qz effektiv ',qxeff,qyeff,qzeff
        end do
        write(*,'(A,3F18.9)') ' Bx By Bz effektiv ',bxeff,byeff,bzeff
      end do

      return
      end

C----------------------------------------------------------------------
      subroutine pra(a,da,ifit,ipr,ierr,b)
C     print list of parameters for b or all b's if b=0
C     ipr = 5 print all parameters
C     ipr = 4 print selected parameters
C     ipr = 3 print parameters not equal zero
C     ipr = 2 print parameters which are included in the fit
C             plus derived parameters.
C     ipr = 1 print parameters which are included in the fit
C     ipr = 0 print nothing
C     ierr =1 print errors
      implicit none
      include 'iam.fi'
      real*8 a(DIMPAR,DIMVB),da(DIMPAR,DIMVB)
      integer ifit(DIMPAR,DIMVB),ipr,ierr,b
      integer i,itop,ip,ix,il,b1,b2,ib,adj,b3,b4
      logical lprint,lx,allb,lcom
      character*10 pstr,pistr
      character*30 xstr
      character*1  c
      integer myand,len_c
      external myand,len_c
      include 'iamdata.fi'
            
      if (myand(ipr,7).eq.0) return
      if (b.eq.0) then
        b1=1
        b2=size(S_NB)
      else
      end if
C     write(*,'(/,5X,A)') 'Parameters and Errors'
      do i=1, DIMPRR
        allb=.false.
        if (b.eq.0) then
          b1=1
          b2=size(S_NB)
          allb=.true.
          do ib=b1,b2-1
            if (a(i,ib).ne.a(i,ib+1)) allb=.false.
          end do
          if (allb) then 
            b1=1
            b2=1
          end if
        else
          b1=b
          b2=b
        end if

        do ib=b1,b2
          lprint=.false.
          lcom=.false.
          c=' ' 
          adj=ctlnb(CB_ADJ,ib)
          if ((i.eq.P_FF).and.(myand(adj,2).ne.0)) lcom=.true.
          if (myand(ipr,7).eq.1) then
            if (ifit(i,ib).ne.0) lprint=.true.
          end if
          if (myand(ipr,7).eq.2) then
            if (ifit(i,ib).ne.0) lprint=.true.
            if (lcom) lprint=.true.
          end if
          if (myand(ipr,7).eq.3) then
            if (a(i,ib).ne.0.0d0) lprint=.true.
          end if
          if (myand(ipr,7).eq.4) then
            lprint=.true.
            if  (parstr(i)(1:1).eq.'_') lprint=.false.
            if ((parstr(i)(1:1).eq.'C').and.
     $           (ctlint(C_SPIN).eq.0)) lprint=.false.
            if ((parstr(i)(1:1).eq.'c').and.
     $           (ctlint(C_SPIN).eq.0)) lprint=.false.
            if (((parstr(i)(1:1).eq.'H').or.(parstr(i)(1:1).eq.'h'))
     $           .and.(a(i,ib).eq.0.0d0)) lprint=.false.
            if (((parstr(i)(1:1).eq.'M').or.(parstr(i)(1:1).eq.'m'))
     $           .and.(a(i,ib).eq.0.0d0)) lprint=.false.
            if (((parstr(i)(1:1).eq.'V').or.(parstr(i)(1:1).eq.'v'))
     $           .and.(ctlint(C_NTOP).eq.0)) lprint=.false.
            if ((parstr(i)(1:1).eq.'F')
     $           .and.(ctlint(C_NTOP).eq.0)) lprint=.false.
            if ((parstr(i)(1:3).eq.'F12')
     $           .and.(ctlint(C_NTOP).le.1)) lprint=.false.
          end if
          if (myand(ipr,7).ge.5) then
            lprint=.true.
          end if
          call fillsp(pstr)
          if (allb) then
            pstr=parstr(i)
          else
            write(pstr,'(2A,I1,A)') parstr(i)(1:len_c(parstr(i)))
     $           ,'(',ib,')'
          end if
          if (lcom) c='\\'
          if (lprint) then
            write(*,'(A,A10,$)') c,pstr!print of BJ BK B-
            call pra_out(a(i,ib),da(i,ib),ifit(i,ib),ierr,lcom,xstr)
            write(*,'(A)') xstr
          end if
        end do
      end do

      if (ctlint(C_NTOP).eq.0) return!exception for no rotor.
      do ix=1 , DIMPIR
        ip=DIMPRR+ix
        allb=.false.
        if (b.eq.0) then
          b1=1
          b2=size(S_NB)
          allb=.true.
          do ib=b1,b2-1
            if (a(ip,ib).ne.a(ip,ib+1)) allb=.false.
          end do
          if (allb) then 
            b1=1
            b2=1
          end if
        else
          b1=b
          b2=b
        end if

        do ib=b1,b2
          adj=ctlnb(CB_ADJ,ib)
          lprint=.false.
          lcom=.false.
          c=' ' 
          if ((ix.eq.PI_F).and.(myand(adj,1).ne.0)) lcom=.true.
          if ((ix.eq.PI_F0).and.(myand(adj,8).eq.0)) lcom=.true.
          if ((ix.eq.PI_RHO).and.(myand(adj,8).ne.0)) lcom=.true.
          if ((ix.eq.PI_BETA).and.(myand(adj,16).ne.0)) lcom=.true.
          if ((ix.eq.PI_GAMA).and.(myand(adj,16).ne.0)) lcom=.true.
          if ((ix.eq.PI_ANGZ).and.(myand(adj,16).eq.0)) lcom=.true.
          if ((ix.eq.PI_ANGX).and.(myand(adj,16).eq.0)) lcom=.true.
          if (myand(ipr,7).eq.1) then!going through various adjf cases
            do itop=1, ctlint(C_NTOP)
              i=ip+(itop-1)*DIMPIR
              if (ifit(i,ib).ne.0) lprint=.true.
            end do
          end if
          if (myand(ipr,7).eq.2) then
            do itop=1, ctlint(C_NTOP)
              i=ip+(itop-1)*DIMPIR
              if (ifit(i,ib).ne.0) lprint=.true.
            end do
            if (lcom) lprint=.true.
          end if
          if (myand(ipr,7).eq.3) then
            do itop=1, ctlint(C_NTOP)
              i=ip+(itop-1)*DIMPIR
              if (a(i,ib).ne.0.0d0) lprint=.true.
            end do
          end if
          if (myand(ipr,7).eq.4) then
            lprint=.true.
            if (parstr(ip)(1:1).eq.'_') lprint=.false.
            lx=.false.
            do itop=1, ctlint(C_NTOP)
              i=ip+(itop-1)*DIMPIR
              if (((parstr(ip)(1:1).eq.'D').or.(parstr(ip)(1:1).eq.'d'))
     $             .and.(a(i,ib).eq.0.0d0)) lx=.true.
            end do
            if (lx) lprint=.false.
          end if
          if (myand(ipr,7).ge.5) then
            lprint=.true.
          end if
          
          if (lprint) then
            call fillsp(pistr)
            do il=1,len(parstr(ip))
              if (parstr(ip)(il:il).eq.'_') goto 10
              pistr(il:il)=parstr(ip)(il:il)
            end do
 10         continue
            if (.not.allb) then
              write(pstr,'(2A,I1,A)') pistr(1:len_c(pistr)),'(',ib,')'
            else
              pstr=pistr
            end if
            if (lcom) c='\\'
            write(*,'(A,A10,$)') c,pstr
            do itop=1, ctlint(C_NTOP)
              i=ip+(itop-1)*DIMPIR
              call pra_out(a(i,ib),da(i,ib),ifit(i,ib),ierr,lcom,xstr)
              write(*,'(A,$)') xstr
            end do
            write(*,*)
          end if
        end do
      end do  
      return
      end

C ---------------------------------------------------------------------
      subroutine pra_out(a,d,ifi,ierr,lcom,xstr)
      implicit none
      real*8 a,d
      integer ifi,ierr
      logical lcom
      character*(*) xstr
      character*15 astr
      character*12 dstr
      integer aexp,myand
      external myand

      call fillsp(xstr)
      call pra_f(a,d,astr,dstr)
      if (ifi.eq.0) then
        if (lcom) then
          write(dstr,'(A12)') ' derived'
        else
          write(dstr,'(A12)') ' fixed  '
        end if
      end if
      if (myand(ierr,1).ne.0) then
        write(xstr,110) astr,'{',dstr,'}'
      else
        write(xstr,130) astr,'            '
      end if

 110  format(A15,1X,A1,A12,A1)
 130  format(A15,A5)

      return
      end

C ---------------------------------------------------------------------
      subroutine pra_f(a,d,astr,dstr)
      implicit none
      real*8 a,d
      character*(*) astr
      character*(*) dstr
      integer aexp

      if (dabs(a).gt.0.0d0) then
        aexp=int(log10(dabs(a)))
      else
        aexp=-9
      end if
      if (aexp.ge.3) then
        write(astr,'(F15.6)') a
        write(dstr,'(F12.6)') d
      end if
      if ((aexp.lt.3).and.(aexp.ge.-3)) then
        write(astr,'(F15.9)') a
        write(dstr,'(F12.9)') d
      end if
      if ((aexp.lt.-3).and.(aexp.ge.-6)) then !Sven2023
        write(astr,'(F12.6,A)') a*1.d6,'E-6' !Sven2023
        write(dstr,'(F9.6,A)')  d*1.d6,'E-6' !Sven2023
      end if                                 !Sven2023
      if ((aexp.lt.-6).and.(aexp.ge.-9)) then !Sven2023
        write(astr,'(F12.6,A)') a*1.d9,'E-9' !Sven2023
        write(dstr,'(F9.6,A)')  d*1.d9,'E-9' !Sven2023
      end if                                   !Sven2023
      if (aexp.lt.-9) then                    !Sven2023
        write(astr,'(F12.6,A)') a*1.d9,'E-9' !Sven2023
        write(dstr,'(F9.6,A)')  d*1.d9,'E-9' !Sven2023
      end if
      return
      end

C ---------------------------------------------------------------------
      subroutine ltxpr(rofit)
      implicit none
      real*8  rofit
      include 'iam.fi'
      integer i
      integer ii,ref,avg,fp4x
      real*8  reffreq,refcalc,fcalc,ffreq
      character*1 bsl,com
      character*40 fp1str
      character*40 fp2str
      character*40 fp3str
      integer myand
      external myand
      include 'iamdata.fi'

      fp4x=0
      do i=1, size(S_G)
        do ii=1, len(qnostr(1)) 
          if (qnostr(MAXQC+i)(ii:ii).ne.' ') then
            if (fp4x.lt.ii) fp4x=ii
          end if
        end do
      end do

      bsl='\\'
      
      write(*,*) bsl,'tt'
      write(*,*) bsl
     $     ,'begin{longtable}{r@{',bsl,' }r@{',bsl,' }r c r@{'
     $     ,bsl,' }r@{',bsl,' }r l r r r}' 
      write(*,*) bsl,'multicolumn{3}{c}{$J''$ $K_-''$ $K_+''$} &'
      write(*,*) ' - &'
      write(*,*) bsl,'multicolumn{3}{c}{$J$ $K_-$ $K_+$} &'
      write(*,*)
     $     ' Sym. &'
     $     ,' $',bsl,'nu$(obs.) &'
c     $     ,' $',bsl,'Delta',bsl,'nu$(obs.) &'
     $     ,' Err.  &'
     $     ,' $',bsl,'delta$ & ',bsl,bsl
      do i=1, ctlint(C_NDATA)
        com=' '
        if (qlin(i,Q_REF,Q_UP).eq.i) then
          write(fp1str,'(1X,3(I3,A),A,3(I3,A))') 
     $      qlin(i,Q_J,Q_UP),' &'
     $     ,qlin(i,Q_TJ,Q_UP)/2,' &'
     $     ,(2*qlin(i,Q_J,Q_UP)+2-qlin(i,Q_TJ,Q_UP))/2,' &'
     $     ,' - &'  
     $     ,qlin(i,Q_J,Q_LO),' &'
     $     ,qlin(i,Q_TJ,Q_LO)/2,' &'
     $     ,(2*qlin(i,Q_J,Q_LO)+2-qlin(i,Q_TJ,Q_LO))/2,' &'
        else
          write(fp1str,'(1X,3(3X,A),A,3(3X,A))') 
     $         ' &',' &',' &',' - &',' &',' &',' &'
        end if

        avg=qlin(i,Q_AVG,Q_UP)
        fcalc=(dnv(avg,NV_ENG,Q_UP)-dnv(avg,NV_ENG,Q_LO)
     $        +dnv(i,NV_ENG,Q_UP)-dnv(i,NV_ENG,Q_LO))/2.0d0
        ffreq=dln(i,LN_FREQ)
        
        ref=qlin(i,Q_REF,Q_UP)
        refcalc=0.0d0
        reffreq=0.0d0
        if ((ref.ne.0).and.(ref.ne.i)) then
          refcalc=dnv(ref,NV_ENG,Q_UP)-dnv(ref,NV_ENG,Q_LO)
          reffreq=dln(ref,LN_FREQ)
        end if

        if (fp4x.ge.1) then
          call fillsp(fp3str)
          if ((avg.ne.0).and.(avg.ne.i)) then
            write(fp2str,'(A,A,A,A)')
     $           qnostr(MAXQC+qlin(avg,Q_S,Q_UP))(2:fp4x),'+'
     $           ,qnostr(MAXQC+qlin(i,Q_S,Q_UP))(2:fp4x),' &'
          else
            write(fp2str,'(A,A,A,A)')
     $           qnostr(MAXQC+qlin(i,Q_S,Q_UP))(2:fp4x)
     $           ,' ',fp3str(2:fp4x),' &'
          end if
        else
          if ((avg.ne.0).and.(avg.ne.i)) then
            write(fp2str,'(I2,A2,I2,A)')
     $           qlin(avg,Q_S,Q_UP),' +',qlin(i,Q_S,Q_UP),' &'
          else
            write(fp2str,'(I6,A)')
     $           qlin(i,Q_S,Q_UP),' &'
          end if
        end if

        if ((ref.ne.0).and.(ref.ne.i)) then
c          write(fp3str,'(F12.3,A,F8.3,A,F8.3,A,2A)')
c     $         ffreq*1000.0d0,' &'
c     $         ,(ffreq-reffreq)*1000.0d0,' &'
c     $         ,(ffreq-reffreq
c     $           -fcalc+refcalc)*1000.0d3,' &',bsl,bsl
          write(fp3str,'(F12.3,A,F6.3,A,F8.3,A,2A)')
     $         (ffreq-reffreq)*1000.0d0,' &'
     $         ,dln(i,LN_ERR)*1000.0d0,' &'
     $         ,(ffreq-reffreq
     $           -fcalc+refcalc)*1000.0d3,' &',bsl,bsl
        else
c          write(fp3str,'(F12.3,A,8X,A,F8.3,A,2A)')
c     $         ffreq*1000.0d0,' &',' &'
c     $         ,(ffreq-fcalc)*1000.0d3,' &',bsl,bsl
          write(fp3str,'(F12.3,A,F6.3,A,F8.3,A,2A)')
     $         ffreq*1000.0d0,' &'
     $         ,dln(i,LN_ERR)*1000.0d0,' &'
     $         ,(ffreq-fcalc)*1000.0d3,' &',bsl,bsl
        end if
        
        com=' '
        if (dln(i,LN_ERR).eq.NOFIT) com='%'        
        write(*,'(4A)') com,fp1str(1:36),fp2str(1:14),fp3str(1:40)
      end do
      write(*,'(A)')'\\end{longtable}'
      
      return
      end 

C----------------------------------------------------------------------
      subroutine prderiv(dfit,nfit)
C     copy of the routine funcs
      implicit none
      include 'iam.fi'
      integer nfit
      integer dfit(DIMFIT)
      real*8  df,dfor,dfcr,dfobs,dfcal
      real*8  dfdao(DIMFIT),dfda(DIMFIT)
      real*8  dupda(DIMFIT),dloda(DIMFIT)
      integer avg,i,ix,ref
      external myand,myor
      integer myand,myor

      write(21,'(A33,$)')' No.   calc  obs-calc '
      do i=1, nfit
        if (dfit(i).gt.0) write(21,'(A16,$)') '-'
        if (dfit(i).lt.0) write(21,'(A16,$)') 'diff-quotient'
      end do
      write(21,*)

      do ix=1, ctlint(C_NDATA)
        dfor=0.0d0
        dfcr=0.0d0
        do  i=1, nfit
          dfdao(i)=0.0d0
          dfda(i)=0.0d0
        end do
        
        dfobs=dln(ix,LN_FREQ)

C     average of the data, normaly avg=ix: this has no effect 
        avg=qlin(ix,Q_AVG,Q_UP)
        if (avg.eq.0) avg = ix
        dfcal=((dnv(avg,NV_ENG,Q_UP)-dnv(avg,NV_ENG,Q_LO))
     $       +(dnv(ix,NV_ENG,Q_UP)-dnv(ix,NV_ENG,Q_LO)))/2.0d0
        do i=1, nfit
          dupda(i)=(dnv(avg,i+NV_DEF,Q_UP)+dnv(ix,i+NV_DEF,Q_UP))/2.0d0
          dloda(i)=(dnv(avg,i+NV_DEF,Q_LO)+dnv(ix,i+NV_DEF,Q_LO))/2.0d0
        end do
     
C     check if a difference is fitted 
        ref=qlin(ix,Q_REF,Q_UP)
        if ((ref.ne.0).and.(ref.ne.ix)) then
          dfor=dln(ref,LN_FREQ)
          dfcr=dnv(ref,NV_ENG,Q_UP)-dnv(ref,NV_ENG,Q_LO)
          do i=1, nfit
            dfdao(i)=-(dnv(ref,i+NV_DEF,Q_UP)-dnv(ref,i+NV_DEF,Q_LO))
          end do
        end if  
      
        df=dfobs-dfor-dfcal+dfcr
        do i=1, nfit
          dfda(i)=dfdao(i)+dupda(i)-dloda(i)
        end do

        write(21,'(I3,2F15.7,$)') ix,dfcal-dfcr,df
        do i=1, nfit
          if (dfit(i).ne.0) write(21,'(F16.8,$)') dfda(i)
        end do
        write(21,*)
      end do

      return
      end
      
