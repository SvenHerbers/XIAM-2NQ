# Makefile XIAM  (H.Hartwig  Mai 1996)
       F77 = gfortran
  F77FLAGS = -O # -O2 # # -g -C # -O2 -funroll-loops -m486 -fexpensive-optimizations -fstrength-reduce # 
      SRCS = iam.f iamm.f iamv.f iamv2.f iamio.f iamint.f iamfit.f iamadj.f iamsys.f 
      OBJS = iam.o iamm.o iamv.o iamv2.o iamio.o iamint.o iamfit.o iamadj.o iamsys.o
    LIBOBJ = mgetx.o iamlib.o 
    LIBSRC = mgetx.f iamlib.f 
   EXENAME = XiamNQ2

       LOCAL_LIBS = -ldiv

  LOCAL_LIBS_PATH = -L../../lib -L../lib -L./lib

iam:     $(OBJS) $(LIBOBJ)
	$(F77) -static -o  $(EXENAME) $(OBJS) $(LIBOBJ) 

#iam:     $(OBJS) 
#	$(F77) -o $(EXENAME) $(OBJS) $(LOCAL_LIBS) $(LOCAL_LIBS_PATH)

iam.o :   iam.f iam.fi iamdata.fi  

iamio.o :  iamio.f iam.fi iamdata.fi 

iamint.o :  iamint.f iam.fi  

iamadj.o :  iamadj.f iam.fi  

iamm.o :  iamm.f iam.fi

iamv.o :  iamv.f iam.fi

iamv2.o : iamv2.f iam.fi

iamfit.o : iamfit.f iam.fi

mgetx.o : mgetx.f mgetx.fi


.f.o: 
	$(F77) $(F77FLAGS) -c $<

#    for SGI   
#iamv.o : iamv.f iam.fi
#	gfortran -c -SWP:=ON iamv.f


