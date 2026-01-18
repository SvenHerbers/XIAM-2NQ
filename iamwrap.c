
#include <stdio.h>
#include <time.h>
#include <signal.h>

#ifdef __GNUC__

void sig_func__(int * );

void thiss(int i) {
  sig_func__(&i);
}

c_signal__() {
  signal(SIGINT,thiss);
}

dateput_() {
  time_t now;
  now=time(NULL);
  puts(ctime(&now));
}

#elif __sgi 

void sig_func_(int * );

void thiss(int i) {
  sig_func_(&i);
}

void c_signal_() {
  signal(SIGINT,thiss);
}

void dateput_() {
  time_t now;
  now=time(NULL);
  puts(ctime(&now));
  fflush(stdout);
}

#else

void sig_func(int * );

void thiss(int i) {
  sig_func(&i);
}

void c_signal() {
  signal(SIGINT,thiss);
}

void dateput() {
  time_t now;
  now=time(NULL);
  puts(ctime(&now));
  fflush(stdout);
}
 
#endif
