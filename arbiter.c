#ifndef ARBITER_C
#define ARBITER_C

#include <pthread.h>
#include <unistd.h>
#include "arbiter.h"

static volatile char arbiter_frame_limit;

static void arbiter_eexit(char *s)
{
  printf("%s\n",s);
  exit(1);
}

static void *arbiter_frame_limit_thread(void *arg)
{
  while(1) {
    sleep(ARBITER_PROGRESS_RATE);
    arbiter_frame_limit = 0;
  }

  return NULL;
}

void arbiter_init()
{
  static int init;  pthread_attr_t a;  pthread_t t;  

  if( init ) {
    init = 1;
    return;
  }

  if( pthread_attr_init(&a) )                                    arbiter_eexit("arbiter: pthread_attr_init()");
  if( pthread_attr_setscope(&a,PTHREAD_SCOPE_SYSTEM) )           arbiter_eexit("arbiter: pthread_attr_setscope()");
  if( pthread_attr_setdetachstate(&a,PTHREAD_CREATE_DETACHED) )  arbiter_eexit("arbiter: pthread_attr_setdetachstate()");
  if( pthread_create(&t, &a, arbiter_frame_limit_thread, NULL) ) arbiter_eexit("arbiter: pthread_create()");
}

int arbiter_progress(int progress)
{
  char buf[5];

  // Check for frame lock
  if( arbiter_frame_limit )
    return 1;
  
  // Send progress (percentage)
  if( progress > 100 ) progress = 100;
  if( progress < 0   ) progress = 0;
  sprintf(buf,"p%d\n",progress);
  if( write(FD_IN, buf, strlen(buf)) != strlen(buf) )
    return 0;

  // Lock the frame
  arbiter_frame_limit = 1;
  return 1;
}

int arbiter_custom(char *c)
{
  // Send custom string
  if( strlen(c)+2 > MAX_CMD )                   return 0;
  if( write(FD_IN, "c", 1) != 1 )               return 0;
  if( write(FD_IN, c, strlen(c)) != strlen(c) ) return 0;

  return 1;
}

int arbiter_exception(char *e)
{
  // Send fatal exception
  if( strlen(e)+2 > MAX_CMD )                   return 0;
  if( write(FD_IN, "e", 1) != 1 )               return 0;
  if( write(FD_IN, e, strlen(e)) != strlen(e) ) return 0;

  return 1;
}

int arbiter_checkpoint(char *f)
{
  // Send checkpoint
  if( strlen(f)+2 > MAX_CMD )                   return 0;
  if( write(FD_IN, "k", 1) != 1 )               return 0;
  if( write(FD_IN, f, strlen(f)) != strlen(f) ) return 0;
  if( write(FD_IN, "\n", 1) != 1 )              return 0;

  return 1;
}

int arbiter_finished()
{
  // send progress over pipe
  if( write(FD_IN, "f\n", 2) != 2 ) return 0;
  return 1;
}

#endif // #ifndef ARBITER_C
