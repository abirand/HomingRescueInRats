#ifndef ARBITER_H
#define ARBITER_H

// Global
#define FD_IN       137
#define CLIENT_TTL  300
#define MAX_CMD     128
#define ARBITER_PROGRESS_RATE 60

// Interface
#ifndef ARBITER_C
extern void arbiter_init();
extern int  arbiter_progress(int progress);
extern int  arbiter_custom(char *c);
extern int  arbiter_exception(char *e);
extern int  arbiter_checkpoint(char *f);
extern int  arbiter_finished();
#endif

#endif
