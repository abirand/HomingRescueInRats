PROGRAMS = niche view

# Basic build types
all:
	@echo "!!"
	@echo "!! Building \"niche\" for `uname` on `arch`."
	@echo "!!"
	@make Darwin

debug:
	@echo "!!"
	@echo "!! Building \"niche\" for `uname` on `arch` in DEBUG mode."
	@echo "!!"
	@make clean
	@make Darwin_DEBUG
	@touch makefile

profile:
	@echo "!!"
	@echo "!! Building \"niche\" for `uname` on `arch` in PROFILE mode."
	@echo "!!"
	@make clean
	@make Darwin
	@touch makefile

# GNU/Linux settings
Linux:
	@\
	__NICHE_FLAGS="-O3 -Wall"  \
	__NICHE_LIBS="-lm -lpthread"         \
	__VIEW_FLAGS="-O3 -Wall"   \
	__VIEW_LIBS="-lm -lGL -lGLU"         \
	make `arch`

Linux_DEBUG:
	@\
	__NICHE_FLAGS="-gdwarf-2 -g3 -Wall"  \
	__NICHE_LIBS="-lm -lpthread"         \
	__VIEW_FLAGS="-gdwarf-2 -g3 -Wall"   \
	__VIEW_LIBS="-lm -lGL -lGLU"         \
	 make `arch`

Linux_PROFILE:
	@\
	__NICHE_FLAGS="-O3 -fno-inline -Wall"  \
	__NICHE_LIBS="-lm -lpthread"           \
	__VIEW_FLAGS="-w"                      \
	__VIEW_LIBS="-lm -lGL -lGLU"           \
	make `arch`

# SunOS/Unix settings
SunOS:
	@\
	__NICHE_FLAGS="-O3 -fmove-all-movables -fschedule-insns2 -w -includeieeefp.h"  \
	__NICHE_LIBS="-lm -lpthread -lrt"                                              \
	__VIEW_FLAGS="-w -I/usr/openwin/include/ -includeieeefp.h"                     \
	__VIEW_LIBS="-L/usr/openwin/lib/ -lm -lGL -lGLU -lX11"                         \
	make `arch`

SunOS_DEBUG:
	@\
	__NICHE_FLAGS="-gdwarf-2 -g3 -Wall -includeieeefp.h"         \
	__NICHE_LIBS="-lrt -lm -lpthread"                            \
	__VIEW_FLAGS="-gdwarf-2 -g3 -Wall -I/usr/openwin/include/"   \
	__VIEW_LIBS="-L/usr/openwin/lib/ -lm -lGL -lGLU -lX11"       \
	make `arch`

#  Mac/Apple/Darwin settings
Darwin:
	@\
	__NICHE_FLAGS="-O3 -fschedule-insns2 -fdelete-null-pointer-checks -finline-limit=100000 -w"  \
	__NICHE_LIBS="-lm -lpthread"                                                                                     \
	__VIEW_FLAGS="-w -I/opt/X11/include/ -I/opt/homebrew/Cellar/gsl/2.7.1/include/"                                 \
	__VIEW_LIBS="-L/opt/X11/lib/ -lm -lGL -lGLU -lX11"                                                               \
	make `arch`


Darwin_DEBUG:
	@\
	__NICHE_FLAGS="-g3 -Wall"                             \
	__NICHE_LIBS="-lm -lpthread"                          \
	__VIEW_FLAGS="-g3 -Wall  -I/opt/X11/include/ "        \
	__VIEW_LIBS="-L/opt/X11/lib/ -lm -lGL -lGLU -lX11"  \
	make `arch`


# Arch types (Just in case there is somthing specific that needs to go here)
i386:   
	@__ENDIAN="-DENDIAN=NICHE_LITTLE_ENDIAN" make compile
i486:   
	@__ENDIAN="-DENDIAN=NICHE_LITTLE_ENDIAN" make compile
i586:   
	@__ENDIAN="-DENDIAN=NICHE_LITTLE_ENDIAN" make compile
i686:   
	@__ENDIAN="-DENDIAN=NICHE_LITTLE_ENDIAN" make compile
x86_64: 
	@__ENDIAN="-DENDIAN=NICHE_LITTLE_ENDIAN" make compile
arm64:
	@__ENDIAN="-DENDIAN=NICHE_LITTLE_ENDIAN" make compile
sun4:
	@__ENDIAN="-DENDIAN=NICHE_BIG_ENDIAN"    make compile
ppc:    
	@__ENDIAN="-DENDIAN=NICHE_BIG_ENDIAN"    make compile

# Actually compile...
compile: $(PROGRAMS)

# parser for _CPUS
get_CPUS: get_CPUS.c parameters.h
	gcc -o get_CPUS get_CPUS.c

# niche
threads.o: get_CPUS thread.c shared.h random.h parameters.h makefile
	@echo gcc $$__NICHE_FLAGS -c thread.c $$__ENDIAN -DTHREAD=Thread_0 -o t.o;
	@gcc $$__NICHE_FLAGS -c thread.c $$__ENDIAN -DTHREAD=Thread_0 -o threads.o
	@a=`./get_CPUS`;                                                              \
	i=1;                                                                          \
	while [ "$$i" != "$$a" ] ; do                                                 \
	  echo gcc $$__NICHE_FLAGS -c thread.c $$__ENDIAN -DTHREAD=Thread_$$i -o t.o; \
	  gcc $$__NICHE_FLAGS -c thread.c $$__ENDIAN -DTHREAD=Thread_$$i -o t.o;      \
	  ld -r threads.o t.o -o tt.o;                                                \
	  rm -f threads.o;                                                            \
	  mv tt.o threads.o;                                                          \
	  i=`expr $$i + 1`;                                                           \
	done

niche.o: niche.c niche.h shared.h random.h parameters.h makefile
	@echo gcc $$__NICHE_FLAGS -DVERSION=\"`head -n 1 Changelog`\" -c niche.c
	@gcc $$__NICHE_FLAGS -DVERSION=\"`head -n 1 Changelog`\" -c niche.c

niche: niche.o threads.o random.o makefile
	@echo gcc $$__NICHE_FLAGS -o niche niche.o threads.o random.o $$__NICHE_LIBS
	@gcc $$__NICHE_FLAGS -o niche niche.o threads.o random.o $$__NICHE_LIBS
	@cp niche simulation

random.o: random.c random.h types.h makefile
	@echo gcc $$__VIEW_FLAGS -c random.c
	@gcc $$__VIEW_FLAGS -c random.c

# niche_view
view.o: view.c view.h ps.h makefile
	@echo gcc $$__VIEW_FLAGS -c view.c
	@gcc $$__VIEW_FLAGS -c view.c

view: view.o makefile
	@echo gcc $$__VIEW_FLAGS -o view view.o $$__LIBS $$__VIEW_LIBS
	@gcc $$__VIEW_FLAGS -o view view.o $$__LIBS $$__VIEW_LIBS

# maintenence
clean: 
	rm -f $(PROGRAMS) get_CPUS *.o

strip:
	rm -f $(PROGRAMS) get_CPUS *.o *~ *.par *.sys *.deme *.ind *.png *.eps \#*
