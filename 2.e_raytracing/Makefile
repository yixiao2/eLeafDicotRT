all: main inner_check trace_recal trace_recalRGB trace_recalRGB_v3

main: main_all.o Global.o precalculate.o Trace_dicot.o
	gcc -g -Wall main_all.o Global.o precalculate.o Trace_dicot.o -o rt4abprofiletraining -lm

main_all.o: main_all.c Defs.h Global.h Trace.h precalculate.h
	gcc -g -Wall -c main_all.c

inner_check: inner_check.o Global.o precalculate.o Trace_dicot.o
	gcc -g -Wall inner_check.o Global.o precalculate.o Trace_dicot.o -o inner_check -lm

inner_check.o: inner_check.c Defs.h Global.h Trace.h precalculate.h
	gcc -g -Wall -c inner_check.c

trace_recal: trace_recal.o Global.o precalculate.o
	gcc -g -Wall trace_recal.o Global.o precalculate.o -o trace_recal -lm

trace_recal.o: trace_recal.c Defs.h Global.h precalculate.h
	gcc -g -Wall -c trace_recal.c

trace_recalRGB: trace_recalRGB.o Global.o precalculate.o
	gcc -g -Wall trace_recalRGB.o Global.o precalculate.o -o trace_recalRGB -lm

trace_recalRGB.o: trace_recalRGB.c Defs.h Global.h precalculate.h
	gcc -g -Wall -c trace_recalRGB.c

trace_recalRGB_v3: trace_recalRGB_v3.o Global.o precalculate.o
	gcc -g -Wall trace_recalRGB_v3.o Global.o precalculate.o -o trace_recalRGB_v3 -lm

trace_recalRGB_v3.o: trace_recalRGB_v3.c Defs.h Global.h precalculate.h
	gcc -g -Wall -c trace_recalRGB_v3.c	

Global.o: Global.c Defs.h
	gcc -g -Wall -c Global.c

precalculate.o: precalculate.c Defs.h Global.h Trace.h
	gcc -g -Wall -c precalculate.c

Trace_dicot.o: Trace_dicot.c Defs.h Global.h
	gcc -g -Wall -c Trace_dicot.c

clean:
	@echo "cleanning project"
	-rm *.o
	@echo "complete"
