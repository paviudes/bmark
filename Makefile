CFLAGS = -Wall
CC = cc
LIBS = -lm

benchmarking:obj/main.o obj/tiling.o obj/arrays.o obj/rgens.o obj/holes.o obj/dual.o obj/correlations.o obj/simulation.o obj/plot.o obj/memory.o obj/cell.o obj/grid.o obj/performance.o obj/load.o obj/save.o obj/regulartilings.o obj/homological.o obj/weights.o obj/notations.o obj/remote.o obj/noise.o obj/mcresult.o obj/rbim.o obj/submission.o obj/components.o obj/draw.o
	$(CC) $(CFLAGS) -o benchmarking $(LIBS) obj/main.o obj/tiling.o obj/arrays.o obj/rgens.o obj/holes.o obj/dual.o obj/correlations.o obj/simulation.o obj/plot.o obj/memory.o obj/cell.o obj/grid.o obj/performance.o obj/load.o obj/save.o obj/regulartilings.o obj/homological.o obj/weights.o obj/notations.o obj/remote.o obj/noise.o obj/mcresult.o obj/rbim.o obj/submission.o obj/components.o obj/draw.o

obj/main.o: source/main.c Makefile
	$(CC) $(CFLAGS) -c source/main.c -o obj/main.o

obj/tiling.o: source/tiling.c source/tiling.h Makefile
	$(CC) $(CFLAGS) -c source/tiling.c -o obj/tiling.o

obj/arrays.o: source/arrays.c source/arrays.h Makefile
	$(CC) $(CFLAGS) -c source/arrays.c -o obj/arrays.o

obj/draw.o: source/draw.c source/draw.h Makefile
	$(CC) $(CFLAGS) -c source/draw.c -o obj/draw.o

obj/rgens.o: source/rgens.c source/rgens.h Makefile
	$(CC) $(CFLAGS) -c source/rgens.c -o obj/rgens.o

obj/holes.o: source/holes.c source/holes.h Makefile
	$(CC) $(CFLAGS) -c source/holes.c -o obj/holes.o

obj/dual.o: source/dual.c source/dual.h Makefile
	$(CC) $(CFLAGS) -c source/dual.c -o obj/dual.o

obj/correlations.o: source/correlations.c source/correlations.h Makefile
	$(CC) $(CFLAGS) -c source/correlations.c -o obj/correlations.o

obj/simulation.o: source/simulation.c source/simulation.h Makefile
	$(CC) $(CFLAGS) -c source/simulation.c -o obj/simulation.o

obj/plot.o: source/plot.c source/plot.h Makefile
	$(CC) $(CFLAGS) -c source/plot.c -o obj/plot.o

obj/components.o: source/components.c source/components.h Makefile
	$(CC) $(CFLAGS) -c source/components.c -o obj/components.o

obj/memory.o: source/memory.c source/memory.h Makefile
	$(CC) $(CFLAGS) -c source/memory.c -o obj/memory.o

obj/cell.o: source/cell.c source/cell.h Makefile
	$(CC) $(CFLAGS) -c source/cell.c -o obj/cell.o

obj/grid.o: source/grid.c source/grid.h Makefile
	$(CC) $(CFLAGS) -c source/grid.c -o obj/grid.o

obj/noise.o: source/noise.c source/noise.h Makefile
	$(CC) $(CFLAGS) -c source/noise.c -o obj/noise.o

obj/mcresult.o: source/mcresult.c source/mcresult.h Makefile
	$(CC) $(CFLAGS) -c source/mcresult.c -o obj/mcresult.o

obj/rbim.o: source/rbim.c source/rbim.h Makefile
	$(CC) $(CFLAGS) -c source/rbim.c -o obj/rbim.o

obj/performance.o: source/performance.c source/performance.h Makefile
	$(CC) $(CFLAGS) -c source/performance.c -o obj/performance.o

obj/load.o: source/load.c source/load.h Makefile
	$(CC) $(CFLAGS) -c source/load.c -o obj/load.o

obj/save.o: source/save.c source/save.h Makefile
	$(CC) $(CFLAGS) -c source/save.c -o obj/save.o

obj/regulartilings.o: source/regulartilings.c source/regulartilings.h Makefile
	$(CC) $(CFLAGS) -c source/regulartilings.c -o obj/regulartilings.o

obj/homological.o: source/homological.c source/homological.h Makefile
	$(CC) $(CFLAGS) -c source/homological.c -o obj/homological.o

obj/weights.o: source/weights.c source/weights.h Makefile
	$(CC) $(CFLAGS) -c source/weights.c -o obj/weights.o

obj/notations.o: source/notations.c source/notations.h Makefile
	$(CC) $(CFLAGS) -c source/notations.c -o obj/notations.o

obj/remote.o: source/remote.c source/remote.h Makefile
	$(CC) $(CFLAGS) -c source/remote.c -o obj/remote.o

obj/submission.o: source/submission.c source/submission.h Makefile
	$(CC) $(CFLAGS) -c source/submission.c -o obj/submission.o

clean:
	$(RM) obj/main.o obj/tiling.o obj/arrays.o obj/rgens.o obj/holes.o obj/dual.o obj/correlations.o obj/simulation.o obj/plot.o obj/memory.o obj/cell.o obj/grid.o obj/performance.o obj/load.o obj/save.o obj/regulartilings.o obj/homological.o obj/weights.o obj/notations.o obj/remote.o obj/noise.o obj/mcresult.o obj/rbim.o obj/submission.o obj/components.o obj/draw.o
