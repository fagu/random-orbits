#ifndef PROGRESS_H
#define PROGRESS_H

#include <cstdio>
#include <ctime>

inline void clear_progress() {
	fprintf(stderr, "\033[2K\033[1G");
	fflush(stderr);
}

inline void show_progress(long long done, long long todo) {
	static long long firstclock;
	static long long firstdone = -1;
	static long long lastclock;
	if (firstdone == -1 || clock()-lastclock > CLOCKS_PER_SEC/5) {
		fprintf(stderr, "\033[2K\033[1G");
		fflush(stderr);
		fprintf(stderr, "%3d%% done, ETA", (int)((100ll*done)/todo));
		if (firstdone != -1 && done > firstdone) {
			double invspeed = (double)(clock()-firstclock)/CLOCKS_PER_SEC/(done-firstdone);
			int remsecs = invspeed*(todo-done);
			if (remsecs >= 3600)
				fprintf(stderr, " %dh", remsecs/3600);
			if (remsecs >= 60)
				fprintf(stderr, " %dm", (remsecs%3600)/60);
			fprintf(stderr, " %2ds", remsecs%60);
		} else {
			firstclock = clock();
			firstdone = done;
		}
		fflush(stderr);
		lastclock = clock();
	}
}

#endif
