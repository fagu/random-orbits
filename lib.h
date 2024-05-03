#include <chrono>
using namespace std;

// Measures the time between construction and destruction and adds it to the duration variable passed during construction.
class mytimer {
	chrono::duration<double> &delta;
	chrono::time_point<chrono::steady_clock> start;
public:
	mytimer(chrono::duration<double> &delta_) : delta(delta_) {
		start = chrono::steady_clock::now();
	}
	~mytimer() {
		chrono::time_point<chrono::steady_clock> end = chrono::steady_clock::now();
		delta += end - start;
	}
};
