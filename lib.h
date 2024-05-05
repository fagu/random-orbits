#include <chrono>
using namespace std::chrono;

// Measures the time between construction and destruction and adds it to the duration variable passed during construction.
class mytimer {
	duration<double> &delta;
	time_point<steady_clock> start;
public:
	mytimer(duration<double> &delta_) : delta(delta_) {
		start = steady_clock::now();
	}
	~mytimer() {
		time_point<steady_clock> end = steady_clock::now();
		delta += end - start;
	}
};
