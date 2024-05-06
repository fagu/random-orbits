#include <cstdio>
#include <cstring>
#include <optional>
#include <random>
#include <unistd.h>
#include <getopt.h>
#include "lib.h"
#include "arbxx.h"
#include "progress.h"
using namespace std;

slong max_precision = 0;
duration<double> time_random(0);
duration<double> time_reduce(0);
long long fail_smax = 0;
long long fail_smin = 0;
long long fail_prob = 0;
long long fail_azero = 0;
long long fail_signature = 0;
long long fail_abs_disc = 0;
long long fail_in_ball = 0;
long long fail_irreducible = 0;
long long fail_primitive = 0;
long long fail_maximal = 0;

struct quadratic_form {
	// The quadratic form a X^2 + b X Y + c Y^2.
	fmpzxx a, b, c;
	// The discriminant of the quadratic form.
	fmpzxx disc() const {
		return b*b - 4 * a * c;
	}
	// The polynomial f(X, 1) = a X^2 + b X + c.
	fmpz_polyxx polynomial() const {
		fmpz_polyxx f;
		f.set_coeff(2, a);
		f.set_coeff(1, b);
		f.set_coeff(0, c);
		return f;
	}
	
	// Whether the given nonzero quadratic form is primitive.
	bool is_primitive() const {
		fmpzxx content;
		fmpz_gcd3(content.inner, a.inner, b.inner, c.inner);
		return content == 1;
	}
	
	// Assume disc(f) != 0.
	// This function returns whether the corresponding quadratic ring is a maximal order.
	bool is_maximal() const {
		// The ring can be nonmaximal only at the primes p such that p^2 divides the discriminant of f.
		for (const auto& pe : disc().factor_abs()) {
			if (pe.second >= 2) {
				// The ring is maximal at p with p^2 | disc(f) if and only if p = 2 and disc(f)/4 = 2 or 3 mod 4.
				if (pe.first != 2) {
					return false;
				} else {
					fmpzxx r = remainder(divexact(disc(), 4), 4);
					if (r == 0 || r == 1)
						return false;
				}
			}
		}
		return true;
	}
	
	// - f(Y, X)
	quadratic_form swap_X_Y() const {
		return {-c, -b, -a};
	}
	
	// f(X + k Y, Y)
	quadratic_form add_kY_to_X(const fmpzxx& k) const {
		return {a, 2*a*k + b, a*k*k + b*k + c};
	}
	
	// - f(-X, Y)
	quadratic_form flip_X() const {
		return {-a, b, -c};
	}
	
	void do_step() {
		fmpzxx D = disc();
		*this = swap_X_Y();
		if (a < 0)
			*this = flip_X();
		if (a*a > D) {
			// printf("A\n");
			while(b > a)
				*this = add_kY_to_X(-1);
			while(b <= -a)
				*this = add_kY_to_X(+1);
		} else {
			// printf("B\n");
			while(b < 0)
				*this = add_kY_to_X(+1);
			// printf("B2\n");
			while((2*a + b) * (2*a + b) < D)
				*this = add_kY_to_X(+1);
			// printf("B3\n");
			while(b > 0 && b*b >= D)
				*this = add_kY_to_X(-1);
			// printf("B4\n");
		}
	}
	
	friend bool operator!=(const quadratic_form& f, const quadratic_form& g) {
		return !(f.a == g.a && f.b == g.b && f.c == g.c);
	}
	
	// Assume f is irreducible over Q.
	// Replaces f by the unique reduced quadratic form in the GL_2(Z)-orbit of f.
	void reduce() {
		if (disc() > 0) {
			fmpzxx D = disc();
			if (a < 0)
				*this = flip_X();
			while(c > 0 || b < (a+c).abs())
				do_step();
			quadratic_form F = *this, G = *this;
			do {
				if (make_tuple(F.a, F.b, F.c) < make_tuple(a, b, c))
					*this = F;
				F.do_step();
				G.do_step();
				G.do_step();
			} while(F != G);
		} else {
			while(true) {
				if (a < 0) {
					*this = flip_X();
				} else if (a > c || (a == c && b < 0)) {
					*this = swap_X_Y();
				} else if (b > a) {
					*this = add_kY_to_X(-1);
				} else if (b <= -a) {
					*this = add_kY_to_X(+1);
				} else {
					break;
				}
			}
		}
	}
};

struct parameters {
	int nr_real_embeddings;
	fmpzxx T; // discriminant bound
	bool only_maximal;
	bool reduce;
};

template<class Generator>
optional<quadratic_form> try_generate(const parameters& params, Generator& gen) {
	rand_real prob, randt, rands, randk1, randk2, randk3;
	for (slong prec = 20; ; prec *= 2) {
		if (max_precision < prec)
			max_precision = prec;
		// fprintf(stderr, "Precision = %ld\n", prec);
		{
			mytimer tim(time_random);
			prob.refine(prec, gen);
			randt.refine(prec, gen);
			rands.refine(prec, gen);
			randk1.refine(prec, gen);
			randk2.refine(prec, gen);
			randk3.refine(prec, gen);
		}
		arbxx sqrt3 = SQRT(3);
		arbxx sqrt3over8 = SQRT(DIV(3,8));
		arbxx sqrt3over2 = SQRT(DIV(3,2));
		// R = 5 / 4 for real fields
		// R = 2 for imaginary fields
		arbxx R = params.nr_real_embeddings == 2 ? DIV(5, 4) : DIV(2, 1);
		// lambda = R * T^(1/2)
		arbxx lambda = MUL(R, ROOT(params.T, 2));
		// smin = sqrt(sqrt(3) / 2)
		arbxx smin = SQRT(DIV(sqrt3, 2));
		// smax = (lambda * sqrt(3/8))^(1/2)
		arbxx smax = ROOT(MUL(lambda, sqrt3over8), 2);
		// fprintf(stderr, "smin = %s, smax = %s\n", to_str(smin, prec).c_str(), to_str(smax, prec).c_str());
		// L1p = smax^2 + sqrt(3/2) * lambda
		arbxx L1p = ADD(POW_SI(smax, 2), MUL(sqrt3over2, lambda));
		// L2p = 1 + 2 * lambda
		arbxx L2p = ADD(1, MUL(2, lambda));
		// L3p = smin^-2 + sqrt(3/2) * lambda
		arbxx L3p = ADD(POW_SI(smin, -2), MUL(sqrt3over2, lambda));
		// Lprod = L1p*L2p*L3p
		arbxx Lprod = MUL(MUL(L1p, L2p), L3p);
		try {
			// t = random number in [-1/2,1/2)
			arbxx t = SUB(randt.get(), DIV(1,2));
			// s = random number in [smin, infinity) with probability density proportional to s^-2 d*s.
			arbxx s = DIV(smin, SQRT(rands.get()));
			// If s >= smax:
			if (!lt(s, smax)) {
				// Start over.
				fail_smax++;
				return nullopt;
			}
			// If 1 - t^2 >= s^4:
			if (!lt(SUB(1, POW_UI(t, 2)), POW_UI(s, 4))) {
				// Start over.
				fail_smin++;
				return nullopt;
			}
			// l1p = floor(1 + sqrt(3/2) * lambda * s^-2)
			fmpzxx l1p = to_int(FLOOR(ADD(1, MUL(MUL(sqrt3over2, lambda), POW_SI(s, -2)))));
			// l2p = floor(1 + 2 * lambda)
			fmpzxx l2p = to_int(FLOOR(ADD(1, MUL(2, lambda))));
			// l3p = floor(1 + sqrt(3/2) * lambda * s^2)
			fmpzxx l3p = to_int(FLOOR(ADD(1, MUL(MUL(sqrt3over2, lambda), POW_SI(s, 2)))));
			// lprod = l1p*l2p*l3p
			fmpzxx lprod(l1p*l2p*l3p);
			// We always have 0 <= lprod <= Lprod.
			assert(!certainly_negative(lprod));
			assert(!certainly_positive(SUB(lprod, Lprod)));
			// With probability 1 - lprod / Lprod:
			if (!lt(prob.get(), DIV(lprod, Lprod))) {
				// Start over.
				fail_prob++;
				return nullopt;
			}
			// Random integers 0 <= ki < lip for i = 1,...,3
			fmpzxx k1 = to_int(FLOOR(MUL(randk1.get(), l1p)));
			fmpzxx k2 = to_int(FLOOR(MUL(randk2.get(), l2p)));
			fmpzxx k3 = to_int(FLOOR(MUL(randk3.get(), l3p)));
			// amin = ceil(- l1p / 2)
			fmpzxx amin(cdiv_q(-l1p, 2));
			// a = amin + k1
			fmpzxx a(amin + k1);
			// bmin_float = - l2p / 2 + 2 * t * a
			arbxx bmin_float = neg(DIV(l2p, 2));
			bmin_float = ADD(bmin_float, MUL(MUL(2, t), a));
			// bmin = ceil(bmin_float)
			fmpzxx bmin = to_int(CEIL(bmin_float));
			// b = bmin + k2
			fmpzxx b(bmin + k2);
			// cmin_float = - l3p / 2 - t^2 * a + t * b
			arbxx cmin_float = neg(DIV(l3p, 2));
			cmin_float = SUB(cmin_float, MUL(POW_UI(t, 2), a));
			cmin_float = ADD(cmin_float, MUL(t, b));
			// cmin = ceil(cmin_float)
			fmpzxx cmin = to_int(CEIL(cmin_float));
			// c = cmin + k3
			fmpzxx c(cmin + k3);
			// f = a*X^2 + b*X*Y + c*Y^2
			quadratic_form f = {a, b, c};
			if (a == 0) {
				// Start over.
				fail_azero++;
				return nullopt;
			}
			fmpzxx disc = f.disc();
			// fmpz_poly_discriminant(disc._fmpz(), f._poly());
			assert(disc == b*b - 4*a*c);
			int disc_sgn = params.nr_real_embeddings == 2 ? 1 : -1;
			if (disc.sgn() != disc_sgn) {
				// Wrong signature.
				// Start over.
				fail_signature++;
				return nullopt;
			}
			fmpzxx absdisc(disc.abs());
			if (absdisc > params.T) {
				// Start over.
				fail_abs_disc++;
				return nullopt;
			}
			// ap = a * s^2
			arbxx ap = a;
			ap = MUL(ap, POW_UI(s, 2));
			// bp = b - 2 * t * a
			arbxx bp = b;
			bp = SUB(bp, MUL(MUL(2, t), a));
			// cp = (c + t^2 * a - t * b) * s^-2
			arbxx cp = c;
			cp = ADD(cp, MUL(POW_UI(t, 2), a));
			cp = SUB(cp, MUL(t, b));
			cp = DIV(cp, POW_UI(s, 2));
			// q = 3 * ap^2 + bp^2 + 3 * cp^2 + 2 * ap * cp
			arbxx q;
			q = MUL(3, POW_UI(ap, 2));
			q = ADD(q, POW_UI(bp, 2));
			q = ADD(q, MUL(3, POW_UI(cp, 2)));
			q = ADD(q, MUL(2, MUL(ap, cp)));
			// If q >= |disc| * R^2:
			if (!lt(q, MUL(absdisc, POW_UI(R, 2)))) {
				// Start over.
				fail_in_ball++;
				return nullopt;
			}
			// f is irreducible over Q if and only if its discriminant is not a square.
			if (disc.is_square()) {
				// Start over.
				fail_irreducible++;
				return nullopt;
			}
			// If the polynomial is not primitive:
			if (!f.is_primitive()) {
				// Start over.
				fail_primitive++;
				return nullopt;
			}
			if (params.only_maximal && !f.is_maximal()) {
				fail_maximal++;
				return nullopt;
			}
			if (params.reduce) {
				mytimer tim(time_reduce);
				f.reduce();
			}
			return f;
		} catch (insufficient_precision) {
			// Increase precision.
			// fprintf(stderr, "Increasing precision.\n");
		}
	}
}

template<class Generator>
quadratic_form generate(const parameters& params, Generator& gen) {
	while(true) {
		optional<quadratic_form> f = try_generate(params, gen);
		if (f.has_value())
			return f.value();
	}
}

void show_help(const char* program_name) {
	fprintf(stderr,
		"Usage: %s [options] N r T\n"
		"\n"
		"This program generates N random pairs (R,C), where R is a quadratic\n"
		"integral domain with |disc(R)| <= T and C is an invertible ideal class\n"
		"of R.\n"
		"\n"
		"The probability of generating a particular pair (R,C) is proportional to\n"
		"Reg(R).\n"
		"\n"
		"Output format:\n"
		"  The output (stdout) contains one line\n"
		"    a b c\n"
		"  per pair (R,C), consisting of the integer coefficients of a binary\n"
		"  quadratic form a X^2 + b X Y + c Y^2 corresponding to (R,C) in the Gauss\n"
		"  parameterization.\n"
		"\n"
		"Parameters:\n"
		"  N  Number of pairs (R,C) to generate.\n"
		"  r  Number of real embeddings (either 0 or 2).\n"
		"  T  Upper bound on the absolute value of the discriminant.\n"
		"     Instead of an integer, you can pass '-' as an argument.\n"
		"     In that case, pass the bound T to the program through stdin.\n"
		"     (This can be helpful if you want to use a very large number T.)\n"
		"\n"
		"Options:\n"
		"  --only-maximal   Only generate maximal orders.\n"
		"                   Note: This can be slow because the discriminant needs\n"
		"                         to be factored.\n"
		"  --reduce         Reduce the quadratic forms. Always show a canonical\n"
		"                   quadratic form for each pair (R,C).\n"
		"                   Note: This can be slow for r = 2 (real quadratic\n"
		"                         number fields).\n"
		"  --verbose        Print extra information to stderr.\n"
		"  --seed [SEED]    Unsigned 32 bit integer to use as a seed for the\n"
		"                   random number generator.\n"
		"  --progress       Show progress and ETA.\n"
		"  -h               Print this help message.\n",
		program_name);
}

[[noreturn]] void err_help(const char* program_name) {
	show_help(program_name);
	exit(1);
}

long long nr_orbits;
parameters params;
int verbose;
unsigned int seed;
int progress;

void parse_args(int argc, char **argv) {
	const char* program_name = argc > 0 ? argv[0] : "random-quadratic-ideal-class";
	int only_maximal = 0;
	int reduce = 0;
	verbose = 0;
	seed = 471932630;
	progress = 0;
	static option long_options[] = {
		{"only-maximal", no_argument, &only_maximal, 1},
		{"reduce", no_argument, &reduce, 1},
		{"verbose", no_argument, &verbose, 1},
		{"seed", required_argument, 0, 's'},
		{"progress", no_argument, &progress, 1},
		{0, 0, 0, 0}
	};
	int c;
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
		switch (c) {
		case 0:
			break;
		case 's': {
			if (sscanf(optarg, "%u", &seed) != 1)
				err_help(program_name);
			break;
		}
		case 'h': {
			show_help(program_name);
			exit(0);
			break;
		}
		case '?': {
			err_help(program_name);
			break;
		}
		default:
			abort();
		}
	}
	params.only_maximal = only_maximal;
	params.reduce = reduce;
	if (optind >= argc || sscanf(argv[optind], "%lld", &nr_orbits) != 1 || !(nr_orbits >= 0))
		err_help(program_name);
	optind++;
	if (optind >= argc || sscanf(argv[optind], "%d", &params.nr_real_embeddings) != 1 || !(params.nr_real_embeddings == 0 || params.nr_real_embeddings == 2))
		err_help(program_name);
	optind++;
	if (optind >= argc)
		err_help(program_name);
	if (strcmp(argv[optind], "-") == 0) {
		if (fmpz_fread(stdin, params.T._fmpz()) <= 0)
			err_help(program_name);
	} else {
		if (fmpz_set_str(params.T._fmpz(), argv[optind], 10) != 0)
			err_help(program_name);
	}
	optind++;
	if (params.nr_real_embeddings == 2 && params.T < 5) {
		fprintf(stderr, "There is no irreducible orbit with 0<disc<5.\n");
		exit(1);
	}
	if (params.nr_real_embeddings == 0 && params.T < 3) {
		fprintf(stderr, "There is no irreducible orbit with 0<-disc<3.\n");
		exit(1);
	}
	if (optind != argc)
		err_help(program_name);
}

int main(int argc, char **argv) {
	parse_args(argc, argv);
	
	if (verbose) {
		fprintf(stderr, "N = %lld\n", nr_orbits);
		fprintf(stderr, "r = %d\n", params.nr_real_embeddings);
		fprintf(stderr, "T = ");
		fmpz_fprint(stderr, params.T.inner);
		fprintf(stderr, "\n");
		fprintf(stderr, "seed = %u\n", seed);
		if (params.only_maximal)
			fprintf(stderr, "Only generating maximal orders.\n");
		if (params.reduce)
			fprintf(stderr, "Reducing the quadratic forms.\n");
	}
	
	mt19937 gen(seed);
	
	for(long long roun = 0; roun < nr_orbits; roun++) {
		if (progress)
			show_progress(roun, nr_orbits);
		quadratic_form f = generate(params, gen);
		fmpz_print(f.a.inner);
		printf(" ");
		fmpz_print(f.b.inner);
		printf(" ");
		fmpz_print(f.c.inner);
		printf("\n");
	}
	if (progress)
		clear_progress();
	if (verbose) {
		fprintf(stderr, "Maximum precision used: %ld\n", max_precision);
		fprintf(stderr, "Running times:\n");
		fprintf(stderr, "  generate random numbers: %lfs\n", time_random.count());
		fprintf(stderr, "  reduce quadratic form: %lfs\n", time_reduce.count());
		fprintf(stderr, "Number of failures due to:\n");
		fprintf(stderr, "  s > smax: %lld\n", fail_smax);
		fprintf(stderr, "  s < (1-t^2)^(1/4): %lld\n", fail_smin);
		fprintf(stderr, "  probability: %lld\n", fail_prob);
		fprintf(stderr, "  a = 0: %lld\n", fail_azero);
		fprintf(stderr, "  wrong signature: %lld\n", fail_signature);
		fprintf(stderr, "  |disc| > T: %lld\n", fail_abs_disc);
		fprintf(stderr, "  f lying outside transformed ball: %lld\n", fail_in_ball);
		fprintf(stderr, "  reducibility: %lld\n", fail_irreducible);
		fprintf(stderr, "  maximality: %lld\n", fail_maximal);
	}
	return 0;
}
