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
chrono::duration<double> time_random(0);
chrono::duration<double> time_irreducible(0);
chrono::duration<double> time_reduce(0);
long long fail_smax = 0;
long long fail_smin = 0;
long long fail_prob = 0;
long long fail_azero = 0;
long long fail_signature = 0;
long long fail_abs_disc = 0;
long long fail_in_ball = 0;
long long fail_irreducible = 0;
long long fail_aut = 0;
long long fail_uniform = 0;
long long fail_maximal = 0;

struct cubic_form {
	// The cubic form a X^3 + ... + d Y^3.
	fmpzxx a, b, c, d;
	// The discriminant of the cubic form.
	fmpzxx disc() const {
		return b*b*c*c - 4*a*c*c*c - 4*b*b*b*d - 27*a*a*d*d + 18*a*b*c*d;
	}
	// The polynomial f(X, 1) = a X^3 + ... + d.
	fmpz_polyxx polynomial() const {
		fmpz_polyxx f;
		f.set_coeff(3, a);
		f.set_coeff(2, b);
		f.set_coeff(1, c);
		f.set_coeff(0, d);
		return f;
	}
	// Assume f is irreducible over Q.
	// The automorphism group of the corresponding cubic ring R is either C_3 or trivial.
	// This function returns whether the automorphism group of R is trivial.
	bool is_triv_aut() const {
		// Write f = a X^3 + ... + d Y^3.
		// The ring has nontrivial automorphism group if and only if all of the following hold:
		//  - disc(f) = s^2 for some integer s.
		//  - s divides 3ac - b^2 and 3bd - c^2.
		fmpzxx di = disc();
		if (!di.is_square())
			return true;
		fmpzxx s = di.sqrt();
		if (!((3*a*c - b*b).divisible_by(s)))
			return true;
		if (!((3*b*d - c*c).divisible_by(s)))
			return true;
		return false;
	}
	
	// Assume a != 0 and disc(f) != 0.
	// This function returns whether the corresponding cubic ring is a maximal order.
	bool is_maximal() const {
		// The ring can be nonmaximal only at the primes p such that p^2 divides the discriminant of f.
		for (const auto& pe : disc().factor_abs()) {
			if (pe.second >= 2) {
				fmpzxx p = pe.first;
				// Write f = a X^3 + ... + d Y^3.
				// The ring is nonmaximal at p if and only if one of the following holds:
				//  - p divides f
				//  - p^2 divides a and p divides b.
				//  - There is a root [r] of the polynomial f(X,1) modulo p such that
				//    p^2 divides f(r,1) and f(r+p,1).
				//    This does not depend on the choice of the integer representative r of
				//    the residue class [r].
				if (a.divisible_by(p) && b.divisible_by(p) && c.divisible_by(p) && d.divisible_by(p))
					return false;
				if (a.divisible_by(p * p) && b.divisible_by(p))
					return false;
				fmpz_polyxx pol = polynomial();
				for (const auto& re : pol.roots_mod(p)) {
					fmpzxx r = re.first;
					assert(pol(r).divisible_by(p));
					if (pol(r).divisible_by(p * p) && pol(r + p).divisible_by(p * p))
						return false;
				}
			}
		}
		return true;
	}
	
	// f(Y, X)
	cubic_form swap_X_Y() const {
		return {d, c, b, a};
	}
	
	// f(X + k Y, Y)
	cubic_form add_kY_to_X(const fmpzxx& k) const {
		return {a, 3*a*k + b, 3*a*k*k + 2*b*k + c, a*k*k*k + b*k*k + c*k + d};
	}
	
	// f(-X, Y)
	cubic_form flip_X() const {
		return {-a, b, -c, d};
	}
	
	// f(X, -Y)
	cubic_form flip_Y() const {
		return {a, -b, c, -d};
	}
	
	// Assume f is irreducible over Q.
	// Replaces f by the unique reduced cubic form in the GL_2(Z)-orbit of f.
	void reduce() {
		// See Belabas, A fast algorithm to compute cubic fields.
		if (disc() > 0) {
			while(true) {
				fmpzxx P = b*b - 3*a*c;
				fmpzxx Q = b*c - 9*a*d;
				fmpzxx R = c*c - 3*b*d;
				if (a < 0) {
					*this = flip_X();
				} else if (b < 0 || (b == 0 && d < 0)) {
					*this = flip_Y();
				} else if (P > R || (P == R && (a > d.abs() || (a == d.abs() && b >= c.abs())))) {
					*this = swap_X_Y();
				} else if (Q > P) {
					*this = add_kY_to_X(-1);
				} else if (-Q > P) {
					*this = add_kY_to_X(+1);
				} else if (Q == 0 && d >= 0) {
					// In Definition 3.2, Belabas demands that if Q = 0, then d < 0.
					// This condition is redundant.
					abort();
				} else if (P == Q && b >= (3*a-b).abs()) {
					*this = add_kY_to_X(-1);
				} else {
					break;
				}
			}
		} else {
			while(true) {
				if (d*d - a*a + a*c - b*d <= 0) {
					*this = swap_X_Y();
				} else if (a < 0) {
					*this = flip_X();
				} else if (b < 0 || (b == 0 && d < 0)) {
					*this = flip_Y();
				} else if (a*d - b*c >= (a+b)*(a+b) + a*c) {
					*this = add_kY_to_X(+1);
				} else if (a*d - b*c <= -(a-b)*(a-b) - a*c) {
					*this = add_kY_to_X(-1);
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
	bool only_triv_aut;
	bool uniform;
	bool reduce;
};

template<class Generator>
optional<cubic_form> try_generate(const parameters& params, Generator& gen) {
	rand_real prob, randt, rands, randk1, randk2, randk3, randk4;
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
			randk4.refine(prec, gen);
		}
		arbxx sqrt3 = SQRT(3);
		arbxx sqrt5 = SQRT(5);
		// R = 5 / 4 for totally real fields
		// R = 7 / 4 for not totally real fields
		arbxx R = params.nr_real_embeddings == 3 ? DIV(5, 4) : DIV(7, 4);
		// lambda = R * T^(1/4)
		arbxx lambda = MUL(R, ROOT(params.T, 4));
		// smin = sqrt(sqrt(3) / 2)
		arbxx smin = SQRT(DIV(sqrt3, 2));
		// smax = (lambda / 2)^(1/3)
		arbxx smax = ROOT(DIV(lambda, 2), 3);
		// fprintf(stderr, "smin = %s, smax = %s\n", to_str(smin, prec).c_str(), to_str(smax, prec).c_str());
		// L1p = smax^3 + lambda
		arbxx L1p = ADD(POW_SI(smax, 3), lambda);
		// L2p = smax + sqrt(5) * lambda
		arbxx L2p = ADD(smax, MUL(sqrt5, lambda));
		// L3p = smin^-1 + sqrt(5) * lambda
		arbxx L3p = ADD(POW_SI(smin, -1), MUL(sqrt5, lambda));
		// L4p = smin^-3 + lambda
		arbxx L4p = ADD(POW_SI(smin, -3), lambda);
		// Lprod = L1p*L2p*L3p*L4p
		arbxx Lprod = MUL(MUL(L1p, L2p), MUL(L3p, L4p));
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
			// l1p = floor(1 + lambda * s^-3)
			fmpzxx l1p = to_int(FLOOR(ADD(1, MUL(lambda, POW_SI(s, -3)))));
			// l2p = floor(1 + sqrt(5) * lambda * s^-1)
			fmpzxx l2p = to_int(FLOOR(ADD(1, MUL(sqrt5, MUL(lambda, POW_SI(s, -1))))));
			// l3p = floor(1 + sqrt(5) * lambda * s^1)
			fmpzxx l3p = to_int(FLOOR(ADD(1, MUL(sqrt5, MUL(lambda, s)))));
			// l4p = floor(1 + lambda * s^3)
			fmpzxx l4p = to_int(FLOOR(ADD(1, MUL(lambda, POW_UI(s, 3)))));
			// lprod = l1p*l2p*l3p*l4p
			fmpzxx lprod(l1p*l2p*l3p*l4p);
			// We always have 0 <= lprod <= Lprod.
			assert(!certainly_negative(lprod));
			assert(!certainly_positive(SUB(lprod, Lprod)));
			// With probability 1 - lprod / Lprod:
			if (!lt(prob.get(), DIV(lprod, Lprod))) {
				// Start over.
				fail_prob++;
				return nullopt;
			}
			// Random integers 0 <= ki < lip for all i
			fmpzxx k1 = to_int(FLOOR(MUL(randk1.get(), l1p)));
			fmpzxx k2 = to_int(FLOOR(MUL(randk2.get(), l2p)));
			fmpzxx k3 = to_int(FLOOR(MUL(randk3.get(), l3p)));
			fmpzxx k4 = to_int(FLOOR(MUL(randk4.get(), l4p)));
			// amin = ceil(- l1p / 2)
			fmpzxx amin(cdiv_q(-l1p, 2));
			// a = amin + k1
			fmpzxx a(amin + k1);
			// bmin_float = - l2p / 2 + 3 * t * a
			arbxx bmin_float = neg(DIV(l2p, 2));
			bmin_float = ADD(bmin_float, MUL(MUL(3, t), a));
			// bmin = ceil(bmin_float)
			fmpzxx bmin = to_int(CEIL(bmin_float));
			// b = bmin + k2
			fmpzxx b(bmin + k2);
			// cmin_float = - l3p / 2 - 3 * t^2 * a + 2 * t * b
			arbxx cmin_float = neg(DIV(l3p, 2));
			cmin_float = SUB(cmin_float, MUL(MUL(3, POW_UI(t, 2)), a));
			cmin_float = ADD(cmin_float, MUL(MUL(2, t), b));
			// cmin = ceil(cmin_float)
			fmpzxx cmin = to_int(CEIL(cmin_float));
			// c = cmin + k3
			fmpzxx c(cmin + k3);
			// dmin_float = - l4p / 2 + t^3 * a - t^2 * b + t * c
			arbxx dmin_float = neg(DIV(l4p, 2));
			dmin_float = ADD(dmin_float, MUL(POW_UI(t, 3), a));
			dmin_float = SUB(dmin_float, MUL(POW_UI(t, 2), b));
			dmin_float = ADD(dmin_float, MUL(t, c));
			// dmin = ceil(dmin_float)
			fmpzxx dmin = to_int(CEIL(dmin_float));
			// d = dmin + k4
			fmpzxx d(dmin + k4);
			// f = a X^3 + ... + d Y^3
			cubic_form f = {a, b, c, d};
			if (a == 0) {
				// Start over.
				fail_azero++;
				return nullopt;
			}
			fmpzxx disc = f.disc();
			//fmpz_poly_discriminant(disc._fmpz(), f._poly());
			int disc_sgn = params.nr_real_embeddings == 3 ? 1 : -1;
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
			// ap = a * s^3
			arbxx ap = a;
			ap = MUL(ap, POW_UI(s, 3));
			// bp = (b - 3 * t * a) * s
			arbxx bp = b;
			bp = SUB(bp, MUL(MUL(3, t), a));
			bp = MUL(bp, s);
			// cp = (c + 3 * t^2 * a - 2 * t * b) * s^-1
			arbxx cp = c;
			cp = ADD(cp, MUL(MUL(3, POW_UI(t, 2)), a));
			cp = SUB(cp, MUL(MUL(2, t), b));
			cp = DIV(cp, s);
			// dp = (d - t^3 * a + t^2 * b - t * c) * s^-3
			arbxx dp = d;
			dp = SUB(dp, MUL(POW_UI(t, 3), a));
			dp = ADD(dp, MUL(POW_UI(t, 2), b));
			dp = SUB(dp, MUL(t, c));
			dp = DIV(dp, POW_UI(s, 3));
			// q = 5 * ap^2 + bp^2 + cp^2 + 5 * dp^2 + 2 * ap * cp + 2 * bp * dp
			arbxx q;
			q = MUL(5, POW_UI(ap, 2));
			q = ADD(q, POW_UI(bp, 2));
			q = ADD(q, POW_UI(cp, 2));
			q = ADD(q, MUL(5, POW_UI(dp, 2)));
			q = ADD(q, MUL(2, MUL(ap, cp)));
			q = ADD(q, MUL(2, MUL(bp, dp)));
			// If q >= sqrt(|disc|) * R^2:
			if (!lt(q, MUL(SQRT(absdisc), POW_UI(R, 2)))) {
				// Start over.
				fail_in_ball++;
				return nullopt;
			}
			assert(lt(ap, DIV(lambda, 2)));
			assert(lt(bp, DIV(MUL(sqrt5, lambda), 2)));
			assert(lt(cp, DIV(MUL(sqrt5, lambda), 2)));
			assert(lt(dp, DIV(lambda, 2)));
			// Note: The cubic form aX^3 + ... + dY^3 is irreducible if and only if a != 0 and the polynomial f(X,1) = aX^3 + ... + d is irreducible.
			if (!f.polynomial().irreducible_over_Q()) {
				// Start over.
				fail_irreducible++;
				return nullopt;
			}
			if (params.only_triv_aut && !f.is_triv_aut()) {
				fail_aut++;
				return nullopt;
			}
			if (params.uniform && !params.only_triv_aut && f.is_triv_aut() && uniform_int_distribution<int>(0,2)(gen) != 0) {
				// Fail with probability 2/3 if Aut(R) = 1.
				fail_uniform++;
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
cubic_form generate(const parameters& params, Generator& gen) {
	while(true) {
		optional<cubic_form> f = try_generate(params, gen);
		if (f.has_value())
			return f.value();
	}
}

void show_help(const char* program_name) {
	fprintf(stderr,
		"Usage: %s [options] N r T\n"
		"\n"
		"This program generates N random cubic integral domains R with a given\n"
		"signature r and |disc(R)| <= T.\n"
		"\n"
		"The probability of generating a particular ring R is proportional to\n"
		"1 / #Aut(R).\n"
		"\n"
		"Output format:\n"
		"  The output (stdout) contains one line\n"
		"    a b c d\n"
		"  per ring R, consisting of the integer coefficients of a binary cubic\n"
		"  form a X^3 + b X^2 Y + c X Y^2 + d Y^3 corresponding to R in the Levi\n"
		"  parameterization.\n"
		"\n"
		"Parameters:\n"
		"  N  Number of cubic rings to generate.\n"
		"  r  Number of real embeddings (either 1 or 3).\n"
		"  T  Upper bound on the absolute value of the discriminant.\n"
		"     Instead of an integer, you can pass '-' as an argument.\n"
		"     In that case, pass the bound T to the program through stdin.\n"
		"     (This can be helpful if you want to use a very large number T.)\n"
		"\n"
		"Options:\n"
		"  --only-maximal   Only generate maximal orders.\n"
		"                   Note: This can be slow because the discriminant needs\n"
		"                         to be factored.\n"
		"  --only-triv-aut  Only generate rings with trivial automorphism group.\n"
		"  --uniform        Generate all rings with the same probability, instead\n"
		"                   of with probability proportional to 1 / #Aut(R).\n"
		"  --reduce         Reduce the cubic forms. Always show the same cubic\n"
		"                   form for a cubic ring.\n"
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
	const char* program_name = argc > 0 ? argv[0] : "random-cubic";
	int only_maximal = 0;
	int only_triv_aut = 0;
	int uniform = 0;
	int reduce = 0;
	verbose = 0;
	seed = 471932630;
	progress = 0;
	static option long_options[] = {
		{"only-maximal", no_argument, &only_maximal, 1},
		{"only-triv-aut", no_argument, &only_triv_aut, 1},
		{"uniform", no_argument, &uniform, 1},
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
	params.only_triv_aut = only_triv_aut;
	params.uniform = uniform;
	params.reduce = reduce;
	if (optind >= argc || sscanf(argv[optind], "%lld", &nr_orbits) != 1 || !(nr_orbits >= 0))
		err_help(program_name);
	optind++;
	if (optind >= argc || sscanf(argv[optind], "%d", &params.nr_real_embeddings) != 1 || !(params.nr_real_embeddings == 1 || params.nr_real_embeddings == 3))
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
	if (params.nr_real_embeddings == 3 && params.T < 49) {
		fprintf(stderr, "There is no irreducible orbit with 0<disc<49.\n");
		exit(1);
	}
	if (params.nr_real_embeddings == 1 && params.T < 23) {
		fprintf(stderr, "There is no irreducible orbit with 0<-disc<23.\n");
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
		if (params.only_triv_aut)
			fprintf(stderr, "Only generating rings with trivial automorphism group.\n");
		if (params.reduce)
			fprintf(stderr, "Reducing the cubic forms.\n");
	}
	
	mt19937 gen(seed);
	
	for(long long roun = 0; roun < nr_orbits; roun++) {
		if (progress)
			show_progress(roun, nr_orbits);
		cubic_form f = generate(params, gen);
		fmpz_print(f.a.inner);
		printf(" ");
		fmpz_print(f.b.inner);
		printf(" ");
		fmpz_print(f.c.inner);
		printf(" ");
		fmpz_print(f.d.inner);
		printf("\n");
	}
	if (progress)
		clear_progress();
	if (verbose) {
		fprintf(stderr, "Maximum precision used: %ld\n", max_precision);
		fprintf(stderr, "Running times:\n");
		fprintf(stderr, "  generate random numbers: %lfs\n", time_random.count());
		fprintf(stderr, "  check irreducibility: %lfs\n", time_irreducible.count());
		fprintf(stderr, "  reduce cubic form: %lfs\n", time_reduce.count());
		fprintf(stderr, "Number of failures due to:\n");
		fprintf(stderr, "  s > smax: %lld\n", fail_smax);
		fprintf(stderr, "  s < (1-t^2)^(1/4): %lld\n", fail_smin);
		fprintf(stderr, "  probability: %lld\n", fail_prob);
		fprintf(stderr, "  a = 0: %lld\n", fail_azero);
		fprintf(stderr, "  wrong signature: %lld\n", fail_signature);
		fprintf(stderr, "  |disc| > T: %lld\n", fail_abs_disc);
		fprintf(stderr, "  f lying outside transformed ball: %lld\n", fail_in_ball);
		fprintf(stderr, "  reducibility: %lld\n", fail_irreducible);
		fprintf(stderr, "  automorphism group: %lld\n", fail_aut);
		fprintf(stderr, "  uniformness: %lld\n", fail_uniform);
		fprintf(stderr, "  maximality: %lld\n", fail_maximal);
	}
	return 0;
}
