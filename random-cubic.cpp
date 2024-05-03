#include <cstdio>
#include <optional>
#include <random>
#include "lib.h"
#include "arbxx.h"
using namespace std;

#define ADD(a,b) add(a, b, prec)
#define SUB(a,b) sub(a, b, prec)
#define MUL(a,b) mul(a, b, prec)
#define DIV(a,b) div(a, b, prec)
#define POW_UI(a,b) pow_ui(a, b, prec)
#define POW_SI(a,b) pow_si(a, b, prec)
#define ROOT(a,b) root(a, b, prec)
#define INV(a) inv(a, prec)
#define SQRT(a) sqrt(a, prec)
#define LOG(a) log(a, prec)
#define FLOOR(a) floor(a, prec)
#define CEIL(a) ceil(a, prec)

slong max_precision = 0;
chrono::duration<double> time_random(0);
chrono::duration<double> time_irreducible(0);
long long fail_smax = 0;
long long fail_prob = 0;
long long fail_azero = 0;
long long fail_signature = 0;
long long fail_abs_disc = 0;
long long fail_in_ball = 0;
long long fail_irreducible = 0;

template<class Generator>
optional<fmpz_polyxx> try_generate(fmpzxx T, int nr_real_embeddings, Generator& gen) {
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
		arbxx R = nr_real_embeddings == 3 ? DIV(5, 4) : DIV(7, 4);
		// lambda = R * T^(1/4)
		arbxx lambda = MUL(R, ROOT(T, 4));
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
			// alpha = 1 / multiplicity in the Siegel set
			arbxx alpha = 1;
			// If s^4 - 2 / sqrt(3) * s^2 + t^2 < 0:
			if (lt(ADD(SUB(POW_UI(s, 4), MUL(DIV(2, sqrt3), POW_UI(s, 2))), POW_UI(t, 2)), 0)) {
				// We are in the non-unique part of the Siegel set.
				alpha = DIV(1, 2);
			}
			// With probability 1 - alpha * lprod / Lprod:
			if (!lt(prob.get(), MUL(alpha, DIV(lprod, Lprod)))) {
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
			// f = a*X^3 + b*X^2 + c*X + d
			fmpz_polyxx f;
			f.set_coeff(3, a);
			f.set_coeff(2, b);
			f.set_coeff(1, c);
			f.set_coeff(0, d);
			if (a == 0) {
				// Start over.
				fail_azero++;
				return nullopt;
			}
			fmpzxx disc = f.disc();
			//fmpz_poly_discriminant(disc._fmpz(), f._poly());
			int disc_sgn = nr_real_embeddings == 3 ? 1 : -1;
			if (disc.sgn() != disc_sgn) {
				// Wrong signature.
				// Start over.
				fail_signature++;
				return nullopt;
			}
			fmpzxx absdisc(disc.abs());
			if (absdisc > T) {
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
			// Note: The cubic form aX^3 + ... + dY^3 is irreducible if and only if a != 0 and the polynomial f = aX^3 + ... + d is irreducible.
			if (!f.irreducible_over_Q()) {
				// Start over.
				fail_irreducible++;
				return nullopt;
			}
			return f;
		} catch (insufficient_precision) {
			// Increase precision.
			// fprintf(stderr, "Increasing precision.\n");
		}
	}
}

void help(int argc, char** argv) {
	fprintf(stderr, "Usage: %s N r T\n", argc > 0 ? argv[0] : "enumerate");
	fprintf(stderr, "This will generate N independent random binary cubic forms f with integer coefficients with exactly r real roots which are irreducible over Q and satisfy |disc(f)| <= T, with probability proportional to 1/#Stab(f).\n");
	exit(1);
}

int main(int argc, char **argv) {
	if (argc != 4)
		help(argc, argv);
	long long nr_orbits;
	if (sscanf(argv[1], "%lld", &nr_orbits) != 1)
		help(argc, argv);
	if (nr_orbits < 0)
		help(argc, argv);
	int nr_real_embeddings;
	if (sscanf(argv[2], "%d", &nr_real_embeddings) != 1)
		help(argc, argv);
	if (nr_real_embeddings != 3 && nr_real_embeddings != 1)
		help(argc, argv);
	fmpzxx T;
	if (fmpz_set_str(T._fmpz(), argv[3], 10) != 0)
		help(argc, argv);
	if (T <= 0)
		help(argc, argv);
	if (nr_real_embeddings == 3 && T < 49) {
		fprintf(stderr, "There is no irreducible orbit with 0<disc<49.\n");
		exit(1);
	}
	if (nr_real_embeddings == 1 && T < 23) {
		fprintf(stderr, "There is no irreducible orbit with 0<-disc<23.\n");
		exit(1);
	}

	mt19937 gen(42); // TODO
	
	// map<fmpzxx,long long> count;
	for(long long roun = 0; roun < nr_orbits; roun++) {
		// for (int att = 1; ; att++) {
		while(true) {
			optional<fmpz_polyxx> f = try_generate(T, nr_real_embeddings, gen);
			if (f.has_value()) {
				//cout << f.value().pretty("x") << endl;
				for (int i = 0; i < 4; i++) {
					if (i)
						printf(" ");
					fmpz_print(f.value().coeff(i).inner);
				}
				printf("\n");
				// cout << " after " << att << " attempts." << endl;
				// fmpzxx disc;
				// fmpz_poly_discriminant(disc._fmpz(), f.value()._poly());
				// count[disc]++;
				break;
			}
		}
	}
	fprintf(stderr, "Maximum precision used: %ld\n", max_precision);
	fprintf(stderr, "Running times:\n");
	fprintf(stderr, "  generate random numbers: %lfs\n", time_random.count());
	fprintf(stderr, "  check irreducibility: %lfs\n", time_irreducible.count());
	fprintf(stderr, "Number of failures due to:\n");
	fprintf(stderr, "  s > smax: %lld\n", fail_smax);
	fprintf(stderr, "  probability: %lld\n", fail_prob);
	fprintf(stderr, "  a = 0: %lld\n", fail_azero);
	fprintf(stderr, "  wrong signature: %lld\n", fail_signature);
	fprintf(stderr, "  |disc| > T: %lld\n", fail_abs_disc);
	fprintf(stderr, "  f lying outside transformed ball: %lld\n", fail_in_ball);
	fprintf(stderr, "  reducibility: %lld\n", fail_irreducible);
	// fprintf(stderr, "Statistics:\n");
	// for (auto [d,c] : count) {
	// 	fprintf(stderr, " disc = %s: %lld times\n", d.to_string().c_str(), c);
	// }
	return 0;
}
