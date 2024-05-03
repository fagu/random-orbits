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
long long fail_smax = 0;
long long fail_prob = 0;
long long fail_azero = 0;
long long fail_signature = 0;
long long fail_abs_disc = 0;
long long fail_in_ball = 0;
long long fail_irreducible = 0;
long long fail_primitive = 0;

template<class Generator>
optional<fmpz_polyxx> try_generate(fmpzxx T, int nr_real_embeddings, Generator& gen) {
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
		arbxx R = nr_real_embeddings == 2 ? DIV(5, 4) : DIV(2, 1);
		// lambda = R * T^(1/2)
		arbxx lambda = MUL(R, ROOT(T, 2));
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
			// f = a*X^2 + b*X + c
			fmpz_polyxx f;
			f.set_coeff(2, a);
			f.set_coeff(1, b);
			f.set_coeff(0, c);
			if (a == 0) {
				// Start over.
				fail_azero++;
				return nullopt;
			}
			fmpzxx disc = f.disc();
			// fmpz_poly_discriminant(disc._fmpz(), f._poly());
			assert(disc == b*b - 4*a*c);
			int disc_sgn = nr_real_embeddings == 2 ? 1 : -1;
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
			if (f.content() != 1) {
				// Start over.
				fail_primitive++;
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
	fprintf(stderr, "This will generate N independent random binary quadratic forms f with integer coefficients with exactly r real roots which are irreducible over Q and satisfy |disc(f)| <= T, with probability proportional to Reg(f).\n");
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
	if (nr_real_embeddings != 2 && nr_real_embeddings != 0)
		help(argc, argv);
	fmpzxx T;
	if (fmpz_set_str(T._fmpz(), argv[3], 10) != 0)
		help(argc, argv);
	if (T <= 0)
		help(argc, argv);
	mt19937 gen(42); // TODO
	
	// map<fmpzxx,long long> count;
	for(long long roun = 0; roun < nr_orbits; roun++) {
		// for (int att = 1; ; att++) {
		while(true) {
			optional<fmpz_polyxx> f = try_generate(T, nr_real_embeddings, gen);
			if (f.has_value()) {
				// cout << f.value().pretty("x") << endl;
				for (int i = 0; i < 3; i++) {
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
	fprintf(stderr, "Number of failures due to:\n");
	fprintf(stderr, "  s > smax: %lld\n", fail_smax);
	fprintf(stderr, "  probability: %lld\n", fail_prob);
	fprintf(stderr, "  a = 0: %lld\n", fail_azero);
	fprintf(stderr, "  wrong signature: %lld\n", fail_signature);
	fprintf(stderr, "  |disc| > T: %lld\n", fail_abs_disc);
	fprintf(stderr, "  f lying outside transformed ball: %lld\n", fail_in_ball);
	fprintf(stderr, "  reducibiliy: %lld\n", fail_irreducible);
	fprintf(stderr, "  nonprimitivity: %lld\n", fail_primitive);
	// fprintf(stderr, "Statistics:\n");
	// for (auto [d,c] : count) {
	// 	fprintf(stderr, " disc = %s: %lld times\n", d.to_string().c_str(), c);
	// }
	return 0;
}
