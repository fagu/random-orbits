#include <flint/fmpz.h>
#include <flint/fmpz_factorxx.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>

class fmpz_factorxx;

class fmpzxx {
public:
	fmpz_t inner;
	fmpzxx() {
		fmpz_init(inner);
	}
	fmpzxx(slong val) {
		fmpz_init(inner);
		fmpz_set_si(inner, val);
	}
	fmpzxx(const fmpzxx& x) {
		fmpz_init(inner);
		fmpz_set(inner, x.inner);
	}
	fmpzxx& operator=(const fmpzxx& x) {
		fmpz_set(inner, x.inner);
		return *this;
	}
	~fmpzxx() {
		fmpz_clear(inner);
	}
	fmpz_t& _fmpz() {
		return inner;
	}
	// 2^e
	static fmpzxx exp2(ulong e) {
		fmpzxx res = 1;
		fmpz_mul_2exp(res.inner, res.inner, e);
		return res;
	}
	friend fmpzxx& operator*=(fmpzxx& a, const fmpzxx& b) {
		fmpz_mul(a.inner, a.inner, b.inner);
		return a;
	}
	void set_bit(ulong i) {
		fmpz_setbit(inner, i);
	}
	friend fmpzxx operator+(const fmpzxx& a, const fmpzxx& b) {
		fmpzxx res;
		fmpz_add(res.inner, a.inner, b.inner);
		return res;
	}
	friend fmpzxx operator-(const fmpzxx& a, const fmpzxx& b) {
		fmpzxx res;
		fmpz_sub(res.inner, a.inner, b.inner);
		return res;
	}
	friend fmpzxx operator*(const fmpzxx& a, const fmpzxx& b) {
		fmpzxx res;
		fmpz_mul(res.inner, a.inner, b.inner);
		return res;
	}
	friend fmpzxx operator-(const fmpzxx& a) {
		fmpzxx res;
		fmpz_neg(res.inner, a.inner);
		return res;
	}
	friend fmpzxx cdiv_q(const fmpzxx& a, const fmpzxx& b) {
		fmpzxx res;
		fmpz_cdiv_q(res.inner, a.inner, b.inner);
		return res;
	}
	fmpzxx abs() const {
		fmpzxx res;
		fmpz_abs(res.inner, inner);
		return res;
	}
	friend bool operator==(const fmpzxx& a, const fmpzxx& b) {
		return fmpz_equal(a.inner, b.inner);
	}
	friend bool operator!=(const fmpzxx& a, const fmpzxx& b) {
		return !fmpz_equal(a.inner, b.inner);
	}
	friend bool operator>(const fmpzxx& a, const fmpzxx& b) {
		return fmpz_cmp(a.inner, b.inner) > 0;
	}
	friend bool operator<=(const fmpzxx& a, const fmpzxx& b) {
		return fmpz_cmp(a.inner, b.inner) <= 0;
	}
	friend bool operator<(const fmpzxx& a, const fmpzxx& b) {
		return fmpz_cmp(a.inner, b.inner) < 0;
	}
	int sgn() const {
		return fmpz_sgn(inner);
	}
	bool is_square() const {
		return fmpz_is_square(inner);
	}
	// Returns the (nonnegative) square root of this number, assuming it has a square root in Z.
	fmpzxx sqrt() const {
		fmpzxx res;
		fmpz_sqrt(res.inner, inner);
		return res;
	}
	bool divisible_by(const fmpzxx& b) const {
		return fmpz_divisible(inner, b.inner);
	}
	// Factorization of the absolute value of this number, which is assumed to be nonzero.
	std::vector<std::pair<fmpzxx,ulong>> factor_abs() const {
		fmpz_factor_t fa;
		fmpz_factor_init(fa);
		fmpz_factor(fa, inner);
		std::vector<std::pair<fmpzxx,ulong>> res;
		res.reserve(fa->num);
		for (slong i = 0; i < fa->num; i++) {
			fmpzxx p;
			fmpz_set(p.inner, fa->p + i);
			res.emplace_back(p, fa->exp[i]);
		}
		fmpz_factor_clear(fa);
		return res;
	}
};

class fmpz_polyxx {
public:
	fmpz_poly_t inner;
	fmpz_polyxx() {
		fmpz_poly_init(inner);
	}
	fmpz_polyxx(const fmpz_polyxx& x) {
		fmpz_poly_init(inner);
		fmpz_poly_set(inner, x.inner);
	}
	fmpz_polyxx& operator=(const fmpz_polyxx& x) {
		fmpz_poly_set(inner, x.inner);
		return *this;
	}
	~fmpz_polyxx() {
		fmpz_poly_clear(inner);
	}
	slong length() const {
		return fmpz_poly_length(inner);
	}
	void set_coeff(slong n, const fmpzxx& c) {
		fmpz_poly_set_coeff_fmpz(inner, n, c.inner);
	}
	fmpzxx coeff(slong n) const {
		fmpzxx res;
		fmpz_poly_get_coeff_fmpz(res.inner, inner, n);
		return res;
	}
	fmpzxx disc() const {
		fmpzxx res;
		fmpz_poly_discriminant(res.inner, inner);
		return res;
	}
	fmpzxx content() const {
		fmpzxx res;
		fmpz_poly_content(res.inner, inner);
		return res;
	}
	fmpzxx operator()(const fmpzxx& x) const {
		fmpzxx res;
		fmpz_poly_evaluate_fmpz(res.inner, inner, x.inner);
		return res;
	}
	// Returns whether a given nonconstant polynomial is irreducible over Q.
	bool irreducible_over_Q() const {
		fmpz_poly_factor_t pfac;
		fmpz_poly_factor_init(pfac);
		fmpz_poly_factor(pfac, inner);
		// FLINT removes the content of the polynomial before factorization, so this is
		// effectively a factorization over Q.
		// Hence, f is irreducible over Q if and only if the number of factors is 1.
		bool irred = pfac->num == 1;
		fmpz_poly_factor_clear(pfac);
		return irred;
	}
	// Let p be a prime number and assume that the polynomial f is not divisible by p.
	// Returns the roots of f modulo p together with their multiplicities.
	std::vector<std::pair<fmpzxx, slong>> roots_mod(const fmpzxx& p) const {
		fmpz_mod_ctx_t ctx;
		fmpz_mod_ctx_init(ctx, p.inner);
		fmpz_mod_poly_t fm;
		fmpz_mod_poly_init(fm, ctx);
		for (slong i = 0; i < length(); i++)
			fmpz_mod_poly_set_coeff_fmpz(fm, i, coeff(i).inner, ctx);
		fmpz_mod_poly_factor_t fac;
		fmpz_mod_poly_factor_init(fac, ctx);
		fmpz_mod_poly_roots(fac, fm, 1, ctx);
		std::vector<std::pair<fmpzxx, slong>> res;
		res.reserve(fac->num);
		for (slong i = 0; i < fac->num; i++) {
			fmpzxx r;
			fmpz_mod_poly_get_coeff_fmpz(r.inner, fac->poly + i, 0, ctx);
			res.emplace_back(-r, fac->exp[i]);
		}
		fmpz_mod_poly_factor_clear(fac, ctx);
		fmpz_mod_poly_clear(fm, ctx);
		fmpz_mod_ctx_clear(ctx);
		return res;
	}
};
