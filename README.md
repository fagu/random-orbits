This command line application generates cubic integral domains (or invertible ideal classes in quadratic integral domains) with bounded discriminant and given signature uniformly at random.

# Installation

## Dependencies

- Linux (although it should be not too difficult to port this to other operating systems)
- A recent C++ compiler.
- CMake
- pkgconf
- GMP library, version >= 6.2 (https://gmplib.org/)

Moreover, you either need:

- FLINT library, version >= 3.0 (https://flintlib.org/index.html)

Or:

- FLINT library, version >= 2.8 (https://flintlib.org/index.html)
- Arb library, version >= 2.22 (https://arblib.org/)

### Arch Linux

On Arch Linux, install the following packages:

```
base-devel cmake flint
```

### Ubuntu

On Ubuntu >= 22.04, install the following packages:

```
g++ cmake pkgconf libflint-dev libflint-arb-dev
```
(Once FLINT 3 appears in the Ubuntu package repositories, `libflint-arb-dev` will probably disappear and become unnecessary.)

## Building

```
mkdir build
cd build
cmake ..
make
```

# Usage

In the build directory, run `./random-cubic-ring` or `./random-quadratic-ideal-class` for help.

```
Usage: ./random-cubic-ring [options] N r T

This program generates N random cubic integral domains R with a given
signature r and |disc(R)| <= T.

The probability of generating a particular ring R is proportional to
1 / #Aut(R).

Output format:
  The output (stdout) contains one line
    a b c d
  per ring R, consisting of the integer coefficients of a binary cubic
  form a X^3 + b X^2 Y + c X Y^2 + d Y^3 corresponding to R in the Levi
  parameterization.

Parameters:
  N  Number of cubic rings to generate.
  r  Number of real embeddings (either 1 or 3).
  T  Upper bound on the absolute value of the discriminant.
     Instead of an integer, you can pass '-' as an argument.
     In that case, pass the bound T to the program through stdin.
     (This can be helpful if you want to use a very large number T.)

Options:
  --only-maximal   Only generate maximal orders.
                   Note: This can be slow because the discriminant needs
                         to be factored.
  --only-triv-aut  Only generate rings with trivial automorphism group.
  --uniform        Generate all rings with the same probability, instead
                   of with probability proportional to 1 / #Aut(R).
  --reduce         Reduce the cubic forms. Always show the same cubic
                   form for a cubic ring.
  --verbose        Print extra information to stderr.
  --seed [SEED]    Unsigned 32 bit integer to use as a seed for the
                   random number generator.
  --progress       Show progress and ETA.
  -h               Print this help message.
```

# Examples

Generate 100000 random cubic integral domains with discriminant in the interval [1,10^30]:
```
% ./random-cubic-ring 100000 3 1000000000000000000000000000000 --progress > out.txt
```

Generate 100000 random cubic integral domains with discriminant in the interval [-10^30,-1]:
```
% ./random-cubic-ring 100000 1 1000000000000000000000000000000 --progress > out.txt
```

Generate 100 random cubic integral domains with discriminant in the interval [1,10^100000]:
```
% python -c 'print("1"+"0"*100000)' | ./random-cubic-ring 100 3 - --progress > out.txt
```

The following examples are sanity checks recovering all fields of small discriminant.

Generate 10000 random cubic number fields with discriminant in the interval [1,200]. Each line corresponds to a cubic number field. The first number is the number of times this field was generated. The other four numbers are the coefficients of the cubic form. Note that the number field with trivial automorphism group is found about 3 times as often as the three extensions with automorphism group C3:
```
% ./random-cubic-ring 10000 3 200 --only-maximal --reduce | sort | uniq -c
   1627 1 0 -3 1
   1695 1 1 -2 -1
   4921 1 1 -3 -1
   1757 1 1 -4 1
```

The same, but with discriminants in the interval [-100,-1]. All such number fields have trivial automorphism group:
```
% ./random-cubic-ring 10000 1 100 --only-maximal --reduce | sort | uniq -c
   1366 1 0 1 1
   1461 1 0 2 1
   1417 1 1 1 2
   1438 1 1 2 1
   1422 1 1 3 1
   1454 1 2 2 2
   1442 1 2 3 3
```

Find all 4804 cubic number fields with discriminant in the interval [1,100000]:
```
% ./random-cubic-ring 80000 3 100000 --only-maximal --reduce | sort | uniq | wc -l
4804
```

Find all 4753 cubic number fields with discriminant in the interval [1,100000] and trivial automorphism group:
```
% ./random-cubic-ring 50000 3 100000 --only-maximal --only-triv-aut --reduce | sort | uniq | wc -l
4753
```
