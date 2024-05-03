This command line application generates cubic rings with bounded discriminant and given signature uniformly at random.

# Installation

## Dependencies

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

In the build directory, run `./random-cubic` or `./random-quadratic` for help.
