# How to compile and link Libint

## 1. Introduction

The HFX modules in cp2k interface the C++ library Libint in order to calculate
the four center electron repulsion integrals. I you want to use that part of
the code, you have to compile and link in these libraries (`libint.a` and
`libderiv.a`).

## 2. How to compile Libint

Compilation and installation of Libint is relatively straightforward.
You can download the library from

http://sourceforge.net/projects/libint/files/v1-releases/libint-1.1.4.tar.gz/download

(Note: Do not use Libint 1.1.3. All what follows has been tested for
       version 1.1.4.)

Once you have unzipped the tarball, just execute the following commands:

```console
> ./configure --prefix='your_install_directory'
> make
> make install
```

In that case, the library is compiled with default settings for the angular
momentum (i.e. g-functions for energies, f-functions for derivatives).

If you want to use basis functions with higher angular momenta, you have to
provide this information to the configure script. For example:

```console
> ./configure --prefix='your_install_directory' --with-libint-max-am=5 --with-libderiv-max-am1=4
```

for h- and g-functions respectively.

In that case you have to specify two additional DFLAGS that tell cp2k
what is the maximum angular momentum that the libint libraries can deal
with. For h-functions energies and g-functions derivatives for example,
that is

-  `-D__LIBINT_MAX_AM=6`
-  `-D__LIBDERIV_MAX_AM1=5`

Note: This is `max_am + 1` as specified in the configure script. These values can
      also be found in `libint.h` and `lideriv.h` respectively.


## 3. How to link to Libint

If you want to link against the Libint library, you have to add `-D__LIBINT`
to the `DFLAGS` in the cp2k architecture file. If you try to run cp2k using
the Hartree-Fock module without specifying this flag, you will get a
run-time error message and cp2k will abort.


You only have to link in the Libint libraries in the correct
order via the arch file:

```
LIBS        = [...] \
              'path_to_libint_libs'/lib/libderiv.a \
              'path_to_libint_libs'/lib/libint.a \
              -lstdc++
```

`lstdc++` is needed if you use the GNU C++ compiler.
