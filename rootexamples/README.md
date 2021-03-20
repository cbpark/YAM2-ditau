In YAM2, run

```
make lib
```

to create the shared library, `libYAM2.so`. In macOS, it is `libYAM2.dylib`.

Modify the below lines in [`calcm2.C`](./calcm2.C) to adjust the paths to the shared libraries and the headers:

```
R__LOAD_LIBRARY(/usr/local/lib/libYAM2.so)
R__LOAD_LIBRARY(/usr/lib/libnlopt.so)
```

and

``` c++
    // For the include path of YAM2, where YAM2/yam2.h exists.
    gSystem->AddIncludePath("/usr/local/include");
```

In macOS, the shared libraries would have the `.dylib` suffix, not `.so`.

```
$ root -b -q calcm2.C
   ------------------------------------------------------------------
  | Welcome to ROOT 6.22/06                        https://root.cern |
  | (c) 1995-2020, The ROOT Team; conception: R. Brun, F. Rademakers |
  | Built for linuxx8664gcc on Nov 27 2020, 15:14:08                 |
  | From tags/v6-22-06@v6-22-06                                      |
  | Try '.help', '.demo', '.license', '.credits', '.quit'/'.q'       |
   ------------------------------------------------------------------


Processing calcm2.C...
M2 = 0.847415
where
  k1x: -1.79922, k1y: 0.30767, k1z: 0.560096
  k2x: 4.14272, k2y: 0.20543, k2z: 0.23148
```
