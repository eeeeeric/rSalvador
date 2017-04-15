rSalvador
=========

rSalvador: An R Tool for the Luria-Delbruck Fluctuation Assay

Installing from source
----------------------
```
R CMD check rsalvador
R CMD INSTALL rsalvador
```

Building on Windows - Environment Preparation
---------------------------------------------
Download and install R, rtools, Strawberry Perl, and MiKTeX. Strawberry Perl,
MiKTeX, and rtools will automatically add themselves to the PATH, but you will
need to add R manually. Note that MiKTeX is only necessary if you wish to build
documentation,

Open a R prompt, and install the hypergeo and gdata packages.
```
install.packages('hypergeo')
install.packages('gdata')
```

Now you are ready to build (see below).

Building an installable archive (all platforms)
-----------------------------------------------
This will check the source files (requires LaTeX), okay to skip
```
R CMD check rsalvador
```

This will produce the archive
```
R CMD INSTALL --build rsalvador
```

Author
------
Qi Zheng
