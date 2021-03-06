rSalvador
=========

rSalvador: An R Tool for the Luria-Delbruck Fluctuation Assay

📣 For those who prefer the convenience of a web tool, webSalvador is a new app powered by rSalvador: https://websalvador.eeeeeric.com/

Installing from GitHub
----------------------

- Download and install R 3.6.1 or later.
- Install Rtools 3.5 by downloading and executing the file Rtools35.exe.
- Launch R.
- Install the R package devtools by executing the following command from within R.
```install.packages('devtools')```
- Execute the following devtools command from within R.
```devtools::install_github("eeeeeric/rSalvador", subdir = "rsalvador")```

Installing from prebuilt packages
---------------------------------
`R CMD INSTALL <package>` where `<package>` is the appropriate download from
the [release page](https://github.com/eeeeeric/rSalvador/releases).

These releases are built using [GitHub Actions](.github/workflows/build.yml).

Installing from source
----------------------
```
R CMD check rsalvador
R CMD INSTALL rsalvador
```

For Python 3 users
------------------
Some frequently used rSalvador functions are also available in pySalvador for Python users ([source](./pysalvador/pysalvador.py), [Jupyter manual](./pysalvador/userManual.ipynb)).

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
