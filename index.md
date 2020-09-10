rSalvador
=========

rSalvador: An R Tool for the Luria-Delbruck Fluctuation Assay

Installation and Usage
----------------------
Follow the [guide](https://github.com/eeeeeric/rSalvador/raw/master/docs/guide.pdf) for installation instructions and a basic tutorial.

Download builds for 
- [Windows](https://github.com/eeeeeric/rSalvador/releases/download/v1.7/rsalvador_1.7.zip)
- [macOS](https://github.com/eeeeeric/rSalvador/releases/download/v1.7/rsalvador_1.7.tgz)
- [Linux](https://github.com/eeeeeric/rSalvador/releases/download/v1.7/rsalvador_1.7.tar)

The [rSalvador manual](https://github.com/eeeeeric/rSalvador/raw/master/docs/rsalvador-manual.pdf) is also available.

The example files mentioned in the guide: [example1.txt](https://github.com/eeeeeric/rSalvador/raw/master/example/example1.txt) and [example2.xlsx](https://github.com/eeeeeric/rSalvador/raw/master/example/example2.txt).

Installing from GitHub
----------------------

- Download and install R 3.6.1.
- Install Rtools 3.5 by downloading and executing the file Rtools35.exe.
- Launch R.
- Install the R package devtools by executing the following command from within R.
```install.packages('devtools')```
- Execute the following devtools command from within R.
```devtools::install_github("eeeeeric/rSalvador",subdir = "rsalvador")```

For Python 3 users
------------------
Some frequently used rSalvador functions are also available in pySalvador for Python users ([source](https://github.com/eeeeeric/rSalvador/blob/master/pysalvador/pysalvador.py), [manual](https://github.com/eeeeeric/rSalvador/blob/master/pysalvador/userManual.html), [Jupyter manual](https://github.com/eeeeeric/rSalvador/blob/master/pysalvador/userManual.ipynb)).

Author
------
Qi Zheng