rSalvador
=========

Technical details are given in [Zheng (2017)](https://github.com/eeeeeric/rSalvador/blob/master/docs/G3Paper.pdf).

webSalvador
-----------
📣 For those who prefer the convenience of a web tool, webSalvador is a new app powered by rSalvador: [➡️ webSalvador](https://websalvador.eeeeeric.com)

Installation and Usage
----------------------
Follow the [guide](https://github.com/eeeeeric/rSalvador/raw/master/docs/guide.pdf) for installation instructions and a basic tutorial.

Download builds for 
- [Windows](https://github.com/eeeeeric/rSalvador/releases/download/v1.8/rsalvador_1.8.zip)
- [macOS](https://github.com/eeeeeric/rSalvador/releases/download/v1.8/rsalvador_1.8.tgz)
- [Linux](https://github.com/eeeeeric/rSalvador/releases/download/v1.8/rsalvador_1.8_R_x86_64-pc-linux-gnu.tar.gz)

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
```devtools::install_github("eeeeeric/rSalvador", subdir = "rsalvador")```

For Python 3 users
------------------
Some frequently used rSalvador functions are also available in pySalvador for Python users ([source](https://github.com/eeeeeric/rSalvador/blob/master/pysalvador/pysalvador.py), [Jupyter manual](https://github.com/eeeeeric/rSalvador/blob/master/pysalvador/userManual.ipynb)).

Estimating rates of non-neutral mutations
-----------------------------------------

- [kessler.py](https://github.com/eeeeeric/rSalvador/blob/master/python-examples/kessler.py)
- [stable_likely.py](https://github.com/eeeeeric/rSalvador/blob/master/python-examples/stable_likely.py)
- [confint.py](https://github.com/eeeeeric/rSalvador/blob/master/python-examples/confint.py)
- [simuKessler.py](https://github.com/eeeeeric/rSalvador/blob/master/python-examples/simuKessler.py)
- [meander.py](https://github.com/eeeeeric/rSalvador/blob/master/python-examples/meander.py)
- [naive.py](https://github.com/eeeeeric/rSalvador/blob/master/python-examples/naive.py)



Author
------
Qi Zheng
