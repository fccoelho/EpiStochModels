# EpiStochModels
Continuous time stochastic epidemic models implemented in D. Python wrapper has been removed because PyD is no longer maintained.

So now models must be written in `D` like exemplified in [`app.d`](source/app.d) and on [`CTMC_example.d`](docs/notebooks/CTMC_example.d).  

D models writen followind the second example above, can be run like this:

```bash
dub run --build=release --single my_model.d
```

Check the [documentation](https://epistochmodels.readthedocs.io/en/latest/) for detailed instructions on how to install and use. You can also have a look at the example notebooks [here](docs/notebooks).


## Development

In order to work on `EpiStochModels` you need to have [DMD](https://dlang.org/download.html#dmd) or [LDC](https://github.com/ldc-developers/ldc#installation) dlang compilers installed. The you can activate the dlang environment in your shell:

```bash
source ~/dlang/ldc-1.27.1/activate.fish
```
the line above, is for fish shell, if you use bash, do this (to use the DMD compiler)

```bash
source ~/dlang/dmd-2.097.2/activate 
```

to deactivate the virtual environment run `deactivate` after you are done.

To test the package:

```bash
dub test
```

