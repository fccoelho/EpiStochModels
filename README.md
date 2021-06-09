# EpiStochModels
Continuous time stochastic epidemic models implemented in D and wrapped in Python.

Check the [documentation](https://epistochmodels.readthedocs.io/en/latest/) for detailed instructions on how to install and use. You can also have a look at the example notebooks [here](docs/notebooks).


## Development

In order to work on `EpiStochModels` you need to have [DMD](https://dlang.org/download.html#dmd) or [LDC](https://github.com/ldc-developers/ldc#installation) dlang compilers installed. The you can activate the dlang environment in your shell:

```bash
source ~/dlang/ldc-1.26.0/activate.fish
```

To test the package:

```bash
dub test
```

