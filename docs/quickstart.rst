Getting started with EpiStochModels
================================================
EpistochModels can be installed from source (if you clone our `Github repository <https://github.com/fccoelho/EpiStochModels>`_).

But first you must have a **D** compiler installed. If your operating system does not have a pre-built package:

.. code-block:: bash

    $ curl https://dlang.org/install.sh | bash -s

Then you can build and install the Python package. 

.. code-block:: bash

    $ python setup.py build
    $ python setup.py install

We will publish the package on PyPI soon so that you will be able to simply pip install it.

Using the Built-in Models
----------------
The models can be used from Python  (see `notebooks <docs/notebooks>`_) or straight from D code:

.. literalinclude:: notebooks/sir_example.d
    :language: D
    :emphasize-lines: 26-27,29
    :linenos:

The *D* code above can be run with the following command:

.. code-block:: bash

    $ dub run --build=release --single sir_example.d

Implementing your own CTMC models
-----------------------------------
If you are willing to write `D` code, EpiStochModels also provides 
a class to simulate any CTMC model you want, as long as you provide 
it with a transition matrix and propensity functions.

The transition matrix encodes the discrete changes in each of the compartments (columns) when a specific transitions (lines) happens.
The propensity array is a 1-dimensional array conatining the probabilities of each transition happening. In EpistochModels, the elements of this array 
are functions of the parameters and state-variables of the model.


The code below is the implementation of a `SIRD model <https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIRD_model>`_.

.. mermaid::
   :caption: SIRD Model transitions:
   graph TD 
        S -- infection --> I
        I -- recovery --> R
        I -- death --> D

.. code-block:: D
    // define the transition matrix
    int[][] tmat = [[-1,1,0,0],
                    [0,-1,1,0], 
                    [0,-1,0,1]];

    // Define the parameter map
    double[string] pars = [
        "beta": 0.3,
        "gam": 0.05,
        "mu": 0.01
    ];

    // Define the propensity array as an array of Lambdas
    alias double delegate(int[], double[string]) dy;
    dy[] props = [
        (v,p)=> p["beta"]*v[0]*v[1]/sum(v[0..$-1]), // infection
        (v,p)=> p["gam"]*v[1], // recovery
        (v,p)=> p["mu"]*v[1] // death
        ];

    // Instantiate the model, initialize it and run.
    CTMC model = new CTMC(tmat, props);
    model.initialize([995,5,0,0],pars);
    auto res = model.run(0, 1000);
    // save the results
    model.save("sird.csv", res);