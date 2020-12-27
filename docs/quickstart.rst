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

Using the models
----------------
The models can be used from Python  (see `notebooks <docs/notebooks>`_) or straight from D code:

.. literalinclude:: notebooks/sir_example.d
    :language: D
    :emphasize-lines: 26-27,29
    :linenos:

The *D* code above can be run with the following command:

.. code-block:: bash

    $ dub run --build=release --single sir_example.d