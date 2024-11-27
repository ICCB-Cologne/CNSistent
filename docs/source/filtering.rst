Filtering
---------

To help with defining cutoff regions for samples, it is possible to use a utility function that finds a knee/elbow point for a dataset.

The primary function is ``find_bends``, which will convert a value (e.g. coverage) into a cumulative distribution between ``min_val`` and ``max_val`` with ``steps`` steps. The function will then find the knee/elbow point in the distribution.

The function will return the point where the slope of the cumulative distribution has the highest convex (knee) or concave (elbow) curvature. To avoid finding local minimum, ``dist`` can be set to consider the angle between the mean of the nearest ``dist`` points. If ``allow_pad`` is set to true, the endpoints are also considered as potential knees/elbows, with the slope at the beginning and the end being 0.

.. image:: files/example_knee.png
    :alt: CN Tracks