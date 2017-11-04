This package was originally developed using SageMath kernel 7.3. This package depends on the sage package being available. This means that importing it into your own project requires that you are on a Mac or Linux machine (last I checked it wasn't available on Windows) and have SageMath installed. Alternatively, you can make a project on Cocalc.com. If you import this package in a Cocalc project, then you'll have to run your code using either a SageMath kernel (stable is recommended) or with the Python 2 (SageMath) kernel.

# Usage

## calculatePartition
This function will compute Kostant's partition function or the q-analog of that function on the given weight in the given Lie algebra.

### parameters
name - string: The lie algebra you are working with (i.e. G2, A4, ...)

weight - array of numbers: The weight you want to partition. It is entered as a linear combination of simple roots or fundamental roots based on the parameter "simple". In the case of simple roots in G2, [1,2] represents the weight $1\alpha_1 + 2\alpha_2$.

positive_roots - array of vectors: This parameter can be changed to restrict which roots will be used to partition, but otherwise it should be left alone. Default value is the entire set of positive roots.

translations - should be left alone

q_analog - boolean: This tells us to give either the q_analog of Kostant's partition function or the usual version. Default value is false.

simple (not implemented) - boolean: This determines if the weight entered is interpreted as a linear combination of simple roots or fundamental weights. Default value is true.

## findAltSet
This function gives the subset of the Weyl group of the given Lie algebra that will have a nonzero contribution to the sum in the multiplicity formula, given parameters $\mu$ and $\lambda$.

### parameters
name - string: The lie algebra you are working with (i.e. G2, A4, ...)

lamb - array of numbers: The highest weight of an irreducible representation entered as linear combination of simple roots or fundamental weights. Default value is the highest positive root.

mu - array of numbers: Any weight entered as a linear combination of fundamental weights or simple roots. Default value is the 0 weight.

simple (not implemented yet) - boolean: A boolean flag telling us to interpret the weights entered as a linear combination of simple roots or fundamental weights. Default value is true.

## calculateMultiplicity
This function computes the multiplicity formula, or it's q-analog, for the given $\lambda$ and $\mu$ in the given Lie algebra.

### parameters
name - string: The lie algebra you are working with (i.e. G2, A4, ...)

lamb - array of numbers: The highest weight of an irreducible representation entered as linear combination of simple roots or fundamental weights. Default value is the highest root.

mu - array of numbers: Any weight entered as a linear combination of fundamental weights or simple roots. Default value is the 0 weight.

q_analog - boolean: This tells us to give either the q_analog of Kostant's partition function or the usual version. Default value is false.

simple - boolean: A boolean flag telling us to interpret the weights entered as a linear combination of simple roots or fundamental weights. The default value is true.

## printPartitions
This function gets all of the possible partitions of the entered weight with the positive roots of the entered Lie algebra.

### parameters
name - string: The lie algebra you are working with (i.e. G2, A4, ...)

weight - array of numbers: The weight you want to partition. It is entered as a linear combination of simple roots or fundamental roots based on the parameter "simple". In the case of simple roots in G2, [1,2] represents the weight $1\alpha_1 + 2\alpha_2$.

tex - string: The name of the file the outputted Latex will be saved to. You must include the ".tex" ending.
