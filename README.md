This package was developed using SageMath kernel 9.1 and depends on the sage package being available locally. Alternatively, you can make a project on Cocalc.com. If you import this package in a Cocalc project, then you'll have to run your code using a SageMath kernel (stable is recommended).

# Usage

## calculatePartition
This function will compute Kostant's partition function or the q-analog of that function on the given weight in the given Lie algebra.

### parameters
name - string: The lie algebra you are working with (i.e. G2, A4, ...)

weight - array of numbers: The weight you want to partition. It is entered as a linear combination of simple roots or fundamental roots based on the parameter "simple". In the case of simple roots in G2, [1,2] represents the weight <img src="https://rawgit.com/antman1935/lie_algebras/master/svgs/7f8aa090d23837855fa5b83b9db06b98.svg?invert_in_darkmode" align=middle width=71.4879pt height=21.18732pt/>.

positive_roots - array of vectors: This parameter can be changed to restrict which roots will be used to partition, but otherwise it should be left alone. Default value is the entire set of positive roots.

translations - should be left alone

q_analog - boolean: This tells us to give either the q_analog of Kostant's partition function or the usual version. Default value is false.

simple - boolean: This determines if the weight entered is interpreted as a linear combination of simple roots or fundamental weights. Default value is true.

## findAltSet
This function gives the subset of the Weyl group of the given Lie algebra that will have a nonzero contribution to the sum in the multiplicity formula, given parameters <img src="https://rawgit.com/antman1935/lie_algebras/master/svgs/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode" align=middle width=9.90495pt height=14.15535pt/> and <img src="https://rawgit.com/antman1935/lie_algebras/master/svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58914pt height=22.83138pt/>.

### parameters
name - string: The lie algebra you are working with (i.e. G2, A4, ...)

lamb - array of numbers: The highest weight of an irreducible representation entered as linear combination of simple roots or fundamental weights. Default value is the highest positive root.

mu - array of numbers: Any weight entered as a linear combination of fundamental weights or simple roots. Default value is the 0 weight.

simple - boolean: A boolean flag telling us to interpret the weights entered as a linear combination of simple roots or fundamental weights. Default value is true.

## calculateMultiplicity
This function computes the multiplicity formula, or it's q-analog, for the given <img src="https://rawgit.com/antman1935/lie_algebras/master/svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58914pt height=22.83138pt/> and <img src="https://rawgit.com/antman1935/lie_algebras/master/svgs/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode" align=middle width=9.90495pt height=14.15535pt/> in the given Lie algebra.

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

weight - array of numbers: The weight you want to partition. It is entered as a linear combination of simple roots or fundamental roots based on the parameter "simple". In the case of simple roots in G2, [1,2] represents the weight <img src="https://rawgit.com/antman1935/lie_algebras/master/svgs/7f8aa090d23837855fa5b83b9db06b98.svg?invert_in_darkmode" align=middle width=71.4879pt height=21.18732pt/>.

tex - string: The name of the file the outputted Latex will be saved to. You must include the ".tex" ending.

simple - boolean: Tells us whether to interpret the input weight as a linear combination of simple roots or fundamental weights. The also determines whether we print the positive roots with fundamental weights or simple roots.
