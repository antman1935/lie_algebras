import __future__
from Weight import Weight
from PartitionTree import PartitionTree
import util as util
from sage.all import vector
from sage.all import RootSystem, WeylGroup

def Eval(string):
    return eval(compile(str(string), '<string>', 'eval', __future__.division.compiler_flag))

def convertWeightParameters(name, weight, simple):
    weyl_group = WeylGroup(name, prefix = "s")
    standard_to_simple = getBasisChange(name)
    fund_to_simple = getFundamentalToSimple(name)

    while not len(weight) == standard_to_simple.ncols():
        weight.append(0)
    if not simple:
        weight = standard_to_simple.inverse() * (fund_to_simple * vector(weight))
    else:
        weight = standard_to_simple.inverse() * vector(weight)

    return weyl_group.domain()(list(eval(str(weight))))

def changeFundToSimple(name, weight):
    basis_change = getFundamentalToSimple(name)
    weight = list(weight)
    while not len(weight) == basis_change.ncols():
            weight.append(0)
    return basis_change * vector(weight)

def getFundamentalToSimple(name):
    simples = RootSystem(name).ambient_space().simple_roots()
    simple_basis = ([Eval((str(x).replace("(", "[").replace(")", "]"))) for x in simples])
    for i in range(len(simple_basis), len(simple_basis[0])):
        simple_basis.append([1 if j==i else 0 for j in range(0, len(simple_basis[0]))])

    fundamentals = RootSystem(name).ambient_space().fundamental_weights()
    fund_basis = ([Eval((str(x).replace("(", "[").replace(")", "]"))) for x in fundamentals])
    for i in range(len(fund_basis), len(fund_basis[0])):
        fund_basis.append([1 if j==i else 0 for j in range(0, len(fund_basis[0]))])

    return matrix(simple_basis).transpose().inverse() * matrix(fund_basis).transpose()

def getBasisChange(name):
    simples = RootSystem(name).ambient_space().simple_roots()
    basis_change = ([Eval((str(x).replace("(", "[").replace(")", "]"))) for x in simples])
    for i in range(len(basis_change), len(basis_change[0])):
        basis_change.append([1 if j==i else 0 for j in range(0, len(basis_change[0]))])
    return matrix(basis_change).transpose().inverse()

"""
Calculate the Partition (or q-analog) function for the weight of the given Lie
Algebra.

lie_algebra_name: name of the Lie Algebra (i.e. A2, B3, G2, E7, etc.)
weight: The weight to be partitioned, as a list of coefficients to simple roots.
q_analog: Boolean, defaults to False. If True, compute the q-analog of the
          partition function.
positive_roots: The set of positive weights to use for partitioning. Defaults to
                the entire set of positive roots. Pass in a subset to restrict
                which partitons will be counted.
translations: Dictionary from string to variable objects used in the geometric
              sum. Generate with util.getVariableDictionary or leave empty.
"""
def calculatePartition(lie_algebra_name, weight, q_analog = False, positive_roots = [], translations = {}):
    lie_algebra = RootSystem(lie_algebra_name).ambient_space()

    if positive_roots == []:
        positive_roots = util.getPositiveRoots(lie_algebra_name)

    if translations == {}:
        translations = util.getVariableDictionary(lie_algebra, q_analog)

    termsForSum = [util.geometricSumForPartition(list(positive_root), translations, q_analog) for positive_root in positive_roots]
    answer = 1
    for x in termsForSum:
        answer *= x

    # find the partition (or q-analog) using geometric series expansion. Our
    # answer is the coefficient of the term with variable exponents matching the
    # coefficents of their corresponding simple root.
    for i in range(len(weight)):
        answer = answer.series(translations["A" + str(i+1)], weight[i] + 1).truncate()
        answer = answer.coefficient(translations["A" + str(i+1)], weight[i])

    answer = answer.expand()
    return answer

"""
Return a list containing the alternation set for the given lambda and mu. These
are the elements of the Weyl Group of the Lie Algebra that contribute non-zero
terms to the sum in the muliplicity formula.

lie_algebra_name: name of the Lie Algebra (i.e. A2, B3, G2, E7, etc.)
lamb: The value of lambda in the multiplicity formula, as a list of coefficients
      to simple roots.
mu: The value of mu in the multiplicity formula, as a list of coefficients to
    simple roots.
"""
def findAltSet(lie_algebra_name, lamb = None, mu = None):
    # initialize constants and vector space for the lie algebra
    lie_algebra = RootSystem(lie_algebra_name).ambient_space()
    weyl_group = WeylGroup(lie_algebra_name, prefix = "s")
    simples = weyl_group.gens()

    altset = [weyl_group.one()]

    # used to change the basis from the standard basis of R^n to simple roots
    standard_to_simple_basis_change = util.getStandardToSimpleBasisChange(lie_algebra_name)

    lamb_weight = util.getLambda(lie_algebra, standard_to_simple_basis_change, lamb)
    mu_weight = util.getMu(lie_algebra, standard_to_simple_basis_change, mu)
    rho = lie_algebra.rho()

    # check to see if the alt set is the empty set
    init = (lamb_weight + rho) - (rho + mu_weight)
    init = standard_to_simple_basis_change * vector(util.convertWeightToList(init))
    init = Weight(list(init))

    if init.isNegative():
        return []

    # conduct breadth first search of weyl group elements that can be applied
    # without making the weight to partition have any negative or fractional
    # coefficients. We can stop appending them once we get a negative/fractional
    # coefficient because applying another permutation will not remove a
    # fraction/make a weight less negative.
    length = len(altset)
    i = 0
    while i < length:
        for simple in simples:
            elm = altset[i] * simple
            if (altset[i] == simple) or (altset[i] == elm):
                continue

            res = elm.action(lamb_weight + rho) - (rho + mu_weight)
            res = standard_to_simple_basis_change * vector(util.convertWeightToList(res))
            res = Weight(list(res))

            if not res.isNegative() and not res.hasFraction():
                if not elm in altset: # TODO: check if this is this necessary
                    altset.append(altset[i]*simple)
                    length += 1
        i+=1

    return altset

"""
Compute Kostant's multiplicity formula (or the q-analog) for lambda and mu in
the given Lie Algebra.

lie_algebra_name: name of the Lie Algebra (i.e. A2, B3, G2, E7, etc.)
lamb: The value of lambda in the multiplicity formula, as a list of coefficients
      to simple roots.
mu: The value of mu in the multiplicity formula, as a list of coefficients to
    simple roots.
q_analog: Boolean, defaults to False. If True, compute the q-analog of the
          multiplicity formula.
"""
def calculateMultiplicity(lie_algebra_name, lamb = None, mu = None, q_analog = False):
    lie_algebra = RootSystem(lie_algebra_name).ambient_space()
    weyl_group = WeylGroup(lie_algebra_name, prefix = "s")

    # used to change the basis from the standard basis of R^n to simple roots
    standard_to_simple_basis_change = util.getStandardToSimpleBasisChange(lie_algebra_name)

    lamb_weight = util.getLambda(lie_algebra, standard_to_simple_basis_change, lamb)
    mu_weight = util.getMu(lie_algebra, standard_to_simple_basis_change, mu)
    rho = lie_algebra.rho()

    positive_roots = util.getPositiveRoots(lie_algebra_name)
    translations = util.getVariableDictionary(lie_algebra, q_analog)

    mult = 0
    altset = findAltSet(lie_algebra_name, lamb, mu)
    for elm in altset:
        # expression in partition function
        res = elm.action(lamb_weight + rho) - (mu_weight + rho)

        #change basis from standard basis to simple roots
        res = vector(util.convertWeightToList(res))
        res = standard_to_simple_basis_change * res

        term = calculatePartition(lie_algebra_name, list(res), q_analog, positive_roots, translations)

        term *= (-1)**elm.length()
        mult += term

    return mult

def printPartitions(name, weight, tex, simple = True):
    lie_algebra = RootSystem(name).ambient_space()

    # used to change the basis from the standard basis of R^n to simple roots
    change_basis = getBasisChange(name)
    symb = "\\alpha_"

    # if the weight is in terms of fundamentals, change basis to simples for partitioning
    if not simple:
        weight = changeFundToSimple(name, weight)
        symb = "\\omega_"

    positive_roots = [vector(list(Eval(x))) for x in RootSystem(name).ambient_space().positive_roots()]
    positive_roots = [change_basis * x for x in positive_roots]
    weight_positive_roots = [Weight(list(root)) for root in positive_roots]

    weight = Weight(weight)
    tree = PartitionTree(weight, weight_positive_roots, 0, 0)

    # here we change the positive roots to be in terms of the fundamental weights so what is printed matches the input
    if not simple:
        change_basis = getFundamentalToSimple(name).inverse()
        positive_roots = [change_basis * x for x in positive_roots]
    latex_roots = ["".join([str(root[j] if root[j] != 1 else "") + symb + str(j+1)+"+" if root[j] != 0 else "" for j in range(0, len(root))]) for root in positive_roots]
    latex_roots = [root[0:len(root)-1].replace("+-1", "-").replace("-1", "-").replace("+-", "-") for root in latex_roots]

    partitions = []
    tree.getPartitions(partitions)

    output = open(tex, "w")
    output.write("\\documentclass{article}\n\\begin{document}\n")
    for partition in partitions:
        string = ""
        for i in range(len(partition)):
            coeff = partition[i]
            if coeff > 1:
                string += str(coeff) + "(" + latex_roots[i] + ") +"
            elif coeff > 0:
                string += "(" + latex_roots[i] + ")+"
        string = "$"+ string[0:len(string)-1]
        string += "$\\\\\\\\\n\n"
        output.write(string)
    output.write("\\end{document}")
    output.close()
    print "Done"

if __name__ == "__main__":
    #calculateMultiplicity("A3", [3,3,3], [2,2,2], True)
    #calculateMultiplicity("A3", [3,3,3], q_analog=True)
    #calculateMultiplicity("A3", mu=[2,2,2], q_analog=True)
    x = calculatePartition("A4", [9,9,8,9],q_analog=False)
    #printPartitions("C2", [2,2], "output.tex")
    print(x)
    #findAltSet("D4")
