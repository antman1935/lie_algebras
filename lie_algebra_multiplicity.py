import __future__
from Weight import Weight
from PartitionTree import PartitionTree
from sage.all import *

def Eval(string):
    return eval(compile(str(string), '<string>', 'eval', __future__.division.compiler_flag))

def getFundamentalToSimple(name):
    simples = RootSystem(name).ambient_space().simple_roots()
    simpleM = ([Eval((str(x).replace("(", "[").replace(")", "]"))) for x in simples])
    for i in range(len(simpleM), len(simpleM[0])):
        simpleM.append([1 if j==i else 0 for j in range(0, len(simpleM[0]))])

    fundamentals = RootSystem(name).ambient_space().fundamental_weights()
    fundM = ([Eval((str(x).replace("(", "[").replace(")", "]"))) for x in fundamentals])
    for i in range(len(fundM), len(fundM[0])):
        fundM.append([1 if j==i else 0 for j in range(0, len(fundM[0]))])

    return matrix(simpleM).transpose().inverse() * matrix(fundM).transpose()

def getBasisChange(name):
    simples = RootSystem(name).ambient_space().simple_roots()
    bChange = ([Eval((str(x).replace("(", "[").replace(")", "]"))) for x in simples])
    for i in range(len(bChange), len(bChange[0])):
        bChange.append([1 if j==i else 0 for j in range(0, len(bChange[0]))])
    return matrix(bChange).transpose().inverse()

def geometricSumForPartition(positive_root, translations, q_analog):
    x = 1
    for i in range(0, len(positive_root)):
        for j in range(0, positive_root[i]):
            x = x * translations["A" + str(i+1)]
    return 1/(1 - x) if not q_analog else 1/(1 -translations['q']*x)

def calculatePartition(name, weight, positive_roots = [], translations = {}, q_analog = False):
    if positive_roots == []:
        bChange = getBasisChange(name)
        positive_roots = [vector(list(Eval(x))) for x in RootSystem(name).ambient_space().positive_roots()]
        positive_roots = [bChange * x for x in positive_roots]

    if translations == {}:
        s = ''
        if q_analog:
            s = 'q, '
        s += 'A1'
        for i in range(1, len(weight)):
            s += ', A' + str(i + 1)
        variables = var(s)
        for i in range(0, len(weight)):
            translations["A" + str(i+1)] = eval("A" + str(i+1))

        if q_analog:
            translations['q'] = q

    termsForSum = [geometricSumForPartition(list(x), translations, q_analog) for x in positive_roots]
    answer = 1
    for x in termsForSum:
        answer *= x

    for i in range(0, len(weight)):
        answer = answer.series(translations["A" + str(i+1)], weight[i] + 1).truncate()
        answer = answer.coefficient(translations["A" + str(i+1)], weight[i])

    answer = answer.expand()
    return answer

def findAltSet(name, lamb = None, mu = None):
    # initialize constants and vector space for the lie algebra
    lie_algebra = RootSystem(name).ambient_space()
    weyl_group = WeylGroup(name, prefix = "s")
    simples = weyl_group.gens()

    altset = [weyl_group.one()]

    # used to change the basis from the standard basis of R^n to simple roots
    changeBasis = getBasisChange(name)

    # if lambda is not specified, the highest root is used
    if lamb == None:
        lamb = lie_algebra.highest_root()

    # if mu is not specified, 0 vector is used
    if mu == None:
        mu = weyl_group.domain()([0 for i in range(0, len(lie_algebra.simple_roots()))])

    # check to see if the alt set is the empty set
    init = (lamb + mu)
    init = changeBasis * vector(list(Eval(init)))
    init = Weight(init)

    if init.isNegative():
        return []

    rho = lie_algebra.rho()
    length = len(altset)
    i=0
    while i < length:
        for simple in simples:
            if ((altset[i] == simple)or (altset[i] == altset[i] * simple)):
                continue
            res = (altset[i]*simple).action(lamb + rho) - (rho + mu)
            res = changeBasis * vector(list(Eval(res)))
            res = Weight(res)

            if not (res.isNegative() or res.hasFraction()):
                if not (altset[i]*simple in altset):
                    altset.append(altset[i]*simple)
                    length += 1
        i+=1

    return altset


def calculateMultiplicity(name, lamb = None, mu = None, q_analog = False, simple=True):
    mult = 0
    lie_algebra = RootSystem(name).ambient_space()
    weyl_group = WeylGroup(name, prefix = "s")

    # used to change the basis from the standard basis of R^n to simple roots
    standard_to_simple = getBasisChange(name)
    fund_to_simple = getFundamentalToSimple(name)

    positive_roots = [vector(list(Eval(x))) for x in RootSystem(name).ambient_space().positive_roots()]
    positive_roots = [standard_to_simple * x for x in positive_roots]

    # if lambda is not specified, the highest root is used
    if lamb == None:
        lamb = lie_algebra.highest_root()
    else:
        while not len(lamb) == standard_to_simple.ncols():
            lamb.append(0)
        if not simple:
            lamb = standard_to_simple.inverse() * (fund_to_simple * vector(lamb))
        else:
            lamb = standard_to_simple.inverse() * vector(lamb)
        lamb = weyl_group.domain()(list(eval(str(lamb))))

    # if mu is not specified, 0 vector is used
    if mu == None:
        mu = weyl_group.domain()([0 for i in range(0, len(lie_algebra.simple_roots()))])
    else:
        while not len(mu) == standard_to_simple.ncols():
            mu.append(0)
        if not simple:
            mu = standard_to_simple.inverse() * (fund_to_simple * vector(mu))
        else:
            mu = standard_to_simple.inverse() * vector(mu)
        mu = weyl_group.domain()(list(eval(str(mu))))

    rho = lie_algebra.rho()
    altset = findAltSet(name, lamb, mu)
    #print(changeBasis * vector(list(Eval(lamb))))
    translations = {}
    for elm in altset:
        # expression in partition function
        res = elm.action(lamb + rho) - (mu + rho)

        #change basis from standard basis to simple roots
        res = vector(list(Eval(res)))
        res = standard_to_simple * res

        term = calculatePartition(name, list(res), positive_roots, translations = translations, q_analog=q_analog)

        term *= (-1)**elm.length()
        mult += term

    return mult

def printPartitions(name, weight, tex):
    lie_algebra = RootSystem(name).ambient_space()

    # used to change the basis from the standard basis of R^n to simple roots
    changeBasis = getBasisChange(name)

    positive_roots = [vector(list(Eval(x))) for x in RootSystem(name).ambient_space().positive_roots()]
    positive_roots = [getBasisChange(name) * x for x in positive_roots]
    latex_roots = ["".join([str(root[j] if root[j] != 1 else "") + "\\alpha_" + str(j+1)+"+" if root[j] != 0 else "" for j in range(0, len(root))]) for root in positive_roots]
    latex_roots = [root[0:len(root)-1] for root in latex_roots]
    weight_positive_roots = [Weight(list(root)) for root in positive_roots]
    weight = Weight(weight)
    tree = PartitionTree(weight, weight_positive_roots, 0, 0)

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
