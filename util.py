from sage.all import vector, matrix
from sage.all import RootSystem
from sage.all import var

def convertWeightToList(weight):
    arr = []
    for index, value in sorted(weight):
        while len(arr) < index:
            arr.append(0)
        arr.append(value)
    while len(arr) < weight.parent().dimension():
        arr.append(0)
    return arr

def getStandardToSimpleBasisChange(name):
    simples = RootSystem(name).ambient_space().simple_roots()
    bChange = [convertWeightToList(x) for x in simples]

    # insert standard e_i for i > len(simple roots) to make a square vector
    for i in range(len(bChange), len(bChange[0])):
        bChange.append([1 if j==i else 0 for j in range(0, len(bChange[0]))])

    # transpose gets simple to standard basis change, inverse reverses
    return matrix(bChange).transpose().inverse()

def getPositiveRoots(name):
    standard_to_simple_basis_change = getStandardToSimpleBasisChange(name)
    positive_roots = [vector(convertWeightToList(x)) for x in RootSystem(name).ambient_space().positive_roots()]
    return [standard_to_simple_basis_change * x for x in positive_roots]

def getVariableDictionary(lie_algebra, q_analog):
    s = ""
    if q_analog:
        s = "q, "
    s += ", ".join([f"A{i + 1}" for i in range(lie_algebra.dimension())])
    variables = var(s)
    return {str(variable): variable for variable in variables}

def geometricSumForPartition(positive_root, translations, q_analog):
    x = 1 if not q_analog else translations["q"]
    for i in range(len(positive_root)):
            x = x * (translations["A" + str(i+1)] ** positive_root[i])
    print(positive_root, " becomes ", x)
    return 1/(1 - x)

def getLambda(lie_algebra, standard_to_simple_basis_change, lamb):
    # if lamb is not specified, highest root is used
    if lamb == None:
        return lie_algebra.highest_root()
    elif type(lamb) is list:
        while not len(lamb) == lie_algebra.dimension():
            lamb.append(0)
        temp = standard_to_simple_basis_change.inverse() * vector(lamb)
        return lie_algebra(list(temp))

    return lamb

def getMu(lie_algebra, standard_to_simple_basis_change, mu):
    # if mu is not specified, 0 vector is used
    if mu == None:
        return lie_algebra(0)
    elif type(mu) is list:
        while not len(mu) == lie_algebra.dimension():
            mu.append(0)
        temp = standard_to_simple_basis_change.inverse() * vector(mu)
        return lie_algebra(list(temp))

    return mu
