from Weight import Weight
from copy import copy

TerminalNode = '0'

class PartitionTree:

    def __init__(self, r, rs, index, terms):
        self.weight = r;
        self.children = [];
        self.terms = terms
        
        # if all coefficients are zero, then we've found a successful partition
        if (r.isZero()):
            self.children.append(TerminalNode);
            
        else:
            if (index < len(rs)):
                child = PartitionTree(Weight(list(r.coefficients)), rs, index + 1, 0)
                child.cleanTree()
                self.children.append(child); # exclude current root
                
                # start looking at combinations of current root
                newChild = r - rs[index]
                terms = 1
                
                # continue to subtract current root until a coefficient is below zero
                while (newChild.isPositive()):
                    child = PartitionTree(newChild, rs, index + 1, terms)
                    child.cleanTree()
                    self.children.append(child);
                    newChild = newChild - rs[index]
                    terms += 1
                    
    # deletes branches that do not end in the terminal node
    def cleanTree(self):
        if (len(self.children) == 0):
            del(self)
            return True
        
        delete = True
        for child in self.children:
            if child == TerminalNode: # if a branch ends in TerminalNode, don't delete
                return False
            
            # if branch has multiple branches following, check to see if at least one ends in TerminalNode
            delete = delete and child.cleanTree()
            
        if (delete):
            del(self)
            return True
        else:
            return False
                    
    def getPartitions(self, array):
        for x in self.children:
            x.getPartitionHelper([], array)
        
    def getPartitionHelper(self, ar, array):
        ar.append(self.terms)
        for x in self.children:
            if (x == TerminalNode):
                array.append(ar)
            else:
                x.getPartitionHelper(copy(ar),array)
    
    def countPartitions(self):
        parts = 0;
        for child in self.children:
            if (child == TerminalNode):
                return 1;
            else:
                parts += child.countPartitions();
        return parts
        
    def generatePq(self, e):
        for child in self.children:
            child.generatePqHelper(self.weight, 0, e)
        
    def generatePqHelper(self, prevRoot, rootsUsed, equationCoefficients):
        rootsUsed += self.terms
            
        for child in self.children:
            if (TerminalNode == child):
                equationCoefficients[rootsUsed] = equationCoefficients[rootsUsed] + 1
                return
            else:
                child.generatePqHelper(self.weight, rootsUsed, equationCoefficients)
    

                
def partitionInteger(i):
    roots = [Weight([i]) for i in range(1, i + 1)]
    partitionOfI = PartitionTree(roots[len(roots) - 1], roots, 0)
    partitionOfI.cleanTree()
    return partitionOfI

def main():
    p5 = partitionInteger(10)
    num = p5.countPartitions()
    print("10 has " + str(num) + " partitions")
    m_q = [0 for i in range(0, 11)]
    p5.generatePq(m_q)
    print(m_q)
    