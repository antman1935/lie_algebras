from copy import copy

class Weight:
    def __init__(self, cs):
        self.coefficients = cs
        
    def height(self):
        ht = 0
        for x in self.coefficients:
            ht += x
        return ht
    
    def __sub__(self, other):
        c = copy(self.coefficients);
        for i in range(0, len(c)):
            c[i] -= other.coefficients[i]
        return Weight(c)
            
    def __add__(self, other):
        c = copy(self.coefficients);
        for i in range(0, len(c)):
            c[i] += other.coefficients[i]
        return Weight(c)
        
    def __str__(self):
        return str(self.coefficients)
    
    def equals(self, other):
        for i in range(0, len(self.coefficients)):
            if not (self.coefficients[i] == other.coefficients[i]):
                return False
        return True
    
    def __ne__(self, other):
        return not self == other
        
    def isZero(self):
        for coeff in self.coefficients:
            if coeff != 0:
                return False
        return True
        
    def isPositive(self):
        for coeff in self.coefficients:
            if coeff < 0:
                return False
        return True
        
    def isNegative(self):
        for coeff in self.coefficients:
            if coeff < 0:
                return True
        return False
    
    def hasFraction(self):
        for coeff in self.coefficients:
             if ((coeff - int(coeff)) != 0):
                return True
        return False
    
    
            
    