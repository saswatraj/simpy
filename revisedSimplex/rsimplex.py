from numpy import *


class RevisedSimplex:

    def __init__(self, simplexElement):
        self.objFunc = simplexElement.objFunc
        self.Aj = simplexElement.Aj
        self.b = simplexElement.b
        self.XbInverse = simplexElement.XbInverse
        self.BasicIndexes = simplexElement.BasicIndexes
        self.CoeffBasic = simplexElement.CoeffBasic
        self.Vars = simplexElement.Vars
        self.MAX_VAR = 99999
        self.optimizedVal = None

    def display(self):
        """Displays each step information"""
        print("-------------------------------------------------------------")
        print("Objective Func:")
        print(self.objFunc)
        print("\nAj:")
        print(self.Aj)
        print("\nb:")
        print(self.b)
        print("\nXbInverse:")
        print(self.XbInverse)
        print("\nBasicIndexes:")
        print(self.BasicIndexes)
        print("\nCoeffBasic:")
        print(self.CoeffBasic)
        print("\nVars:")
        print(self.Vars)
        print("-------------------------------------------------------------")

    def solve(self):
        """Solves the optimization problem
            returns : optimized value and last state table"""
        self.display()
        smallest, enteringVar = self.getEnteringVar()
        #print(smallest,enteringVar)
        if smallest >= 0:
            Xb = self.XbInverse * self.b
            val = self.CoeffBasic.getT() * Xb
            self.optimizedVal = float(val)
            return
        leaveVar, leaveRow, Yent = self.getLeavingVar(enteringVar)
        self.BasicIndexes[leaveRow] = enteringVar
        self.Vars[leaveVar][0] = 'nb'
        self.Vars[enteringVar][0] = 'b'
        self.changeCoeffOfBasic()
        eMatrix = self.getEMatrix(leaveRow, Yent)
        self.XbInverse = eMatrix * self.XbInverse
        self.solve()

    def getEnteringVar(self):
        """Finds the entering variable and the smallest value"""
        smallest = self.MAX_VAR
        enteringVar = None
        for i in range(0, len(self.Vars)):
            if self.Vars[i][0] == 'nb':
                temp = self.XbInverse * self.Aj[i]
                temp = self.CoeffBasic.getH() * temp
                temp = temp - self.Vars[i][1]
                if temp <= smallest:
                    smallest = float(temp)
                    enteringVar = i
        return (smallest, enteringVar)

    def getLeavingVar(self, enteringVar):
        """Finds the leaving variable and the row of the leaving variable"""
        Xb = self.XbInverse * self.b
        Yent = self.XbInverse * self.Aj[enteringVar]
        leaveRow = None
        lsmall = self.MAX_VAR
        for j in range(0, size(self.b)):
            temp = Xb[j] / Yent[j]
            if temp > 0 and temp <= lsmall:
                lsmall = float(temp)
                leaveRow = j
        leaveVar = self.BasicIndexes[leaveRow]
        return (leaveVar, leaveRow, Yent)

    def changeCoeffOfBasic(self):
        """Changes the coeff of basic variables"""
        string = ''
        for j in range(0, size(self.BasicIndexes)):
            string += str(self.Vars[self.BasicIndexes[j]][1])
            if j < size(self.CoeffBasic) - 1:
                string += ';'
        self.CoeffBasic = matrix(string)

    def getEMatrix(self, leaveRow, Yent):
        """Finds the new B Inverse from the old"""
        eMatrix = []
        dim = len(self.BasicIndexes)
        for i in range(0, dim):
            temp = [0] * dim
            temp[i] = 1
            eMatrix.append(temp)
        #replace the column
        for i in range(0, len(self.BasicIndexes)):
            if i == leaveRow:
                eMatrix[i][leaveRow] = float(1 / Yent[leaveRow])
            else:
                eMatrix[i][leaveRow] = float(
                    -Yent[i] / Yent[leaveRow])
        string = ''
        for i in range(0, len(eMatrix)):
            st = ''
            for j in range(0, len(eMatrix)):
                st += str(eMatrix[i][j]) + ' '
            if i < len(eMatrix) - 1:
                string += st + ';'
            else:
                string += st
        eMatrix = matrix(string)
        return eMatrix
