from numpy import *


class BadArgumentError(Exception):

    def __init__(self):
        self.message = "Arguments given not valid"

    def __str__(self):
        return self.message


class FileParseError(Exception):

    def __init__(self):
        self.message = "File not valid to be parsed"

    def __str__(self):
        return self.message


class SimplexElement:

    def __init__(self):
        self.objFunc = []
        self.Aj = []
        self.b = None
        self.XbInverse = None
        self.BasicIndexes = []
        self.CoeffBasic = None
        self.Vars = []
        self.varConst = []  # for the greater and less than sign

    def readFromFile(self, filename):
        if filename is None:
            raise BadArgumentError
        else:
            validLines = self.getValidLines(filename)
            validLines = self.getObjectiveFunc(validLines)
            validLines = self.parseConstraints(validLines)
            self.getVarsFilled()
            self.fillBasicCoeffMatrix()
            self.getXbInverse()

    def getValidLines(self, filename):
        FILE = open(filename, "r")
        lines = []
        for line in FILE:
            if line.startswith("//"):
                pass
            elif line.startswith("#"):
                line = line[1:len(line) - 1]
                lines.append(line)
            else:
                raise FileParseError
        return lines

    def getObjectiveFunc(self, vlines):
        if vlines[0].startswith("[obj]"):
            vlines[1] = vlines[1].rstrip("]")
            vlines[1] = vlines[1].lstrip("[")
            func = vlines[1].split(",")
            self.objFunc = list(map(int, func))
        else:
            raise FileParseError
        vlines = vlines[3:]
        return vlines

    def parseConstraints(self, vlines):
        if vlines[0].startswith("[const]"):
            i = 1
            AMatrix = []
            bMatrix = []
            while not vlines[i].startswith("[!const]"):
                vlines[i] = vlines[i].split(" ")
                vlines[i][0] = vlines[i][0][1:len(vlines[i][0]) - 1]
                vlines[i][0] = vlines[i][0].split(",")
                vlines[i][0] = list(map(int, vlines[i][0]))
                AMatrix.append(vlines[i][0])
                self.varConst.append(vlines[i][1])
                bMatrix.append(int(vlines[i][2]))
                i += 1
            vlines = vlines[i + 1:]
            AMatrix = self.addSlackVariables(AMatrix)
            self.convertToAj(AMatrix)
            self.convertToB(bMatrix)
            return vlines
        else:
            raise FileParseError

    def addSlackVariables(self, AMatrix):
        for i in range(0, len(self.varConst)):
            if self.varConst[i] == 'lt':
                for j in range(0, len(AMatrix)):
                    if j == i:
                        AMatrix[j].append(1)
                        self.BasicIndexes.append(len(AMatrix[j]) - 1)
                        self.objFunc.append(0)
                    else:
                        AMatrix[j].append(0)
            elif self.varConst[i] == 'gt':
                for j in range(0, len(AMatrix)):
                    if j == i:
                        AMatrix.append(-1)
                        AMatrix.append(1)
                        self.objFunc.append(0)
                        self.objFunc.append(0)
                        self.BasicIndexes.append(len(AMatrix[j]))
                    else:
                        AMatrix.append(0)
                        AMatrix.append(0)
        return AMatrix

    def convertToAj(self, AMatrix):
        dim1 = len(AMatrix)
        dim2 = len(AMatrix[0])
        for j in range(0, dim2):
            string = ''
            for i in range(0, dim1):
                if i < dim1 - 1:
                    string += str(AMatrix[i][j]) + ";"
                else:
                    string += str(AMatrix[i][j])
            self.Aj.append(matrix(string))

    def convertToB(self, bMatrix):
        string = ''
        for i in range(0, len(bMatrix)):
            if i < len(bMatrix) - 1:
                string += str(bMatrix[i]) + ";"
            else:
                string += str(bMatrix[i])
        self.b = matrix(string)

    def getVarsFilled(self):
        for i in range(0, len(self.objFunc)):
            item = []
            if i in self.BasicIndexes:
                item.append('b')
            else:
                item.append('nb')
            item.append(self.objFunc[i])
            self.Vars.append(item)

    def fillBasicCoeffMatrix(self):
        string = ''
        for i in range(0, len(self.Vars)):
            if self.Vars[i][0] == 'b':
                string += str(self.Vars[i][1]) + ";"
        print(string)
        string = string[:len(string) - 1]
        self.CoeffBasic = matrix(string)

    def getXbInverse(self):
        nbv = 0
        for i in range(0, len(self.Vars)):
            if self.Vars[i][0] == 'b':
                nbv += 1
        self.XbInverse = eye(nbv)
