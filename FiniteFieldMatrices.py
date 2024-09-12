from FiniteFieldRepresentation import fieldInt

class fieldIntMatrix:
    def __init__(self, intMatrix : list[list[int]], prime : int, power : int, prodPol = None):
        self.dim = (len(intMatrix), len(intMatrix[0]))
        self.primeBase = prime
        self.primePower = power
        self.fieldSize = prime ** power
        if prodPol is None:
            prodPol = fieldInt(0, prime, power).prodPol
        self.prodPol = prodPol
        
        self.matrix = [[fieldInt(intMatrix[i][j], prime, power, prodPol) for j in range(self.dim[1])] for i in range(self.dim[0])]

    @staticmethod
    def matrixEye(size : int, prime : int, power : int, prodPol = None):
        """
        Returns the size x size identity matrix of
        the finite field of size prime^power.
        """
        
        I = [[1 if (i == j) else 0 for j in range(size)] for i in range(size)]
        return fieldIntMatrix(I, prime, power, prodPol)

        
    
    def scalarProd(self, scalar : int):
        """
        Returns the matrix multiplied by the scalar value
        """
        n, m = self.dim
        s = fieldInt(scalar, self.primeBase, self.primePower)
        result = [[(s * self.matrix[i][j]).value for j in range(m)] for i in range(n)]
        return fieldIntMatrix(result, self.primeBase, self.primePower, self.prodPol)

    def T(self):
        """
        Returns the transpose of the matrix.
        """
        n, m = self.dim
        selfT = [[self.matrix[i][j].value for i in range(n)] for j in range(m)]
        return fieldIntMatrix(selfT, self.primeBase, self.primePower, self.prodPol)
        
    def GaussElim(self):
        """
        Returns a Gaussian Elimination of the matrix.
        """

        n, m = self.dim
        G = fieldIntMatrix.matrixEye(n, self.primeBase, self.primePower) * self
        steps = 0
        for i in range(min(n, m)):
            # finds the row that has a non zero value
            # at the needed column
            pivot = i
            allZero = False
            while G.matrix[pivot][i].value == 0:
                pivot += 1
                if pivot >= n:
                    allZero = True
                    break
            if allZero:
                continue

            # T is the matrix that exchanges the rows "i" and "pivot"
            T = matrixRowExchange(i, pivot, n, self.primeBase, self.primePower)
            G = T * G
            

            steps += 1
            #print(steps)
            
            # Now the actual multiplying and adding rows
            for j in range(i + 1, n):
                if int(G.matrix[j][i]) == 0:
                    continue

                # T2 changes the row "j" by adding a * the row "i"
                a = - (G.matrix[j][i] // G.matrix[i][i])
                T2 = matrixRowMulSum( a, j, i, n, self.primeBase, self.primePower)
                G = T2 * G
                
                steps += 1
                #print(steps)
                #print(G)
        return G

    def GaussElimTotal(self):
        """
        Returns a Gaussian Elimination of the Matrix
        whose upper-left most submatrix is the identity matrix
        """

        n, m = self.dim
        G = self.GaussElim()
        one = fieldInt(1, self.primeBase, self.primePower)
        for i in range(min(n, m)):
            if G.matrix[i][i].value == 0:
                break
            T = fieldIntMatrix.matrixEye(n, self.primeBase, self.primePower)
            T.matrix[i][i] = one // G.matrix[i][i]
            G = T * G
            
            for j in range(i):
                if int(G.matrix[j][i]) == 0:
                    continue
                T2 = matrixRowMulSum(- G.matrix[j][i], j, i, n, self.primeBase, self.primePower)
                G = T2 * G
                #print(G)
        return G

    def weight(self):
        """
        Returns the amount of non-zero entries
        of the matrix.
        """

        boolMatrix = self.toBoolMatrix()
        w = 0 # weight
        for l in boolMatrix:
            w += len(list(filter(lambda x: x, l)))
    
        return w
        
    def toIntMatrix(self):
        """
        Returns the same matrix but whose entries
        are on the integers.
        """
        n, m = self.dim
        intMatrix = [[(self.matrix[i][j].value) for j in range(m)] for i in range(n)]
        return intMatrix

    def toBoolMatrix(self):
        """
        Returns a matrix of boolean values which are
        true where the matrix has nonzero values and
        false otherwise.
        """
        n, m = self.dim
        boolMatrix = [[bool(self.matrix[i][j]) for j in range(m)] for i in range(n)]
        return boolMatrix

    
    def __str__(self):
        t = "Matriz del campo F" + str(self.fieldSize) + "\n"
        for l in self.matrix:
            for i in l:
                t += str(i.value) + ' '
            t += '\n'
        return t
    
    def __add__(self, other):
        # Resulting Matrix
        n, m = self.dim
        sumMatrix = [[0 for __ in range(m)] for _ in range(n)]
        if self.dim != other.dim:
            raise NotImplementedError("Error, matrices de distintos tamaños no se pueden sumar")
        for i in range(n):
            for j in range(m):
                sumMatrix[i][j] += (self.matrix[i][j] + other.matrix[i][j]).value
        return fieldIntMatrix(sumMatrix, self.primeBase, self.primePower, self.prodPol)

    def __neg__(self):
        n, m = self.dim
        negSelf = [[(- a).value for a in li] for li in self.matrix]
        return fieldIntMatrix(negSelf, self.primeBase, self.primePower, self.prodPol)

    def __sub__(self, other):
        # Resulting Matrix
        n, m = self.dim
        subMatrix = [[0 for __ in range(m)] for _ in range(n)]
        if self.dim != other.dim:
            return -1
        
        return self + (- other)

    def __mul__(self, other):
        n, p = self.dim
        p1, m = other.dim
        if p != p1:
            return -1

        # Resulting Matrix
        C = [[fieldInt(0, self.primeBase, self.primePower, self.prodPol) for __ in range(m)] for _ in range(n)]
        for i in range(n):
            for j in range(m):
                for k in range(p):
                    C[i][j] += self.matrix[i][k] * other.matrix[k][j]
                C[i][j] = C[i][j].value
        return fieldIntMatrix(C, self.primeBase, self.primePower, self.prodPol)

    def __invert__(self):
        n, m = self.dim
        I = fieldIntMatrix.matrixEye(n, self.primeBase, self.primePower, self.prodPol)
        if n != m:
            return I
        G = (self | I).GaussElimTotal()
        for i in range(n):
            if G.matrix[i][i].value != 1:
                break

        invSelf = [[G.matrix[i][n + j].value for j in range(n)] for i in range(n)]
        return fieldIntMatrix(invSelf, self.primeBase, self.primePower)
    
    def __or__(self, other):
        # Concatenation of matrices
        n, m1 = self.dim
        nn, m2 = other.dim

        if n != nn:
            return -1

        C = [[self.matrix[i][j].value if (j < m1) else other.matrix[i][j - m1].value for j in range(m1 + m2)] for i in range(n)]

        return fieldIntMatrix(C, self.primeBase, self.primePower, self.prodPol)


def matrixRowExchange(row1 : int, row2: int, size: int, prime : int, power: int):
    """
    Returns a matrix that, when multiplied by the right with some
    other matrix A gives as a result the matrix A with the rows
    row1 and row2 exchanged.
    """

    one = fieldInt(1, prime, power)
    Itemp = fieldIntMatrix.matrixEye(size, prime, power)
    Itemp.matrix[row1][row1] -= one
    Itemp.matrix[row2][row2] -= one
    Itemp.matrix[row1][row2] += one
    Itemp.matrix[row2][row1] += one

    return Itemp

def matrixRowMulSum(a, row1 : int, row2 : int, size : int, prime : int, power : int):
    """
    Returns a matrix that, when multiplied by the right with some
    other matrix A gives as a result the matrix A with the row
    row1 replaced by row1 + a * row2. That is:
    row1 <- row1 + a * row2
    """
    one = fieldInt(1, prime, power)
    Itemp = fieldIntMatrix.matrixEye(size, prime, power)
    Itemp.matrix[row1][row2] += a

    return Itemp

def main():
    A = [[1, 1, 1],
         [0, 1, 1],
         [1, 0, 1]]

    A = fieldIntMatrix(A, 3, 1)
    
    I = fieldIntMatrix.matrixEye(3, 3, 1)
    
    z = [[0, 0, 1, 0, 1, 1]]
    z = fieldIntMatrix(z, 3, 1)
    
    C = (A | I).T()
    D = A * A
    E = - A
    X = (~ A)
    print(A)
    print(C, D, E)
    print(A | A.T())
    print(X)
    print(A * X)

if __name__ == "__main__":
    main()
