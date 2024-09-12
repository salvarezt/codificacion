from random import randint, choice
from FiniteFieldRepresentation import fieldInt as fint
from FiniteFieldMatrices import fieldIntMatrix as fmatrix

def bitFlippingAlgorithm(A : fmatrix, x : fmatrix) -> fmatrix:
    """
    By flipping bits from the vector x this
    outputs a solution vector for the equation
    Ax = 0.
    """
    
    n, _ = x.dim
    # degree[i] will be the number of instances of
    # the variable x_i in the equations
    # represented by the matrix A.
    degree = [0 for _ in range(n)]

    Abool = A.toBoolMatrix()
    m, _ = A.dim # A is an m x n matrix.

    for i in range(m):
        for j in range(n):
            if Abool[i][j]:
                # If variable x_j appears in equation (i) then:
                degree[j] += 1 # The amount of times x_j appears gets increased.
    
    s = A * x # syndrome (column vector)
    sbool = s.toBoolMatrix()

    iszero = False

    while not iszero: # While the syndrome isn't 0

        # fail[i] is the amount of wrong equations
        # where the variable x_i appears.
        fail = [0 for _ in range(n)]
        
        iszero = True # Assume that the syndrome was 0.
        
        for i in range(m):
            """
            If the syndrome has a non-zero
            value. Then the boolean version
            of that value will be True. Since
    
            0 -> False
            otherwise -> True
            """
            if not sbool[i][0]: # If s[i] is 0, then:
                continue # Check the next coordinate

            # To get here, at least one coordinate
            # of the syndrome was non-zero
            iszero = False
            
            for j in range(n): # For each variable x_j
                if Abool[i][j]: # If the variable appears on equation (i)
                    fail[j] += 1 # Then, equation (i) failed.
                # This is because the i-th value of the syndrome is non-zero.
                # And we want the syndrome to be zero.

        if iszero: # If the syndrome is zero, then we're done.
            break
        
        failureRate = [fail[i] / degree[i] for i in range(n)]
        maxFailureRate = max(failureRate)

        if(maxFailureRate <= 0):
            print("Error :c")

        etemp = [0 for _ in range(n)]
        for i in range(n):
            if failureRate[i] == maxFailureRate:
                etemp[i] += randint(1, x.fieldSize - 1)
                break

        etemp = fmatrix([etemp], x.primeBase, x.primePower)
        x = x + etemp.T()
        
        s = A * x # Recalculate the syndrome (column vector)
        sbool = s.toBoolMatrix()
        
    return x

def solveForXi(A : fmatrix, x : fmatrix, i : int) -> tuple[fint, int]:
    """
    Given the matrix A and the vector x,
    this function solves for the component x_i in the
    equations given by Ax = 0.

    It returns the solution that corrects
    the greatest amount of equations.
    """

    n, _ = x.dim
    prime, power = x.primeBase, x.primePower
    
    if i >= n:
        return fint(0, prime, power, x.prodPol)

    # Current value of the i-th component
    currentValue = x.matrix[i][0]
    
    s = A * x # sindrome
    numberOfEquations, _ = s.dim

    d = 0
    connectedEquations = []
    
    for equation in range(numberOfEquations):
        c = A.matrix[equation][i]
        if int(c) != 0:
            d += 1
            connectedEquations.append((equation, c))

    if d == 0:
        # Changing the current value won't change
        # anything in the syndrome.
        return (currentValue, d) 
    
    potencialValues = []
    tmax = 0
    
    for equation, coefficient in connectedEquations:
        potencialValue = s.matrix[equation][0] - coefficient * currentValue
        potencialValue = - potencialValue // coefficient
        
        x2 = fmatrix(x.toIntMatrix(), prime, power, x.prodPol)
        x2.matrix[i][0] = potencialValue

        s2 = A * x2
        t = 0 # Satisfied equations by the potencialValue
        
        for eq, _ in connectedEquations:
            if int(s2.matrix[eq][0]) == 0:
                t += 1

        potencialValues.append((t, potencialValue))

        if t > tmax:
            tmax = t

    goodValues = list(filter(lambda x: x[0] == tmax, potencialValues))
    goodValues = list(map(lambda x: x[1], goodValues))
    
    if currentValue in goodValues:
        return (currentValue, tmax)
    
    return (choice(goodValues), tmax)

    
def extendedBitFlippingAlgorithm(A : fmatrix, x : fmatrix) -> fmatrix:
    """
    By flipping digits from the vector x this
    outputs a solution vector for the equation
    Ax = 0 in the field Fq.
    """

    n, _ = x.dim
    # degree[i] will be the number of instances of
    # the variable x_i in the equations
    # represented by the matrix A.
    degree = [0 for _ in range(n)]

    Abool = A.toBoolMatrix()
    m, _ = A.dim # A is an m x n matrix.

    for i in range(m):
        for j in range(n):
            if Abool[i][j]:
                # If variable x_j appears in equation (i) then:
                degree[j] += 1 # The amount of times x_j appears gets increased.
    
    s = A * x # syndrome (column vector)

    while s.weight() != 0: # While the syndrome isn't 0

        potencialChoices = []
        for i in range(n):
            potencialChoices.append(solveForXi(A, x, i))

        successRates = [- potencialChoices[i][1] / degree[i] if potencialChoices[i][0] == x.matrix[i][0] else potencialChoices[i][1] / degree[i] for i in range(n)]
        maxSuccessRate = max(successRates)

        # print(maxSuccessRate, end=' ')
        
        if maxSuccessRate > 0:
            for i in range(n):
                if successRates[i] == maxSuccessRate:
                    x.matrix[i][0] = potencialChoices[i][0]
                    break
        else:
            pass
        
        s = A * x # Recalculate the syndrome (column vector)
        
    return x

def main():

    # F11 (enviar numeros de tel, cedulas) 0..9 " "
    # F27 (alfabeto usual) 0..25 " "
    # F2, F16, F256 (Imagenes)
    # F2 -> 24-tuplas
    # F16 -> 6-tuplas 16a + b, 0..15
    # F256 -> 3-tuplas

    """
    # An [12, 10, 3] code in F11:
    P = [[10, 10],
         [10, 9],
         [10, 8],
         [10, 7],
         [10, 6],
         [10, 5],
         [10, 4],
         [10, 3],
         [10, 2],
         [10, 1]]
    """

    """
    # A [4, 2, 3] code in F3:
    P = [[2, 2],
        [2, 1]]
    """

    P = [[1, 1 + i] for i in range(31)]
    
    prime = 2
    power = 8
    
    P = fmatrix(P, prime, power, [1, 0, 1, 1, 1, 0, 0, 0, 1])

    n, m = P.dim
    I = fmatrix.matrixEye(n, prime, power, [1, 0, 1, 1, 1, 0, 0, 0, 1])
    G = (I | P)

    #print(G)

    I = fmatrix.matrixEye(m, prime, power, [1, 0, 1, 1, 1, 0, 0, 0, 1])
    H = (- P.T() | I)
    
    print(H) # La version booleana de la matriz H.
    
    #x = [[5, 0, 5, 1, 0, 10, 2, 7, 5, 2]]
    x = [[0 for i in range(31)]]
    x = fmatrix(x, prime, power, [1, 0, 1, 1, 1, 0, 0, 0, 1])
    xG = x * G
    #print(x) # El mensaje original y el extendido.
    
    #e = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2]]
    e = [[0 for i in range(33)]]
    e[0][20] = 23
    e = fmatrix(e, prime, power, [1, 0, 1, 1, 1, 0, 0, 0, 1])
    
    print(xG, xG + e)
    #print(bitFlippingAlgorithm(H, (xG + e).T()).T())

    print((extendedBitFlippingAlgorithm(H, (xG + e).T())).T())
    
    return 0

if __name__ == '__main__':
    main()
