from FiniteFieldRepresentation import fieldInt as fint
from FiniteFieldMatrices import fieldIntMatrix as fmatrix
from GallagersMethod import extendedBitFlippingAlgorithm as gallager

def encode(message : list[int], P : fmatrix) -> list[int]:
    """
    Given a message represented with numbers from 0 to q - 1
    for a fixed q, and a codification matrix P from the field Fq.

    Returns the encoded message represented with another list of numbers.
    """

    n, m = P.dim
    prime, power = P.primeBase, P.primePower

    In = fmatrix.matrixEye(n, prime, power)
    G = (In | P)
    
    sizeMessage = len(message)
    simbolsMissing = sizeMessage % n

    for _ in range((n - simbolsMissing) % n):
        message.append(0)

    encodedMessage = []

    #print("Encoding In process")
    
    count = 0
    for i in range(len(message) // n):
        x = [message[i * n : (i + 1) * n]]
        x = fmatrix(x, prime, power)
        
        xG = x * G # encoded Message
        
        encodedMessage += xG.toIntMatrix()[0]
        if i > count * 0.1 * len(message) // n:
            #print(f"{count}0%")
            count += 1


    return encodedMessage


def decode(messageWithErrors : list[int], P : fmatrix) -> list[int]:
    """
    Given a message that may or may not have errors,
    and the codification matrix P.

    Returns the original message provided that there aren't enough errors
    in the codification.

    """

    n, m = P.dim
    prime, power = P.primeBase, P.primePower
    
    Im = fmatrix.matrixEye(m, prime, power)
    H = (- P.T() | Im)

    decodedMessage = []

    print("Decoding In process")
    
    count = 0
    for i in range(len(messageWithErrors) // (n + m)):
        x = [messageWithErrors[i * (n + m) : (i + 1) * (n + m)]]
        x = fmatrix(x, prime, power, P.prodPol)
        
        xNoErrors = gallager(H, x.T()).T()
        
        decodedMessage += xNoErrors.toIntMatrix()[0][:n]

        if i > count * 0.1 * len(messageWithErrors) // (n + m):
            print(f"{count}0%")
            count += 1

    return decodedMessage

def main():
    prime, power = 11, 1
    
    P = [[10, 10],
         [10, 9],
         [10, 8],
         [10, 7],
         [10, 6]]
    P = fmatrix(P, prime, power)
    

    msg = [1, 0, 2, 7, 5, 2, 2, 5, 6, 4, 3, 1, 1, 2, 0, 2, 5, 1, 3, 0]
    code = encode(msg, P)
    code2 = code.copy()
    code2[1] += 5
    msg2 = decode(code2, P)

    print(msg)
    print(code)
    print(code2)
    print(msg2)

    return 0

if __name__ == "__main__":
    main()
