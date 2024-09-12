from FiniteFieldCodingDecoding import encode, decode
from FiniteFieldMatrices import fieldIntMatrix as fmatrix
from PIL import Image
from numpy import array

def textToMessage(text : str) -> list[int]:
    """
    Given a text, returns the same message represented
    by numbers from 0 to 26 with the following encoding:
    1 -> A, 2 -> B, ..., 26 -> Z and 0 represents anything else.
    """

    t = text.upper()
    msg = [ord(c) - ord('A') + 1 if 0 <= ord(c) - ord('A') < 26 else 0 for c in t]
    return msg

def messageToText(msg : list[int]) -> str:
    """
    Given a list of numbers from 0 to 26, it returns a text
    inverting the encoding given by the function textToMessage.
    """

    text = [chr(i + ord('A') - 1) if 1 <= i <= 26 else ' ' for i in msg]
    return ''.join(text)

def numberToMessage(b : str) -> list[int]:
    """
    Given a number b in base 10, we store
    each of its digits in a list.
    """
    
    numbers = [ord('0') + i for i in range(10)]
    msg = [int(d) if ord(d) in numbers else 10 for d in b]
    return msg

def messageToNumber(msg : list[int]) -> str:
    """
    Given a list of numbers from 0 to 9, we
    return the number in base 10 that they
    represent.
    """

    num = ''.join([str(i) if i != 10 else ' ' for i in msg])
    
    return num

def imageToMessage(img : Image) -> list[int]:
    """
    Given an image, it returns a list whose values are given
    by the following idea:

    If the pixel values are R, G, B and A then:
    msg = [R(0, 0), G(0, 0), B(0, 0), A(0, 0), R(1, 0),
          ..., B(width - 1, 0), A(width - 1), R(0, 1), ...]
    
    """

    rgba = img.convert('RGBA')
    width, height = img.size
    rgba_matrix = array(rgba)
    
    message = []

    for j in range(height):
        for i in range(width):
            message += list(rgba_matrix[i][j])

    return message

def messageToImage(msg : list[int], width : int, height : int, name = 'new') -> Image:
    """
    Given a message encoded by the imageToMessage function,
    it returns the corresponding image.
    """

    size = len([str(i) for i in msg])

    #print(size)
    
    if size != width * height * 4:
        print("Error! The dimentions of the image given aren't correct")

    pixels = [[[0, 0, 0, 0] for _ in range(height)] for __ in range(width)]
    for i in range(size // 4):
        pixels[i % width][i // width] = array(msg[i * 4 : (i + 1) * 4]).astype('uint8')

    for i in range(width):
        pixels[i] = array(pixels[i])

    image = Image.fromarray(array(pixels), mode="RGBA")
    image.save(name + '.png')

    return image

def main():

    P = [[255, 255 - i] for i in range(62)]
    P = fmatrix(P, 2, 8)
    
    a = "Hola Mundo"
    ta = textToMessage(a)
    mta = messageToText(ta)
    print(a)
    print(ta)
    print(mta, end='\n\n')

    b = '10275225641001098141'
    nb = numberToMessage(b)
    mnb = messageToNumber(nb)
    print(b)
    print(nb)
    print(mnb, end='\n\n')

    c = Image.open("Alice_geez.png")
    ic = imageToMessage(c)
    encodeic = encode(ic, P)
    mencodeic = messageToImage(encodeic, 64, 62, 'encoded')
    decodeic = decode(encodeic, P)
    mic = messageToImage(decodeic, 62, 62, 'decoded')
    
    
    #mic = messageToImage(ic, *c.size)
    #print(c.format, c.size, c.mode)

    
    
    
    return 0

if __name__ == "__main__":
    main()
