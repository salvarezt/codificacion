from FiniteFieldCodingDecoding import encode, decode
from FiniteFieldMatrices import fieldIntMatrix as fmatrix
from NoisyChannel import Canal
import EncodingAndDecoding as code

from PIL import Image

def main():

    Canal11 = Canal(11, 0.1)
    Canal27 = Canal(27, 0.1)
    Canal256 = Canal(256, 0.01)
    
    P11 = [[10, 10 - i] for i in range(8)] # Numbers
    P27 = [[26, 26 - i] for i in range(8)] # Text
    P256 = [[255, 255 - i] for i in range(25)] # Images

    P11 = fmatrix(P11, 11, 1) # 11 = 11 ** 1
    P27 = fmatrix(P27, 3, 3) # 27 = 3 ** 3
    P256 = fmatrix(P256, 2, 8, [1, 0, 1, 1, 1, 0, 0, 0, 1]) # 256 = 2 ** 8

    number = '1027522564 3112025130'
    text = "Este es un mensaje bastante importante asi que debe enviarse sin errores"
    
    print("Numero a enviar:", number)
    print("Texto a enviar:", text, '\n')
    
    encodeNumber = encode(code.numberToMessage(number), P11)
    encodeText = encode(code.textToMessage(text), P27)
    
    print("Numero codificado:", code.messageToNumber(encodeNumber))
    print("Texto codificado:", code.messageToText(encodeText), '\n')
    
    recievedNumber = Canal11.sendMessage(encodeNumber)
    recievedText = Canal27.sendMessage(encodeText)

    print("Numero recibido:",code.messageToNumber(recievedNumber))
    print("Texto recibido:",code.messageToText(recievedText), '\n')

    decodeNumber = code.messageToNumber(decode(recievedNumber, P11))
    decodeText = code.messageToText(decode(recievedText, P27))
    
    print("Numero recuperado:", decodeNumber)
    print("Texto recuperado:", decodeText)

    print("\nAhora a enviar la imagen")
    
    img = Image.open("Alice_geez2.png")
    sizeOfImg = img.size
    
    print("Tamaño de la imagen en pixeles:", sizeOfImg)

    firstEncodingImage = code.imageToMessage(img)
    encodeImage = encode(firstEncodingImage, P256)
    
    # Genera la imagen codificada (opcional). El 8 solo funciona si la imagen es 100x100 px
    code.messageToImage(encodeImage, sizeOfImg[0] + 8, sizeOfImg[1], 'A2encoding')
    print("¡Imagen codificada!")
        
    recievedImage = Canal256.sendMessage(encodeImage)
    
    # Genera la imagen recibida (opcional). El 8 solo funciona si la imagen es 100x100 px
    code.messageToImage(recievedImage, sizeOfImg[0] + 8, sizeOfImg[1], 'A2recieved')
    print("¡Imagen recibida!")
    
    decodedImage = decode(recievedImage, P256)

    # Genera la imagen decodificada
    code.messageToImage(decodedImage, sizeOfImg[0], sizeOfImg[1], 'A2decoding')
    print("¡La imagen recibida ha sido decodificada!")
    
    return 0

if __name__ == "__main__":
    main()
