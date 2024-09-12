# Proyecto de Codificación

En este proyecto se implementa la teoría de códigos lineales sobre campos finitos además de F2.

Como ejemplo, se codifica texto en F27, datos numéricos en F11 e imágenes en F256.

## Explicación de cada archivo:

EncodingAndDecoding: Se implementa la logica para pasar de texto, imagenes y digitos numericos a datos que se puedan codificar correctamente por el codigo escrito en FiniteFieldCodingDecoding. Y se implementa la forma de decodificar estos datos de vuelta a su formato original.
FiniteFieldRepresentation: Se implementa la logica para trabajar en un campo finito arbitrario. (OJO: Si el campo es Fq, con q = p ** n, entonces el algoritmo actual genera correctamente un campo finito solo si n es 1 o es primo. Asi que para F256 hay que suministrar un polinomio irreducible de grado 8 que genera la regla con la cual se multiplican los elementos).
FiniteFieldMatrices: Se implementa la logica para trabajar matrices cuyas entradas sean elementos de un campo finito arbitrario.
FiniteFieldCodingDecoding: Se implementa la codificacion y decodificacion explicada en codigos lineales. Dados un mensaje x, y una matriz de paridad P, el mensaje codificado que se envia a través del canal ruidoso es (x, xP).
GallagersMethod: Se implementa el metodo de Gallager (extendedBitFlippingAlgorithm) para campos finitos arbitrarios.
NoisyChannel: Implementacion del canal ruidoso. Se le suministra al canal la probabilidad de que un digito cambie (probabilidad de error).
main: Archivo principal, aqui se pone en accion el algoritmo mostrando el proceso completo para numeros, texto e imagenes.
