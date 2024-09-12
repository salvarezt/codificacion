import numpy as np
import random

class Canal:
    def __init__(self, field, probability):
        self.field = list(range(field))
        self.probability = probability

    def sendMessage(self, msg):
        npMsg = np.array(msg)
        randomMatrix = np.random.rand(*npMsg.shape)
        mask = randomMatrix < self.probability
        replacementValues = np.random.choice(self.field, size=npMsg.shape)
        resultMatrix = np.where(mask, replacementValues, npMsg)

        return [int(i) for i in resultMatrix]




canal = Canal(10,0.1)
