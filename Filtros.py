import matplotlib.pyplot as plt
import numpy as np

def lowpass(data,Delta_t,Fc):
    alpha = (2*np.pi*Delta_t*Fc)/(2*np.pi*Delta_t*Fc+1)
    y = np.zeros(len(data))
    y[0] = data[0]
    for i in range(1,len(y)):
        y[i] = alpha*data[i] + (1-alpha)*y[i-1]
    return y

def Media_Movel(x_n):
    y=[]

    for i in range(0,len(x_n)):
        resultado = (x_n[i]-x_n[i-1])/2
        y.append(resultado)
    return y

def highpass(data,Delta_t,Fc):
    alpha = (1)/(2*np.pi*Delta_t*Fc+1)
    y = np.zeros(len(data))
    y[0] = data[0]
    for i in range(1,len(y)):
        y[i] = alpha*y[i-1] + alpha*(data[i]-data[i-1])
    return y
    