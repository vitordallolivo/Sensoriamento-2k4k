import matplotlib.pyplot as plt

def lowpass(x_n,dt,RC,IGNORE):
    y=[]
    alpha = dt/(RC+dt)
    y.append(alpha*x_n[0])

    for i in range(1,len(x_n)):
        resultado=alpha*x_n[i]+(1-alpha)*y[i-1]
        if (resultado>IGNORE[0] and resultado<IGNORE[1]):
            resultado=0
        y.append(resultado)
    return y

def media(x_n,dt,RC,IGNORE):
    y=[]

    for i in range(1,len(x_n)):
        resultado = (x_n[i]-x_n[i-1])/2
        y.append(resultado)
    return y

def highpass(x_n,dt,RC,IGNORE):
    y=[]
    alpha = RC/(RC+dt)
    y.append(x_n[0])

    for i in range(1,len(x_n)):
        resultado=alpha*y[i-1]+alpha*(x_n[i]-x_n[i-1])
        if (resultado>IGNORE[0] and resultado<IGNORE[1]):
            resultado=0
        y.append(resultado)

    return y