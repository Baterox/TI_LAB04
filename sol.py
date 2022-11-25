# Librerías necesarias
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

def ajuste_lineal(distancia, medicion):
    sum = [0, 0, 0, 0]
    sigma = 0
    n = len(distancia)
    for i in range(n):
        sum[1] += medicion[i]
        sum[2] += distancia[i] * medicion[i]
        sum[0] += distancia[i]
        sum[3] += distancia[i]**2
    
    mu = np.round(sum[1]/n, decimals=2)

    for i in range(len(medicion)):
        sigma += (medicion[i]-mu)**2

    sigma = np.sqrt(sigma/n)

    m = (sum[2]-sum[0]*sum[1]/n)/(sum[3]-(sum[0]**2)/n)
    b = (sum[1]/n) - m*(sum[0]/n)
    return ([m*d + b for d in distancia], np.round(m, decimals=2), np.round(b, decimals=2), mu, sigma)

def PL(sigma, frecuencia, distancia, mediciones):
    # Calculo de lambda
    lamb = (3*np.exp(8)) / frecuencia
    # Contante N
    Ns = 2
    # d0
    d0 = distancia[0]
    # PL(d0)
    PL_Sd0 = 10 * Ns * np.log10((4*np.pi*d0) / lamb)
    # Ruido Gaussiano X
    X = sigma * np.random.randn(len(mediciones))
    # Perdida de potencia log-shadow
    return [ PL_Sd0 + 10 * Ns * np.log10(distancia[i] / d0) + X[i] for i in range(len(X))]

def PR(Gt, Gr, L, Pt, PL):
    Z = Pt + Gt + Gr - L
    return [Z - pl for pl in PL]

def CDF(sigma, mu, N):
    data = sigma*np.random.randn(N) + mu
    count, bins_count = np.histogram(data, bins=N)
    pdf = count / sum(count)
    return (bins_count[1:], np.cumsum(pdf))

def graficar(distancia, curvas, apariencia, titulo, legendas, xlabel, ylabel):
    for i in range(len(curvas)):
        plt.plot(distancia, curvas[i], apariencia[i])
    plt.title(titulo)
    plt.legend(legendas)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    plt.show()

if __name__ == "__main__":
    # Datos
    frecuencia03 = 2640*np.exp(6)
    frecuencia04 = 2650*np.exp(6)
    # Extracción de la data contenida en el dataset (mediciones lab 03 y 04)
    file = open('./datalab.csv')
    file.readline()

    mediciones = [row.replace("\n", "").split(",") for row in file]
    distancialog = [float(i[1]) for i in mediciones]
    medicioneslab03 = [float(i[2]) for i in mediciones]
    medicioneslab04 = [float(i[3]) for i in mediciones]

    # Ajuste lineal de las mediciones del lab 03 y 04
    ajustelab03, ml3, bl3, medial3, sigmal3 = ajuste_lineal(distancialog, medicioneslab03)
    graficar(distancialog,
            [medicioneslab03, ajustelab03],
            ["bo", "r"],
            "Ajuste Lineal - Laboratorio 03",
            ["Mediciones Laboratorio 03", "Ajuste Lineal {}d{}{}".format(ml3, '-' if bl3 < 0 else '+', np.abs(bl3))],
            "d [10log(m)]",
            "Potencia [dBm]")

    ajustelab04, ml4, bl4, medial4, sigmal4 = ajuste_lineal(distancialog, medicioneslab04)
    graficar(distancialog,
            [medicioneslab04, ajustelab04],
            ["go", "r"],
            "Ajuste Lineal - Laboratorio 04",
            ["Mediciones Laboratorio 04", "Ajuste Lineal {}d{}{}".format(ml4, '-' if bl4 < 0 else '+', np.abs(bl4))],
            "d [10log(m)]",
            "Potencia [dBm]")

    graficar(distancialog,
            [medicioneslab03, ajustelab03, medicioneslab04, ajustelab04],
            ["bo", "r", "go", "r"],
            "Ajuste Lineal - Laboratorio 03 y 04",
            ["Mediciones Laboratorio 03", "Ajuste Lineal {}d{}{}".format(ml3, '-' if bl3 < 0 else '+', np.abs(bl3)), "Mediciones Laboratorio 04", "Ajuste Lineal {}d{}{}".format(ml4, '-' if bl4 < 0 else '+', np.abs(bl4))],
            "d [10log(m)]",
            "Potencia [dBm]")

    distancialog_lab34 = []
    mediciones_lab34 = []
    mediciones_lab3 = []
    mediciones_lab4 = []
    for i in range(len(distancialog)):
        if medicioneslab03[i] != None:
            distancialog_lab34.append(distancialog[i])
            mediciones_lab34.append(medicioneslab03[i])
            mediciones_lab3.append(medicioneslab03[i])
            mediciones_lab4.append(None)
        if medicioneslab04[i] != None:
            distancialog_lab34.append(distancialog[i])
            mediciones_lab34.append(medicioneslab04[i])
            mediciones_lab4.append(medicioneslab04[i])
            mediciones_lab3.append(None)

    ajustelab34, ml34, bl34, medial34, sigmal34 = ajuste_lineal(distancialog_lab34, mediciones_lab34)
    graficar(distancialog_lab34,
            [mediciones_lab3, mediciones_lab4, ajustelab34],
            ["bo", "go", "r"],
            "Ajuste Lineal - Laboratorio 03 + 04",
            ["Mediciones Laboratorio 03", "Mediciones Laboratorio 04", "Ajuste Lineal {}d{}{}".format(ml34, '-' if bl34 < 0 else '+', np.abs(bl34))],
            "d [10log(m)]",
            "Potencia [dBm]")

    print("Curva\tPendiente\tMedia\tDesviación Estándar")
    print("Lab03\t{}\t\t{}\t{}".format(ml3, medial3, sigmal3))
    print("Lab04\t{}\t\t{}\t{}".format(ml4, medial4, sigmal4))
    print("Lab3+4\t{}\t\t{}\t{}".format(ml34, medial34, sigmal34))

    # Path loss de cada curva
    PL03 = PL(sigmal3, frecuencia03, distancialog, medicioneslab03)
    PL04 = PL(sigmal4, frecuencia04, distancialog, medicioneslab04)

    graficar(distancialog,
            [PL03, PL04],
            ["bo", "go"],
            "Path Loss - Laboratorio 03 y 04",
            ["Path Loss Laboratorio 03", "Path Loss Laboratorio 04"],
            "d [10log(m)]",
            "Potencia [dBm]")

    xcdf03, ycdf03 = CDF(sigmal3, medial3, len(distancialog))
    xcdf04, ycdf04 = CDF(sigmal4, medial4, len(distancialog))
    xcdf34, ycdf34 = CDF(sigmal34, medial34, len(distancialog))
    
    graficar(xcdf03,
            [ycdf03],
            ["b"],
            "CDF - Lab 03",
            ["CDF"],
            "Potencia [dBm]",
            "Probabilidad")

    graficar(xcdf04,
            [ycdf04],
            ["g"],
            "CDF - Lab 04",
            ["CDF"],
            "Potencia [dBm]",
            "Probabilidad")

    graficar(xcdf34,
            [ycdf34],
            ["r"],
            "CDF - Lab 3+4",
            ["CDF"],
            "Potencia [dBm]",
            "Probabilidad")

    plt.plot(xcdf03, ycdf03)
    plt.plot(xcdf04, ycdf04)
    plt.plot(xcdf34, ycdf34)
    plt.legend(["CDF Lab 03", "CDF Lab 04", "CDF Lab 3+4"])
    plt.title("CDF")
    plt.xlabel("Potencia [dBm]")
    plt.ylabel("Probabilidad")
    plt.grid()
    plt.show()