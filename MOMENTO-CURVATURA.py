import numpy as np
from MODELO_MANDER import (funcion_mander, curva_mander, curva_no_confinada, Ec_concreto)
def generar_barras_rectangular(b, h, c, num_barras, As_barra):
    barras = []
    x_izq = c
    x_der = b - c
    y_inf = c
    y_sup = h - c

    # Caso mínimo
    if num_barras < 4:
        raise ValueError("Se requieren al menos 4 barras")
    # Esquinas
    barras.append((x_izq, y_inf, As_barra))
    barras.append((x_der, y_inf, As_barra))
    barras.append((x_der, y_sup, As_barra))
    barras.append((x_izq, y_sup, As_barra))

    restantes = num_barras - 4
    if restantes > 0:
        n_superior = restantes // 2
        n_inferior = restantes - n_superior
        
        if n_superior > 0:
            xs = np.linspace(x_izq, x_der, n_superior + 2)[1:-1]

            for x in xs:
                barras.append((x, y_sup, As_barra))

        if n_inferior > 0:

            xs = np.linspace(x_izq, x_der, n_inferior + 2)[1:-1]

            for x in xs:
                barras.append((x, y_inf, As_barra))
    return barras
 
def acero_bilineal(eps, fy, Es):
  ey = fy / Es

    if eps >= ey:
        return fy

    elif eps <= -ey:
        return -fy

    else:
        return Es * eps

def deformacion_fibra(phi, c_neutral, y):
    return phi * (c_neutral - y)

def discretizar_seccion_rectangular(b, h, c, n_fibras=100):
  
    fibras = []
    dy = h / n_fibras
    for i in range(n_fibras):

        y = (i + 0.5) * dy

        dentro_nucleo = (y >= c and y <= (h - c))

        if dentro_nucleo:
            ancho = b - 2*c
            fibras.append({"y": y, "area": ancho * dy, "tipo": "confinado"})
        else:
            fibras.append({"y": y, "area": b * dy, "tipo": "no_confinado"})

    return fibras
  def esfuerzo_concreto(eps, tipo, fco, eco, Ec, fcc, ecc):
    if eps <= 0:
        return 0.0
    if tipo == "confinado":

        Esec = fcc / ecc
        r = Ec / (Ec - Esec)
        x = eps / ecc
        return (fcc * x * r) / (r - 1 + x**r)

    else:
        ecu1 = 2 * eco
        ecu = 0.005

        Esec = fco / eco
        r = Ec / (Ec - Esec)

        if eps <= ecu1:
            x = eps / eco
            return (fco * x * r) / (r - 1 + x**r)

        elif eps <= ecu:

            x1 = ecu1 / eco
            f1 = (fco * x1 * r) / (r - 1 + x1**r)

            return f1 * (ecu - eps) / (ecu - ecu1)

        return 0.0
