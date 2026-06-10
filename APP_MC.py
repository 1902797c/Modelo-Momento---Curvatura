# ====================================================================
#  ANÁLISIS MOMENTO – CURVATURA
#  Rectangular : Método Analítico  α – γ   (igual que MATLAB compañero)
#  Circular    : Método de Fibras           (igual que MATLAB compañero)
#
#  Marina Fierros Marcelo — MIAE
# ====================================================================

import numpy as np
from MODELO_MANDER import funcion_mander, Ec_concreto, eco as ECO_REF


# ────────────────────────────────────────────────────────────────────
#  CURVA DE MANDER  (vectorizada)
# ────────────────────────────────────────────────────────────────────

def _curva_mander_vec(eps_vec, fco, fcc, eco, ecc, Ec):
    """
    Devuelve (fc_confinado, fc_no_confinado) para un vector de
    deformaciones positivas.  Misma lógica que MATLAB del compañero.
    """
    esp = 0.005

    Esec_c = fcc / ecc
    r_c    = Ec / (Ec - Esec_c)
    Esec_u = fco / eco
    r_u    = Ec / (Ec - Esec_u)

    eps_s  = np.maximum(eps_vec, 1e-12)
    fc_c   = (fcc * (eps_s / ecc) * r_c) / (r_c - 1.0 + (eps_s / ecc) ** r_c)

    fc_u   = np.zeros_like(eps_vec)
    for i, e in enumerate(eps_vec):
        if e <= 2.0 * eco:
            fc_u[i] = (fco * (e / eco) * r_u) / (r_u - 1.0 + (e / eco) ** r_u)
        elif e <= esp:
            f2e     = (fco * 2.0 * r_u) / (r_u - 1.0 + 2.0 ** r_u)
            fc_u[i] = f2e * (1.0 - (e - 2.0 * eco) / (esp - 2.0 * eco))

    return fc_c, fc_u


# ────────────────────────────────────────────────────────────────────
#  POSICIÓN DEL ACERO
# ────────────────────────────────────────────────────────────────────

def _barras_rectangular(h, c, num_vars):
    """
    Lechos equidistantes entre c y h-c, 2 barras por lecho.
    Igual que MATLAB: y_lechos = linspace(dp, h-dp, num_lechos).
    y medida desde la fibra extrema de compresión.
    """
    dp        = c
    restantes = num_vars - 4
    n_interm  = restantes // 2
    n_lechos  = 2 + n_interm
    y_lechos  = np.linspace(dp, h - dp, n_lechos)

    y = []
    for yi in y_lechos:
        y.extend([yi, yi])
    if restantes % 2:
        y.append(h - dp)
    return np.array(y)


def _barras_circular(D, c, num_vars):
    """
    Barras distribuidas angularmente en el radio D/2 - c.
    y medida desde la fibra extrema de compresión (cima).
    Igual que MATLAB.
    """
    radio   = D / 2.0 - c
    angulos = np.linspace(0.0, 2.0 * np.pi, num_vars, endpoint=False)
    return D / 2.0 - radio * np.cos(angulos)


# ────────────────────────────────────────────────────────────────────
#  RECTANGULAR  —  método α-γ  (igual que calcular_M_phi_analitico MATLAB)
# ────────────────────────────────────────────────────────────────────

def _mc_rectangular(b, h, c, As_var, num_vars,
                    Pu, fco, fcc, eco, ecc, ecu, Ec,
                    fy, Es, n_puntos):

    y_acero = _barras_rectangular(h, c, num_vars)
    A_acero = np.full(len(y_acero), As_var)
    y_cent  = h / 2.0
    ecm_max = min(ecu, 0.025)
    ecm_vec = np.linspace(0.0001, ecm_max, n_puntos)

    Phi = np.zeros(n_puntos)
    M   = np.zeros(n_puntos)

    for i, ecm in enumerate(ecm_vec):

        # α y γ por integración numérica de la curva de Mander
        eps_i  = np.linspace(0.0, ecm, 200)
        fc_i, _ = _curva_mander_vec(eps_i, fco, fcc, eco, ecc, Ec)
        A_fc   = np.trapezoid(fc_i, eps_i)
        Mo_fc  = np.trapezoid(eps_i * fc_i, eps_i)

        if A_fc > 0:
            alpha = A_fc / (fcc * ecm)
            gamma = 1.0 - Mo_fc / (ecm * A_fc)
        else:
            alpha, gamma = 0.0, 0.5

        # Bisección para kd  (igual que MATLAB: c_min=1e-3, c_max=h*5, tol=100)
        c_min, c_max = 1e-3, h * 5.0
        for _ in range(100):
            c_g   = (c_min + c_max) / 2.0
            Cc    = alpha * fcc * b * c_g
            eps_s = ecm * (c_g - y_acero) / c_g
            fsi   = np.clip(Es * eps_s, -fy, fy)
            P_int = Cc + np.sum(fsi * A_acero)
            if abs(P_int - Pu) < 100.0: break
            elif P_int > Pu: c_max = c_g
            else:            c_min = c_g

        M[i]   = Cc * (y_cent - gamma * c_g) + np.sum(fsi * A_acero * (y_cent - y_acero))
        Phi[i] = ecm / c_g

    return np.concatenate([[0.0], Phi]), np.concatenate([[0.0], M])


# ────────────────────────────────────────────────────────────────────
#  CIRCULAR  —  método de fibras  (igual que calcular_M_phi_fibras MATLAB)
#
#  NOTA CLAVE: el MATLAB del compañero aplica la MISMA curva de Mander
#  a toda la sección (no separa recubrimiento), porque usa
#  fc_c para A_confinado y fc_u para A_no_confinado pero con la misma
#  función calcular_curva_mander que recibe fcc del núcleo.
#  Separar el recubrimiento con una curva que cae a 0 en eps=0.005
#  produce una discontinuidad artificial (la "panza").
#  Solución: toda la sección usa la curva confinada, truncada en ecu.
# ────────────────────────────────────────────────────────────────────

def _mc_circular(D, c, As_var, num_vars,
                 Pu, fco, fcc, eco, ecc, ecu, Ec,
                 fy, Es, n_puntos, n_capas=100):

    # Geometría de fibras — sección circular completa
    dy    = D / n_capas
    y_fib = np.linspace(dy / 2.0, D - dy / 2.0, n_capas)
    R     = D / 2.0

    A_fib = np.zeros(n_capas)
    for k in range(n_capas):
        dist       = abs(R - y_fib[k])
        A_fib[k]   = 2.0 * np.sqrt(max(R**2 - dist**2, 0.0)) * dy

    # Acero
    y_acero = _barras_circular(D, c, num_vars)
    A_acero = np.full(num_vars, As_var)
    y_cent  = R

    ecm_max = min(ecu, 0.025)
    ecm_vec = np.linspace(0.0001, ecm_max, n_puntos)

    Phi = np.zeros(n_puntos)
    M   = np.zeros(n_puntos)

    for i, ecm in enumerate(ecm_vec):

        # Bisección  (c_min=1e-3, c_max=D*5, tol=100 — igual que MATLAB)
        c_min, c_max = 1e-3, D * 5.0

        for _ in range(100):
            c_g     = (c_min + c_max) / 2.0
            phi     = ecm / c_g
            eps_fib = phi * (c_g - y_fib)

            # Curva confinada para toda la sección, truncada en ecu
            fc_fib  = np.zeros(n_capas)
            idx_c   = eps_fib > 0
            if np.any(idx_c):
                ev            = eps_fib[idx_c]
                fc_tmp, _     = _curva_mander_vec(ev, fco, fcc, eco, ecc, Ec)
                fc_tmp[ev > ecu] = 0.0          # falla del núcleo
                fc_fib[idx_c] = fc_tmp

            F_conc  = np.sum(fc_fib * A_fib)
            eps_s   = phi * (c_g - y_acero)
            fsi     = np.clip(Es * eps_s, -fy, fy)
            F_acero = np.sum(fsi * A_acero)

            P_int = F_conc + F_acero
            if abs(P_int - Pu) < 100.0: break
            elif P_int > Pu: c_max = c_g
            else:            c_min = c_g

        M[i]   = (np.sum(fc_fib * A_fib   * (y_cent - y_fib)) +
                  np.sum(fsi    * A_acero  * (y_cent - y_acero)))
        Phi[i] = ecm / c_g

    return np.concatenate([[0.0], Phi]), np.concatenate([[0.0], M])


# ────────────────────────────────────────────────────────────────────
#  FUNCIÓN PRINCIPAL  —  interfaz hacia APP_MC.py
# ────────────────────────────────────────────────────────────────────

def calcular_momento_curvatura(
        tipo_seccion: str,
        b: float = None, h_sec: float = None,
        D: float = None,
        c: float = 5.0,
        Pu: float = 0.0,
        fco: float = 250.0,
        fy: float = 4200.0,
        Es: float = 2_000_000.0,
        num_barras: int = 8,
        As_barra: float = 2.87,
        fyh: float = 4200.0,
        rho_s: float = 0.012,
        s: float = 8.0,
        ds: float = 40.0,
        esm: float = 0.12,
        tipo_confinamiento: str = "Circular con espiral",
        Asx: float = None, Asy: float = None,
        wi: np.ndarray = None, s_prima: float = None,
        n_fibras: int = 100,
        n_puntos: int = 60,
):
    """Devuelve (phi [rad/cm], M [kg·cm], info_mander)."""
    eco = ECO_REF
    Ec  = Ec_concreto(fco)
    ds_mander = (D - 2 * c) if (tipo_seccion == "circular" and D is not None) else ds

    _, _, _, _, info_est, _ = funcion_mander(
        fco, fyh, rho_s, s, ds_mander, esm, tipo_confinamiento,
        b=b, h=h_sec, c=c, Asx=Asx, Asy=Asy, wi=wi, s_prima=s_prima,
    )
    fcc = info_est["fcc"]
    ecc = info_est["ecc"]
    ecu = info_est["ecu"]

    if tipo_seccion == "rectangular":
        phi_arr, M_arr = _mc_rectangular(
            b, h_sec, c, As_barra, num_barras,
            Pu, fco, fcc, eco, ecc, ecu, Ec,
            fy, Es, n_puntos)
    else:
        phi_arr, M_arr = _mc_circular(
            D, c, As_barra, num_barras,
            Pu, fco, fcc, eco, ecc, ecu, Ec,
            fy, Es, n_puntos, n_capas=n_fibras)

    return phi_arr, np.abs(M_arr), info_est


# ────────────────────────────────────────────────────────────────────
#  PUNTOS CLAVE  (agrietamiento, fluencia, último)
# ────────────────────────────────────────────────────────────────────

def puntos_clave(phi_arr, M_arr, fco, Ec, fy, Es):
    if len(M_arr) < 4:
        return {}

    dphi = np.diff(phi_arr)
    dM   = np.diff(M_arr)
    K    = np.where(dphi > 0, dM / dphi, np.inf)

    idx_agr, idx_y = 0, len(M_arr) - 1
    for i in range(1, len(K)):
        if K[i] < 0.70 * K[0]:
            idx_agr = i
            break
    for i in range(1, len(K)):
        if K[i] < 0.05 * K[0]:
            idx_y = i
            break
    idx_u = int(np.argmax(M_arr))

    return {
        "agr": (phi_arr[idx_agr], M_arr[idx_agr]),
        "y":   (phi_arr[idx_y],   M_arr[idx_y]),
        "u":   (phi_arr[idx_u],   M_arr[idx_u]),
    }
