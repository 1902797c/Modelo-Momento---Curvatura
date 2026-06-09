# ====================================================================
#  ANÁLISIS MOMENTO – CURVATURA POR FIBRAS
#  Marina Fierros Marcelo — MIAE
# ====================================================================
 
import numpy as np
#from scipy.optimize import brentq
from MODELO_MANDER import (funcion_mander, Ec_concreto, eco as ECO_REF)
 
# ────────────────────────────────────────────────────────────────────
#    Modeli elasto-plástico perfecto
# ────────────────────────────────────────────────────────────────────
 
def acero_bilineal(eps: float, fy: float, Es: float) -> float:
    ey = fy / Es
    if eps >= ey:
        return fy
    elif eps <= -ey:
        return -fy
    else:
        return Es * eps
 
 
# ────────────────────────────────────────────────────────────────────
# COMPATIBILIDAD DE DEFORMACIONES
# ────────────────────────────────────────────────────────────────────
 
def deformacion_fibra(phi: float, kd: float, y: float) -> float:
    return phi * (kd - y)
 
 
# ────────────────────────────────────────────────────────────────────
#  CONCRETO
# ────────────────────────────────────────────────────────────────────
 
def esfuerzo_concreto(eps: float, tipo: str,
                      fco: float, eco: float, Ec: float,
                      fcc: float, ecc: float, ecu: float) -> float:
    
    if eps <= 0.0:          # fibra en tensión → concreto no sirve
        return 0.0
 
    if tipo == "confinado":
        if eps > ecu:       # falla del núcleo confinado
            return 0.0
        Esec = fcc / ecc
        r    = Ec / (Ec - Esec)
        x    = eps / ecc
        return (fcc * x * r) / (r - 1.0 + x ** r)
 
    else:   
        ecu1 = 2.0 * eco    #  0.004  
        ecu_nc = 0.005      
 
        if eps > ecu_nc:
            return 0.0
 
        Esec = fco / eco
        r    = Ec / (Ec - Esec)
 
        if eps <= ecu1:
            x = eps / eco
            return (fco * x * r) / (r - 1.0 + x ** r)
        else:
            x1 = ecu1 / eco
            f1 = (fco * x1 * r) / (r - 1.0 + x1 ** r)
            return f1 * (ecu_nc - eps) / (ecu_nc - ecu1)
 
 
# ────────────────────────────────────────────────────────────────────
# DISCRETIZACIÓN 
# ────────────────────────────────────────────────────────────────────
 
def discretizar_rectangular(b: float, h: float, c: float,
                             n_fibras: int = 200):
    dy     = h / n_fibras
    fibras = []
    for i in range(n_fibras):
        y_centro = (i + 0.5) * dy

        en_nucleo = (y_centro >= c) and (y_centro <= h - c)
 
        if en_nucleo:
            ancho_conf   = b - 2.0 * c
            ancho_noconf = 2.0 * c         
        else:
            ancho_conf   = 0.0
            ancho_noconf = b                
 
        fibras.append({
            "y":            y_centro,
            "dy":           dy,
            "ancho_conf":   ancho_conf,
            "ancho_noconf": ancho_noconf,})
    return fibras
 
def discretizar_circular(D:float,c:float,n_fibras:int = 200):
    R  = D / 2.0
    Rc = R - c          # radio del núcleo confinado
    dy = D / n_fibras
    fibras = []
 
    for i in range(n_fibras):
        y_centro = (i + 0.5) * dy
        y_local  = y_centro - R   # coordenada desde el centro del círculo
 
        arg_total = R ** 2 - y_local ** 2
        if arg_total <= 0:
            b_total = 0.0
        else:
            b_total = 2.0 * np.sqrt(arg_total)
 
        arg_nucleo = Rc ** 2 - y_local ** 2
        if arg_nucleo <= 0:
            b_nucleo = 0.0
        else:
            b_nucleo = 2.0 * np.sqrt(arg_nucleo)
 
        ancho_noconf = max(b_total - b_nucleo, 0.0)
 
        fibras.append({
            "y":            y_centro,
            "dy":           dy,
            "ancho_conf":   b_nucleo,
            "ancho_noconf": ancho_noconf,
        })
    return fibras
 
 
# ────────────────────────────────────────────────────────────────────
# POSICIÓN ACERO LONGITUDINAL
# ────────────────────────────────────────────────────────────────────
 
def generar_barras_rectangular(b: float, h: float, c: float,
                                num_barras: int, As_barra: float):
  
    if num_barras < 4:
        raise ValueError("Se requieren al menos 4 barras longitudinales.")
 
    y_sup = c           # fibra de compresión (y = 0 es la cara superior)
    y_inf = h - c       # fibra de tensión
 
    barras = []
 
    # Esquinas
    for _ in range(2):
        barras.append((y_sup, As_barra))   # capa de compresión
    for _ in range(2):
        barras.append((y_inf, As_barra))   # capa de tensión
 
    restantes = num_barras - 4
    if restantes > 0:
        n_sup = restantes // 2
        n_inf = restantes - n_sup
 
        # Barras intermedias en capa de compresión
        for _ in range(n_sup):
            barras.append((y_sup, As_barra))
        # Barras intermedias en capa de tensión
        for _ in range(n_inf):
            barras.append((y_inf, As_barra))
 
    return barras   # [(y_cm, As_cm²), ...]
 
 
def generar_barras_circular(D: float, c: float, num_barras: int, As_barra: float):
  
    R_barras = D / 2.0 - c   # radio del círculo de barras
    angulos  = np.linspace(0, 2 * np.pi, num_barras, endpoint=False)
    barras   = []
    for theta in angulos:
        # y_local = -R_barras·cos(θ)  →  +R desplazamiento hasta y=0 en cima
        y = D / 2.0 - R_barras * np.cos(theta)
        barras.append((y, As_barra))
    return barras   # [(y_cm, As_cm²), ...]
 
 
# ────────────────────────────────────────────────────────────────────
#  EQUILIBRIO
# ────────────────────────────────────────────────────────────────────
 
def fuerzas_internas(kd: float, phi: float,
                     fibras: list, barras: list,
                     fco: float, eco: float, Ec: float,
                     fcc: float, ecc: float, ecu: float,
                     fy: float, Es: float,
                     h: float):
    N = 0.0
    M = 0.0
 
    # --- Contribución del concreto ---
    for f in fibras:
        y    = f["y"]
        dy   = f["dy"]
        eps  = deformacion_fibra(phi, kd, y)   # compresión > 0
 
        # Núcleo confinado
        if f["ancho_conf"] > 0:
            fc   = esfuerzo_concreto(eps, "confinado",
                                     fco, eco, Ec, fcc, ecc, ecu)
            dA   = f["ancho_conf"] * dy
            dF   = fc * dA
            N   += dF
            M   += dF * (h / 2.0 - y)
 
        # Recubrimiento no confinado
        if f["ancho_noconf"] > 0:
            fc   = esfuerzo_concreto(eps, "no_confinado",
                                     fco, eco, Ec, fcc, ecc, ecu)
            dA   = f["ancho_noconf"] * dy
            dF   = fc * dA
            N   += dF
            M   += dF * (h / 2.0 - y)
 
    # --- Contribución del acero ---
    for (y_b, As_b) in barras:
        eps_s = deformacion_fibra(phi, kd, y_b)
        fs    = acero_bilineal(eps_s, fy, Es)
        dF    = fs * As_b
        N    += dF
        M    += dF * (h / 2.0 - y_b)
 
    return N, M
 
 
# ────────────────────────────────────────────────────────────────────
#   curva M–φ
# ────────────────────────────────────────────────────────────────────
 
def calcular_momento_curvatura(

        tipo_seccion: str,          
        # Rectangular
        b: float = None, h_sec: float = None,
        # Circular
        D: float = None,
        # Recubrimiento 
        c: float = 5.0,
        # Carga axial
        Pu: float = 0.0,            # kg  (positivo = compresión)
        # Materiales — concreto
        fco: float = 250.0,         # kg/cm²
        # Materiales — acero longitudinal
        fy: float = 4200.0,         # kg/cm²
        Es: float = 2_000_000.0,    # kg/cm²
        num_barras: int = 8,
        As_barra: float = 2.87,     # cm²  
        fyh: float = 4200.0,
        rho_s: float = 0.012,
        s: float = 8.0,
        ds: float = 40.0,
        esm: float = 0.12,
        tipo_confinamiento: str = "Circular con espiral",
        Asx: float = None, Asy: float = None,
        wi: np.ndarray = None, s_prima: float = None,
        n_fibras: int = 200,
        n_puntos: int = 80):
    
 
    # ── Propiedades del concreto confinado (Mander) ──────────
    eco = ECO_REF      # 0.002
    Ec  = Ec_concreto(fco)
 
    # diámetro del núcleo)
    if tipo_seccion == "circular" and D is not None:
        ds_mander = D - 2 * c   
    else:
        ds_mander = ds
 
    _, _, _, _, info_est, _ = funcion_mander(
        fco, fyh, rho_s, s, ds_mander, esm,
        tipo_confinamiento,
        b=b, h=h_sec, c=c,
        Asx=Asx, Asy=Asy, wi=wi, s_prima=s_prima)
 
    fcc = info_est["fcc"]
    ecc = info_est["ecc"]
    ecu = info_est["ecu"]
 
    # ── Dimensión de altura para momentos ───────────────────────────
    if tipo_seccion == "rectangular":
        h = h_sec
    else:
        h = D
 
    # ── Discretización de la sección ────────────────────────────────
    if tipo_seccion == "rectangular":
        fibras = discretizar_rectangular(b, h, c, n_fibras)
        barras = generar_barras_rectangular(b, h, c, num_barras, As_barra)
    else:
        fibras = discretizar_circular(D, c, n_fibras)
        barras = generar_barras_circular(D, c, num_barras, As_barra)
 
    # ── Loop principal sobre εcm ─────────────────────────────────────
 
    phi_arr = []
    M_arr   = []
 
    ecm_max  = min(ecu, 0.025)
    ecm_grid = np.linspace(1e-4, ecm_max, n_puntos)
 
    # Rango de kd: extendido hasta 5·h para cubrir alto Pu o flexión pura
    KD_SEARCH = np.linspace(0.01 * h, 5.0 * h, 150)
 
    for ecm in ecm_grid:
 
        # Función de desequilibrio: N_interno(kd) − Pu = 0
        def desequilibrio(kd_trial, _ecm=ecm):
            if kd_trial <= 0:
                return -abs(Pu) - 1.0
            phi_trial = _ecm / kd_trial
            N, _ = fuerzas_internas(
                kd_trial, phi_trial,
                fibras, barras,
                fco, eco, Ec, fcc, ecc, ecu,
                fy, Es, h,
            )
            return N - Pu
 
        # Buscar el intervalo [kd_a, kd_b] donde hay cambio de signo
        f_vals = np.array([desequilibrio(k) for k in KD_SEARCH])
 
        kd_a = kd_b = None
        for i in range(len(f_vals) - 1):
            if f_vals[i] * f_vals[i + 1] <= 0:
                kd_a = KD_SEARCH[i]
                kd_b = KD_SEARCH[i + 1]
                break
 
        if kd_a is None:
            # Sin cruce para este εcm — omitir (continue, no break)
            continue
 
       # try:
            # kd_sol = brentq(desequilibrio, kd_a, kd_b, xtol=1e-7, maxiter=100)
        except ValueError:
            continue
 
        phi_sol = ecm / kd_sol
        _, M_sol = fuerzas_internas(
            kd_sol, phi_sol,
            fibras, barras,
            fco, eco, Ec, fcc, ecc, ecu,
            fy, Es, h)
 
        phi_arr.append(phi_sol)
        M_arr.append(abs(M_sol))
 
    return np.array(phi_arr), np.array(M_arr), info_est
 
 
# ────────────────────────────────────────────────────────────────────
#  Agrietamiento, fluencia y último
# ────────────────────────────────────────────────────────────────────
 
def puntos_clave(phi_arr: np.ndarray, M_arr: np.ndarray,
                 fco: float, Ec: float,
                 fy: float, Es: float) -> dict:

    if len(M_arr) < 4:
        return {}
 
    # Rigidez secante en cada segmento
    dphi = np.diff(phi_arr)
    dM   = np.diff(M_arr)
    K    = np.where(dphi > 0, dM / dphi, np.inf)
 
    # Agrietamiento
    idx_agr = 0
    for i in range(1, len(K)):
        if K[i] < 0.70 * K[0]:
            idx_agr = i
            break
 
    # Fluencia
    idx_y = len(M_arr) - 1
    for i in range(1, len(K)):
        if K[i] < 0.05 * K[0]:
            idx_y = i
            break
 
    # ── Último: máximo de M 
    idx_u = int(np.argmax(M_arr))
 
    return {
        "agr": (phi_arr[idx_agr], M_arr[idx_agr]),
        "y":   (phi_arr[idx_y],   M_arr[idx_y]),
        "u":   (phi_arr[idx_u],   M_arr[idx_u])}
