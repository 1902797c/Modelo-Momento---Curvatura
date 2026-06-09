import numpy as np

# ====================================================================
#  MODELO DE MANDER (1988) — CON Y SIN EFECTO DE TASA DE DEFORMACIÓN
#		MARINA FIERROS MARCELO - MIAE
# ====================================================================


# ----------------------------------------------------------------
#  			CONSTANTES DE CONVERSIÓN
# ----------------------------------------------------------------
KG_CM2_A_MPA = 1.0 / 10.197   # 1 MPa = 10.197 kg/cm²
MPA_A_KG_CM2 = 10.197

# ----------------------------------------------------------------
#  		    PROPIEDADES DE LOS MATERIALES
# ----------------------------------------------------------------
eco = 0.002

def Ec_concreto(fco):
    
    fco_MPa = fco * KG_CM2_A_MPA
    Ec_MPa  = 5000 * np.sqrt(fco_MPa)
    return Ec_MPa * MPA_A_KG_CM2   # kg/cm²


rho_cc = 0.02   # Cuantía longitudinal típica (Mander 1988)

# ----------------------------------------------------------------
#  			PARÁMETROS DINÁMICOS 
# ----------------------------------------------------------------

def parametros_dinamicos(fco_kgcm2, eco, Ec_kgcm2, strain_rate):
    fco_MPa = fco_kgcm2 * KG_CM2_A_MPA

    e_ref = 0.00001   # Tasa de referencia cuasi-estática 
    base_f  = 0.035 * (fco_MPa ** 2)
    Df = ((1 + (strain_rate  / base_f) ** (1 / 6)) / 
	(1 + (e_ref / base_f) ** (1 / 6)))

    base_E  = 0.035 * (fco_MPa ** 3)
    DE = ((1 + (strain_rate  / base_E) ** (1 / 6)) /
        (1 + (e_ref / base_E) ** (1 / 6)))

    De = (1 / (3 * Df)) * (1 + np.sqrt(1.0 + (3.0 * Df ** 2) / DE))

    fco_dyn = Df * fco_kgcm2
    eco_dyn = De * eco
    Ec_dyn  = DE * Ec_kgcm2

    return fco_dyn, eco_dyn, Ec_dyn, Df, DE, De


# ----------------------------------------------------------------
#  		   COEFICIENTE DE EFICACIA (ke) 
# ----------------------------------------------------------------

def ke_estribos_espirales(s, ds):
    return (1 - (s / (2 * ds))) / (1 - rho_cc)

def ke_estribos_circulares(s, ds):
    return ((1 - (s / (2 * ds))) ** 2) / (1 - rho_cc)

# ----------------------------------------------------------------
#          PRESIÓN Y RESISTENCIA DE CONFINAMIENTO NOMINAL 
# ----------------------------------------------------------------

def presion_confinamiento_nominal(rho_s, fyh):
    return 0.5 * rho_s * fyh


def fcc_mander(fco, fl):
    return fco * (-1.254 + 2.254*np.sqrt(1 + 7.94 * (fl / fco)) - 2*(fl / fco))

def ecc_mander(eco, fcc, fco):
    return eco * (1 + 5 * (fcc / fco - 1))


# ----------------------------------------------------------------
#  		   CURVA DE MANDER (CONFINADA) 
# ----------------------------------------------------------------

def curva_mander(eps, fcc, ecc, Ec):
    Esec = fcc / ecc
    r = Ec / (Ec - Esec)
    x = eps / ecc
    return (fcc * x * r) / (r - 1 + x ** r)

# ----------------------------------------------------------------
# 			CURVA NO CONFINADA 
# ----------------------------------------------------------------

def curva_no_confinada(eps, fco, eco, Ec):
    
    fc   = np.zeros_like(eps)
    ecu1 = 2 * eco       # 0.004
    ecu  = 0.005         # deformación de falla total

    Esec = fco / eco
    r = Ec / (Ec - Esec)

    for i in range(len(eps)):
        if eps[i] <= ecu1:
            x     = eps[i] / eco
            fc[i] = (fco * x * r) / (r - 1.0 + x ** r)
        elif eps[i] <= ecu:
            x1    = ecu1 / eco
            f1    = (fco * x1 * r) / (r - 1.0 + x1 ** r)
            fc[i] = f1 * (ecu - eps[i]) / (ecu - ecu1)
        else:
            fc[i] = 0

    return fc
# ----------------------------------------------------------------
#			  Selección de ke
# ----------------------------------------------------------------

def _calcular_fl_ke(tipo_col, rho_s, fyh, s, ds, b, h, c, Asx, Asy, wi, s_prima):
 
    if tipo_col == "Circular con espiral":
        ke = ke_estribos_espirales(s, ds)
        fl = ke * presion_confinamiento_nominal(rho_s, fyh)

    elif tipo_col == "Circular con estribos circulares":
        ke = ke_estribos_circulares(s, ds)
        fl = ke * presion_confinamiento_nominal(rho_s, fyh)

    elif tipo_col == "Rectangular con estribos":
        bc = b - 2 * c
        dc = h - 2 * c

        if bc <= 0 or dc <= 0:
            raise ValueError("Dimensiones del núcleo inválidas (bc o dc ≤ 0)")

        if s_prima is None:
            s_prima = s
        if wi is None:
            wi = np.array([bc / 3, bc / 3, bc / 3])

        sum_wi = np.sum(wi ** 2)
        ke = ((1 - sum_wi / (6 *bc*dc)) * (1 - s_prima / (2*bc)) *(1 - s_prima / (2*dc)) / (1 - rho_cc))

        if Asx is None or Asy is None:
            raise ValueError("Faltan Asx o Asy para sección rectangular")

        flx = (Asx * fyh) / (s * dc)
        fly = (Asy * fyh) / (s * bc)
        fl  = np.sqrt((ke * flx) * (ke * fly))

    else:
        raise ValueError(f"Tipo no válido: '{tipo}'")

    return fl, ke

# ----------------------------------------------------------------
#  			 FUNCIÓN PRINCIPAL
# =----------------------------------------------------------------

def funcion_mander(fco, fyh, rho_s, s, ds, esm, tipo, b=None, h=None, c=None,
                   Asx=None, Asy=None, wi=None, s_prima=None, strain_rate=0.00001):
    
    Ec = Ec_concreto(fco)
    fco_dyn, eco_dyn, Ec_dyn, Df, DE, De = parametros_dinamicos(fco, eco, Ec, strain_rate)
    fl, ke = _calcular_fl_ke(tipo, rho_s, fyh, s, ds, b, h, c, Asx, Asy, wi, s_prima)


    # CONFINADO
    fcc_est = fcc_mander(fco, fl)
    ecc_est = ecc_mander(eco, fcc_est, fco)
    ecu_est = 0.004 + (1.4 * rho_s * fyh * esm) / fcc_est

    # CONFINADO DINAMICO (Tasa de deformación)
    fcc_din = fcc_mander(fco_dyn, fl)
    ecc_din = ecc_mander(eco_dyn, fcc_din, fco_dyn)
    ecu_din = 0.004 + (1.4 * rho_s * fyh * esm) / fcc_din

# ----------------------------------------------------------------
#  			 DUCTILIDAD
# =----------------------------------------------------------------
    mu_est = ecu_est / ecc_est
    mu_din = ecu_din / ecc_din

    eps        = np.linspace(0.0, 0.02, 500)
    fc_no_conf = curva_no_confinada(eps, fco, eco, Ec)

    fc_conf_est = curva_mander(eps, fcc_est, ecc_est, Ec)
    fc_conf_est[eps > ecu_est] = 0.0

    fc_conf_din = curva_mander(eps, fcc_din, ecc_din, Ec_dyn)
    fc_conf_din[eps > ecu_din] = 0.0

    info_est = {"fco":  fco,  "eco":  eco,  "Ec":  Ec, "fcc":  fcc_est, "ecc": ecc_est, "ecu": ecu_est,
        "fl":   fl,   "ke":   ke, "mu": mu_est,}

    info_din = {"Df": Df, "DE": DE, "De": De,
        "fco_est": fco,     "eco_est": eco,     "Ec_est": Ec,
        "fco_dyn": fco_dyn, "eco_dyn": eco_dyn, "Ec_dyn": Ec_dyn,
        "fcc": fcc_din, "ecc": ecc_din, "ecu": ecu_din,
        "fl": fl, "ke": ke,
        "strain_rate": strain_rate,"mu": mu_din,}

    return eps, fc_no_conf, fc_conf_est, fc_conf_din, info_est, info_din
