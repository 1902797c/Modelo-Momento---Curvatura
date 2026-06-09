
# ====================================================================
#  APP — MOMENTO CURVATURA  (M – φ)
#  Marina Fierros Marcelo — MIAE
# ====================================================================

import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
 
from MODELO_MANDER import Ec_concreto, eco as ECO_REF
from MOMENTO_CURVATURA import (calcular_momento_curvatura, puntos_clave)
st.set_page_config(page_title="Momento – Curvatura",layout="wide",initial_sidebar_state="expanded",)
 
st.title("Curva Momento – Curvatura   |   M – φ")
st.caption("Marina Fierros Marcelo · MIAE")
 
st.markdown("""<style>
h1 { font-size: 30px !important; }
[data-testid="stMetricValue"] { font-size: 18px; }
[data-testid="stMetricLabel"] { font-size: 13px; }
</style>
""", unsafe_allow_html=True)
# ─────────────────────────────────────────────────────────────────
#  MATERIAL
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Parámetros del Material")
 
fco  = st.sidebar.number_input("f'co — Resistencia del concreto (kg/cm²)",
                                value=250.0, min_value=100.0, step=10.0)
fyh  = st.sidebar.number_input("fyh — Fluencia acero transversal (kg/cm²)",
                                value=4200.0, step=100.0)
esm  = st.sidebar.number_input("εsm — Deformación máx. acero transversal",
                                value=0.12, format="%.3f")
# ─────────────────────────────────────────────────────────────────
# Acero longitudinal
# ───────────────────────────────────────────────────────────────── 

st.sidebar.header("Acero longitudinal")
 
AREAS_VARILLA = {
    "#3": 0.71, "#4": 1.27, "#5": 1.99,
    "#6": 2.87, "#8": 5.07, "#10": 7.94, "#12": 11.40,}
 
fy_long  = st.sidebar.number_input("fy longitudinal (kg/cm²)", value=4200.0, step=100.0)
Es_long  = st.sidebar.number_input("Es longitudinal (kg/cm²)", value=2_000_000.0, step=100_000.0)
var_long = st.sidebar.selectbox("Varilla longitudinal", list(AREAS_VARILLA.keys()), index=3)
n_barras = st.sidebar.number_input("Número de barras longitudinales", min_value=4, value=8, step=1)
As_barra = AREAS_VARILLA[var_long]

# ─────────────────────────────────────────────────────────────────
#  TASA DE DEFORMACIÓN
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Tasa de Deformación")
strain_rate = st.sidebar.number_input(
    "Tasa de deformación ε̇ (s⁻¹)", value=0.05000, format="%.5f",
    help="0.000001 = cuasi-estático · 0.01–0.10 = sísmico típico")

# ────────────────────────────────────────────────────────────────────
# TIPO DE SECCIÓN Y CONFINAMIENTO
# ────────────────────────────────────────────────────────────────────
st.sidebar.header("Geometría y confinamiento")
 
tipo_col = st.sidebar.selectbox(
    "Tipo de columna",
    ["Circular con espiral",
     "Circular con estribos circulares",
     "Rectangular con estribos"])
 
AREAS_VARILLA_T = AREAS_VARILLA   
 
# Valores paea evitar error
b = h_sec = D = None
c = rho_s = s = ds = None
Asx = Asy = wi = s_prima = None


# ────────────────────────────────────────────────────────────────────
#                                 CIRCULAR
# ────────────────────────────────────────────────────────────────────
if tipo_col in ["Circular con espiral", "Circular con estribos circulares"]:
    D  = st.sidebar.number_input("Diámetro exterior D (cm)", value=50.0, step=5.0)
    c  = st.sidebar.number_input("Recubrimiento c (cm)", value=5.0, step=0.5)
    s  = st.sidebar.number_input("Espaciamiento s (cm)", value=8.0, step=1.0)
 
    var_t    = st.sidebar.selectbox("Varilla transversal", list(AREAS_VARILLA_T.keys()), index=2)
    Ash_barra = AREAS_VARILLA_T[var_t]
 
    if tipo_col == "Circular con espiral":
        Ash = Ash_barra
    else:
        ramas = st.sidebar.number_input("Número de ramas", value=2, min_value=1, step=1)
        Ash   = ramas * Ash_unit
 
    db_t = np.sqrt(4 * Ash_barra / np.pi)
    ds   = D - 2 * c - db_t
    rho_s = (4.0 * Ash) / (ds * s)
    
    # =========================
    # CUANTÍA LONGITUSINAL 
    # =========================
    Ag   = np.pi * D ** 2 / 4.0
    tipo_seccion = "circular"
    tipo_mander  = tipo_col

    st.sidebar.caption(f"ds = {ds:.2f} cm  |  ρs = {rho_s*100:.2f}%")
    
else:    # Rectangular
    b     = st.sidebar.number_input("Ancho b (cm)",    value=30.0, step=5.0)
    h_sec = st.sidebar.number_input("Peralte h (cm)",  value=65.0, step=5.0)
    c     = st.sidebar.number_input("Recubrimiento c (cm)", value=5.0, step=0.5)
    s     = st.sidebar.number_input("Espaciamiento s (cm)", value=8.0, step=1.0)
 
    var_t    = st.sidebar.selectbox("Varilla transversal (estribo)", list(AREAS_VARILLA_T.keys()), index=2)
    area_t   = AREAS_VARILLA_T[var_t]
 
    ramas_x = st.sidebar.number_input("Ramas en X", min_value=2, value=2, step=1)
    ramas_y = st.sidebar.number_input("Ramas en Y", min_value=2, value=2, step=1)
 
    Asx = ramas_x * area_t
    Asy = ramas_y * area_t
 
    s_prima = st.sidebar.number_input("Espaciamiento libre s' (cm)", value=6.0)
 
    bc    = b - 2.0 * c
    dc    = h_sec - 2.0 * c
    ds    = min(bc, dc)
    n_x   = max(int(ramas_x) - 1, 1)
    n_y   = max(int(ramas_y) - 1, 1)
    wi_x  = np.full(n_x, bc / n_x)
    wi_y  = np.full(n_y, dc / n_y)
    wi    = np.concatenate([wi_x, wi_y])
    rho_s = (2.0 * (Asx + Asy)) / (s * (bc + dc))
 
    Ag   = b * h_sec
    tipo_seccion = "rectangular"
    tipo_mander  = "Rectangular con estribos"
    
    st.sidebar.caption(f"ρs = {rho_s*100:.2f}%")
 
As_total = n_barras * As_barra
rho_g    = As_total / Ag
 
# ────────────────────────────────────────────────────────────────────
#  SIDEBAR — CARGA AXIAL (DOS NIVELES)
# ────────────────────────────────────────────────────────────────────
st.sidebar.header("Carga axial")
 
Pu1_frac = st.sidebar.slider("Nivel 1  —  Pu / (f'c · Ag)",
                              min_value=0.0, max_value=0.60,
                              value=0.10, step=0.05,
                              help="0 = flexión pura")
Pu2_frac = st.sidebar.slider("Nivel 2  —  Pu / (f'c · Ag)",
                              min_value=0.0, max_value=0.60,
                              value=0.30, step=0.05)
 
Pu1 = Pu1_frac * fco * Ag   # kg
Pu2 = Pu2_frac * fco * Ag   # kg
 
st.sidebar.write(f"Pu₁ = {Pu1:,.0f} kg  |  Pu₂ = {Pu2:,.0f} kg")
 
# ────────────────────────────────────────────────────────────────────
#  BOTÓN DE CÁLCULO
# ────────────────────────────────────────────────────────────────────
calcular = st.sidebar.button("▶  Calcular M–φ", use_container_width=True)
 
# ────────────────────────────────────────────────────────────────────
#  RESULTADOS
# ────────────────────────────────────────────────────────────────────
if calcular:
 
    Ec  = Ec_concreto(fco)
    eco = ECO_REF
 
    kwargs_comunes = dict(
        tipo_seccion      = tipo_seccion,
        b                 = b,
        h_sec             = h_sec,
        D                 = D,
        c                 = c,
        fco               = fco,
        fy                = fy_long,
        Es                = Es_long,
        num_barras        = int(n_barras),
        As_barra          = As_barra,
        fyh               = fyh,
        rho_s             = rho_s,
        s                 = s,
        ds                = ds,
        esm               = esm,
        tipo_confinamiento= tipo_mander,
        Asx               = Asx,
        Asy               = Asy,
        wi                = wi,
        s_prima           = s_prima,
        n_fibras          = 200,
        n_puntos          = 80)
 
    with st.spinner("Calculando curvas…"):
        try:
            phi1, M1, info1 = calcular_momento_curvatura(Pu=Pu1, **kwargs_comunes)
            phi2, M2, info2 = calcular_momento_curvatura(Pu=Pu2, **kwargs_comunes)
        except Exception as e:
            st.error(f"Error en el cálculo: {e}")
            st.stop()
 
    # ── Puntos clave ─────────────────────────────────────────────────
    pk1 = puntos_clave(phi1, M1, fco, Ec, fy_long, Es_long)
    pk2 = puntos_clave(phi2, M2, fco, Ec, fy_long, Es_long)
 
    # ── MÉTRICAS ─────────────────────────────────────────────────────
    st.markdown("---")
    st.markdown("### Propiedades del concreto confinado — Mander (1988)")
 
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("f'co (kg/cm²)",  f"{fco:.1f}")
    c2.metric("f'cc (kg/cm²)",  f"{info1['fcc']:.2f}",
              delta=f"+{info1['fcc']-fco:.2f}")
    c3.metric("εcc",            f"{info1['ecc']:.4f}",
              delta=f"+{info1['ecc']-eco:.4f}")
    c4.metric("εcu",            f"{info1['ecu']:.4f}")
    c5.metric("ke",             f"{info1['ke']:.3f}")
 
    st.markdown("### Sección transversal")
    sc1, sc2, sc3, sc4 = st.columns(4)
    sc1.metric("Ag (cm²)",             f"{Ag:.2f}")
    sc2.metric("As total (cm²)",        f"{As_total:.2f}")
    sc3.metric("ρg (%)",                f"{rho_g*100:.2f}")
    sc4.metric("ρs (%)",                f"{rho_s*100:.2f}")
 
    # ── DUCTILIDAD ───────────────────────────────────────────────────
    if pk1 and pk2:
        st.markdown("### Ductilidad de curvatura  μ = φu / φy")
        d1, d2 = st.columns(2)
 
        def _mu(pk):
            if "y" in pk and "u" in pk and pk["y"][0] > 0:
                return pk["u"][0] / pk["y"][0]
            return float("nan")
 
        with d1:
            st.markdown(f"**Nivel 1 — Pu = {Pu1_frac:.0%} f'c·Ag**")
            da, db_col, dc = st.columns(3)
            da.metric("φy (rad/cm)", f"{pk1['y'][0]:.2e}")
            db_col.metric("φu (rad/cm)", f"{pk1['u'][0]:.2e}")
            dc.metric("μφ", f"{_mu(pk1):.2f}")
 
        with d2:
            st.markdown(f"**Nivel 2 — Pu = {Pu2_frac:.0%} f'c·Ag**")
            da2, db2, dc2 = st.columns(3)
            da2.metric("φy (rad/cm)", f"{pk2['y'][0]:.2e}")
            db2.metric("φu (rad/cm)", f"{pk2['u'][0]:.2e}")
            dc2.metric("μφ", f"{_mu(pk2):.2f}")
 
    # ── GRÁFICA ──────────────────────────────────────────────────────
    st.markdown("### Curvas M – φ")
 
    fig, ax = plt.subplots(figsize=(11, 5.5))
 
    color1, color2 = "black", "steelblue"
 
    # Nivel 1
    ax.plot(phi1 * 1e4, M1 / 1e5,
            color=color1, lw=2.2,
            label=f"Pu₁ = {Pu1_frac:.0%} f'c·Ag  —  Pu₁ = {Pu1:,.0f} kg")
 
    # Nivel 2
    ax.plot(phi2 * 1e4, M2 / 1e5,
            color=color2, lw=2.2, ls="--",
            label=f"Pu₂ = {Pu2_frac:.0%} f'c·Ag  —  Pu₂ = {Pu2:,.0f} kg")
 
    # Puntos clave
    marcadores = {"agr": ("o", "Agrietamiento"),
                  "y":   ("s", "Fluencia"),
                  "u":   ("^", "Último")}
 
    for llave, (marca, nombre) in marcadores.items():
        for pk, col in [(pk1, color1), (pk2, color2)]:
            if llave in pk:
                phi_pt, M_pt = pk[llave]
                ax.plot(phi_pt * 1e4, M_pt / 1e5,
                        marker=marca, color=col,
                        ms=8, zorder=5)
 
    # Leyenda adicional de marcadores
    from matplotlib.lines import Line2D
    leyenda_extra = [
        Line2D([0], [0], marker="o", color="gray", ls="none", ms=8, label="Agrietamiento"),
        Line2D([0], [0], marker="s", color="gray", ls="none", ms=8, label="Fluencia"),
        Line2D([0], [0], marker="^", color="gray", ls="none", ms=8, label="Último"),]
 
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + leyenda_extra,
              loc="upper center", bbox_to_anchor=(0.5, -0.14),
              ncol=3, fontsize=8.5, frameon=False)
 
    titulo = (f"Curva Momento – Curvatura · {tipo_col}\n"
              f"f'co = {fco:.0f} kg/cm²  |  fy = {fy_long:.0f} kg/cm²  |  "
              f"f'cc = {info1['fcc']:.1f} kg/cm²")
 
    ax.set_xlabel("Curvatura φ  (×10⁻⁴ rad/cm)", fontsize=10)
    ax.set_ylabel("Momento M  (×10⁵ kg·cm)", fontsize=10)
    ax.set_title(titulo, fontsize=9.5)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    plt.tight_layout()
 
    st.pyplot(fig)
 
    # ── TABLA DE RESULTADOS ──────────────────────────────────────────
    with st.expander("Ver tabla de resultados"):
        col_t1, col_t2 = st.columns(2)
 
        df1 = pd.DataFrame({
            "φ (rad/cm)":    phi1,
            "M (kg·cm)":     M1,
            "M (t·m)":       M1 / 100_000,
        })
        df2 = pd.DataFrame({
            "φ (rad/cm)":    phi2,
            "M (kg·cm)":     M2,
            "M (t·m)":       M2 / 100_000,
        })
 
        with col_t1:
            st.markdown(f"**Nivel 1 — Pu = {Pu1_frac:.0%}**")
            st.dataframe(df1.style.format("{:.6f}"),
                         use_container_width=True, height=250)
 
        with col_t2:
            st.markdown(f"**Nivel 2 — Pu = {Pu2_frac:.0%}**")
            st.dataframe(df2.style.format("{:.6f}"),
                         use_container_width=True, height=250)
 
    # ── DESCARGA ─────────────────────────────────────────────────────
    df_export = pd.DataFrame({
        "phi_nivel1_rad_cm": phi1,
        "M_nivel1_kgcm":     M1,
        "phi_nivel2_rad_cm": phi2,
        "M_nivel2_kgcm":     M2,
    })
 
    st.download_button(
        "📥  Descargar resultados (CSV)",
        df_export.to_csv(index=False),
        file_name="curva_M_phi.csv",
    )
 
else:
    st.info("Configura los parámetros en el panel izquierdo y presiona **▶ Calcular M–φ**.")
 


