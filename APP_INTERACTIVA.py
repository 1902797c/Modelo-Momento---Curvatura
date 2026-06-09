import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
from MODELO_MANDER import funcion_mander

st.set_page_config(page_title="Modelo de Mander", layout="wide")

# ─────────────────────────────────────────────────────────────────
st.title("Modelo de Mander — Curva Esfuerzo-Deformación")
st.caption("Marina Fierros Marcelo")
st.markdown("""
<style>
/* Título principal */
h1 {
    font-size: 35px !important;
}
/* Número grande */
[data-testid="stMetricValue"] {font-size: 20px;}

/* Texto arriba */
[data-testid="stMetricLabel"] {font-size: 14px;}

/* Delta */
[data-testid="stMetricDelta"] {font-size: 13px;}

</style>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────
#  MATERIAL
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Parámetros del Material")
fco = st.sidebar.number_input("f'co (kg/cm²)", value=280.0)
fyh = st.sidebar.number_input("fyh (kg/cm²)", value=4200.0)
esm = st.sidebar.number_input("εsm (deformación máxima del acero)", value=0.12,
                               help="Ductilidad máxima del acero de confinamiento")

# ─────────────────────────────────────────────────────────────────
#  TASA DE DEFORMACIÓN
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Tasa de Deformación")
strain_rate = st.sidebar.number_input(
    "Tasa de deformación ε̇ (s⁻¹)", value=0.05000, format="%.5f",
    help="0.000001 = cuasi-estático · 0.01–0.10 = sísmico típico")

# ─────────────────────────────────────────────────────────────────
#  GEOMETRÍA Y CONFINAMIENTO
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Geometría y Confinamiento")

AREAS_VARILLA = {
    "#3": 0.71, "#4": 1.27, "#5": 1.99, "#6": 2.87,
    "#8": 5.07, "#10": 7.94, "#12": 11.40}

tipo = st.sidebar.selectbox(
    "Tipo de columna",
    ["Selecciona una opción",
     "Circular con espiral",
     "Circular con estribos circulares",
     "Rectangular con estribos"]
)
# Inicialización para evitar NameError
rho_s = s = ds = None
b = h = c = Asx = Asy = wi = s_prima = None

RHO_MIN = 0.01
RHO_MAX = 0.025

# ── Secciones circulares ──────────────────────────────────────────
if tipo in ["Circular con espiral", "Circular con estribos circulares"]:

    s  = st.sidebar.number_input("Espaciamiento s (cm)", value=8.0)
    D = st.sidebar.number_input("Diámetro exterior D (cm)", value=50.0)
    c = st.sidebar.number_input("Recubrimiento c (cm)", value=5.0)

    varilla   = st.sidebar.selectbox("Varilla", list(AREAS_VARILLA.keys()))
    Ash_barra = AREAS_VARILLA[varilla]

    if tipo == "Circular con espiral":
        Ash = Ash_barra  
    else:
        ramas = st.sidebar.number_input("Número de ramas", value=2, min_value=1, step=1)
        Ash   = ramas * Ash_barra

# Diámetro de la varilla transversal (a partir del área)
    db = np.sqrt(4 * Ash_barra / np.pi)

# Diámetro del núcleo confinado (Mander)
    ds = D - 2*c - db

# Mostrar al usuario
    st.sidebar.caption(f"Diámetro núcleo ds = {ds:.2f} cm")
    # =========================
    # CUANTÍA ρs (CIRCULAR)
    # =========================
    rho_s = (4.0 * Ash) / (ds * s)

    st.sidebar.write(f"Ash = {Ash:.2f} cm²")
    st.sidebar.write("ρs = ", f"{rho_s*100:.2f}%")

    # =========================
    # VALIDACIÓN INTELIGENTE
    # =========================
    if rho_s < RHO_MIN:
        st.sidebar.error("⚠️ ρs < 1% (BAJO)")

        if tipo == "Circular con espiral":
            st.sidebar.markdown("🔧 **Recomendación:**")
            st.sidebar.write("- Reducir espaciamiento *s*")
            st.sidebar.write("- Usar varilla de mayor diámetro")

        else:
            st.sidebar.markdown("🔧 **Recomendación:**")
            st.sidebar.write("- Aumentar número de ramas")
            st.sidebar.write("- Usar varilla de mayor diámetro")
            st.sidebar.write("- Reducir espaciamiento *s*")

    elif rho_s > RHO_MAX:
        st.sidebar.warning("⚠️ ρs > 2.5% (ALTO)")

        if tipo == "Circular con espiral":
            st.sidebar.markdown("🔧 **Recomendación:**")
            st.sidebar.write("- Aumentar espaciamiento *s*")
            st.sidebar.write("- Usar varilla de menor diámetro")

        else:
            st.sidebar.markdown("🔧 **Recomendación:**")
            st.sidebar.write("- Reducir número de ramas")
            st.sidebar.write("- Usar varilla de menor diámetro")
            st.sidebar.write("- Aumentar espaciamiento *s*")

    else:
        st.sidebar.success("✅ ρs dentro del rango (1%–2.5%)")

# ── Sección rectangular ───────────────────────────────────────────
if tipo == "Rectangular con estribos":

    b = st.sidebar.number_input("Ancho b (cm)", value=40.0)
    h = st.sidebar.number_input("Peralte h (cm)", value=60.0)
    c = st.sidebar.number_input("Recubrimiento c (cm)", value=5.0)
    s = st.sidebar.number_input("Espaciamiento s (cm)", value=8.0)

    varilla    = st.sidebar.selectbox("Varilla (estribo)", list(AREAS_VARILLA.keys()))
    area_barra = AREAS_VARILLA[varilla]

    ramas_x = st.sidebar.number_input("Número de ramas en X", min_value=1.0, max_value=10.0, value=2.0, step=0.1, format="%.2f")
    ramas_y = st.sidebar.number_input("Número de ramas en Y", min_value=1.0, max_value=10.0, value=2.0, step=0.1, format="%.2f")

    Asx = ramas_x * area_barra
    Asy = ramas_y * area_barra
    st.sidebar.write(f"Asx = {Asx:.2f} cm²  |  Asy = {Asy:.2f} cm²")

    s_prima = st.sidebar.number_input("Espaciamiento libre s' (cm)", value=6.0)

    bc = b - 2.0 * c
    dc = h - 2.0 * c
    ds = min(bc, dc)

    # wi automáticos: Espaciado uniforme entre ramas
    # X → ramas_x ramas crean (ramas_x - 1) espacios de bc/(ramas_x-1)
    # Y → ramas_y ramas crean (ramas_y - 1) espacios de dc/(ramas_y-1)
    n_esp_x = max(int(ramas_x) - 1, 1)
    n_esp_y = max(int(ramas_y) - 1, 1)
    wi_x    = np.full(n_esp_x, bc / n_esp_x)
    wi_y    = np.full(n_esp_y, dc / n_esp_y)
    wi      = np.concatenate([wi_x, wi_y])

    # Cuantía automática: ρs = 2(Asx + Asy) / (s·(bc + dc))
    rho_s = (2.0 * (Asx + Asy)) / (s * (bc + dc))
    st.sidebar.write("ρs = ", f"{rho_s*100:.2f}%")
    st.sidebar.caption(
        f"wi_x: {n_esp_x} espacios de {bc/n_esp_x:.2f} cm  |  "
        f"wi_y: {n_esp_y} espacios de {dc/n_esp_y:.2f} cm")

    if rho_s < RHO_MIN:
        st.sidebar.error("⚠️ ρs < 1% (BAJO)")
        st.sidebar.markdown("🔧 **Recomendación:**")
        st.sidebar.write("- Aumentar número de ramas en X o Y")
        st.sidebar.write("- Usar varilla de mayor diámetro")
        st.sidebar.write("- Reducir espaciamiento s")

    elif rho_s > RHO_MAX:
        st.sidebar.warning("⚠️ ρs > 2.5% (ALTO)")
        st.sidebar.markdown("🔧 **Recomendación:**")
        st.sidebar.write("- Reducir número de ramas")
        st.sidebar.write("- Usar varilla de menor diámetro")
        st.sidebar.write("- Aumentar espaciamiento s")

    else:
        st.sidebar.success("✅ ρs dentro del rango (1%–2.5%)")
# ─────────────────────────────────────────────────────────────────
#  				BOTÓN 
# ─────────────────────────────────────────────────────────────────
if st.sidebar.button("Generar curva"):

    if rho_s is None or s is None:
        st.error("Completa los parámetros de confinamiento antes de generar la curva.")
        st.stop()

    try:
        eps, fc_no_conf, fc_conf_est, fc_conf_din, info_est, info_din = funcion_mander(
            fco, fyh, rho_s, s, ds, esm, tipo,
            b=b, h=h, c=c, Asx=Asx, Asy=Asy, wi=wi, s_prima=s_prima,
            strain_rate=strain_rate)
    except Exception as e:
        st.error(f"Error en el modelo: {e}")
        st.stop()

    # ────────────────── Métricas ──────────────────────
    st.markdown("#### Propiedades del concreto confinado")
    col_e, col_d = st.columns(2)

    with col_e:
        st.markdown("**Estático**")
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("f'cc (kg/cm²)", f"{info_est['fcc']:.2f}")
        m2.metric("εcc",            f"{info_est['ecc']:.3f}")
        m3.metric("εcu",            f"{info_est['ecu']:.3f}")
        m4.metric("μ",              f"{info_est['mu']:.2f}")

    with col_d:
        st.markdown(f"**Dinámico** — ε̇ = {strain_rate:.2e} s⁻¹")
        m4, m5, m6, m7 = st.columns(4)
        m4.metric("f'cc (kg/cm²)", f"{info_din['fcc']:.2f}",
                  delta=f"{info_din['fcc'] - info_est['fcc']:+.2f}")
        m5.metric("εcc",            f"{info_din['ecc']:.3f}",
                  delta=f"{info_din['ecc'] - info_est['ecc']:+.3f}")
        m6.metric("εcu",            f"{info_din['ecu']:.3f}",
                  delta=f"{info_din['ecu'] - info_est['ecu']:+.3f}")
        m7.metric("μ", f"{info_din['mu']:.2f}",
          delta=f"{info_din['mu'] - info_est['mu']:+.2f}")

    # ── Factores DIF ─────────────────────────────────────────────
    if strain_rate != 0.00001:

        st.info(
            f"f'co: {info_din['fco_est']:.2f} → **{info_din['fco_dyn']:.2f}** kg/cm²  |  "
            f"εco: {info_din['eco_est']:.3f} → **{info_din['eco_dyn']:.3f}**  |  "
            f"Ec: {info_din['Ec_est']:.0f} → **{info_din['Ec_dyn']:.0f}** kg/cm²")

    # ── Gráfica ────────────────────────────────────────
    mask_nc = eps <= 0.005

    fig, ax = plt.subplots(figsize=(10, 5.5))

    # No confinado
    ax.plot(eps[mask_nc], fc_no_conf[mask_nc],
            color="#555555", lw=2.0, ls="--", label="No confinado")

    # Confinado estático
    ax.plot(eps, fc_conf_est,
            color="Black", lw=2.2,
            label=f"Confinado estático  f'cc = {info_est['fcc']:.1f} kg/cm²")

    # Confinado dinámico
    ax.plot(eps, fc_conf_din,
            color="darkblue", lw=2.5,
            label=f"Confinado dinámico  ε̇ = {strain_rate:.1e} s⁻¹  f'cc = {info_din['fcc']:.1f} kg/cm²")

    # Áreas
    ax.fill_between(eps, fc_no_conf, fc_conf_est,
                    where=(fc_conf_est > fc_no_conf),
                    color="black", alpha=0.08, label="Aporte confinamiento")
    ax.fill_between(eps, fc_conf_est, fc_conf_din,
                    where=(fc_conf_din > fc_conf_est),
                    color="blue", alpha=0.10, label="Incremento por DIF")

    # Líneas de referencia

    

    titulo = f"Curva esfuerzo-deformación — Modelo de Mander (1988) — {tipo}"
    ax.set_xlabel("Deformación unitaria", fontsize=10)
    ax.set_ylabel("Esfuerzo (kg/cm²)", fontsize=10)
    ax.set_title(titulo, fontsize=10)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13),
              ncol=2, fontsize=8.5, frameon=False)
    ax.grid(alpha=0.35)
    plt.tight_layout()

    st.pyplot(fig)
  
    df = pd.DataFrame({
        "eps": eps,
        "No_confinado": fc_no_conf,
        "Confinado_est": fc_conf_est,
        "Confinado_din": fc_conf_din})
    
    st.dataframe(df, use_container_width=True, height=215)
    st.download_button(
        "📥 Descargar resultados (CSV)",
        df.to_csv(index=False),
        file_name="curva_mander.csv")
