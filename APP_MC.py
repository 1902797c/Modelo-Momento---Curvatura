# ====================================================================
#  APP — MOMENTO CURVATURA  (M – φ)
#  Marina Fierros Marcelo — MIAE
# ====================================================================

import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd

from MODELO_MANDER import Ec_concreto, eco as ECO_REF
from MOMENTO_CURVATURA import calcular_momento_curvatura, puntos_clave,coordenadas

st.set_page_config(page_title="Momento – Curvatura", layout="wide", initial_sidebar_state="expanded")

st.title("Curva Momento – Curvatura   |   M – φ")
st.caption("Marina Fierros Marcelo · MIAE")

st.markdown("""<style>
h1 { font-size: 30px !important; }
[data-testid="stMetricValue"] { font-size: 18px; }
[data-testid="stMetricLabel"] { font-size: 13px; }
</style>""", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────
#  MATERIAL
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Parámetros del Material")
fco  = st.sidebar.number_input("f'co — Resistencia del concreto (kg/cm²)", value=250.0, min_value=100.0, step=10.0)
fyh  = st.sidebar.number_input("fyh — Fluencia acero (kg/cm²)", value=4200.0, step=100.0)
esm  = st.sidebar.number_input("εsm — Deformación máx. acero transversal", value=0.12, format="%.3f")

# ─────────────────────────────────────────────────────────────────
#  ACERO LONGITUDINAL
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Acero longitudinal")

AREAS_VARILLA = {
    "#3": 0.71, "#4": 1.27, "#5": 1.99,
    "#6": 2.87, "#8": 5.07, "#10": 7.94, "#12": 11.40,
}

Es_long  = st.sidebar.number_input("Es longitudinal (kg/cm²)", value=2_000_000.0, step=100_000.0)
var_long = st.sidebar.selectbox("Varilla longitudinal", list(AREAS_VARILLA.keys()), index=3)
n_vars   = st.sidebar.number_input("Número de varillas longitudinales", min_value=4, value=8, step=1)
As_var   = AREAS_VARILLA[var_long]

# ─────────────────────────────────────────────────────────────────
#  GEOMETRÍA Y CONFINAMIENTO
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Geometría y confinamiento")

tipo_col = st.sidebar.selectbox(
    "Tipo de columna",
    ["Circular con espiral",
     "Circular con estribos circulares",
     "Rectangular con estribos"])

AREAS_VARILLA_T = AREAS_VARILLA

# Inicialización para evitar errores
b = h_sec = D = None
c = rho_s = s = ds = None
Asx = Asy = wi = s_prima = None

RHO_MIN = 0.01
RHO_MAX = 0.025
# ─────────────────────────────────────────────────────────────────
#  CIRCULAR
# ─────────────────────────────────────────────────────────────────
if tipo_col in ["Circular con espiral", "Circular con estribos circulares"]:
    D  = st.sidebar.number_input("Diámetro exterior D (cm)", value=50.0, step=5.0)
    c  = st.sidebar.number_input("Recubrimiento c (cm)", value=5.0, step=0.5)
    s  = st.sidebar.number_input("Espaciamiento s (cm)", value=8.0, step=1.0)

    var_t   = st.sidebar.selectbox("Varilla transversal", list(AREAS_VARILLA_T.keys()), index=0)
    Ash_var = AREAS_VARILLA_T[var_t]

    if tipo_col == "Circular con espiral":
        Ash = Ash_var
    else:
        ramas = st.sidebar.number_input("Número de ramas", value=2, min_value=1, step=1)
        Ash   = ramas * Ash_var

    db_t  = np.sqrt(4 * Ash_var / np.pi)
    ds    = D - 2 * c - db_t
    rho_s = (4.0 * Ash) / (ds * s)
    Ag    = np.pi * D ** 2 / 4.0

    tipo_seccion = "circular"
    tipo_mander  = tipo_col
    st.sidebar.caption(f"ds = {ds:.2f} cm  |  ρs = {rho_s*100:.2f}%")

    if rho_s < RHO_MIN:
        st.sidebar.error("⚠️ ρs < 1% (BAJO)")

        if tipo_col == "Circular con espiral":
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

        if tipo_col == "Circular con espiral":
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
# ─────────────────────────────────────────────────────────────────
#  RECTANGULAR
# ─────────────────────────────────────────────────────────────────
else:
    b     = st.sidebar.number_input("Ancho b (cm)",   value=50.0, step=5.0)
    h_sec = st.sidebar.number_input("Peralte h (cm)", value=50.0, step=5.0)
    c     = st.sidebar.number_input("Recubrimiento c (cm)", value=5.0, step=0.5)
    s     = st.sidebar.number_input("Espaciamiento s (cm)", value=8.0, step=1.0)

    var_t  = st.sidebar.selectbox("Varilla transversal (estribo)", list(AREAS_VARILLA_T.keys()), index=2)
    area_t = AREAS_VARILLA_T[var_t]

    ramas_x = st.sidebar.number_input("Ramas en X", min_value=2, value=2, step=1)
    ramas_y = st.sidebar.number_input("Ramas en Y", min_value=2, value=2, step=1)

    Asx     = ramas_x * area_t
    Asy     = ramas_y * area_t
    s_prima = st.sidebar.number_input("Espaciamiento libre s' (cm)", value=6.0)

    bc   = b - 2.0 * c
    dc   = h_sec - 2.0 * c
    ds   = min(bc, dc)
    n_x  = max(int(ramas_x) - 1, 1)
    n_y  = max(int(ramas_y) - 1, 1)
    wi_x = np.full(n_x, bc / n_x)
    wi_y = np.full(n_y, dc / n_y)
    wi   = np.concatenate([wi_x, wi_y])
    rho_s = (2.0 * (Asx + Asy)) / (s * (bc + dc))
    Ag    = b * h_sec

    tipo_seccion = "rectangular"
    tipo_mander  = "Rectangular con estribos"
    st.sidebar.caption(f"ρs = {rho_s*100:.2f}%")

As_total = n_vars * As_var
rho_g    = As_total / Ag

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
#  CARGA AXIAL
# ─────────────────────────────────────────────────────────────────
st.sidebar.header("Carga axial")

Pu1_frac = st.sidebar.slider("Nivel 1  —  Pu / (f'c · Ag)",
                              min_value=0.0, max_value=0.50, value=0.05, step=0.05)
Pu2_frac = st.sidebar.slider("Nivel 2  —  Pu / (f'c · Ag)",
                              min_value=0.0, max_value=0.50, value=0.10, step=0.05)

Pu1 = Pu1_frac * fco * Ag
Pu2 = Pu2_frac * fco * Ag
st.sidebar.write(f"Pu₁ = {Pu1:,.0f} kg  |  Pu₂ = {Pu2:,.0f} kg")

# ─────────────────────────────────────────────────────────────────
#  BOTÓN
# ─────────────────────────────────────────────────────────────────
calcular = st.sidebar.button("▶  Calcular M–φ", use_container_width=True)

# ─────────────────────────────────────────────────────────────────
#  CÁLCULO Y RESULTADOS
# ─────────────────────────────────────────────────────────────────
if calcular:

    Ec  = Ec_concreto(fco)
    eco = ECO_REF

    # ── kwargs limpios — nombres exactos de calcular_momento_curvatura ──
    kwargs_comunes = dict(
        tipo_seccion       = tipo_seccion,
        b                  = b,
        h_sec              = h_sec,
        D                  = D,
        c                  = c,
        fco                = fco,
        fy                 = fyh,
        Es                 = Es_long,
        num_barras         = int(n_vars),   # nombre correcto del módulo
        As_barra           = As_var,        # nombre correcto del módulo
        fyh                = fyh,
        rho_s              = rho_s,
        s                  = s,
        ds                 = ds,
        esm                = esm,
        tipo_confinamiento = tipo_mander,
        Asx                = Asx,
        Asy                = Asy,
        wi                 = wi,
        s_prima            = s_prima,
        n_fibras           = 100,
        n_puntos           = 60,
    )

    with st.spinner("Calculando curvas…"):
        try:
            phi1, M1, info1 = calcular_momento_curvatura(Pu=Pu1, **kwargs_comunes)
            phi2, M2, info2 = calcular_momento_curvatura(Pu=Pu2, **kwargs_comunes)
        except Exception as e:
            st.error(f"Error en el cálculo: {e}")
            st.stop()

    pk1 = puntos_clave(phi1, M1, fco, Ec, fyh, Es_long)
    pk2 = puntos_clave(phi2, M2, fco, Ec, fyh, Es_long)


    # ── Gráfica ──────────────────────────────────────────────────
    st.markdown("### Curvas M – φ")
    fig, ax = plt.subplots(figsize=(11, 5.5))

    ax.plot(phi1 * 1e4, M1 / 1e5, color="black",     lw=2.2,
            label=f"Pu₁ = {Pu1_frac:.0%} f'c·Ag  —  Pu₁ = {Pu1:,.0f} kg")
    ax.plot(phi2 * 1e4, M2 / 1e5, color="steelblue", lw=2.2, ls="--",
            label=f"Pu₂ = {Pu2_frac:.0%} f'c·Ag  —  Pu₂ = {Pu2:,.0f} kg")

    marcadores = {"agr": ("o", "Agrietamiento"), "y": ("s", "Fluencia"), "u": ("^", "Último")}
    for llave, (marca, _) in marcadores.items():
        for pk, col in [(pk1, "black"), (pk2, "steelblue")]:
            if llave in pk:
                ax.plot(pk[llave][0] * 1e4, pk[llave][1] / 1e5,
                        marker=marca, color=col, ms=8, zorder=5)

    from matplotlib.lines import Line2D
    extra = [
        Line2D([0],[0], marker="o", color="gray", ls="none", ms=8, label="Agrietamiento"),
        Line2D([0],[0], marker="s", color="gray", ls="none", ms=8, label="Fluencia"),
        Line2D([0],[0], marker="^", color="gray", ls="none", ms=8, label="Último"),
    ]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + extra, loc="upper center",
              bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=8.5, frameon=False)

    ax.set_title(
        f"Curva Momento – Curvatura · {tipo_col}\n"
        f"f'co = {fco:.0f} kg/cm²  |  fy = {fyh:.0f} kg/cm²  |  f'cc = {info1['fcc']:.1f} kg/cm²",
        fontsize=9.5)
    ax.set_xlabel("Curvatura φ  (×10⁻⁴ rad/cm)", fontsize=10)
    ax.set_ylabel("Momento M  (×10⁵ kg·cm)", fontsize=10)
    ax.set_xlim(left=0); ax.set_ylim(bottom=0); ax.grid(alpha=0.3)
    plt.tight_layout()
    st.pyplot(fig)

   
    st.subheader("Esquema de la Sección Transversal y Refuerzo longitudinal")

    fig_sec, ax_sec = plt.subplots(figsize=(4, 5))

    if tipo_seccion == "rectangular":

        ax_sec.add_patch(plt.Rectangle((0, 0), b, h_sec, edgecolor='black', facecolor='#E0E0E0', lw=2, label='Recubrimiento'))
        # 2. Dibujar la trayectoria del estribo perimetral
        ax_sec.add_patch(plt.Rectangle((c, c), b - 2*c, h_sec - 2*c, edgecolor='grey', facecolor='white', lw=1, ls='--', label='Núcleo Confinado'))
        # 3. Obtener y graficar los puntos de acero
        x_s, y_s = coordenadas("rectangular", b, h_sec, c, n_vars)
        ax_sec.scatter(x_s, y_s, color='black', s=130, zorder=5, label='Varillas Long.')
        ax_sec.set_xlim(-5, b + 5)
        ax_sec.set_ylim(-5, h_sec + 5)
        ax_sec.set_xlabel("Ancho b (cm)")
        ax_sec.set_ylabel("Peralte h (cm)")
    else:
    # Caso circular
        circ_ext = plt.Circle((D/2.0, D/2.0), D/2.0, edgecolor='black', facecolor='#E0E0E0', lw=2, label='Recubrimiento')
        circ_int = plt.Circle((D/2.0, D/2.0), D/2.0 - c, edgecolor='grey', facecolor='white', lw=1, ls='--', label='Núcleo Confinado')
        ax_sec.add_patch(circ_ext)
        ax_sec.add_patch(circ_int)
        x_s, y_s = coordenadas("circular", D, D, c, n_vars)
        ax_sec.scatter(x_s, y_s, color='black', s=130, zorder=5, label='Varillas Long.')
        ax_sec.set_xlim(-5, D + 5)
        ax_sec.set_ylim(-5, D + 5)
        ax_sec.set_xlabel("X (cm)")
        ax_sec.set_ylabel("Y (cm)")

    ax_sec.set_aspect('equal')
    ax_sec.grid(True, linestyle=':', alpha=0.5)
    ax_sec.invert_yaxis()  
    ax_sec.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))

    st.pyplot(fig_sec)
    # ── Métricas ─────────────────────────────────────────────────
    st.markdown("---")
    st.markdown("### Propiedades del concreto confinado — Mander (1988)")
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("f'co (kg/cm²)", f"{fco:.1f}")
    c2.metric("f'cc (kg/cm²)", f"{info1['fcc']:.2f}", delta=f"+{info1['fcc']-fco:.2f}")
    c3.metric("εcc", f"{info1['ecc']:.4f}", delta=f"+{info1['ecc']-eco:.4f}")
    c4.metric("εcu", f"{info1['ecu']:.4f}")
    c5.metric("ke",  f"{info1['ke']:.3f}")

    st.markdown("### Sección transversal")
    sc1, sc2, sc3, sc4 = st.columns(4)
    sc1.metric("Ag (cm²)",      f"{Ag:.2f}")
    sc2.metric("As total (cm²)", f"{As_total:.2f}")
    sc3.metric("ρg (%)",         f"{rho_g*100:.2f}")
    sc4.metric("ρs (%)",         f"{rho_s*100:.2f}")

    # ── Ductilidad ───────────────────────────────────────────────
    def _mu(pk):
        if "y" in pk and "u" in pk and pk["y"][0] > 0:
            return pk["u"][0] / pk["y"][0]
        return float("nan")

    if pk1 and pk2:
        st.markdown("### Ductilidad de curvatura  μ = φu / φy")
        d1, d2 = st.columns(2)
        with d1:
            st.markdown(f"**Nivel 1 — Pu = {Pu1_frac:.0%} f'c·Ag**")
            a, b_, cc = st.columns(3)
            a.metric("φy (rad/cm)",  f"{pk1['y'][0]:.2e}")
            b_.metric("φu (rad/cm)", f"{pk1['u'][0]:.2e}")
            cc.metric("μφ",          f"{_mu(pk1):.2f}")
        with d2:
            st.markdown(f"**Nivel 2 — Pu = {Pu2_frac:.0%} f'c·Ag**")
            a2, b2, c2_ = st.columns(3)
            a2.metric("φy (rad/cm)",  f"{pk2['y'][0]:.2e}")
            b2.metric("φu (rad/cm)",  f"{pk2['u'][0]:.2e}")
            c2_.metric("μφ",          f"{_mu(pk2):.2f}")

    # ── Tabla ────────────────────────────────────────────────────
    with st.expander("Ver tabla de resultados"):
        t1, t2 = st.columns(2)
        with t1:
            st.markdown(f"**Nivel 1 — Pu = {Pu1_frac:.0%}**")
            st.dataframe(pd.DataFrame({
                "φ (rad/cm)": phi1, "M (kg·cm)": M1, "M (t·m)": M1/100_000
            }).style.format("{:.6f}"), use_container_width=True, height=250)
        with t2:
            st.markdown(f"**Nivel 2 — Pu = {Pu2_frac:.0%}**")
            st.dataframe(pd.DataFrame({
                "φ (rad/cm)": phi2, "M (kg·cm)": M2, "M (t·m)": M2/100_000
            }).style.format("{:.6f}"), use_container_width=True, height=250)

    # ── Descarga ─────────────────────────────────────────────────
    st.download_button(
        "📥  Descargar resultados (CSV)",
        pd.DataFrame({
            "phi_nivel1_rad_cm": phi1, "M_nivel1_kgcm": M1,
            "phi_nivel2_rad_cm": phi2, "M_nivel2_kgcm": M2,
        }).to_csv(index=False),
        file_name="curva_M_phi.csv",
    )

else:
    st.info("Configura los parámetros en el panel izquierdo y presiona **▶ Calcular M–φ**.")
