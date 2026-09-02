import streamlit as st
import numpy as np
from propy import WageningenBPropeller

st.set_page_config(page_title="Cátedra Propulsión - Optimización", layout="wide", page_icon="⚓")

# Estilos CSS personalizados para agrandar los números y mejorar la estética de la tabla
st.markdown("""
    <style>
    .big-table-header {
        font-size: 22px !important;
        font-weight: bold !important;
        color: #1A365D;
        border-bottom: 2px solid #1A365D;
        padding-bottom: 5px;
        margin-top: 20px;
    }
    .metric-row {
        background-color: #F7FAFC;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 8px;
        border-left: 5px solid #2B6CB0;
        display: flex;
        justify-content: space-between;
        align-items: center;
    }
    .metric-label {
        font-size: 18px !important;
        font-weight: 500;
        color: #2D3748;
    }
    .metric-value {
        font-size: 26px !important;
        font-weight: bold;
        color: #2B6CB0;
        font-family: 'Courier New', Courier, monospace;
    }
    </style>
""", unsafe_allow_html=True)

st.title("⚓ Sistema de Optimización Propulsiva - Método de la Cátedra")
st.markdown("Automatización estricta del **Anexo 1: Procedimiento analítico de diseño y optimización polinómica**.")

# =============================================================================
# PANEL LATERAL: ENTRADAS DEL BUQUE Y ENTORNO
# =============================================================================
st.sidebar.header("📊 1. Datos del Buque y Ensayo")
V_nudos = st.sidebar.slider("Velocidad del buque (V en Nudos)", 5.0, 30.0, 12.0, 0.5)
R_T = st.sidebar.number_input("Resistencia al avance (R_T en Newtons)", 1000, 1000000, 35000, 1000)

st.sidebar.subheader("🔗 Coeficientes de Interacción")
w = st.sidebar.slider("Fracción de estela (w)", 0.05, 0.40, 0.20, 0.01)
t = st.sidebar.slider("Deducción de empuje (t)", 0.05, 0.35, 0.15, 0.01)
t_0 = st.sidebar.slider("Deducción de empuje en bita (t_0)", min_value=0.01, max_value=0.10, value=0.04, step=0.01)
num_ejes = st.sidebar.radio("Líneas de eje", options=[1, 2], index=0)

st.sidebar.subheader("🧼 Verificación de Cavitación (Keller)")
h_eje = st.sidebar.slider("Inmersión del eje (h en metros)", 0.5, 10.0, 2.0, 0.1)
Z = st.sidebar.selectbox("Número de palas (z)", options=[3, 4, 5, 6], index=1)

st.sidebar.subheader("📐 Restricción Física del Codaste")
D = st.sidebar.number_input("Diámetro disponible de la hélice (D en metros)", 0.5, 8.0, 1.40, 0.05)
rho = st.sidebar.number_input("Densidad del agua (kg/m³)", value=1025.0)
eficiencia_ejes = st.sidebar.slider("Eficiencia de línea de ejes (η_s)", min_value=0.90, max_value=1.0, value=0.97, step=0.01)

# =============================================================================
# DESARROLLO MATEMÁTICO DEL ANEXO 1
# =============================================================================
V_A = (V_nudos * 0.514444) * (1 - w)
T_total = R_T / (1 - t)
T_est_por_helice = T_total / num_ejes  

P_atm = 101000.0  
P_0 = P_atm + (rho * 9.81 * h_eje)
P_v = 1700.0
k_keller = 0.2 if num_ejes == 1 else 0.1
Fa_F_min = ((1.3 + 0.3 * Z) * T_est_por_helice) / ((P_0 - P_v) * (D**2)) + k_keller

AE_A0 = max(0.35, min(1.0, float(np.ceil(Fa_F_min * 20) / 20)))
C_star = T_est_por_helice / (rho * (D**2) * (V_A**2))

# Barrido analítico para intersectar KT = C* * J^2
pd_rango = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])
lista_J_optimos = []
lista_eta_optimos = []
pares_validos_pd = []

for pd_test in pd_rango:
    prop_eval = WageningenBPropeller(diameter=D, blades=Z, pd_ratio=pd_test, area_ratio=AE_A0)
    j_interseccion = None
    min_dif = float('inf')
    
    for j_candidato in np.linspace(0.01, 1.3, 500):
        kt_hélice = prop_eval.kt(j_candidato)
        kt_parabola = C_star * (j_candidato**2)
        dif = abs(kt_hélice - kt_parabola)
        if dif < min_dif:
            min_dif = dif
            j_interseccion = j_candidato
            
    if min_dif < 0.01:
        eta = prop_eval.eta(j_interseccion)
        lista_J_optimos.append(j_interseccion)
        lista_eta_optimos.append(eta)
        pares_validos_pd.append(pd_test)

# =============================================================================
# BLOQUE GENERAL DE CONTROL Y RENDERIZADO
# =============================================================================
col_main1, col_main2 = st.columns([3, 2])

with col_main1:
    st.subheader("PASO 1: Optimización de la Hélice (Aguas Abiertas)")
    ejecutar_opt = st.button("🚀 Calcular Hélice Óptima", type="primary")

if ejecutar_opt or 'opt_realizada' in st.session_state:
    st.session_state['opt_realizada'] = True
    
    if len(pares_validos_pd) >= 3:
        coefs = np.polyfit(pares_validos_pd, lista_eta_optimos, 2)
        pd_optimo_analitico = -coefs[1] / (2 * coefs[0])
        
        if pd_optimo_analitico < 0.5 or pd_optimo_analitico > 1.4:
            pd_optimo_analitico = pares_validos_pd[np.argmax(lista_eta_optimos)]
            
        j_optimo_analitico = np.interp(pd_optimo_analitico, pares_validos_pd, lista_J_optimos)
        
        prop_optima = WageningenBPropeller(diameter=D, blades=Z, pd_ratio=pd_optimo_analitico, area_ratio=AE_A0)
        eta_max = prop_optima.eta(j_optimo_analitico)
        kt_opt = prop_optima.kt(j_optimo_analitico)
        kq_opt = prop_optima.kq(j_optimo_analitico)
        
        n_hélice_rps = V_A / (D * j_optimo_analitico)
        Q_nm = kq_opt * rho * (n_hélice_rps**2) * (D**5)
        
        P_entregada_kW = (2 * np.pi * n_hélice_rps * Q_nm) / 1000.0
        DHP = P_entregada_kW / 0.735499  # cv
        BHP_minimo = DHP / eficiencia_ejes
        
        # Guardamos en variables de estado para que persistan al usar el Paso 2
        st.session_state['n_hélice_rps'] = n_hélice_rps
        st.session_state['prop_optima'] = prop_optima
        
        # RENDERIZADO VISUAL EXCLUSIVO: TABLA CON "BIG NUMBERS"
        with col_main1:
            st.markdown('<div class="big-table-header">Tabla de Resultados Finales</div>', unsafe_allow_html=True)
            
            rows = [
                ("η_MÁX (Eficiencia máxima en aguas abiertas)", f"{eta_max:.4f}"),
                ("P/D (ηMÁX) (Relación paso/diámetro óptima)", f"{pd_optimo_analitico:.4f}"),
                ("J (ηMÁX) (Coeficiente de avance óptimo)", f"{j_optimo_analitico:.4f}"),
                ("nh(mínimas) [rps] (Velocidad de rotación hélice)", f"{n_hélice_rps:.2f} rps"),
                ("Kq (ηMÁX) (Coeficiente de torque óptimo)", f"{kq_opt:.5f}"),
                ("Kt (ηMÁX) (Coeficiente de empuje óptimo)", f"{kt_opt:.4f}"),
                ("Q(mínima) [N.m] (Torque absorbido mínimo)", f"{Q_nm:.1f} N.m"),
                ("DHP mínimo [cv] (Potencia entregada a hélice)", f"{DHP:.1f} cv"),
                ("BHP mínimo [cv] (Potencia al freno mínima)", f"{BHP_minimo:.1f} cv")
            ]
            
            for label, val in rows:
                st.markdown(f"""
                <div class="metric-row">
                    <span class="metric-label">{label}</span>
                    <span class="metric-value">{val}</span>
                </div>
                """, unsafe_allow_html=True)

        # =============================================================================
        # PASO 2: INTERFAZ DINÁMICA DE SELECCIÓN DE MOTOR POR PARTE DEL USUARIO
        # =============================================================================
        with col_main2:
            st.subheader("PASO 2: Selección y Entrada de Maquinaria")
            st.info("Utiliza los datos hidrodinámicos obtenidos a la izquierda para elegir un motor de catálogo comercial y cargar sus datos reales abajo.")
            
            # Entradas de texto numérico para el motor elegido por el alumno
            motor_hp = st.number_input("Potencia nominal del motor seleccionado (HP)", min_value=10.0, max_value=20000.0, value=float(np.round(BHP_minimo * 1.10, -1)))
            motor_rpm = st.number_input("Revoluciones de diseño del motor seleccionado (RPM)", min_value=100.0, max_value=6000.0, value=1800.0, step=100.0)
            
            # Verificaciones cinemáticas y relación de reducción
            n_helice_rpm = st.session_state['n_hélice_rps'] * 60.0
            r_max_calculo = motor_rpm / n_helice_rpm
            
            st.markdown("---")
            st.markdown(f"### ⚙️ Relación de Reducción Máxima de Cálculo:")
            st.markdown(f"<div style='font-size: 38px; font-weight: bold; color: #E53E3E; text-align: center; background-color: #FFF5F5; padding: 15px; border-radius: 8px; border: 2px dashed #E53E3E;'>r máx = {r_max_calculo:.2f} : 1</div>", unsafe_allow_html=True)
            st.markdown(f"*Nota de cátedra: Debes buscar una caja reductora comercial cuyo valor de reducción sea el más próximo e inferior a este límite teórico para evitar subcargar o sobrecargar la hélice.*")
            
            # Entrada de la reducción comercial definitiva elegida por el alumno
            i_comercial = st.number_input("Caja comercial final adoptada por el alumno (i)", min_value=0.5, max_value=15.0, value=float(np.floor(r_max_calculo * 2) / 2), step=0.1)
            
            # Verificación del Tiro a punto fijo basado en la combinación cargada
            n_real_bita_rps = (motor_rpm / i_comercial) / 60.0
            p_opt = st.session_state['prop_optima']
            KT0 = p_opt.kt(0.0)
            T0 = KT0 * rho * (n_real_bita_rps**2) * (D**4)
            T_bita_total = T0 * (1 - t_0) * num_ejes
            
            st.markdown("---")
            st.markdown(f"### 💪 Tiro de la Bita Verificado (J=0):")
            st.markdown(f"<div style='font-size: 32px; font-weight: bold; color: #2F855A; text-align: center; background-color: #F0FFF4; padding: 10px; border-radius: 8px;'>{T_bita_total/9810:.1f} Toneladas ({T_bita_total/1000:.1f} kN)</div>", unsafe_allow_html=True)
    else:
        st.error("La constante C* arroja valores fuera de los diagramas operativos de la Serie B. Modifica los parámetros del codaste en el panel izquierdo.")
