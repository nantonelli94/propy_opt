import streamlit as st
import numpy as np
from propy import WageningenBPropeller

st.set_page_config(page_title="Cátedra Propulsión - Optimización", layout="wide", page_icon="⚓")

st.title("⚓ Sistema de Optimización Propulsiva - Método de la Cátedra")
st.markdown("Automatización estricta del **Anexo 1: Procedimiento analítico de diseño y optimización polinómica**.")

# =============================================================================
# PANEL LATERAL: ENTRADAS DEL USUARIO
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
Z = st.sidebar.selectbox("Número de palas (z)", options=[3, 4, 5, 6, 7], index=1)

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

# Barrido analítico de las curvas para intersectar KT = C* * J^2
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
# PROCESAMIENTO E INTERFAZ DE RESULTADOS
# =============================================================================
if st.button("🚀 Ejecutar Optimización Global", type="primary"):
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
        DHP = P_entregada_kW / 0.735499  # kW a CV (caballos de vapor)
        BHP_minimo = DHP / eficiencia_ejes
        BHP_diseño = BHP_minimo * 1.10  
        
        # Base de datos de motores comerciales para determinar la r_máx de cálculo
        CATALOGO_MOTORES = [
            {"modelo": "Motor Naval Heavy-Duty A", "HP": 350.0, "RPM": 1200.0},
            {"modelo": "Motor Naval Medium-Duty B", "HP": 500.0, "RPM": 1800.0},
            {"modelo": "Motor Naval High-Speed C", "HP": 750.0, "RPM": 2100.0},
            {"modelo": "Motor Naval High-Speed D", "HP": 1200.0, "RPM": 2300.0},
            {"modelo": "Motor Naval Potencia Max E", "HP": 2500.0, "RPM": 1800.0}
        ]
        
        motor_valido = None
        for mot in CATALOGO_MOTORES:
            if mot["HP"] >= BHP_diseño:
                motor_valido = mot
                break
        
        # Relación de reducción máxima teórica de cálculo
        if motor_valido:
            r_max_calculo = motor_valido["RPM"] / (n_hélice_rps * 60.0)
            nombre_motor = motor_valido["modelo"]
        else:
            r_max_calculo = 1800.0 / (n_hélice_rps * 60.0) # RPM estándar por defecto si excede catálogo
            nombre_motor = "Motor fuera de catálogo estándar (>2500 HP)"

        # =============================================================================
        # CONSTRUCCIÓN DE LA TABLA DE RESULTADOS FINALES SOLICITADA
        # =============================================================================
        st.subheader("📋 Tabla de resultados finales")
        
        tabla_datos = {
            "Parámetro / Variable": [
                "η_MÁX (Eficiencia máxima en aguas abiertas)",
                "P/D (ηMÁX) (Relación paso/diámetro óptima)",
                "J (ηMÁX) (Coeficiente de avance óptimo)",
                "nh(mínimas) [rps] (Velocidad de rotación de la hélice)",
                "Kq (ηMÁX) (Coeficiente de torque óptimo)",
                "Kt (ηMÁX) (Coeficiente de empuje óptimo)",
                "Q(mínima) [N.m] (Torque absorbido mínimo)",
                "DHP mínimo [cv] (Potencia entregada a la hélice)",
                "BHP mínimo [cv] (Potencia al freno mínima del motor)",
                "r máx (relación de reducción máxima de cálculo)"
            ],
            "Valor Obtenido": [
                f"{eta_max:.4f}",
                f"{pd_optimo_analitico:.4f}",
                f"{j_optimo_analitico:.4f}",
                f"{n_hélice_rps:.2f} rps",
                f"{kq_opt:.5f}",
                f"{kt_opt:.4f}",
                f"{Q_nm:.1f} N.m",
                f"{DHP:.1f} cv",
                f"{BHP_minimo:.1f} cv",
                f"{r_max_calculo:.2f}:1"
            ]
        }
        
        # Mostrar tabla interactiva de ancho completo en Streamlit
        st.table(tabla_datos)
        
        # Bloque informativo complementario de maquinaria
        st.success(f"⚙️ **Acoplamiento Hidrodinámico:** Para cubrir este requerimiento se sugiere el **{nombre_motor}** con un BHP de diseño corregido (margen 10%) de **{BHP_diseño:.1f} cv**.")
        
    else:
        st.error("La constante C* arroja valores fuera de los diagramas operativos de la Serie B. Intenta cambiar el diámetro o la resistencia del casco.")
