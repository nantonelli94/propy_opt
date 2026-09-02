import streamlit as st
import numpy as np
from propy import WageningenBPropeller

st.set_page_config(page_title="Cátedra Propulsión - Optimización", layout="wide", page_icon="⚓")

st.title("⚓ Sistema de Optimización Propulsiva - Método de la Cátedra")
st.markdown("Automatización estricta del **Anexo 1: Procedimiento analítico de diseño y optimización polinómica**.")

# =============================================================================
# PANEL LATERAL: ENTRADAS DEL USUARIO (MÉTODO DE LA CÁTEDRA)
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
T_est_por_helice = T_total / num_ejes  # Tiro estimado por hélice (T*)

# 2. Área mínima de Keller
P_atm = 101000.0  # Conforme a tu apunte (101000 N/m²)
P_0 = P_atm + (rho * 9.81 * h_eje)
P_v = 1700.0
k_keller = 0.2 if num_ejes == 1 else 0.1
Fa_F_min = ((1.3 + 0.3 * Z) * T_est_por_helice) / ((P_0 - P_v) * (D**2)) + k_keller

# Forzamos un área comercial válida redondeando hacia arriba del límite de Keller
AE_A0 = max(0.35, min(1.0, float(np.ceil(Fa_F_min * 20) / 20)))

# 3. Cálculo de la constante del método C*
C_star = T_est_por_helice / (rho * (D**2) * (V_A**2))

# =============================================================================
# INTERFAZ DE RESULTADOS E INTERSECCIONES
# =============================================================================
st.subheader("📝 Evaluación Hidrodinámica e Intersecciones")
col_inf1, col_inf2, col_inf3 = st.columns(3)
col_inf1.metric("Tiro Necesario por Hélice (T*)", f"{T_est_por_helice/1000:.2f} kN")
col_inf2.metric("Área Mínima de Keller (Fa/F)", f"{Fa_F_min:.3f}", f"Hélice AE/A0 = {AE_A0:.2f}")
col_inf3.metric("Constante C* del Método", f"{C_star:.4f}")

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
# 4. OPTIMIZACIÓN ANALÍTICA: LÍNEA DE TENDENCIA POLINÓMICA
# =============================================================================
if len(pares_validos_pd) >= 3:
    coefs = np.polyfit(pares_validos_pd, lista_eta_optimos, 2)
    
    # Derivamos: 2*a*(P/D) + b = 0
    pd_optimo_analitico = -coefs[1] / (2 * coefs[0])
    
    if pd_optimo_analitico < 0.5 or pd_optimo_analitico > 1.4:
        pd_optimo_analitico = pares_validos_pd[np.argmax(lista_eta_optimos)]
        
    j_optimo_analitico = np.interp(pd_optimo_analitico, pares_validos_pd, lista_J_optimos)
    
    prop_optima = WageningenBPropeller(diameter=D, blades=Z, pd_ratio=pd_optimo_analitico, area_ratio=AE_A0)
    eta_max = prop_optima.eta(j_optimo_analitico)
    kt_opt = prop_optima.kt(j_optimo_analitico)
    kq_opt = prop_optima.kq(j_optimo_analitico)
    
    n_hélice_rps = V_A / (D * j_optimo_analitico)
    n_hélice_rpm = n_hélice_rps * 60.0
    
    # Potencia entregada (DHP)
    Q_nm = kq_opt * rho * (n_hélice_rps**2) * (D**5)
    P_entregada_kW = (2 * np.pi * n_hélice_rps * Q_nm) / 1000.0
    DHP = P_entregada_kW / 0.735499  
    
    # RESOLUCIÓN DEL ERROR DE VARIABLES: Usamos la variable explícita de eficiencia
    BHP_minimo = DHP / eficiencia_ejes
    BHP_diseño = BHP_minimo * 1.10  
    
    st.subheader("📈 Resultado de la Optimización Polinómica")
    
    col_g1, col_g2, col_g3, col_g4 = st.columns(4)
    col_g1.metric("P/D Óptimo Analítico", f"{pd_optimo_analitico:.3f}")
    col_g2.metric("Coeficiente J* Óptimo", f"{j_optimo_analitico:.3f}")
    col_g3.metric("Rendimiento Máximo (η)", f"{eta_max:.2%}")
    col_g4.metric("Giro Requerido Hélice", f"{n_hélice_rpm:.1f} RPM")
    
    st.info(f"**Requisito de Potencia del Motor Calculado:** El motor seleccionado debe entregar al menos **{BHP_diseño:.1f} HP**.")
    
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
            
    if motor_valido:
        i_teorica = motor_valido["RPM"] / n_hélice_rpm
        cajas_comerciales = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
        i_comercial = cajas_comerciales[np.argmin([abs(c - i_teorica) for c in cajas_comerciales])]
        
        st.success("⚙️ Acoplamiento de Maquinaria Sugerido:")
        col_m1, col_m2, col_m3 = st.columns(3)
        col_m1.metric("Motor Comercial Elegido", motor_valido["modelo"], f"{motor_valido['HP']} HP @ {motor_valido['RPM']} RPM")
        col_m2.metric("Relación de Caja Calculada (i)", f"{i_teorica:.2f}:1")
        col_m3.metric("Caja Comercial más Próxima", f"{i_comercial:.1f}:1")
        
        n_real_bita_rps = (motor_valido["RPM"] / i_comercial) / 60.0
        KT0 = prop_optima.kt(0.0)
        T0 = KT0 * rho * (n_real_bita_rps**2) * (D**4)
        T_bita_total = T0 * (1 - t_0) * num_ejes
        
        st.metric("💪 Fuerza Final Verificada en el Tiro de la Bita", f"{T_bita_total/1000:.1f} kN (~{T_bita_total/9810:.1f} Toneladas)")
    else:
        st.warning("⚠️ La potencia requerida supera los límites del catálogo comercial básico.")
else:
    st.error("La constante C* arroja valores fuera de los diagramas operativos de la Serie B. Intenta cambiar el diámetro o la velocidad de diseño.")
