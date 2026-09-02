import streamlit as st
import numpy as np
from propy import WageningenBPropeller

st.set_page_config(page_title="Optimizador Propulsivo Naval", layout="wide", page_icon="⚓")
st.title("⚓ Diseñador de Hélices con Verificación de Cavitación (Keller)")
st.subheader("Cálculo automático de Motor, Caja, Hélice Óptima y Tiro sobre la Bita")

# PANEL LATERAL: ENTRADAS
st.sidebar.header("📊 Parámetros del Casco y Operación")
V_nudos = st.sidebar.slider("Velocidad de diseño (Nudos)", min_value=5.0, max_value=35.0, value=12.0, step=0.5)
R_T = st.sidebar.number_input("Resistencia al avance (Newtons)", min_value=1000, max_value=2000000, value=35000, step=5000)

st.sidebar.subheader("🔗 Coeficientes de Interacción")
w = st.sidebar.slider("Fracción de estela (w)", min_value=0.05, max_value=0.40, value=0.20, step=0.01)
t = st.sidebar.slider("Deducción de empuje (t)", min_value=0.05, max_value=0.35, value=0.15, step=0.01)
t_0 = st.sidebar.slider("Deducción de empuje en bita (t_0)", min_value=0.01, max_value=0.10, value=0.04, step=0.01)

st.sidebar.subheader("⚙️ Configuración de Propulsión")
num_ejes = st.sidebar.radio("Cantidad de líneas de eje (Hélices)", options=[1, 2], index=0)

st.sidebar.subheader("🧼 Parámetros de Cavitación (Keller)")
h_eje = st.sidebar.slider("Inmersión del eje de la hélice (h en metros)", min_value=0.5, max_value=15.0, value=2.0, step=0.1)
AE_A0_propuesto = st.sidebar.slider("Relación de área expandida propuesta (AE/A0 o Fa/F)", min_value=0.30, max_value=1.0, value=0.65, step=0.05)

st.sidebar.subheader("📐 Geometría y Restricciones")
D_max = st.sidebar.slider("Diámetro máximo permitido (m)", min_value=0.5, max_value=6.0, value=1.40, step=0.05)
Z = st.sidebar.selectbox("Número de palas (Z)", options=[3, 4, 5, 6, 7], index=1)
rho = st.sidebar.number_input("Densidad del agua (kg/m³)", value=1025.0)
eta_s = st.sidebar.slider("Eficiencia de línea de ejes", min_value=0.90, max_value=1.0, value=0.97, step=0.01)

# Cálculo de presiones para Keller
P_atm = 101325.0  
P_0 = P_atm + (rho * 9.81 * h_eje)  
P_v = 1700.0  
k_keller = 0.2 if num_ejes == 1 else 0.1

V = V_nudos * 0.514444          
V_A = V * (1 - w)
T_req_total = R_T / (1 - t)
T_req_por_helice = T_req_total / num_ejes  

# =============================================================================
# GENERACIÓN AUTOMÁTICA DEL CATÁLOGO MASIVO DE MOTORES (Hasta 5000 HP y de 800 a 3000 RPM)
# =============================================================================
CATALOGO_MOTORES = []

# Rango de potencias: Motores chicos/medianos de 100 en 100 HP, motores grandes de 500 en 500 HP
rangos_hp = list(range(200, 1000, 100)) + list(range(1000, 5500, 500))

for hp in rangos_hp:
    # Para cada potencia, barremos desde 800 hasta 3000 RPM con pasos de 200 RPM
    for rpm in range(800, 3200, 200):
        # Clasificación pedagógica del motor según sus vueltas
        if rpm <= 1200:
            tipo = "Heavy Duty"
        elif rpm <= 2200:
            tipo = "Medium Duty"
        else:
            tipo = "High Speed"
            
        CATALOGO_MOTORES.append({
            "modelo": f"Catálogo {hp}HP @ {rpm}RPM ({tipo})",
            "HP": float(hp),
            "RPM": float(rpm)
        })

CATALOGO_CAJAS = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]

# =============================================================================
# EJECUCIÓN DEL ALGORITMO
# =============================================================================
if st.button("🚀 Ejecutar Optimización Global", type="primary"):
    resultados = []
    
    # Contenedor para mostrar progreso ya que la base de datos ahora es muy grande
    with st.spinner("Escaneando base de datos masiva de motores, reducciones y hélices..."):
        for motor in CATALOGO_MOTORES:
            P_b_max = motor["HP"] * 745.7
            RPM_motor = motor["RPM"]
            n_motor = RPM_motor / 60.0
            
            for i in CATALOGO_CAJAS:
                n_prop = n_motor / i
                
                for D in np.arange(0.5, D_max + 0.05, 0.05):
                    # Filtro de Cavitación de Keller
                    keller_numerador = (1.3 + 0.3 * Z) * T_req_por_helice
                    keller_denominador = (P_0 - P_v) * (D**2)
                    Fa_F_min = (keller_numerador / keller_denominador) + k_keller
                    
                    if AE_A0_propuesto < Fa_F_min:
                        continue  
                    
                    J = V_A / (n_prop * D)
                    if J <= 0 or J > 1.5: continue
                    
                    KT_req = T_req_por_helice / (rho * (n_prop**2) * (D**4))
                    best_pd = None
                    min_error = float('inf')
                    
                    for pd_test in np.linspace(0.5, 1.4, 60): # Reducido a 60 pasos para mantener fluida la web
                        prop_test = WageningenBPropeller(diameter=D, blades=Z, pd_ratio=pd_test, area_ratio=AE_A0_propuesto)
                        try:
                            error = abs(prop_test.kt(J) - KT_req)
                            if error < min_error:
                                min_error = error
                                best_pd = pd_test
                        except:
                            continue
                    
                    if min_error < 0.015:
                        prop_optima = WageningenBPropeller(diameter=D, blades=Z, pd_ratio=best_pd, area_ratio=AE_A0_propuesto)
                        KQ = prop_optima.kq(J)
                        P_d = (2 * np.pi * rho * (n_prop**3) * (D**5) * KQ) / 1000.0
                        P_disponible_eje = (P_b_max * eta_s) / 1000.0
                        
                        # Filtro estricto: El motor debe cubrir la demanda pero no estar sobredimensionado en exceso (Margen de carga > 65%)
                        if P_d <= P_disponible_eje and P_d >= (P_disponible_eje * 0.65):
                            eta_O = prop_optima.eta(J)
                            KT0 = prop_optima.kt(0.0)
                            KQ0 = prop_optima.kq(0.0)
                            
                            Torque_max_eje = (P_b_max / (2 * np.pi * n_motor)) * eta_s
                            n0_prop = np.sqrt(Torque_max_eje / (2 * np.pi * rho * (D**5) * KQ0))
                            if n0_prop > n_prop: n0_prop = n_prop
                                
                            T0 = KT0 * rho * (n0_prop**2) * (D**4)
                            T_bita_helice = (T0 * (1 - t_0)) / 1000.0
                            T_bita_total_kN = T_bita_helice * num_ejes
                            
                            resultados.append({
                                'Motor Propuesto': motor['modelo'], 
                                'Potencia (HP)': motor['HP'], 
                                'RPM Motor': motor['RPM'],
                                'Relación Caja': f"{i}:1",
                                'Diámetro Hélice (m)': round(D, 2), 
                                'Paso P/D': round(best_pd, 2),
                                'Keller Mín (Fa/F)': round(Fa_F_min, 3),
                                'Eficiencia 𝜂_O': round(eta_O * 100, 1), 
                                'Pot. Demandada/Eje (kW)': round(P_d, 1),
                                'Tiro Total Bita (kN)': round(T_bita_total_kN, 1), 
                                'Tiro Total Bita (Ton)': round(T_bita_total_kN / 9.81, 1),
                                '_raw_eta': eta_O
                            })

    if resultados:
        ganador = max(resultados, key=lambda x: x['_raw_eta'])
        st.success(f"🎯 ¡Combinación óptima encontrada entre {len(CATALOGO_MOTORES)} opciones de motorización!")
        
        col1, col2, col3, col4, col5 = st.columns(5)
        col1.metric("Motor Seleccionado", f"{ganador['Potencia (HP)']} HP", ganador['Motor Propuesto'])
        col2.metric("Reducción Caja", ganador['Relación Caja'])
        col3.metric("Hélice (Diám x Paso)", f"{ganador['Diámetro Hélice (m)']}m x {ganador['Paso P/D']}")
        col4.metric("Keller Mínimo", ganador['Keller Mín (Fa/F)'])
        col5.metric("Tiro Total Bita", f"{ganador['Tiro Total Bita (Ton)']} Ton")
        
        st.subheader("📈 Cuadro de Alternativas Viables Encontradas")
        st.markdown("Usa los filtros de las cabeceras para ordenar por potencia, eficiencia o tiro en la bita:")
        for r in resultados: r.pop('_raw_eta', None)
        st.dataframe(resultados, use_container_width=True)
    else:
        st.error("❌ No se encontraron combinaciones viables. La resistencia del casco exige una potencia o un diámetro de hélice que supera tus límites actuales, o el área propuesta no supera el criterio de Keller.")
