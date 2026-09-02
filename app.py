import streamlit as st
import numpy as np
from propy import WageningenBPropeller

# Configuración de la interfaz web
st.set_page_config(page_title="Optimizador Propulsivo Naval", layout="wide", page_icon="⚓")

st.title("⚓ Diseñador de Hélices y Match de Propulsión")
st.subheader("Cálculo automático de Motor, Caja, Hélice Óptima y Tiro sobre la Bita")
st.markdown("Cátedra de Propulsión. Simulación basada en los polinomios oficiales de la **Serie B de Wageningen**.")

# =============================================================================
# PANEL LATERAL: ENTRADAS DEL USUARIO (Parámetros de diseño del casco)
# =============================================================================
st.sidebar.header("📊 Parámetros del Casco y Operación")
V_nudos = st.sidebar.slider("Velocidad de diseño (Nudos)", min_value=5.0, max_value=25.0, value=12.0, step=0.5)
R_T = st.sidebar.number_input("Resistencia al avance (Newtons)", min_value=1000, max_value=500000, value=35000, step=1000)

st.sidebar.subheader("🔗 Coeficientes de Interacción Casco-Hélice")
w = st.sidebar.slider("Fracción de estela (w)", min_value=0.05, max_value=0.40, value=0.20, step=0.01)
t = st.sidebar.slider("Deducción de empuje (t)", min_value=0.05, max_value=0.35, value=0.15, step=0.01)
t_0 = st.sidebar.slider("Deducción de empuje en bita (t_0)", min_value=0.01, max_value=0.10, value=0.04, step=0.01)

st.sidebar.subheader("📐 Geometría y Restricciones de Espacio")
D_max = st.sidebar.slider("Diámetro máximo permitido (m)", min_value=0.5, max_value=4.0, value=1.40, step=0.05)
Z = st.sidebar.selectbox("Número de palas (Z)", options=[3, 4, 5, 6, 7], index=1)
rho = st.sidebar.number_input("Densidad del agua (kg/m³)", value=1025.0)
eta_s = st.sidebar.slider("Eficiencia de línea de ejes (η_s)", min_value=0.90, max_value=1.0, value=0.97, step=0.01)

# =============================================================================
# CATÁLOGOS COMERCIALES DE MOTORES Y CAJAS REDUCTORAS
# =============================================================================
V = V_nudos * 0.514444          # Conversión de nudos a m/s
V_A = V * (1 - w)               # Velocidad de avance
T_req = R_T / (1 - t)           # Empuje requerido total

CATALOGO_MOTORES = [
    {"modelo": "Motor Comercial A (Ligero)", "HP": 300.0, "RPM": 2100.0},
    {"modelo": "Motor Comercial B (Medio)",  "HP": 400.0, "RPM": 1900.0},
    {"modelo": "Motor Comercial C (Pesado)", "HP": 500.0, "RPM": 1800.0},
    {"modelo": "Motor Comercial D (MaxPot)",  "HP": 600.0, "RPM": 1800.0},
    {"modelo": "Motor Comercial E (Mega)",    "HP": 800.0, "RPM": 1600.0}
]
CATALOGO_CAJAS = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

# =============================================================================
# ALGORITMO PRINCIPAL DE OPTIMIZACIÓN Y CÁLCULO
# =============================================================================
if st.button("🚀 Ejecutar Optimización Global", type="primary"):
    resultados = []
    
    # Triple bucle: Evaluamos cada Motor, cada Caja y cada Diámetro físicamente viable
    for motor in CATALOGO_MOTORES:
        P_b_max = motor["HP"] * 745.7
        RPM_motor = motor["RPM"]
        n_motor = RPM_motor / 60.0
        
        for i in CATALOGO_CAJAS:
            n_prop = n_motor / i  # Frecuencia de la hélice (Hz)
            
            for D in np.arange(0.5, D_max + 0.05, 0.05):
                J = V_A / (n_prop * D)
                if J <= 0 or J > 1.5: 
                    continue
                    
                # Empuje adimensional requerido (KT_req)
                KT_req = T_req / (rho * (n_prop**2) * (D**4))
                
                best_pd = None
                min_error = float('inf')
                
                # Buscador interno de Paso (P/D) óptimo usando propy
                for pd_test in np.linspace(0.5, 1.4, 80):
                    prop_test = WageningenBPropeller(diameter=D, blades=Z, pd_ratio=pd_test)
                    try:
                        error = abs(prop_test.kt(J) - KT_req)
                        if error < min_error:
                            min_error = error
                            best_pd = pd_test
                    except:
                        continue
                
                # Si el paso cumple con el empuje hidrodinámico requerido del casco
                if min_error < 0.015:
                    prop_optima = WageningenBPropeller(diameter=D, blades=Z, pd_ratio=best_pd)
                    KQ = prop_optima.kq(J)
                    
                    # Potencia demandada por la hélice en navegación libre (kW)
                    P_d = (2 * np.pi * rho * (n_prop**3) * (D**5) * KQ) / 1000.0
                    P_disponible_eje = (P_b_max * eta_s) / 1000.0
                    
                    # Condición de aceptación: El motor debe soportar la carga en viaje
                    if P_d <= P_disponible_eje:
                        eta_O = prop_optima.eta(J)
                        
                        # --- Verificación de Tiro sobre la Bita (Bollard Pull) a J = 0 ---
                        KT0 = prop_optima.kt(0.0)
                        KQ0 = prop_optima.kq(0.0)
                        
                        # El torque frena al motor a velocidad cero. Calculamos las RPM reales sostenidas
                        Torque_max_eje = (P_b_max / (2 * np.pi * n_motor)) * eta_s
                        n0_prop = np.sqrt(Torque_max_eje / (2 * np.pi * rho * (D**5) * KQ0))
                        if n0_prop > n_prop:
                            n0_prop = n_prop
                            
                        # Fuerza resultante en Newtons transferida al cable de amarre
                        T0 = KT0 * rho * (n0_prop**2) * (D**4)
                        T_bita_kN = (T0 * (1 - t_0)) / 1000.0
                        
                        resultados.append({
                            'Motor': motor['modelo'],
                            'HP Motor': motor['HP'],
                            'Relación Caja': f"{i}:1",
                            'Diámetro (m)': round(D, 2),
                            'Paso P/D': round(best_pd, 2),
                            'Eficiencia 𝜂_O': round(eta_O * 100, 1),
                            'Pot. Demandada (kW)': round(P_d, 1),
                            'Tiro Bita (kN)': round(T_bita_kN, 1),
                            'Tiro Bita (Ton)': round(T_bita_kN / 9.81, 1),
                            '_raw_eta': eta_O
                        })

    # PRESENTACIÓN DE RESULTADOS EN LA INTERFAZ GRÁFICA
    if resultados:
        ganador = max(resultados, key=lambda x: x['_raw_eta'])
        
        st.success("🎯 ¡Combinación óptima calculada con éxito para la cátedra!")
        
        # Tarjetas de resumen destacadas (KPIs)
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Motor Recomendado", f"{ganador['HP Motor']} HP", ganador['Motor'])
        col2.metric("Reducción Caja Ideal", ganador['Relación Caja'])
        col3.metric("Hélice (Diám x Paso)", f"{ganador['Diámetro (m)']}m x {ganador['Paso P/D']}")
        col4.metric("Fuerza en la Bita", f"{ganador['Tiro Bita (Ton)']} Ton", f"{ganador['Tiro Bita (kN)']} kN")
        
        # Tabla interactiva con el total de alternativas
        st.subheader("📈 Cuadro de Alternativas Viables Encontradas")
        st.markdown("Los alumnos pueden ordenar y comparar las opciones del catálogo haciendo clic en las columnas:")
        
        for r in resultados: 
            r.pop('_raw_eta', None)
            
        st.dataframe(resultados, use_container_width=True)
        
    else:
        st.error("❌ No se encontraron combinaciones viables en el catálogo. Sugerencia académica: Aumenta el Diámetro Máximo o disminuye la Resistencia al avance para que las hélices de la Serie B entren en rango operativo.")
