import obspy
import numpy as np
import os
import tkinter as tk
from tkinter import filedialog
from obspy import UTCDateTime, Trace, Stream
from datetime import timedelta
from datetime import datetime
from colorama import Fore, Style
import inspect
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import json
import requests
import random
import time
import seaborn as sns
from fpdf import FPDF
from matplotlib.dates import DateFormatter
from matplotlib.ticker import FixedLocator, FixedFormatter

from flask import Blueprint, render_template, request

pagina2_bp = Blueprint('pagina2', __name__)

#*****************************************************************************
# Funcion transforma de Raw_data a gravedad
# Fecha: Jun 21 2023 
# Se utiliza dentro de la función info_proc
#*****************************************************************************
def Raw2grav(tr):
    # tr_p es una copia de tr para evitar modificar datos de la trama
    tr_p = tr.copy()
    # Factor para conversión Raw to gravedad
    factor = 520000 #(cuentas to gravedad)
    factor2 = 53007 #(cuentas to m/s^2)
    factor3 = 53.007 #(cuentas to mm/s^2)
    # Se transforma datos según factor2 m/s^2
    tr_p.data = np.divide(tr_p,factor2)
    return(tr_p)
#*****************************************************************************



#*****************************************************************************
# Función para obtener valores pico y pico_pico en una trama
# Fecha: Jul 14 2023 
#*****************************************************************************
def tr_param(tri):
    p_max = np.max(np.abs(tri.data)) 
    pk_pk = np.abs(np.max(tri.data)-np.min(tri.data))
    return(p_max, pk_pk)
#*****************************************************************************



#*****************************************************************************
# Función para procesar la información (Version 2)
# Fecha: Jul 14 2023 
#*****************************************************************************
def info_proc_2(tri, tri_c):
    # Se importa método peak_ground_motion para calcular parámetros pico de movimiento
    from obspy.signal.freqattributes import peak_ground_motion
    # Se importa método classic_sta_lta para calcular parámetros STA/LTA
    from obspy.signal.trigger import classic_sta_lta
    # Para STA/LTA se sugiere valores iniciales de 1 sec y 30 sec.
    # Se importa scipy para realizar integral acumulativa
    from scipy import integrate
    from obspy.signal.invsim import cosine_taper
    # Nuevo Julio 12 de 2023
    from obspy.signal.trigger import plot_trigger
    from obspy.signal.trigger import recursive_sta_lta 
    import time
    # Convierte Raw a gravedad. Por defecto en m/s2
    # Los datos son copia de las tramas (asi se evita modificar datos en tramas 
    # originales)
    tri_aux = Raw2grav(tri_c) 
    tri = Raw2grav(tri)
    
    # Modificación sugerida en artículo obspy
    # Para SSI 5% para cada end (0.1, 0.2 sería 10%) - En ejemplo aparece 0.05
    tri_aux.data = tri_aux.data * cosine_taper(len(tri_aux.data), 0.1)
    tri_aux_fil = tri_aux.copy()
    # Filtro paso banda. Para SSI utilizar 0.2 y 50. En el ej. valores 0.5 y 10 
    tri_aux_fil.filter('bandpass', freqmin=0.2, freqmax=50)
    
    # Velocidad (primera integral) usando tr.integrate (SUGERIDA)
    tri_aux_v2 = tri_aux_fil.copy()
    tri_aux_v2.integrate(method='cumtrapz')
    # Desplazamiento (segunda integral) (SUGERIDA)
    tri_aux_d2 = tri_aux_v2.copy()
    tri_aux_d2.integrate(method='cumtrapz')
    
    # Calculo de maximal displacement - PGD, velocity PGV, aceceleration PGA_1 
    # y despPGA, PGV PGD
    # Función retorna Peak ground aceceleration, maximal displacement, velocity y 
    # acceleration
    freq = 0.3 # Valores típicos (0.3, 1.0, 3.0 Hz)
    
    # Cálculo de parámetros usando posición corregida
    (PGAc, PGDc, PGVc, Ac) = peak_ground_motion(tri_aux_d2, tri_c.stats.delta, freq)         
    
    # Cálculo de parámetros de aceleración
    (A_p, A_pp) = tr_param(tri_aux)
    # Aceleración pico al cuadrado para posterior cálculo de la resultante
    A_p_2 = np.square(A_p,dtype=np.float64)
    
    # Cálculo de parámetros de velocidad
    (V_p, V_pp) = tr_param(tri_aux_v2)
    # Velocidad pico al cuadrado para posterior cálculo de la resultante
    V_p_2 = np.square(V_p,dtype=np.float64)
     
    # Cálculo de parámetros de desplazamiento
    (D_p, D_pp) = tr_param(tri_aux_d2)
    
    # Calcular STA/LTA Classic
    # Sampling rate
    df = tri_c.stats.sampling_rate  

    # Calcular recursive STA/LTA 
    cft3 = recursive_sta_lta(tri_c.data, int(0.5 * df), int(5 * df))

#*****************************************************************************
# Identificando eventos 
# Conveniente implementar como función
# Sep-30-2023
#*****************************************************************************    
    # Se define umbral para filtrar eventos
    #umbral = 3 #Debe fijarse en 10 (pruebas con 5.5)
    # abril 30 de 2024 se cambia temporalement
    umbral = 5.5
    # Ventana de tiempo superior para el evento (segundos)
    w_sup = 40
    # Ventana de tiempo inferior para el evento (segundos)
    w_inf = 20

    # Se encuentra el índice del los elementos que cumplen ser eventos según
    # el umbral definido previamente   

    if cft3.max() > umbral:
        # Lista de índices con los eventos
        event_list = np.where(cft3 == cft3.max())[0]
        # Primer evento de la lista (Se podría obtener de la variable anterior)
        event_first = np.where(cft3 == cft3.max())[0][0]
        # Indice del evento convertido a escala de tiempo
        event_first_t = np.uint(np.multiply(event_first,tri_c.stats.delta))
        # Se calculan los límites inferior y superior para extraer información
        # del evento en la trama corregida
        idx_sup = np.uint(event_first + (np.divide(1,tri_c.stats.delta) * w_sup)) 
        # Comprobando que el índice superior no sea mayor al tamaño de la trama
        if idx_sup >= tri_c.stats.npts:
            idx_sup = tri_c.stats.npts
        # Límite inferior    
        idx_inf = np.int32(event_first - (np.divide(1,tri_c.stats.delta) * w_inf))
        # Se comprueba que el índice inferior no sea menor que 0
        if idx_inf < 0:
            idx_inf = 0 
        event_idx = [event_first, event_first_t, w_inf, idx_inf, w_sup, idx_sup]
    else:
        event_idx = []

    # Revisar # Starttime actual extraído de la traza donde se detectó el evento
    time = UTCDateTime(tri.stats.starttime)
#*****************************************************************************            
    # Nuevos datos desde Jul 14 de 2023
    data = [tri.data.mean(), tri_c.data.mean(), A_p, A_pp, A_p_2, V_p, V_pp, PGVc, V_p_2, D_p, D_pp, PGDc, cft3.max()]
    # Nuevos datos incluyendo al inicio julday, fecha (y/m/d) y hora (h:m:s)
    # Fecha: 2023-10-13 
    data = [time.julday, time.date.isoformat(), time.time.isoformat(), tri.data.mean(), tri_c.data.mean(), A_p, A_pp, A_p_2, V_p*1000, V_pp, PGVc, V_p_2, D_p, D_pp, PGDc, cft3.max()]    ##V_P MULTIPLICADO POR 1000
                #0              1                       2                   3                   4           5    6     7       8       9      10     11    12   13    14      15
    return(data, event_idx)
#*****************************************************************************



#*****************************************************************************
# Funcion Javier obtener Velocidad en función de tiempo y frecuencia
# Fecha: sep 26 2023 
# Se utilizará en el programa principal
# ¡¡¡¡¡ En construcción !!!!! Funciona OK
#*****************************************************************************
def VelTimeFreq(tr_actual_c):
    # current_function_name=inspect.currentframe().f_code.co_name
    # print(f"\n     INICIO FUNCIÓN {Fore.YELLOW}{current_function_name}{Style.RESET_ALL}")
    import tkinter
    import tkinter.filedialog
    import os
    import numpy as np
    import obspy.core as OC
    from obspy.signal.invsim import cosine_taper
    from scipy.fft import fft, fftfreq, rfftfreq, rfft
    import matplotlib.pyplot as plt
    from scipy import signal
    from datetime import timedelta
#*****************************************************************************
    # Incluye JEA
    # Sep-27-2023

    flag = 0
    fSample = 250
    sizeWin = int(fSample*2)
    hannWin = np.hanning(sizeWin)
    xFrecVelMax = []
    yVelMax = []

    sos = signal.butter(5, [1,50], 'bandpass', fs=250, output='sos')
    
    #Funcion para calcular VelMax y su posición
    def calculate_velocity(tr):
        tr.data = tr.data * cosine_taper(len(tr.data), 0.1)
        tr.data = signal.sosfilt(sos, tr.data)
        tr.integrate(method='cumtrapz')
        max_vel_pos = np.argmax(tr.data)
        max_vel = tr.data[max_vel_pos] * 1000
        return max_vel_pos, max_vel

    #Funcion para calcular velocidad frecuencia MAX según la VMAxPos 
    def calculate_frequency(tr, max_vel_pos_global):

        sup = tr.data.size
        WinIni = max_vel_pos_global - fSample
        WinFin = max_vel_pos_global + fSample
        if WinIni < 0:
            WinIni = 0
            WinFin = min(500, sup)
        if WinFin > sup:
            WinIni = max(sup - 500, 0)
            WinFin = sup
        VelHann = tr.data[WinIni:WinFin] * hannWin

        N = len(VelHann)
        fftVelHannAmp = 2 * np.abs(rfft(VelHann)) / N
        xf = rfftfreq(N, d=1 / fSample)

        max_frec_vel_pos = np.argmax(fftVelHannAmp)
        max_frec_vel = xf[max_frec_vel_pos]
        return max_frec_vel

    #Crear una copia de las tramas recortadas y eliminar la tendencia lineal

    # Convierte Raw a gravedad. Por defecto en m/s2
    traces = Raw2grav(tr_actual_c)
    # print(f'Traces VelTimeFre: {traces}')

    # print('\n'.join(map(str, traces)))

    #Llama la función calculate_velocity para cada trama/canal
    max_vel_positions,max_vel_values = calculate_velocity(traces)

    #Calcula la posicion de la máxima velocidad al comparar los 3 canales
    MaxVelPos_global = max_vel_positions

    #Calcula la frecuencia de cada canal según la posición de velocidad máxima global
    frequencies = calculate_frequency(traces,MaxVelPos_global)
    
    xFrecVelMax.extend([frequencies])
    
    #Calcula las velocidades máximas según la posición de la velocidad máxima global
    MaxVel_1 = traces.data[MaxVelPos_global]*1000

    yVelMax.extend([MaxVel_1])

    return (xFrecVelMax,yVelMax)
#*****************************************************************************



#*****************************************************************************
# Función para calcular inclinación
# Se calcula inclinación (tilt) de los datos. Ene-13-2023
# Según Application note Analog Device AN-1057. Eq. 11, 12 y 13. Pág. 7
# Fecha: Jun 17 2023 
# Verificar el orden y correspondencia de las tramas para los cálculos.
# @author: jeamac
#*****************************************************************************
def tilt(tr1,tr2,tr3):
# tr1, tr2, tre = trace1, trace2, trace3 (miniseed format)
#*****************************************************************************
    # Version 2 para encontrar el tamaño mínimo en Traces (mejora velocidad)    
    min_tr = np.min([np.max(tr1.stats.npts), np.max(tr2.stats.npts), np.max(tr3.stats.npts)])
    # Cálculos completos
    # Se cambia el tipo de datos a int64 para solucionar desborde con np.square
    tilt_1 = np.rad2deg(np.arctan(np.divide(tr1[:min_tr],np.sqrt((np.square(tr2[:min_tr],dtype=np.float64)+np.square(tr3[:min_tr],dtype=np.float64))))))
    tilt_2 = np.rad2deg(np.arctan(np.divide(tr2[:min_tr],np.sqrt((np.square(tr1[:min_tr],dtype=np.float64)+np.square(tr3[:min_tr],dtype=np.float64))))))
    tilt_3 = np.rad2deg(np.arctan(np.divide(np.sqrt((np.square(tr1[:min_tr],dtype=np.float64)+np.square(tr2[:min_tr],dtype=np.float64))),tr3[:min_tr])))
    return(tilt_1, tilt_2, tilt_3, min_tr)
#*****************************************************************************




#*****************************************************************************
#Función para adquirir lo datos de los canales y procesar inf
#@autor:oszu
#*****************************************************************************
def adquire_data(selected_files, last_time=60):

    current_function_name=inspect.currentframe().f_code.co_name
    print(f"\n     INICIO FUNCIÓN {Fore.YELLOW}{current_function_name}{Style.RESET_ALL}")

    st1 = obspy.Stream()
    for file in selected_files:
        st1 += obspy.read(file)

    tr1 = obspy.Stream()
    tr2 = obspy.Stream()
    tr3 = obspy.Stream()
    
    #Asignar las tramas a las variables correctas según el canal
    for trace in st1:
        if trace.stats.channel=="HNE":
            tr1+=trace
        elif trace.stats.channel=="HNN":
            tr2+=trace
        elif trace.stats.channel=="HNZ":
            tr3+=trace

    print(f"\nTrace adquire tr1: {tr1}")
    print(f"\nTrace adquire tr2: {tr2}")
    print(f"\nTrace adquire tr3: {tr3}")

    df = pd.DataFrame(columns=['Channel', 'VP', 'End_Time','FP','V_Max','Start_Time'])

    print(f"\nTotal tramas por procesar: {len(tr1)+len(tr2)+len(tr3)}\n")

    for trace in tqdm(st1, desc="Procesando tramas", ncols=100): # Ajuste de tamaño de la barra de progreso

        start_time= trace.stats.starttime
        tr_actual=trace

        while len(tr_actual.data) > 1:

            end_time = start_time + timedelta(minutes=last_time)
            tr_actual = trace.slice(starttime=start_time, endtime=end_time)
            
            try:
                
                tr_actual_c=obspy.Stream(traces=[tr_actual.copy()])
                tr_actual_c.merge(fill_value='interpolate')
                tr_actual_c.detrend("linear")
                tr_actual_c=tr_actual_c[0]

                pd1, eve_inf1 = info_proc_2(tr_actual,tr_actual_c) 

                FrecMax,VelMax = VelTimeFreq(tr_actual_c)

                row= pd.DataFrame({
                    'Channel': [tr_actual.stats.channel],
                    'VP': [pd1[8]],
                    'End_Time': [end_time.isoformat()],
                    'FP':FrecMax[0],
                    'V_Max':VelMax[0],
                    'Start_Time':[start_time]
                })
                df=pd.concat([df, row],ignore_index=True)
                # print(f"Trace adquire: {tr_actual}")
                    
            except:
                # print(f"\n{Fore.RED}Enter Except {Style.RESET_ALL}- Finish time trace {tr_actual.stats.endtime}")
                
                pass

            start_time = end_time

    print(f"\n DataFrame resultante:")
    print(f"{df}")
    return(df)
#*****************************************************************************



#*****************************************************************************
#Función para crear graficas  ABNT_BNR
#@autor:oszu
#*****************************************************************************
def plot_dataframe_ABNT_BNR(df,save_path):
    # Definir los colores específicos para cada canal
    palette = {'HNE': 'blue', 'HNN': 'green', 'HNZ': 'red'}
    markers = {'HNE': 's', 'HNN': '^', 'HNZ': 'o'}  # cuadrado, triángulo, círculo

    # Configuración del gráfico
    plt.figure(figsize=(10, 6))
    
    # Crear el gráfico de dispersión
    scatter = sns.scatterplot(data=df, x='FP', y='VP', hue='Channel', palette=palette, style='Channel', markers=markers)

    # Configurar la escala del eje Y como logarítmica
    plt.yscale('log')
    plt.ylim(0.1, 100)
    
    # Cambiar el formato de los ticks del eje Y
    ax = plt.gca()
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.1f}'.format(y)))

    # Añadir la línea del umbral
    umbral_x = [0, 2, 4, 6.2, 8.4, 10.6, 12.8, 15, 20, 25, 30, 35, 40, 50, 60, 70]
    umbral_y = [15, 15, 15, 16, 17, 18, 19, 20, 26, 32, 38, 44, 50, 50, 50, 50]
    plt.plot(umbral_x, umbral_y, label='Threshold', color='Black', linestyle='solid')

    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Título y etiquetas de los ejes
    plt.title('Brazilian Standard ABNT NBR 9653: 2018', fontsize=16, fontweight='bold')
    plt.xlabel('Frequency (Hz)', fontsize=14, fontweight='bold')
    plt.ylabel('Velocity (mm/s)', fontsize=14, fontweight='bold')
    
    plt.legend(title='Channel')

    # Mostrar gráfico
    abnt_bnr_image_path = os.path.join(save_path, 'graph_ABNT_BNR.png')
    plt.savefig(abnt_bnr_image_path)
    plt.close()
    return abnt_bnr_image_path
#*****************************************************************************



#*****************************************************************************
#Función para crear graficas 
#@autor:oszu
#*****************************************************************************
def plot_dataframe_CETESB(df, save_path):
    # Definir los colores específicos para cada canal
    palette = {'HNE': 'blue', 'HNN': 'green', 'HNZ': 'red'}
    markers = {'HNE': 's', 'HNN': '^', 'HNZ': 'o'}  # cuadrado, triángulo, círculo
    
    # Configuración del gráfico
    plt.figure(figsize=(10, 6))
    
    # Crear el gráfico de dispersión
    scatter = sns.scatterplot(data=df, x='FP', y='VP', hue='Channel', palette=palette, style='Channel', markers=markers)

    # Configurar la escala del eje Y como logarítmica
    plt.yscale('linear')
    plt.ylim(0, 5)
    
    # Cambiar el formato de los ticks del eje Y
    ax = plt.gca()
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.1f}'.format(y)))

    # Añadir la línea del umbral
    umbral_x = [0, 2, 4, 6.2, 8.4, 10.6, 12.8, 15, 20, 25, 30, 35, 40, 50, 60, 70]
    umbral_y = [4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2]
    plt.plot(umbral_x, umbral_y, label='Threshold', color='Black', linestyle='solid')

    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Título y etiquetas de los ejes
    plt.title('CETESB D7.013/2015', fontsize=16, fontweight='bold')
    plt.xlabel('Frequency (Hz)', fontsize=14, fontweight='bold')
    plt.ylabel('Velocity (mm/s)', fontsize=14, fontweight='bold')
    
    plt.legend(title='Channel')
    # Mostrar gráfico
    cetesb_image_path = os.path.join(save_path, 'graph_CETESB.png') 
    plt.savefig(cetesb_image_path) 
    plt.close()
    return cetesb_image_path
#*****************************************************************************



#*****************************************************************************
#Función para crear graficas 
#@autor:oszu
#*****************************************************************************
def plot_vp_vs_time(df, save_path):
    # Convertir 'End_Time' a cadena para usar en la gráfica
    df['End_Time_str'] = df['End_Time'].astype(str)

    # Convertir 'End_Time' a objeto datetime manejando ambos formatos
    df['End_Time'] = pd.to_datetime(df['End_Time_str'], errors='coerce', format='%Y-%m-%dT%H:%M:%S.%f')
    df['End_Time'] = df['End_Time'].fillna(pd.to_datetime(df['End_Time_str'], errors='coerce', format='%Y-%m-%dT%H:%M:%S'))

    # Redondear los tiempos a segundos para evitar diferencias mínimas
    df['End_Time_rounded'] = df['End_Time'].dt.floor('S').astype(str) 

    # Asegurarnos de que los tiempos estén ordenados correctamente
    df = df.sort_values(by=['End_Time_rounded', 'Channel'])

    # Definir los colores específicos para cada canal
    palette = {'HNE': 'blue', 'HNN': 'green', 'HNZ': 'red'}
    markers = {'HNE': 's', 'HNN': '^', 'HNZ': 'o'}  # cuadrado, triángulo, círculo

    # Configuración del gráfico
    plt.figure(figsize=(12, 8))

    # Crear el gráfico de dispersión
    scatter = sns.scatterplot(data=df, x='End_Time_rounded', y='VP', hue='Channel', palette=palette, style='Channel', markers=markers)

    # Ajustar los ticks del eje X para que solo muestren unos 4 valores
    unique_times = df['End_Time_str'].unique()
    num_ticks = len(unique_times)
    tick_positions = [0, int(num_ticks * 0.33), int(num_ticks * 0.66), num_ticks - 1]
    tick_labels = [unique_times[pos].split('T')[1] for pos in tick_positions]

    ax = plt.gca()
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, fontsize=8)

    # Mantener todos los puntos en el gráfico aunque no haya etiquetas en el eje X
    ax.xaxis.set_tick_params(length=0)

    # Añadir la cuadrícula
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    min_date = df['Start_Time'].min().strftime('%Y-%m-%d %H:%M:%S')
    fecha = datetime.strptime(min_date, "%Y-%m-%d %H:%M:%S")
    fecha_formateada = fecha.strftime("%Y-%m-%d")

    # Título y etiquetas de los ejes
    plt.title(f'Peak Particle Velocity Trend - {fecha_formateada}', fontsize=16, fontweight='bold')
    plt.xlabel('Time UTC', fontsize=14, fontweight='bold')
    plt.ylabel('Velocity (mm/s)', fontsize=14, fontweight='bold')

    plt.subplots_adjust(bottom=0.25)

    # Guardar gráfico como imagen
    vp_time_image_path = os.path.join(save_path, 'graph_VP_vs_Time.png')
    plt.savefig(vp_time_image_path)
    plt.close()

    # Imprimir información de los puntos utilizados
    # for i, row in df.iterrows():
    #     print(f"Plotted - VP: {row['VP']} End_Time_rounded: {row['End_Time_rounded']} Channel: {row['Channel']}")

    return vp_time_image_path
#*****************************************************************************


#*****************************************************************************
class PDF(FPDF): 
    def header(self): 
        self.set_font("Arial", "B", 16)     
        self.cell(0, 10, "Particle velocity analysis", 0, 1, "C") 
    def footer(self): 
        self.set_y(-15) 
        self.set_font("Arial", "I", 8) 
        self.cell(0, 10, f"Page {self.page_no()}", 0, 0, "C")
#*****************************************************************************
        

#*****************************************************************************
#Función para generar pdf
#@autor:oszu
#*****************************************************************************
def create_pdf(df, data_project, save_path, logo_path):
    current_function_name=inspect.currentframe().f_code.co_name
    print(f"\n     INICIO FUNCIÓN {Fore.YELLOW}{current_function_name}{Style.RESET_ALL}")
    # Crear los gráficos primero
    abnt_bnr_image_path = plot_dataframe_ABNT_BNR(df, save_path) 
    cetesb_image_path = plot_dataframe_CETESB(df, save_path)
    vp_time_image_path = plot_vp_vs_time(df, save_path)

    # Crear el PDF
    pdf = PDF()
    pdf.add_page()
    pdf.set_auto_page_break(auto=True, margin=15)

    pdf.image(logo_path, x=170, y=10, w=30)

    # Agregar título y fecha del análisis
    pdf.set_font("Arial", size=16,style='B')
    # pdf.cell(200, 10, txt=f"Particle velocity analysis", ln=True, align='C')
    pdf.cell(200, 10, ln=True)  # Línea vacía
    
    min_date = df['Start_Time'].min().strftime('%Y-%m-%d %H:%M:%S')
    try: 
        max_date = df['End_Time'].max().strftime('%Y-%m-%d %H:%M:%S')
    except:
        max_date = df['End_Time'].max()
    
    pdf.set_font("Arial", size=12)

    channels = ['HNE', 'HNN', 'HNZ'] 
    notes = [] 
    for channel in channels: 
        df_channel = df[df['Channel'] == channel].reset_index(drop=True)
        if not df_channel.empty:
            max_vp_index = df_channel['VP'].idxmax() 

            max_vp = df_channel['VP'].max() 
            max_fp = df_channel['FP'].iloc[max_vp_index] 
            max_end_time = df_channel['End_Time'].iloc[max_vp_index] 
            
            notes.append(f"Max velocity {round(max_vp, 2)} mm/s at {max_fp} Hz in {channel} - {max_end_time}") 

    info_fields=[
        ("Project:", data_project[0]),
        ("Time analyzed:", f"From {min_date} to {max_date} UTC"),
        ("Client:", data_project[1]),
        ("Location:", data_project[2]),
        ("Notes: ",notes[0]),
        ("",notes[1]),
        ("",notes[2])
    ] 

    for field_name, field_value in info_fields:
        pdf.set_font("Arial", size=12, style='B')  
        pdf.cell(40, 6, txt=field_name, ln=False, align='L')  
        pdf.set_font("Arial", size=12)
        pdf.cell(0, 6, txt=field_value, ln=True, align='L')

    # Agregar las imágenes de los gráficos
    pdf.image(abnt_bnr_image_path, x=10, y=71, w=180)
    pdf.image(cetesb_image_path, x=10, y=176, w=180)

    # Agregar una nueva página para la gráfica de VP vs Time 
    pdf.add_page() 
    pdf.image(logo_path, x=170, y=10, w=30)
    pdf.image(vp_time_image_path, x=10, y=42, w=180)

    # Guardar PDF en una ubicación específica
    fecha = datetime.strptime(min_date, "%Y-%m-%d %H:%M:%S")
    fecha_formateada = fecha.strftime("%Y-%m-%d")
    pdf_path = os.path.join(save_path, f'Report_{fecha_formateada}.pdf') 
    pdf.output(pdf_path) 
    print(f"Reporte guardado en: {pdf_path}")

    os.startfile(pdf_path)
#*****************************************************************************



#*****************************************************************************
# Función para seleccionar una carpeta
def seleccionar_carpeta():
    root = tk.Tk()
    root.withdraw()  # Ocultar la ventana principal
    root.attributes("-topmost",True)
    carpeta_seleccionada = filedialog.askdirectory(title="Seleccionar carpeta")
    root.destroy()
    return carpeta_seleccionada
#*****************************************************************************

#*****************************************************************************
# Función para obtener todos los archivos .miniseed de la carpeta seleccionada
def obtener_archivos_carpeta(carpeta):
    archivos = []
    for root, dirs, files in os.walk(carpeta):
        for file in files:
            if file.endswith(".miniseed"):
                archivos.append(os.path.join(root, file))
    return archivos
#*****************************************************************************

#*****************************************************************************
# Función para filtrar archivos por día y canal
#*****************************************************************************
def filtrar_archivos_por_dia_y_canal(archivos, dia):
    archivos_filtrados = []
    canales = ["HNE", "HNN", "HNZ"]
    
    for archivo in archivos:
        # Extraer el día y el canal del nombre del archivo
        nombre_archivo = os.path.basename(archivo)
        partes = nombre_archivo.split('.')
        dia_archivo = int(partes[6])
        canal = partes[3]
        
        # Verificar si el archivo corresponde al día y canal deseados
        if dia_archivo == dia and canal in canales:
            archivos_filtrados.append(archivo)
    
    return archivos_filtrados
#*****************************************************************************

#*****************************************************************************
def obtener_dias_disponibles(archivos):
    dias_disponibles = set()
    for archivo in archivos:
        # Extraer el día del nombre del archivo
        nombre_archivo = os.path.basename(archivo)
        partes = nombre_archivo.split('.')
        try:
            if len(partes) > 6 and partes[6].isdigit():
                dia = int(partes[6])  # Convertir 'partes[6]' a entero si es un número válido
                dias_disponibles.add(dia)
        except (IndexError, ValueError):
            # Ignorar archivos que no tienen el formato esperado
            continue
    return sorted(dias_disponibles)
#*****************************************************************************




@pagina2_bp.route('/pagina2')
def pagina2():
    return render_template('pagina2.html')

@pagina2_bp.route('/suma2', methods=['POST'])
def suma2():
    num1 = request.form.get('num1', type=int)
    num2 = request.form.get('num2', type=int)
    resultado = num1 + num2
    return f'La suma de {num1} y {num2} es {resultado}.'
