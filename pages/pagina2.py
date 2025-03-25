from flask import Blueprint, render_template, request

from flask import Flask, render_template, request, redirect, jsonify
import matplotlib.pyplot as plt
import inspect

from scipy.integrate import cumulative_trapezoid
import numpy as np
from obspy.signal.invsim import cosine_taper
from scipy.fft import fft, fftfreq, rfftfreq, rfft
import matplotlib.pyplot as plt
from datetime import timedelta

from datetime import timedelta, datetime
import obspy
import numpy as np
from scipy.signal import butter, filtfilt

from flask_socketio import SocketIO, emit

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
import seaborn as sns
from fpdf import FPDF
from obspy.signal.invsim import cosine_taper
from scipy.signal import butter, filtfilt
import matplotlib.dates as mdates
from obspy.signal.trigger import recursive_sta_lta 

pagina2_bp = Blueprint('pagina2_bp', __name__)


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



def Raw2grav_g(tr):
    # Convierte la señal cruda a aceleración en gravidades (g)
    tr_p = tr.copy()
    factor = 520000   # cuentas a g
    tr_p.data = np.divide(tr_p.data, factor)
    return tr_p



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
        event_indices = np.where (cft3>umbral)[0]

        event_idx = []
        for ev in event_indices:
            # Tiempo del evento en segundos (según la traza completa)
            ev_time = ev * tri_c.stats.delta
            
            # Convertir las ventanas de tiempo a número de muestras:            
            win_sup_samples = int((1 / tri_c.stats.delta) * w_sup)
            win_inf_samples = int((1 / tri_c.stats.delta) * w_inf)

            # Se calculan los límites inferior y superior para extraer información
            # del evento en la trama corregida
            idx_sup = ev + win_sup_samples
            # Comprobando que el índice superior no sea mayor al tamaño de la trama
            if idx_sup >= tri_c.stats.npts:
                idx_sup = tri_c.stats.npts
            # Límite inferior    
            idx_inf = ev - win_inf_samples
            # Se comprueba que el índice inferior no sea menor que 0
            if idx_inf < 0:
                idx_inf = 0 
            event_idx.append([ev, ev_time, w_inf, idx_inf, w_sup, idx_sup])
    else:
        event_idx = []

    # Revisar # Starttime actual extraído de la traza donde se detectó el evento
    time = UTCDateTime(tri.stats.starttime)

#*****************************************************************************            
    # Nuevos datos desde Jul 14 de 2023
    data = [tri.data.mean(), tri_c.data.mean(), A_p, A_pp, A_p_2, V_p, V_pp, PGVc, V_p_2, D_p, D_pp, PGDc, cft3.max()]
    # Nuevos datos incluyendo al inicio julday, fecha (y/m/d) y hora (h:m:s)
    # Fecha: 2023-10-13 
    V_p = V_p * 1000 #Pasar a mm
    A_p = A_p / 9.81002509 #Pasar a gravedades
    D_p = D_p * 1000 #Pasar a mm
    data = [time.julday, time.date.isoformat(), time.time.isoformat(), tri.data.mean(), tri_c.data.mean(), A_p, A_pp, A_p_2, V_p, V_pp, PGVc, V_p_2, D_p, D_pp, PGDc, cft3.max()]    ##V_P MULTIPLICADO POR 1000
                #0              1                       2                   3                   4           5    6     7       8   9      10     11    12   13    14      15
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
        return max_vel_pos, max_vel, tr

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
    max_vel_positions,max_vel_values,traces = calculate_velocity(traces)

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
def adquire_data(selected_files, flag, last_time):

    current_function_name=inspect.currentframe().f_code.co_name
    print(f"\n     INICIO FUNCIÓN {Fore.YELLOW}{current_function_name}{Style.RESET_ALL} para {Fore.YELLOW}{flag}{Style.RESET_ALL}")

    if flag == "todo":
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

    elif flag =="eventos":
        
        st1 = obspy.Stream(traces=list(selected_files.values()))

    print(st1)
    df = pd.DataFrame(columns=['Channel', 'VP', 'DP','AP','FP','V_Max','Start_Time','End_Time'])

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

                # data = [time.julday, time.date.isoformat(), time.time.isoformat(), tri.data.mean(), tri_c.data.mean(), A_p, A_pp, A_p_2, V_p*1000, V_pp, PGVc, V_p_2, D_p, D_pp, PGDc, cft3.max()]    ##V_P MULTIPLICADO POR 1000
                #               0              1                       2                   3                   4           5    6     7       8       9      10     11    12   13    14      15

                row= pd.DataFrame({
                    'Channel': [tr_actual.stats.channel],
                    'VP': [pd1[8]],
                    'DP': [pd1[12]],
                    'AP': [pd1[5]],                
                    'FP':FrecMax[0],
                    'V_Max':VelMax[0],
                    'Start_Time':[start_time.isoformat()],
                    'End_Time': [end_time.isoformat()]
                })
                df=pd.concat([df, row],ignore_index=True)
                        
            except Exception as e:
                # print(f"\n{Fore.RED}Enter Except {Style.RESET_ALL}- Finish time trace {tr_actual.stats.endtime} {e}")       
                pass

            start_time = end_time

    print(f"\n DataFrame resultante:")
    print(f"{df}")
    # Seleccionar la fila con VP máximo por cada canal
    df_max = df.loc[df.groupby("Channel")["VP"].idxmax()]
    print("Filas con VP máximo por canal:")
    print(df_max)

    return(df)
#*****************************************************************************



#*****************************************************************************
def searchSTA_LTA(tri_c):
    
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
        event_indices = np.where (cft3>umbral)[0]

        event_idx = []
        for ev in event_indices:
            # Tiempo del evento en segundos (según la traza completa)
            ev_time = ev * tri_c.stats.delta
            
            # Convertir las ventanas de tiempo a número de muestras:            
            win_sup_samples = int((1 / tri_c.stats.delta) * w_sup)
            win_inf_samples = int((1 / tri_c.stats.delta) * w_inf)

            # Se calculan los límites inferior y superior para extraer información
            # del evento en la trama corregida
            idx_sup = ev + win_sup_samples
            # Comprobando que el índice superior no sea mayor al tamaño de la trama
            if idx_sup >= tri_c.stats.npts:
                idx_sup = tri_c.stats.npts
            # Límite inferior    
            idx_inf = ev - win_inf_samples
            # Se comprueba que el índice inferior no sea menor que 0
            if idx_inf < 0:
                idx_inf = 0 
            event_idx.append([ev, ev_time, w_inf, idx_inf, w_sup, idx_sup, tri_c.stats.channel])
            #event_idx.append([ev, ev_time, w_inf, idx_inf, w_sup, idx_sup, tri_c.stats.channel])
            #event_idx.append([0,      1,     2,      3,      4,      5,        6])
    else:
        event_idx = []

    return event_idx

def group_events(event_idx_lists, time_threshold=60):
    """
    Agrupa eventos de diferentes canales si sus tiempos están dentro de 'time_threshold' segundos.
    
    event_idx_lists: lista de listas, por ejemplo, [event_idx_HNE, event_idx_HNN, event_idx_HNZ]
    time_threshold: diferencia en segundos para considerar que dos eventos son el mismo.
    
    Retorna una lista de grupos. Cada grupo es una lista de eventos (con su información completa).
    """
    # Combinar todos los eventos en una sola lista
    all_events = []
    for lst in event_idx_lists:
        all_events.extend(lst)
    # Ordenar por ev_time (segundo elemento)
    all_events.sort(key=lambda x: x[1])
    
    groups = []
    current_group = []
    for event in all_events:
        if not current_group:
            current_group.append(event)
        else:
            # Si la diferencia con el último evento del grupo es menor o igual a time_threshold, se agrupa.
            if event[1] - current_group[-1][1] <= time_threshold:
                current_group.append(event)
            else:
                groups.append(current_group)
                current_group = [event]
    if current_group:
        groups.append(current_group)
    return groups

def filter_event_groups(groups):
    """
    Filtra los grupos para conservar únicamente aquellos que tengan eventos de al menos dos canales distintos.
    Cada evento se espera que tenga al menos 7 elementos, siendo el séptimo (índice 6) el canal.
    """
    filtered = []
    for group in groups:
        unique_channels = set()
        for event in group:
            if len(event) >= 7:
                unique_channels.add(event[6])
        if len(unique_channels) >= 2:
            filtered.append(group)
    return filtered


def extract_event_segments(concatenated_traces, event_group):
    """
    Dado un grupo de eventos (con índices similares en al menos dos canales),
    extrae el segmento de la traza correspondiente en cada canal.
    
    concatenated_traces: un Stream que contiene una traza por canal (por ejemplo, HNE, HNN, HNZ).
    event_group: lista de eventos agrupados (cada uno es [ev, ev_time, w_inf, idx_inf, w_sup, idx_sup, channel])
    
    Retorna un diccionario con claves 'HNE', 'HNN' y 'HNZ' y sus correspondientes segmentos.
    """
    # Para el grupo, definimos el límite inferior global y el superior global:
    group_idx_inf = min(event[3] for event in event_group)
    group_idx_sup = max(event[5] for event in event_group)
    
    segments = {}
    for trace in concatenated_traces:
        chan = trace.stats.channel
        # Calcular tiempo de inicio y fin usando la delta de la traza
        start_time = trace.stats.starttime + group_idx_inf * trace.stats.delta
        end_time = trace.stats.starttime + group_idx_sup * trace.stats.delta
        segments[chan] = trace.slice(starttime=start_time, endtime=end_time)
    return segments
#*****************************************************************************



#*****************************************************************************
def detect_event(selected_files, save_path, utc_time, last_time, logo_path, data_project):
    last_time=0.05
    current_function_name=inspect.currentframe().f_code.co_name
    print(f"\n     INICIO FUNCIÓN {Fore.YELLOW}{current_function_name}{Style.RESET_ALL}")

    st1 = obspy.Stream()
    for file in selected_files:
        st1 += obspy.read(file)

    event_idx_HNE=[]
    event_idx_HNN=[]
    event_idx_HNZ=[]

    for trace in tqdm(st1, desc="Procesando tramas - Buscando eventos", ncols=100): # Ajuste de tamaño de la barra de progreso

        tr_actual_c=obspy.Stream(traces=[trace.copy()])
        tr_actual_c.merge(fill_value='interpolate')
        tr_actual_c.detrend("linear")
        tr_actual_c=tr_actual_c[0]

        event_idx = searchSTA_LTA(tr_actual_c)
        #event_idx.append([ev, ev_time, w_inf, idx_inf, w_sup, idx_sup, tri_c.stats.channel])
        #                ([0,      1,     2,      3,      4,      5,        6])

        if event_idx:
            if tr_actual_c.stats.channel == 'HNE':
                event_idx_HNE.extend(event_idx)
            elif tr_actual_c.stats.channel == 'HNN':
                event_idx_HNN.extend(event_idx)
            elif tr_actual_c.stats.channel == 'HNZ':
                event_idx_HNZ.extend(event_idx)            

    concatenated_traces = Stream()
    for channel in ['HNE', 'HNN', 'HNZ']:
        channel_traces = st1.select(channel=channel)
        if channel_traces:
            # Puedes aplicar aquí un resample y merge, similar a tu lógica
            try:
                concatenated_trace = channel_traces.merge(fill_value='interpolate')[0]

                concatenated_trace.detrend('linear')
                concatenated_traces.append(concatenated_trace)
            except Exception as e:
                print(f"Error concatenando trazas para {channel}: {e}")

    if event_idx_HNE or event_idx_HNN or event_idx_HNZ:
        # Agrupar eventos de cada canal
        groups = group_events([event_idx_HNE, event_idx_HNN, event_idx_HNZ], time_threshold=60)
        # Filtrar los grupos para conservar solo aquellos con al menos dos canales
        filtered_groups = filter_event_groups(groups)
        print(f"Grupos de eventos (con al menos 2 canales): {len(filtered_groups)}")

    if not filtered_groups:
        print('No se encontraron eventos en las trazas')
    else:
        # Extraer segmentos para cada grupo y generar reportes
        for idx, group in enumerate(filtered_groups):
            print("Generando reportes de eventos")
            segments = extract_event_segments(concatenated_traces, group)
            print(f"segmentos de eventos: {segments}")
            # Por ejemplo, obtener la fecha y hora del evento a partir del primer evento del grupo:
            event_time_sec = group[0][1]
            # Convertir a datetime, asumiendo que la traza tiene starttime
            event_datetime = concatenated_traces[0].stats.starttime + event_time_sec
            # Ajustar a UTC-3 si es necesario:
            event_datetime_adjusted = event_datetime - timedelta(hours=utc_time)            
            # Generar un nombre de reporte basado en la fecha y hora del evento:
            report_date_str = event_datetime_adjusted.strftime("%Y%b%d")
            report_name = f"Event_Report_{report_date_str}T{event_datetime_adjusted.strftime('%H-%M-%S')}.pdf"
            # Llama a la función create_pdf_event (o similar) pasando los segmentos
            df = adquire_data(segments,"eventos", last_time)

            acceleration_path, velocity_path, displacement_path, fft_path = signals_generator(segments,"eventos", save_path, 
                                                                                                1, utc_time, last_time)
            create_pdf(df, data_project, save_path, logo_path,
                        acceleration_path, velocity_path, displacement_path, fft_path, last_time, utc_time,
                        report_name=report_name)

    return None


#*****************************************************************************
#Función para crear graficas  ABNT_BNR
#@autor:oszu
#*****************************************************************************
def plot_dataframe_ABNT_BNR(df,save_path):
    # Definir los colores específicos para cada canal
    palette = {'HNE': 'blue', 'HNN': 'green', 'HNZ': 'red'}
    markers = {'HNE': 's', 'HNN': '^', 'HNZ': 'o'}  # cuadrado, triángulo, círculo

    # Configuración del gráfico
    plt.figure(figsize=(12, 8))
    
    # Crear el gráfico de dispersión
    scatter = sns.scatterplot(data=df, 
                              x='FP', 
                              y='VP', 
                              hue='Channel', 
                              palette=palette, 
                              style='Channel', 
                              markers=markers,
                              s=150)

    # Configurar la escala del eje Y como logarítmica
    plt.yscale('log')
    if df['VP'].max() < 50:
        plt.ylim(0, 50)
    else:
        plt.ylim(0, df['VP'].max() + 1)
    
    plt.xscale('log',base=2)
    
    ax = plt.gca()
    ax.set_xticks([1, 2, 4, 8, 16, 32, 64, 80])  # Establecer manualmente los ticks
    ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))  # Mostrar sin exponente

    # Cambiar el formato de los ticks del eje Y
    ax = plt.gca()
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.1f}'.format(y)))

    # Añadir la línea del umbral
    umbral_x = [0, 2, 4, 6.2, 8.4, 10.6, 12.8, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80]
    umbral_y = [15, 15, 15, 16, 17, 18, 19, 20, 26, 32, 38, 44, 50, 50, 50, 50, 50]
    plt.plot(umbral_x, umbral_y, label='Threshold', color='Black', linestyle='solid')

    plt.grid(True, which='both', linestyle='--', linewidth=1.5, alpha=0.8)

    # Título y etiquetas de los ejes
    plt.title('ABNT NBR 9653: 2018', fontsize=16, fontweight='bold')
    plt.xlabel('Frequency (Hz)', fontsize=14, fontweight='bold')
    plt.ylabel('Velocity (mm/s)', fontsize=14, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(title='Channel',fontsize=15, title_fontproperties={'weight': 'bold', 'size': 15})

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
    plt.figure(figsize=(12, 8))
    
    # Crear el gráfico de dispersión
    scatter = sns.scatterplot(data=df, 
                              x='FP', 
                              y='VP', 
                              hue='Channel', 
                              palette=palette, 
                              style='Channel', 
                              markers=markers,
                              s=150
                              )

    # Configurar la escala del eje Y como logarítmica
    plt.yscale('linear')
    if df['VP'].max() < 5:
        plt.ylim(0, 5)
    else:
        plt.ylim(0, df['VP'].max() + 1)
    
    plt.xscale('log',base=2)
    
    ax = plt.gca()
    ax.set_xticks([1, 2, 4, 8, 16, 32, 64, 80])  # Establecer manualmente los ticks
    ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))  # Mostrar sin exponente

    # Cambiar el formato de los ticks del eje Y
    ax = plt.gca()
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.1f}'.format(y)))

    # Añadir la línea del umbral
    umbral_x = [0, 2, 4, 6.2, 8.4, 10.6, 12.8, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80]
    umbral_y = [4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2]
    plt.plot(umbral_x, umbral_y, label='Threshold', color='Black', linestyle='solid')

    plt.grid(True, which='both', linestyle='--', linewidth=1.5, alpha=0.8)

    # Título y etiquetas de los ejes
    plt.title('CETESB D7.013/2015', fontsize=16, fontweight='bold')
    plt.xlabel('Frequency (Hz)', fontsize=14, fontweight='bold')
    plt.ylabel('Velocity (mm/s)', fontsize=14, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.legend(title='Channel',fontsize=15, title_fontproperties={'weight': 'bold', 'size': 15})
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
def plot_vp_vs_time(df, save_path, utc_time):
    # Convertir 'Start_Time' a cadena para usar en la gráfica
    df['Start_Time_str'] = df['Start_Time'].astype(str)

    # Convertir 'Start_Time' a objeto datetime manejando ambos formatos
    df['Start_Time'] = pd.to_datetime(df['Start_Time_str'], errors='coerce', format='%Y-%m-%dT%H:%M:%S.%f')
    df['Start_Time'] = df['Start_Time'].fillna(pd.to_datetime(df['Start_Time_str'], errors='coerce', format='%Y-%m-%dT%H:%M:%S'))

    # Redondear los tiempos a segundos para evitar diferencias mínimas
    df['Start_Time_rounded'] = df['Start_Time'] - timedelta(hours=utc_time)
    df['Start_Time_rounded_str'] = df['Start_Time_rounded'].dt.strftime('%H:%M:%S')

    # Asegurarnos de que los tiempos estén ordenados correctamente
    df = df.sort_values(by=['Start_Time_rounded', 'Channel'])

    # Definir los colores específicos para cada canal
    palette = {'HNE': 'blue', 'HNN': 'green', 'HNZ': 'red'}
    markers = {'HNE': 's', 'HNN': '^', 'HNZ': 'o'}  # cuadrado, triángulo, círculo

    # Configuración del gráfico
    plt.figure(figsize=(13, 9))

    # Crear el gráfico de dispersión
    scatter = sns.scatterplot(data=df, 
                              x='Start_Time_rounded_str', 
                              y='VP', 
                              hue='Channel', 
                              palette=palette, 
                              style='Channel', 
                              markers=markers,
                              s=150
                              )

    # Ajustar los ticks del eje X para que solo muestren unos 4 valores
    unique_x_values = df['Start_Time_rounded_str'].unique().tolist()
    num_ticks = len(unique_x_values)
    tick_positions = [0, int(num_ticks * 0.25), int(num_ticks * 0.5), int(num_ticks * 0.75), num_ticks - 1]
    tick_labels = [unique_x_values[pos] for pos in tick_positions]

    ax = plt.gca()
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=0, fontsize=8)

    # Mantener todos los puntos en el gráfico aunque no haya etiquetas en el eje X
    ax.xaxis.set_tick_params(length=0)

    # Añadir la cuadrícula
    plt.grid(True, which='both', linestyle='--', linewidth=1.5, alpha=0.8)
        
    df['Start_Time_title'] = df['Start_Time'] - timedelta(hours=utc_time)
    min_date = df['Start_Time_title'].min().strftime('%Y-%m-%d %H:%M:%S')
    fecha_formateada = datetime.strptime(min_date, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")

    # Título y etiquetas de los ejes
    plt.title(f'Peak particle velocity trend - {fecha_formateada}', fontsize=18, fontweight='bold')
    plt.xlabel(f'Time UTC-{utc_time}', fontsize=16, fontweight='bold')
    plt.ylabel('Velocity (mm/s)', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(title='Channel',fontsize=15, title_fontproperties={'weight': 'bold', 'size': 15})

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
def plot_ap_vs_time(df, save_path, utc_time):
    # Convertir 'Start_Time' a cadena para usar en la gráfica
    df['Start_Time_str'] = df['Start_Time'].astype(str)

    # Convertir 'Start_Time' a objeto datetime manejando ambos formatos
    df['Start_Time'] = pd.to_datetime(df['Start_Time_str'], errors='coerce', format='%Y-%m-%dT%H:%M:%S.%f')
    df['Start_Time'] = df['Start_Time'].fillna(pd.to_datetime(df['Start_Time_str'], errors='coerce', format='%Y-%m-%dT%H:%M:%S'))

    # Redondear los tiempos a segundos para evitar diferencias mínimas
    df['Start_Time_rounded'] = df['Start_Time'] - timedelta(hours=utc_time)
    df['Start_Time_rounded_str'] = df['Start_Time_rounded'].dt.strftime('%H:%M:%S')

    # Asegurarnos de que los tiempos estén ordenados correctamente
    df = df.sort_values(by=['Start_Time_rounded', 'Channel'])

    # Definir los colores específicos para cada canal
    palette = {'HNE': 'blue', 'HNN': 'green', 'HNZ': 'red'}
    markers = {'HNE': 's', 'HNN': '^', 'HNZ': 'o'}  # cuadrado, triángulo, círculo

    # Configuración del gráfico
    plt.figure(figsize=(13, 9))

    # Crear el gráfico de dispersión
    scatter = sns.scatterplot(data=df, 
                              x='Start_Time_rounded_str', 
                              y='AP', 
                              hue='Channel', 
                              palette=palette, 
                              style='Channel', 
                              markers=markers,
                              s=150
                              )
    
    # Ajustar los ticks del eje X para que solo muestren unos 4 valores
    unique_x_values = df['Start_Time_rounded_str'].unique().tolist()
    num_ticks = len(unique_x_values)
    tick_positions = [0, int(num_ticks * 0.25), int(num_ticks * 0.5), int(num_ticks * 0.75), num_ticks - 1]
    tick_labels = [unique_x_values[pos] for pos in tick_positions]

    ax = plt.gca()
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=0, fontsize=8)

    # Mantener todos los puntos en el gráfico aunque no haya etiquetas en el eje X
    ax.xaxis.set_tick_params(length=0)

    # Añadir la cuadrícula
    plt.grid(True, which='both', linestyle='--', linewidth=1.5, alpha=0.8)
        
    df['Start_Time_title'] = df['Start_Time'] - timedelta(hours=utc_time)
    min_date = df['Start_Time_title'].min().strftime('%Y-%m-%d %H:%M:%S')
    fecha_formateada = datetime.strptime(min_date, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")

    # Título y etiquetas de los ejes
    plt.title(f'Peak particle acceleration trend - {fecha_formateada}', fontsize=18, fontweight='bold')
    plt.xlabel(f'Time UTC-{utc_time}', fontsize=16, fontweight='bold')
    plt.ylabel('Acceleration (g)', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(title='Channel',fontsize=15, title_fontproperties={'weight': 'bold', 'size': 15})

    plt.subplots_adjust(bottom=0.25)

    # Guardar gráfico como imagen
    ap_time_image_path = os.path.join(save_path, 'graph_AP_vs_Time.png')
    plt.savefig(ap_time_image_path)
    plt.close()

    # Imprimir información de los puntos utilizados
    # for i, row in df.iterrows():
    #     print(f"Plotted - VP: {row['VP']} End_Time_rounded: {row['End_Time_rounded']} Channel: {row['Channel']}")

    return ap_time_image_path
#*****************************************************************************



#*****************************************************************************
def plot_dp_vs_time(df, save_path, utc_time):
    # Convertir 'Start_Time' a cadena para usar en la gráfica
    df['Start_Time_str'] = df['Start_Time'].astype(str)

    # Convertir 'Start_Time' a objeto datetime manejando ambos formatos
    df['Start_Time'] = pd.to_datetime(df['Start_Time_str'], errors='coerce', format='%Y-%m-%dT%H:%M:%S.%f')
    df['Start_Time'] = df['Start_Time'].fillna(pd.to_datetime(df['Start_Time_str'], errors='coerce', format='%Y-%m-%dT%H:%M:%S'))

    # Redondear los tiempos a segundos para evitar diferencias mínimas
    df['Start_Time_rounded'] = df['Start_Time'] - timedelta(hours=utc_time)
    df['Start_Time_rounded_str'] = df['Start_Time_rounded'].dt.strftime('%H:%M:%S') 

    # Asegurarnos de que los tiempos estén ordenados correctamente
    df = df.sort_values(by=['Start_Time_rounded', 'Channel'])

    # Definir los colores específicos para cada canal
    palette = {'HNE': 'blue', 'HNN': 'green', 'HNZ': 'red'}
    markers = {'HNE': 's', 'HNN': '^', 'HNZ': 'o'}  # cuadrado, triángulo, círculo

    # Configuración del gráfico
    plt.figure(figsize=(13, 9))

    # Crear el gráfico de dispersión
    scatter = sns.scatterplot(data=df, 
                              x='Start_Time_rounded_str', 
                              y='DP', 
                              hue='Channel', 
                              palette=palette, 
                              style='Channel', 
                              markers=markers,
                              s=150
                              )

    # Ajustar los ticks del eje X para que solo muestren unos 4 valores
    unique_x_values = df['Start_Time_rounded_str'].unique().tolist()
    num_ticks = len(unique_x_values)
    tick_positions = [0, int(num_ticks * 0.25), int(num_ticks * 0.5), int(num_ticks * 0.75), num_ticks - 1]
    tick_labels = [unique_x_values[pos] for pos in tick_positions]

    ax = plt.gca()
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=0, fontsize=8)

    # Mantener todos los puntos en el gráfico aunque no haya etiquetas en el eje X
    ax.xaxis.set_tick_params(length=0)

    # Añadir la cuadrícula
    plt.grid(True, which='both', linestyle='--', linewidth=1.5, alpha=0.8)
        
    df['Start_Time_title'] = df['Start_Time'] - timedelta(hours=utc_time)
    min_date = df['Start_Time_title'].min().strftime('%Y-%m-%d %H:%M:%S')
    fecha_formateada = datetime.strptime(min_date, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")

    # Título y etiquetas de los ejes
    plt.title(f'Peak particle displacement trend - {fecha_formateada}', fontsize=18, fontweight='bold')
    plt.xlabel(f'Time UTC-{utc_time}', fontsize=16, fontweight='bold')
    plt.ylabel('Displacement (mm)', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(title='Channel',fontsize=15, title_fontproperties={'weight': 'bold', 'size': 15})

    plt.subplots_adjust(bottom=0.25)

    # Guardar gráfico como imagen
    dp_time_image_path = os.path.join(save_path, 'graph_DP_vs_Time.png')
    plt.savefig(dp_time_image_path)
    plt.close()

    # Imprimir información de los puntos utilizados
    # for i, row in df.iterrows():
    #     print(f"Plotted - VP: {row['VP']} End_Time_rounded: {row['End_Time_rounded']} Channel: {row['Channel']}")

    return dp_time_image_path
#*****************************************************************************


#*****************************************************************************

def resample_traces(traces, target_sampling_rate):
    for trace in traces:
        if trace.stats.sampling_rate != target_sampling_rate:
            trace.resample(target_sampling_rate)
        trace.data = trace.data.astype('float32')
    return traces

def adjust_color_intensity(color, factor):
    return tuple(min(1, c * factor + (1 - factor)) for c in color)

def signals_generator(selected_files,flag_trace, save_path, flag, utc_time, last_time=10):
    print("Generando gráficas de tiempos y frecuencias...")
    if flag_trace=="todo":
        st1 = obspy.Stream()
        valid_files = []
        for file in selected_files:
            try:
                tr = obspy.read(file)
                st1 += tr
                valid_files.append(file)
            except Exception as e:
                print(f"Error leyendo el archivo {file}: {e}")

        if len(st1) == 0:
            print("No se encontraron señales en los archivos seleccionados.")
            return
    elif flag_trace=="eventos":
        st1 = obspy.Stream(traces=list(selected_files.values()))
    # Solo se procesan los canales HNE, HNN y HNZ
    ordered_channels = ['HNE', 'HNN', 'HNZ']
    concatenated_traces = obspy.Stream()
    target_rate = 250  

    for channel in ordered_channels:
        traces = st1.select(channel=channel)
        if len(traces) > 0:
            traces = resample_traces(traces, target_rate)
            try:
                concatenated_trace = traces.merge(fill_value='interpolate')[0]
                concatenated_trace.detrend('linear')  # Quitar la media de cada traza
                concatenated_traces.append(concatenated_trace)
            except Exception as e:
                print(f"Error concatenando trazas para el canal {channel}: {e}")
                continue

    if len(concatenated_traces) == 0:
        print("No se pudieron concatenar las trazas debido a errores o tasas de muestreo incompatibles.")
        return

    # Definir colores y etiquetas para cada canal
    colors = {
        'HNE': (0, 0, 1),      # Azul
        'HNN': (0, 1, 0),      # Verde
        'HNZ': (1, 0, 0)       # Rojo
    }
    labels = {
        'HNE': 'HNE',
        'HNN': 'HNN',
        'HNZ': 'HNZ'
    }
    
    # === Graficar Aceleración (usando Raw2grav) ===
    # Convertimos cada traza cruda a aceleración en m/s² mediante Raw2grav y luego a mm/s² (multiplicando por 1000)
    fig, ax = plt.subplots(len(concatenated_traces), 1, figsize=(12, 8), sharex=True)
    if len(concatenated_traces) == 1:
        ax = [ax]
    
    max_y_value = 0
    acc_traces = []  # Almacenar las trazas de aceleración procesadas
    fft_paths = {}

    for trace in concatenated_traces:
        acc_trace = Raw2grav(trace)  # Convierte de cuentas a m/s²
        acc_trace.data = acc_trace.data * 1000  # Convertir a mm/s²
        acc_traces.append(acc_trace)
        local_max = np.abs(acc_trace.data).max()
        if local_max > max_y_value:
            max_y_value = local_max
    max_y_value = max_y_value+(max_y_value*0.5)


    max_y_value_g=0
    acc_traces_g = []
    for trace in concatenated_traces:
        acc_trace_g = Raw2grav_g(trace)  # Convierte de cuentas a m/s²
        acc_traces_g.append(acc_trace_g)
        local_max = np.abs(acc_trace_g.data).max()
        if local_max > max_y_value_g:
            max_y_value_g = local_max
    max_y_value_g = max_y_value_g+(max_y_value_g*0.5)

    for i, acc_trace_g in enumerate(acc_traces_g):

        # Obtener los tiempos en UTC en formato datetime
        times_list = acc_trace_g.times("utcdatetime")
        # Convertir a un array numpy de enteros (segundos desde el epoch)
        times_np = np.fromiter(
                    (t.timestamp() if callable(getattr(t, "timestamp", None)) else t for t in times_list),
                    dtype=np.float64
                    )

        # Restar 3 horas (3*3600 segundos)
        times_np = times_np - utc_time * 3600
        # Convertir de epoch a la representación numérica de Matplotlib
        times_adjusted = times_np / 86400.0 + 719163.0

        n=min(len(times_adjusted),len(acc_trace_g.data))
        times_adjusted=times_adjusted[:n]
        data_plot = acc_trace_g.data[:n]

        channel = acc_trace_g.stats.channel
        color = colors.get(channel, 'k')
        label = labels.get(channel, f'Channel {channel}')
        ax[i].plot(times_adjusted, data_plot, label=label, color=color)
        ax[i].legend(loc='upper right', fontsize=15, title_fontsize=15)
        ax[i].set_ylabel('Acceleration (g)', fontsize=12,fontweight='bold')
        ax[i].set_ylim([-max_y_value_g, max_y_value_g])
        ax[i].tick_params(axis='both', which='major', labelsize=15)
        ax[i].grid(True, which='both', linestyle='--', linewidth=1.5, alpha=0.8)
        ax[i].xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=8))
        ax[i].xaxis_date()
        ax[i].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))

    ax[-1].set_xlabel(f'Time (UTC-{utc_time})', fontsize=15,fontweight='bold')
    fig.suptitle('Signals over time (acceleration)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    aceleration_path = os.path.join(save_path, 'signal_acceleration_gravity.png')
    plt.savefig(aceleration_path)
    plt.close()
    print(f"Gráfica de aceleración guardada en: {aceleration_path}")

    # === Calcular Velocidad y Desplazamiento aplicando el proceso sugerido ===
    # Para cada traza: convertir a aceleración con Raw2grav, aplicar cosine taper, filtrar y luego integrar.
    velocity_traces = obspy.Stream()
    displacement_traces = obspy.Stream()
    
    for trace in concatenated_traces:
        # Convertir la traza cruda a aceleración en m/s²
        acc_trace = Raw2grav(trace)
        # Aplicar cosine taper (por ejemplo, con porcentaje 0.1)        
        acc_trace.data = acc_trace.data * cosine_taper(len(acc_trace.data), 0.1)
        # Crear una copia y aplicar filtro bandpass
        acc_filt = acc_trace.copy()
        acc_filt.filter('bandpass', freqmin=0.2, freqmax=50)
        
        # Velocidad: primera integral para obtener m/s, luego convertir a mm/s
        vel_trace = acc_filt.copy()
        vel_trace.integrate(method='cumtrapz')
        vel_trace.data = vel_trace.data * 1000  # de m/s a mm/s
        velocity_traces.append(vel_trace)
        
        # Desplazamiento: segunda integral para obtener m, luego convertir a mm
        disp_trace = vel_trace.copy()
        disp_trace.integrate(method='cumtrapz')
        disp_trace.data = disp_trace.data  # de m a mm
        displacement_traces.append(disp_trace)

    # === Graficar Velocidad ===
    fig_v, ax_v = plt.subplots(len(velocity_traces), 1, figsize=(12, 8), sharex=True)
    if len(velocity_traces) == 1:
        ax_v = [ax_v]
    max_v_value = max(np.abs(trace.data).max() for trace in velocity_traces)
    max_v_value= max_v_value+(max_v_value*0.5)
    for i, trace in enumerate(velocity_traces):
        
        # Obtener los tiempos en UTC en formato datetime
        times_list = trace.times("utcdatetime")
        # Convertir a un array numpy de enteros (segundos desde el epoch)
        times_np = np.fromiter(
                    (t.timestamp() if callable(getattr(t, "timestamp", None)) else t for t in times_list),
                    dtype=np.float64
                    )

        # Restar 3 horas (3*3600 segundos)
        times_np = times_np - utc_time * 3600
        # Convertir de epoch a la representación numérica de Matplotlib
        times_adjusted = times_np / 86400.0 + 719163.0

        n = min(len(times_adjusted), len(trace.data))
        times_adjusted = times_adjusted[:n]
        data_plot = trace.data[:n]

        channel = trace.stats.channel
        color = colors.get(channel, 'k')
        label = labels.get(channel, f'Channel {channel}')
        ax_v[i].plot(times_adjusted, data_plot, label=label, color=color)
        ax_v[i].legend(loc='upper right', fontsize=15, title_fontsize=15)
        ax_v[i].set_ylabel('Velocity (mm/s)', fontsize=12,fontweight='bold')
        ax_v[i].set_ylim([-max_v_value, max_v_value])
        ax_v[i].tick_params(axis='both',which='major',labelsize=15)
        ax_v[i].grid(True, which='both', linestyle='--', linewidth=1.5, alpha=0.8)
        ax_v[i].xaxis_date()
        ax_v[i].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    ax_v[-1].set_xlabel(f'Time (UTC-{utc_time})',fontsize=15, fontweight='bold')
    fig_v.suptitle('Signals over time (velocity)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    velocity_path = os.path.join(save_path, 'signal_velocity.png')
    plt.savefig(velocity_path)
    plt.close()
    print(f"Gráfica de velocidad guardada en: {velocity_path}")

    # === Graficar Desplazamiento ===
    fig_d, ax_d = plt.subplots(len(displacement_traces), 1, figsize=(12, 8), sharex=True)
    if len(displacement_traces) == 1:
        ax_d = [ax_d]
    max_d_value = max(np.abs(trace.data).max() for trace in displacement_traces)
    max_d_value = max_d_value + (max_d_value*0.5)
    for i, trace in enumerate(displacement_traces):

        # Obtener los tiempos en UTC en formato datetime
        times_list = trace.times("utcdatetime")
        # Convertir a un array numpy de enteros (segundos desde el epoch)
        times_np = np.fromiter(
                    (t.timestamp() if callable(getattr(t, "timestamp", None)) else t for t in times_list),
                    dtype=np.float64
                    )

        # Restar 3 horas (3*3600 segundos)
        times_np = times_np - utc_time * 3600
        # Convertir de epoch a la representación numérica de Matplotlib
        times_adjusted = times_np / 86400.0 + 719163.0

        n = min(len(times_adjusted), len(trace.data))
        times_adjusted = times_adjusted[:n]
        data_plot = trace.data[:n]

        channel = trace.stats.channel
        color = colors.get(channel, 'k')
        label = labels.get(channel, f'Channel {channel}')
        ax_d[i].plot(times_adjusted, data_plot, label=label, color=color)
        ax_d[i].legend(loc='upper right', fontsize=15, title_fontsize=15)
        ax_d[i].set_ylabel('Displacement (mm)',fontsize=12,fontweight='bold')
        ax_d[i].set_ylim([-max_d_value, max_d_value])
        ax_d[i].tick_params(axis='both',which='major',labelsize=15)
        ax_d[i].grid(True, which='both', linestyle='--', linewidth=1.5,alpha=0.8)
        ax_d[i].xaxis_date()
        ax_d[i].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))

    ax_d[-1].set_xlabel(f'Time (UTC-{utc_time})',fontsize=15,fontweight='bold')
    fig_d.suptitle('Signals over time (displacement)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    displacement_path = os.path.join(save_path, 'signal_displacement.png')
    plt.savefig(displacement_path)
    plt.close()
    print(f"Gráfica de desplazamiento guardada en: {displacement_path}")

    ################
    # Supongamos que tienes definida la función highpass_filter y Raw2grav, y diccionarios de colores y etiquetas:
    def highpass_filter(data, cutoff, fs, order=4):
        nyquist = 0.5 * fs
        normal_cutoff = cutoff / nyquist
        b, a = butter(order, normal_cutoff, btype='high', analog=False)
        filtered_data = filtfilt(b, a, data)
        return filtered_data

    # Por cada traza en concatenated_traces se genera una gráfica FFT
    for trace in concatenated_traces:

        sampling_rate = 250
        # Convertir la traza a aceleración en m/s² usando Raw2grav
        acc_trace = Raw2grav(trace)
        # Si es necesario, convierte a mm/s² multiplicando por 1000 (opcional)
        # acc_trace.data = acc_trace.data * 1000

        # Aplicar filtro highpass a los datos de aceleración
        filtered_data = highpass_filter(acc_trace.data, 0.1, sampling_rate)

        # Calcular la transformada de Fourier
        fft_data = np.fft.fft(filtered_data)
        fft_freq = np.fft.fftfreq(len(fft_data), 1/sampling_rate)

        # Seleccionar sólo la parte positiva del espectro
        pos_mask = fft_freq > 0
        fft_freq_pos = fft_freq[pos_mask]
        fft_data_pos = np.abs(fft_data[pos_mask])

        # Crear la gráfica de Fourier para el canal actual
        plt.figure(figsize=(12, 8))
        channel = acc_trace.stats.channel
        color = colors.get(channel, 'k')
        plt.plot(fft_freq_pos, fft_data_pos, label=f"{channel}", color=color)
        plt.xlabel("Frequency (Hz)", fontsize=15,fontweight='bold')
        plt.ylabel("Amplitude", fontsize=15,fontweight='bold')
        plt.title(f"Fourier Transform - Channel {channel}", fontsize=16, fontweight="bold")
        plt.legend(fontsize=15, title_fontsize=15)
        plt.grid(True)
        plt.tick_params(axis='both',which='major',labelsize=15)

        # Guardar la gráfica con un nombre único para cada canal
        fft_path = os.path.join(save_path, f"fft_{channel}.png")
        plt.savefig(fft_path)
        plt.close()
        print(f"Gráfica de transformadas para {channel} guardada en: {fft_path}")
        fft_paths[channel] = fft_path

    return aceleration_path, velocity_path, displacement_path, fft_paths



#*****************************************************************************



#*****************************************************************************
class PDF(FPDF):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.custom_header = "ANALYSIS REPORT"  # Valor por defecto

    def header(self):
        self.set_font("Arial", "B", 16)
        # Usa el valor de custom_header (si lo modificas, se mostrará el nuevo título)
        self.cell(0, 10, self.custom_header, 0, 1, "C")

    def footer(self):
        self.set_y(-15)
        self.set_font("Arial", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}", 0, 0, "C")

#*****************************************************************************
        

#*****************************************************************************
#Función para generar pdf
#@autor:oszu
#*****************************************************************************
def create_pdf(df, data_project, save_path, logo_path,acceleration_path, 
               velocity_path, displacement_path, fft_path, last_time, utc_time, report_name): 

    current_function_name = inspect.currentframe().f_code.co_name
    print(f"\n     INICIO FUNCIÓN {current_function_name}")
    size_general = 10
    df_a=df.copy()
    df_d=df.copy()
    df_fft=df.copy()
    # Crear los gráficos primero
    print("Creando las gráficas...")

    abnt_bnr_image_path = plot_dataframe_ABNT_BNR(df, save_path)
    cetesb_image_path = plot_dataframe_CETESB(df, save_path)    
    vp_time_image_path = plot_vp_vs_time(df, save_path, utc_time)
    ap_time_image_path= plot_ap_vs_time(df_a,save_path, utc_time)
    dp_time_image_path = plot_dp_vs_time(df_d,save_path, utc_time)

    pdf = PDF()
    pdf.add_page()
    pdf.set_auto_page_break(auto=True, margin=15)

    # Logo en la esquina superior derecha
    pdf.image(logo_path, x=170, y=10, w=25)

    # Medidas aproximadas
    left_x = 10       # Columna izquierda
    right_x = 105     # Columna derecha
    start_y = 30      # Margen superior inicial
    column_width = 95 # Ancho de la segunda columna (210 mm A4 - 2 márgenes - left_x)

    # === Columna Izquierda ===
    pdf.set_xy(left_x, start_y)
    pdf.set_font("Arial", size=size_general, style='B')

    # Convertir a datetime para evitar errores de string
    df['Start_Time'] = pd.to_datetime(df['Start_Time'], errors='coerce', format='%Y-%m-%dT%H:%M:%S.%f')
    df['Start_Time'] = df['Start_Time'].fillna(
        pd.to_datetime(df['Start_Time'], errors='coerce', format='%Y-%m-%dT%H:%M:%S')
    )
    df['End_Time'] = pd.to_datetime(df['End_Time'], errors='coerce', format='%Y-%m-%dT%H:%M:%S.%f')
    df['End_Time'] = df['End_Time'].fillna(
        pd.to_datetime(df['End_Time'], errors='coerce', format='%Y-%m-%dT%H:%M:%S')
    )

    # Fechas y notas
    min_date = (df['Start_Time'].min() - timedelta(hours=utc_time)).strftime('%Y-%m-%d %H:%M:%S')
    try:
        max_date = (df['End_Time'].max() - timedelta(hours=utc_time)).strftime('%Y-%m-%d %H:%M:%S')
    except Exception as e:
        max_date = (df['End_Time'].max() - timedelta(hours=utc_time))
        print(f"Error converting max_date: {e}")

    channels = ['HNE', 'HNN', 'HNZ']
    notes = []
    for channel in channels:
        df_channel = df[df['Channel'] == channel].reset_index(drop=True)
        if not df_channel.empty:
            max_vp_index = df_channel['VP'].idxmax()
            max_vp = df_channel['VP'].max()
            max_fp = df_channel['FP'].iloc[max_vp_index]
            max_end_time = (pd.to_datetime(df_channel['Start_Time'].iloc[max_vp_index]) - timedelta(hours=utc_time)).strftime('%Y-%m-%d %H:%M:%S')
            notes.append(f"Max velocity {round(max_vp, 2)} mm/s at {max_fp} Hz in {channel} - {max_end_time}")
        else:
            notes.append(f"No data for {channel}")

    info_fields = [
        ("Date:", f"{datetime.now().strftime('%Y-%m-%d')} UTC-5"),
        ("Project:", data_project[0]),        
        ("Data interval:", f"{(df['Start_Time'].min() - timedelta(hours=utc_time)).strftime('%Y-%m-%d %H:%M:%S')} to {(df['Start_Time'].max() - timedelta(hours=utc_time)+timedelta(minutes=last_time)).strftime('%Y-%m-%d %H:%M:%S')} UTC-{utc_time}"),
        ("Client:", data_project[1]),
        ("Location:", data_project[2]),
        ("Notes:","")
    ]

    # Imprimir info en la columna izquierda
    for field_name, field_value in info_fields:
        pdf.set_font("Arial", size=size_general, style='B')
        pdf.cell(25, 5, txt=field_name, ln=0, align='L')
        pdf.set_font("Arial", size=size_general)
        pdf.cell(0, 5, txt=str(field_value), ln=1, align='L')

    # === Columna Derecha ===
    pdf.set_xy(right_x, start_y+10)  # Comienza en la misma Y, pero X=105
    pdf.set_font("Arial", size=size_general, style='B')
    # Especificamos un ancho fijo (column_width) para mantenernos en la columna derecha
    pdf.cell(column_width, 5, txt="Continuos Vibration Data", ln=1, align='C')
    pdf.ln(2)

    # Datos para la tabla
    table_data = {}
    for channel in ['HNE', 'HNZ', 'HNN']:
        df_channel = df[df['Channel'] == channel].reset_index(drop=True)
        if not df_channel.empty:
            idx = df_channel['VP'].idxmax()
            ppv = round(df_channel['VP'].max(), 3)
            ppa = round(df_channel['AP'].max(), 3)
            ppd = round(df_channel['DP'].max(), 3)
            freq = df_channel['FP'].loc[idx]
            table_data[channel] = {"PPV": ppv, "PPA": ppa, "PPD": ppd, "FREQ": freq}
        else:
            table_data[channel] = {"PPV": "", "PPA": "", "PPD": "", "FREQ": ""}

    # Definir anchos de celdas dentro de la columna derecha
    col1_width = 30
    col_width = (column_width - col1_width) / 3

    # --- Encabezado de la tabla ---
    pdf.set_x(right_x)  # Aseguramos la posición X en la segunda columna
    pdf.set_font("Arial", size=size_general, style='B')
    pdf.cell(col1_width, 5, txt="", border=0, align='C')
    for header in ['HNE', 'HNN', 'HNZ']:
        pdf.cell(col_width, 5, txt=header, border=0, align='C')
    pdf.ln(6)

    # --- Filas de la tabla ---
    pdf.set_font("Arial", size=size_general)
    parameters = [
        ("PPA (g)", "PPA"),
        ("PPV (mm/s)", "PPV"),        
        ("PPD (mm)", "PPD"),
        ("FREQ (Hz)", "FREQ")
    ]
    for label, key in parameters:
        pdf.set_x(right_x)  # Mover a la segunda columna al inicio de cada fila
        pdf.set_font("Arial", size=size_general, style='B')
        pdf.cell(col1_width, 4, txt=label, border=0, align='L')
        for channel in ['HNE', 'HNN', 'HNZ']:
            value = table_data[channel][key]
            pdf.set_font("Arial", size=size_general)
            pdf.cell(col_width, 4, txt=str(value), border=0, align='C')
        pdf.ln(5)

    # --- Datos adicionales debajo de la tabla ---
    pdf.ln(4)
    # Definir los campos extra en una lista similar a info_fields
    if last_time>=1:
        extra_info_fields = [
            ("Sample Rate:", "250 Hz"),
            ("Window Size:", f"{last_time} min"),
            ("Serial:", "ME 379")
        ]
    else:
        extra_info_fields = [
            ("Sample Rate:", "250 Hz"),
            ("Window Size:", f"{last_time*60} seg"),
            ("Serial:", "ME 379")
        ]

# Imprimir los campos extra en la columna derecha
    for field_name, field_value in extra_info_fields:
        pdf.set_x(right_x)  # Asegura que cada línea comienza en la misma posición (columna derecha)
        pdf.set_font("Arial", size=size_general, style='B')
        # Ajusta el ancho (por ejemplo, 30 mm para la parte en negrilla)
        pdf.cell(30, 5, txt=field_name, ln=0, align='L')
        pdf.set_font("Arial", size=size_general, style='')
        # El ancho restante de la columna para el valor (column_width - 30)
        pdf.cell(column_width - 30, 5, txt=str(field_value), ln=1, align='L')

    # Salto de línea antes de imágenes
    pdf.ln(10)

    # Páginas nuevas para cada gráfica (puedes ajustar a tu gusto)
    pdf.image(abnt_bnr_image_path, x=32, y=91, w=150)
    pdf.image(cetesb_image_path, x=32, y=188, w=150)

    pdf.custom_header = "PARTICLE VELOCITY ANALYSIS"  # Cambia el título
    pdf.add_page()    
    pdf.image(logo_path, x=170, y=10, w=25)
    pdf.image(velocity_path, x=12, y=38, w=180)
    pdf.image(vp_time_image_path, x=12, y=170, w=180)

    pdf.custom_header = "PARTICLE ACCELERATION ANALYSIS"
    pdf.add_page()    
    pdf.image(logo_path, x=170, y=10, w=25)
    pdf.image(acceleration_path, x=12, y=38, w=180)
    pdf.image(ap_time_image_path, x=12, y=170, w=180)

    pdf.custom_header = "PARTICLE DISPLACEMENT ANALYSIS"
    pdf.add_page()    
    pdf.image(logo_path, x=170, y=10, w=25)
    pdf.image(displacement_path, x=12, y=38, w=180)
    pdf.image(dp_time_image_path, x=12, y=170, w=180)

    pdf.custom_header = "FAST FOURIER TRANSFORM (FFT)"
    pdf.add_page()    
    pdf.image(logo_path, x=170, y=10, w=25)
    pdf.image(fft_path['HNE'], x=12, y=38, w=180)
    pdf.image(fft_path['HNN'], x=12, y=170, w=180)

    pdf.custom_header = "FAST FOURIER TRANSFORM (FFT)"
    pdf.add_page()    
    pdf.image(logo_path, x=170, y=10, w=25)
    pdf.image(fft_path['HNZ'], x=12, y=38, w=180)

    # Guardar PDF
    fecha = datetime.strptime(min_date, "%Y-%m-%d %H:%M:%S")
    fecha_formateada = fecha.strftime("%Y-%m-%d")
    pdf_path = os.path.join(save_path, f'{report_name}')
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

import re

def get_hour_from_filename(filename):
    """
    Extrae el valor de la hora del nombre de archivo.
    Se asume que el patrón es algo similar a:
    ...D.2025.059.00.miniseed  (donde "00" es la hora)
    """
    match = re.search(r'\.D\.\d{4}\.\d{3}\.(\d{2})\.miniseed', filename)
    if match:
        return int(match.group(1))
    return None



#*****************************************************************************
# Programa principal
# Encuentra los valores de interés (velocidades, aceleraciones, desplazamientos,
# entre otros) del archivo .mseed encontrado en la carpeta seleccionada y los 
# envía al servidor thingsboard.
# Fecha: 28/11/2024
# autor @oszu
#*****************************************************************************

from flask import Flask, render_template, request, jsonify
import os
from datetime import datetime, timedelta

app = Flask(__name__)

# Asumimos que estas funciones ya están definidas/importadas:
# filtrar_archivos_por_dia_y_canal, get_hour_from_filename, detect_event,
# adquire_data, signals_generator, create_pdf

@app.route('/pagina2', methods=['GET', 'POST'])
def pagina2():
    message = None
    if request.method == 'POST':
        try:
            # Validar y obtener los valores ingresados por el usuario
            hours_input = request.form.get('hours_input')
            julian_days = request.form.get('julian_days')
            if not hours_input or not julian_days:
                raise ValueError("Todos los campos son obligatorios.")
            seg_option = int(hours_input)
            if seg_option not in [1, 12, 24]:
                raise ValueError("El valor de 'Analizar por' debe ser 1, 12 o 24.")
            try:
                dias_deseados = [int(day.strip()) for day in julian_days.split(",")]
            except ValueError:
                raise ValueError("El campo 'Días julianos' debe contener números separados por comas.")

            # Obtener los archivos enviados desde el input oculto
            # Es importante que en el input de la plantilla HTML se agregue: name="folderInput"
            folder_files = request.files.getlist('folderInput')
            if not folder_files or len(folder_files) == 0:
                raise ValueError("No se seleccionaron archivos.")

            # Aquí podrías, por ejemplo, guardar temporalmente los archivos o procesarlos directamente.
            # En este ejemplo, asumiremos que 'folder_files' es una lista de objetos FileStorage
            # y que nuestras funciones trabajan con estos objetos.
            archivos_seleccionados = folder_files

            # Variables fijas para el reporte
            logo_path = r'C:\Users\ssi_s\OneDrive\Escritorio\Oscar\Codigo Reporte Jorge\Reports_code_SSI\LogoSSI.png'
            last_time_hour = 5   # ventana (minutos) para 1 hora
            last_time = 10       # ventana (minutos) para 12 y 24 horas
            utc_time = 3
            project = 'Vibração de Solo - Braskem'
            client = 'AGILE - BRASKEM'
            location = 'Maceió - AL'
            notes = ''
            data_project = [project, client, location, notes]

            # Procesar cada día juliano ingresado
            for dia in dias_deseados:
                # Filtrar archivos correspondientes al día juliano 'dia'
                files_dia = filtrar_archivos_por_dia_y_canal(archivos_seleccionados, dia)
                if not files_dia:
                    print(f"No se encontraron archivos para el día {dia}.")
                    continue

                year = 2025
                date_obj = datetime(year, 1, 1) + timedelta(days=dia - 1)
                date_str = date_obj.strftime("%b%d")
                if seg_option == 1:
                    # Análisis por tramos de 1 hora
                    last_time_used = last_time_hour
                    for hr in range(24):
                        # Para cada hora, filtrar los archivos
                        files_hour = [
                            f for f in files_dia
                            if (get_hour_from_filename(f.filename) is not None and get_hour_from_filename(f.filename) == hr)
                        ]
                        if files_hour:
                            print(f"Analizando día {dia} - hora {hr:02d}:00-{hr:02d}:59 UTC...")
                            selected_files = files_hour
                            # Ejecutar la función de detección de eventos (ajusta según tus necesidades)
                            detect_event(selected_files, os.path.dirname(files_hour[0].filename), utc_time, last_time_used, logo_path, data_project)
                            # Obtener datos y generar gráfica
                            df = adquire_data(selected_files, "todo", last_time_used)
                            acceleration_path, velocity_path, displacement_path, fft_path = signals_generator(
                                selected_files, "todo", os.path.dirname(files_hour[0].filename), 1, utc_time, last_time_used
                            )
                            # Ajustar la hora a UTC-3
                            if hr < 3:
                                report_date = date_obj - timedelta(days=1)
                                adjusted_hr = hr + 24 - 3
                            else:
                                report_date = date_obj
                                adjusted_hr = hr - 3
                            report_date_str = report_date.strftime("%b%d")
                            report_name = f"Report_{report_date_str}_{adjusted_hr:02d}.pdf"
                            create_pdf(
                                df, data_project, os.path.dirname(files_hour[0].filename),
                                logo_path, acceleration_path, velocity_path, displacement_path,
                                fft_path, last_time_used, utc_time, report_name=report_name
                            )
                        else:
                            print(f"No se encontraron archivos para el día {dia} en la hora {hr:02d}:00-{hr:02d}:59.")
                elif seg_option == 12:
                    # Análisis por tramos de 12 horas
                    files_morning = [
                        f for f in files_dia
                        if (get_hour_from_filename(f.filename) is not None and get_hour_from_filename(f.filename) < 12)
                    ]
                    files_afternoon = [
                        f for f in files_dia
                        if (get_hour_from_filename(f.filename) is not None and get_hour_from_filename(f.filename) >= 12)
                    ]
                    if files_morning:
                        print(f"Analizando día {dia} (00:00 - 11:59)...")
                        selected_files = files_morning
                        detect_event(selected_files, os.path.dirname(files_morning[0].filename), utc_time, logo_path, data_project)
                        df = adquire_data(selected_files, "todo", last_time)
                        acceleration_path, velocity_path, displacement_path, fft_path = signals_generator(
                            selected_files, "todo", os.path.dirname(files_morning[0].filename), 1, utc_time, last_time
                        )
                        morning_date_obj = date_obj - timedelta(days=1)
                        morning_date_str = morning_date_obj.strftime("%b%d")
                        create_pdf(
                            df, data_project, os.path.dirname(files_morning[0].filename),
                            logo_path, acceleration_path, velocity_path, displacement_path,
                            fft_path, last_time, utc_time, report_name=f"Report_{morning_date_str}_21-09.pdf"
                        )
                    else:
                        print(f"No se encontraron archivos para el tramo 00-11 en el día {dia}.")
                    if files_afternoon:
                        print(f"Analizando día {dia} (12:00 - 23:59)...")
                        selected_files = files_afternoon
                        detect_event(selected_files, os.path.dirname(files_afternoon[0].filename), utc_time, last_time, logo_path, data_project)
                        df = adquire_data(selected_files, "todo", last_time)
                        acceleration_path, velocity_path, displacement_path, fft_path = signals_generator(
                            selected_files, "todo", os.path.dirname(files_afternoon[0].filename), 1, utc_time, last_time
                        )
                        create_pdf(
                            df, data_project, os.path.dirname(files_afternoon[0].filename),
                            logo_path, acceleration_path, velocity_path, displacement_path,
                            fft_path, last_time, utc_time, report_name=f"Report_{date_str}_09-21.pdf"
                        )
                    else:
                        print(f"No se encontraron archivos para el tramo 12-23 en el día {dia}.")
                elif seg_option == 24:
                    # Análisis por 24 horas: usar todos los archivos del día
                    print(f"Analizando día {dia} (24 horas)...")
                    detect_event(files_dia, os.path.dirname(files_dia[0].filename), utc_time, logo_path, data_project)
                    selected_files = files_dia
                    df = adquire_data(selected_files, "todo", last_time)
                    acceleration_path, velocity_path, displacement_path, fft_path = signals_generator(
                        selected_files, "todo", os.path.dirname(files_dia[0].filename), 1, utc_time, last_time
                    )
                    create_pdf(
                        df, data_project, os.path.dirname(files_dia[0].filename),
                        logo_path, acceleration_path, velocity_path, displacement_path,
                        fft_path, last_time, utc_time, report_name=f"Report_{date_str}_24.pdf"
                    )
            message = "Reportes generados correctamente."
        except Exception as e:
            message = f"Error: {str(e)}"
            print(message)
    return render_template('pagina2.html', message=message)
