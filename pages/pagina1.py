from flask import Blueprint, render_template, request

from flask import Flask, render_template, request, redirect, jsonify
import matplotlib.pyplot as plt
import os
import subprocess
import paramiko
import io
import base64
import threading
import subprocess
import time
import inspect

from scipy.integrate import cumulative_trapezoid

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

from datetime import timedelta, datetime
import obspy
import numpy as np
from scipy.signal import butter, filtfilt

from flask_socketio import SocketIO, emit


def read_data(host, port, username, password, data_name):

    # Comando WinSCP para el archivo results
    remote_path_results = f'/mnt/mmcblk0p1/{data_name}-results'
    local_path_results = f'{data_name}-results'
    command_results = f'winscp.com /command "open scp://{username}:{password}@{host}:{port}" "get {remote_path_results} {local_path_results}" "exit"'

    attempts=0
    max_attempts=1
    # Intentar leer el archivo results infinitamente hasta que se genere
    while attempts<max_attempts:
        try:
            subprocess.run(command_results, shell=True, check=True)

            # Leer y almacenar el contenido del archivo results
            with open(local_path_results, 'r') as file:
                results = file.read()
            break  # Salir del bucle una vez que se ha leído el archivo
        except subprocess.CalledProcessError:
            attempts+=1
            print(f"No se encontró el archivo {data_name}-results. Reintentando en 5 segundos...")
            time.sleep(5)  # Esperar 5 segundos antes de intentar de nuevo

    if attempts == max_attempts:
        print(f"Error: No se pudo encontrar el archivo {data_name}-results después de {max_attempts} intentos.")
        return None, None
    
    # Comando WinSCP para el primer archivo
    remote_path = f'/mnt/mmcblk0p1/{data_name}'
    local_path = data_name
    command = f'winscp.com /command "open scp://{username}:{password}@{host}:{port}" "get {remote_path} {local_path}" "exit"'
    subprocess.run(command, shell=True)

    # Leer y almacenar el contenido del primer archivo
    with open(local_path, 'r') as file:
        content = file.read()

    return content, results
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------
def execute_calibration_command(host, port, username, password, data_name, time_selected):
    time_selected=None
    command = f'./calibration {data_name} {time_selected}'
    
    try:
        # Crear un cliente SSH
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        # Conectarse al dispositivo
        client.connect(hostname=host, port=port, username=username, password=password)
        print("\nConectado al dispositivo.")

        # Ejecutar el comando
        stdin, stdout, stderr = client.exec_command(command)
        
        # Leer la salida del comando
        output = stdout.read().decode()
        error = stderr.read().decode()
        
        # Cerrar la conexión
        client.close()
        
        if output:
            print("Salida del comando:")
            print(output)
        if error:
            print("Error del comando:")
            print(error)

    except Exception as e:
        print(f"Ocurrió un error: {e}")
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------
def parse_results(data_results):
    results_dict = {'EHE': '', 'EHN': '', 'EHZ': ''}
    lines = data_results.strip().split('\n')
    current_channel = None

    for line in lines:
        if 'CHANNEL' in line:
            current_channel = line.split()[1]
        elif current_channel:
            results_dict[current_channel] += line + '\n'
    return results_dict
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
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
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------

def calcular_frequency(trace, sampling_rate=1000):

    N = len(trace)
    yf = fft(trace)
    xf = fftfreq(N, 1 / sampling_rate)
    
    yf = np.abs(yf) / N
    
    pos_indices = np.where(xf > 0)
    xf = xf[pos_indices]
    yf = yf[pos_indices]
    
    idx = np.argmax(yf)
    main_frequency = xf[idx]
    main_magnitude = yf[idx]
    
    return main_frequency, main_magnitude
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------
def verificar_valores(tr,traces,tr_name):
    #tri=Raw2grav(tri)
    tri=tr

    min_value=min(tri)
    max_value=max(tri)
    p2p=np.abs(np.max(tri)-np.min(tri))
    mean_value=np.mean(tri)

    cuadrados=[]
    for i in tri:        
        multiplicacion=i*i
        cuadrados.append(multiplicacion)

    rms = np.sqrt(np.mean(cuadrados))
    std_dev = np.std(tri)  
    Frequency, Magnitude = calcular_frequency(tr)
    
    print(f"\nValores calculados en python para {tr_name}")
    
    print(f"For channel {tr_name} frecuency:{Frequency:.6f},Magnitude {Magnitude:.6f}")

    info_proc_str = (f"Frequency: {Frequency:.2f} Hz\n"
                     f"Magnitude: {Magnitude:.6f}\n"
                     f"Min: {min_value:.2f}\n"
                     f"Max: {max_value:.2f}\n"
                     f"Peak to peak: {p2p:.2f}\n"
                     f"Average: {mean_value:.2f}\n"
                     f"Standard deviation: {std_dev:.2f}\n"
                     f"RMS: {rms:.2f}\n")

    print(info_proc_str)
    return info_proc_str
    
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------
def transformar_ejes_y(EHE, EHN, EHZ, axis_y):
    if axis_y == 'Gravities':
        EHE = [e / 520000 for e in EHE]
        EHN = [e / 520000 for e in EHN]
        EHZ = [e / 520000 for e in EHZ]
    elif axis_y == 'Meters / s²':
        EHE = [e / 53007 for e in EHE]
        EHN = [e / 53007 for e in EHN]
        EHZ = [e / 53007 for e in EHZ]
    elif axis_y == 'Milimeters / s²':
        EHE = [e / 53.007 for e in EHE]
        EHN = [e / 53.007 for e in EHN]
        EHZ = [e / 53.007 for e in EHZ]
    elif axis_y == 'Centimeters / s²':
        EHE = [e / 530.07 for e in EHE]
        EHN = [e / 530.07 for e in EHN]
        EHZ = [e / 530.07 for e in EHZ]
    return EHE, EHN, EHZ
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------
def show_test(data_test, data_results, data_name, axis_y):
    axis_x = 'samples'

    traces=[]

    print("\nContenido de data_results:")
    print(data_results)

    lines = data_test.strip().split('\n')
    EHE, EHN, EHZ = [], [], []

    for line in lines[1:]:
        try:
            if line.strip():
                ehe, ehn, ehz = map(int, line.split(','))
                EHE.append(ehe)
                EHN.append(ehn)
                EHZ.append(ehz)
        except ValueError as e:
            print(f"Error al procesar la línea: {line}. Error: {e}")
    
    EHE, EHN, EHZ = transformar_ejes_y(EHE, EHN, EHZ, axis_y)

    traces={

        'EHE':EHE,
        'EHN':EHN,
        'EHZ':EHZ
    }
    for tr_name, tr in traces.items():
        verificar_valores(tr,traces,tr_name)

    

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8))
    fig.suptitle(f"File read: {data_name}")

    ax1.plot(EHE, color='mediumblue')
    ax1.set_title("Channel EHE")
    ax1.set_xlabel(axis_x, fontsize=14)
    ax1.set_ylabel(axis_y, fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=12)

    ax2.plot(EHN, color='lime')
    ax2.set_title("Channel EHN")
    ax2.set_xlabel(axis_x, fontsize=14)
    ax2.set_ylabel(axis_y, fontsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=12)

    ax3.plot(EHZ, color='red')
    ax3.set_title("Channel EHZ")
    ax3.set_xlabel(axis_x, fontsize=14)
    ax3.set_ylabel(axis_y, fontsize=14)
    ax3.tick_params(axis='both', which='major', labelsize=12)

    results_dict_python={}

    for tr_name, tr in traces.items():    
        results_dict_python[tr_name] = verificar_valores(tr,traces,tr_name)

    results_dict_James = parse_results(data_results)
    
    results_dict=results_dict_python

    print(results_dict)
    ax1.text(1.02, 0.5, results_dict['EHE'], transform=ax1.transAxes, verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5),fontsize=14)
    ax2.text(1.02, 0.5, results_dict['EHN'], transform=ax2.transAxes, verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5),fontsize=14)
    ax3.text(1.02, 0.5, results_dict['EHZ'], transform=ax3.transAxes, verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5),fontsize=14)

    plt.tight_layout()

    # Convertir la figura en una imagen en base64 para mostrar en la página web
    img = io.BytesIO() 
    fig.savefig(img, format='png')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()

    return plot_url



pagina1_bp = Blueprint('pagina1_bp', __name__, template_folder='templates', url_prefix='/pagina1')

host = '192.168.0.1'
port = 22
username = 'root'
password = 'admin'

@pagina1_bp.route('/', methods=['GET', 'POST'])
def pagina1():
    plot_url = None
    message = None

    if request.method == 'POST':
        if 'data_name' in request.form:
            data_name = request.form['data_name']
            axis_y = request.form['y_axis']
            data_test, data_results = read_data(host, port, username, password, data_name)
            if data_test is None:
                print("No se pudo leer el archivo de datos")
                message = f"No se pudo leer el archivo de datos {data_name}."
            else:
                plot_url = show_test(data_test, data_results, data_name, axis_y)
        elif 'data_name_write' in request.form:
            data_name_write = request.form['data_name_write']
            time_selected = request.form['time_selected']
            execute_calibration_command(host, port, username, password, data_name_write, time_selected)
            message = f"Comando ejecutado correctamente. Archivo creado: {data_name_write}"

    return render_template('pagina1.html', plot_url=plot_url, message=message)

@pagina1_bp.route('/get_data')
def get_data():
    data_name = request.args.get('data_name')
    if data_name:
        try:
            data_test, data_results = read_data(host, port, username, password, data_name)
            return jsonify(success=True, data_test=data_test, data_results=data_results)
        except Exception as e:
            return jsonify(success=False, message=str(e))
    return jsonify(success=False, message="Nombre de archivo no proporcionado")

@pagina1_bp.route('/status')
def status():
    return jsonify(status="completado")

