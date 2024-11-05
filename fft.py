#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Programa en Python para experimentar con la Transformada Discreta de Fourier y
# la transformada Rápida de Fourier (FFT)
# Copyright (C) 2024 Guillermo R. Friedrich (UTN-FRBB)
#
# Este programa es software libre; puedes redistribuirlo y/o modificarlo
# bajo los términos de la Licencia Pública General de GNU, publicada por
# la Free Software Foundation; ya sea la versión 3 de la Licencia, o
# (a tu elección) cualquier versión posterior.
#
# Este programa se distribuye con la esperanza de que sea útil, pero SIN
# NINGUNA GARANTÍA; ni siquiera la garantía implícita de COMERCIALIZACIÓN o
# ADECUACIÓN PARA UN PROPÓSITO PARTICULAR. Consulta la Licencia Pública General
# de GNU para más detalles.
#
# Deberías haber recibido una copia de la Licencia Pública General de GNU
# junto con este programa. Si no es así, consulta
# <https://www.gnu.org/licenses/>.
#
# Autor: Guillermo R. Friedrich
# Versión: 1.0
# Fecha: 2024-11-05
#

# Las dos líneas siguientes son necesarias para hacer 
# compatible el interfaz Tkinter con los programas basados 
# en versiones anteriores a la 8.5, con las más recientes. 

from tkinter import *    # Carga módulo tk (widgets estándar)
import PIL
from tkinter import ttk  # Carga ttk (para widgets nuevos 8.5+)

import matplotlib
import matplotlib.pylab as pl
import numpy as np
import math

import time

N = 2048
fm = 8000

n = [0 for i in range(N)]
s = [0 for i in range(N)]
x = [0 for i in range(N)]
z = [0 for i in range(N)]

# Coeficientes para un filtro pasabajos IIR Butterworth de orden 3,
# con frecuencia de corte de 2 kHz (para frec. de muestreo de 8 kHz)

u1 = [0,0]
u2 = [0,0]

a1 = [0.10440778E+01, -.47737628E+00]
a2 = [0.41388566E+00, 0] 

b1 = [2,1]
b2 = [1,0]

g = 0.31745307E-01

#-----------------------------------------------------------------
# Comienzo del programa.
# Comentar o descomentar según las opciones que se desea ensayar.
#-----------------------------------------------------------------
def calculo():
    global n, s, x, z, N, fm

    fm = int(tx_fm.get())
    N = int(tx_N.get())

    n = [0 for i in range(N)]
    s = [0 for i in range(N)]
    x = [0 for i in range(N)]
    z = [0 for i in range(N)]

    pl.close()   # Cierra un gráfico anterior, eventualmente abierto
    
    captura()    # Sintetiza una señal

    if int(cb_filtro.get()) != 0:
        n = filtrar()   # Aplica un pasabajos de Butterworth de orden 3

    t = time.time()     # Para calcular el tiempo de cómputo y comparar DFT vs FFT
    
    if int(rb_fft.get()) != 0:
        fft()
    else:
        dft()

    t = (time.time() - t) * 1000

    tx_tiempo.configure(state="normal")
    tx_tiempo.delete(0, END)
    tx_tiempo.insert(0, round(t, 4))
    tx_tiempo.configure(state="disabled")

    out = [0 for i in range(len(amplitud))]   # Para almacenar la salida lineal
    log = [0 for i in range(len(amplitud))]   # Para almacenar la salida en dB

    for i in range(len(amplitud)):
        try:
            j = int(frec[i].get()) * N // fm
            out[i] = x[j]
        except:
            out[i] = 0
        
    try:
        log = 20 * np.log10(out)
    except:
        pass
        
    for i in range(len(amplitud)):
        amplitud[i].configure(state="normal")
        amplitud[i].delete(0,END)
        if cb_notacion.get() != 0:
            amplitud[i].insert(0, "{:.1E}".format(out[i]))
        else:
            amplitud[i].insert(0, round(out[i],12))
        amplitud[i].configure(state="disabled")
        
        dB[i].configure(state="normal")
        dB[i].delete(0,END)
        dB[i].insert(0, str(round(log[i],2))+ " dB")
        dB[i].configure(state="disabled")

    graficar()

#----------------------------------------------------------
# Genera una señal con distintas componentes de frecuencia
#----------------------------------------------------------
def captura():
    for i in range(N):
        n[i]  = 0.00;
        for j in range(len(frec)):
            try:
                n[i] += np.cos( 2.0 * np.pi * int(frec[j].get()) * i / fm )
            except:
                pass

        if cb_ruido.get() != 0:
            n[i] += 1.00 * np.random.randn()

#-------------------------
# calculo de la DFT
#-------------------------
def dft():
    w = [[0 for i in range(2)] for j in range(N)]
    for i in range(N//2+1):
        w[i][0] = 0.0
        w[i][1] = 0.0
        x[i] = 0.0

        for j in range(N):
            w[i][0] += n[j] * np.cos( 2 * np.pi * i * j / N )
            w[i][1] -= n[j] * np.sin( 2 * np.pi * i * j / N )

        if i == 0 or i==(N//2):
            x[i] = np.sqrt( pow( w[i][0], 2 ) + pow( w[i][1], 2 ) ) / N
        else:
            x[i] = np.sqrt( pow( w[i][0], 2 ) + pow( w[i][1], 2 ) ) / (N//2)

#--------------------------------------------------------------------------
# Operación mariposa.
# Parámetros: x10 : entrada de "arriba" (parte real y parte imaginaria)
#             x11 : entrada de "abajo"  (parte real y parte imaginaria)
#             k : multiplicador del exponente de e^(- j 2 PI / 2^n)
#             divisor : divisor del exponente de e^(- j 2 PI k / 2^n)
# Retorna: dos listas con los valores real e imaginario de ambas salidas
#          de la mariposa. La primera lista contiene el resultado de
#          "arriba" y la segunda lista el resultado de "abajo"
#--------------------------------------------------------------------------
def mariposa(x10, x11, k, divisor):
    tita = - 2 * np.pi * k / divisor   # ángulo de giro 

    seno = np.sin(tita)
    coseno = np.cos(tita)

    re = x11[0] * coseno - x11[1] * seno        # aplica el giro a la "muestra" de "abajo"
    im = x11[1] * coseno + x11[0] * seno

    return [[x10[0]+re, x10[1]+im], [x10[0]-re, x10[1]-im]]

#-------------------------
# calculo de la FFT
#-------------------------
def fft():
    X = [ [ [0,0] for j in range(N) ] for i in range(2) ]
    
    r = round(math.log(N)/math.log(2))   # cantidad de etapas de la FFT
    
    pow2 = [pow(2,i) for i in range(N)]  # lista con los valores de potencias de 2.
                                         # Tiene la finalidad de hacer los cálculos una sola vez.
    entrada = r%2
    for i in range(N):
        j=bit_reversal(i,r)
        X[entrada][i][0] = n[j]
        X[entrada][i][1] = 0

    factor_etapa = 1   # Se usa para dos finalidades: como divisor para (-j 2 PI k) en la exponencial compleja
                       # y como tamaño del bloque de salida (o distancia entre los vértices de la mariposa).
                       # En la primera etapa es 2, luego se multiplica sucesivamente por 2.
        
    for m in range(r-1, -1, -1):
        # Con la variable m va recorriendo las distintas etapas de la FFT; m==0 corresponde a los resultados finales.
        factor_etapa *= 2
        entrada = (m+1)%2    # indice al "banco" de entrada
        salida = m%2         # indice a "banco" de salida

        for i in range(pow2[m]):
            # Con la variable i recorre cada uno de los subbloques de la etapa m.
            # En la primera etapa hay N subbloques con una frecuencia (k=0)
            # En la segunda etapa hay N/2 subbloques, con dos frecuencias (k: 0, 1)
            # En la tercera etapa hay N/4 subbloques, con cuatro frecuencias (k: 0, 1, 2, 3)
            # Etc.
            
            for k in range(pow2[r-m-1]):
                # Con la variable k recorre cada una de las frecuencias para cada subbloque de la etapa m
                h0 = i*factor_etapa+k
                h1 = h0 + factor_etapa//2
                
                # Cálculos de la operación mariposa
                # Implementa la mariposa en una función separada
                    
                [ X[salida][h0], X[salida][h1] ] = mariposa( X[entrada][h0], X[entrada][h1], k, factor_etapa )                   

    for i in range(N):
        x[i] = np.sqrt( pow(X[salida][i][0], 2) + pow(X[salida][i][1], 2)) / N
        if i != 0:
            x[i] = x[i] * 2
            
#-----------------------------------------------------------
# Retorna el valor que resulta de la reversión de bits de n
#-----------------------------------------------------------
def bit_reversal(n, nbits):
    m = 0
    for i in range(nbits):
        m |= ((n & (1<<i))!=0) << (nbits-i-1)
    return m

#-----------------------------------------------------------
# Primer etapa del filtro IIR
#-----------------------------------------------------------
def iir1(x):
    u = x + u1[0] * a1[0] + u1[1] * a1[1]
    y = u + u1[0] * b1[0] + u1[1] * b1[1]
    u1[1] = u1[0]
    u1[0] = u
    return y

#-----------------------------------------------------------
# Segunda etapa del filtro IIR
#-----------------------------------------------------------
def iir2(x):
    u = x + u2[0] * a2[0] + u2[1] * a2[1]
    y = u + u2[0] * b2[0] + u2[1] * b2[1]
    u2[1] = u2[0]
    u2[0] = u
    return y

#-----------------------------------------------------------
# Filtra la señal de entrada
#-----------------------------------------------------------
def filtrar():
    for i in range(N):
        aux = iir1(n[i]*g)
        s[i] = iir2(aux)
    return s

#-----------------------------------------------------------
# Presenta el gráfico del resultado de la FFT (módulo)
#-----------------------------------------------------------
def graficar():
    global x
    #x = [2*i for i in x]
    pl.plot(np.linspace(0, fm/2, int(N/2), endpoint=True),x[0:N//2])
    mng = pl.get_current_fig_manager()
#    mng.full_screen_toggle()
    pl.show()   

#-----------------------------------------------------------
# Muestra en consola los valores del módulo de la FFT
#-----------------------------------------------------------
def mostrar():
    for i in range(N//2):
        if x[i]>0.01:
            print(i,i*8000/N,x[i])
            
#---------------------------------------------------------------------------
# Función que implementa la ventana principal usando tkinter.
#---------------------------------------------------------------------------
def winMain():
    global tx_fm, tx_N, rb_fft, cb_ruido, cb_filtro, cb_notacion, frec, amplitud, dB, tx_tiempo
        
    main = Tk()
    main.geometry("650x400+150+150") # anchura x altura + coord x + coord y
    main.title("Programa para experimentar con DFT y FFT")
            
    # Crea un Button con un texto configurable y parámetro de la función a ejecutar
    # también configurable según parámetros de __init__()
    #b1 = Button(self.main, text="Graficar "+fname, width=14, command=lambda:Grafico(f))
    b1 = Button(main, text="Calcular y graficar", width=14, command=calculo)
    b1.place(x=490, y=340)

    Label(main, text="Frecuencia de muestreo:").place(x=10, y=20)
    tx_fm = Entry(main, width=8)
    tx_fm.place(x=175, y=20)                
    tx_fm.insert(0,str(fm))                                         

    Label(main, text="N (cantidad de muestras):").place(x=10, y=50)
    tx_N = Entry(main, width=8)
    tx_N.place(x=175, y=50)                
    tx_N.insert(0,str(N))                                         

    # variables a ser controladas por los radiobuttons y/o checkbuttons
    rb_fft = IntVar(main, 1)
    cb_filtro = IntVar(main, 0)
    cb_ruido = IntVar(main, 0)
    cb_notacion = IntVar(main, 0)

    # Crea un marco para los radiobuttons
    rb_frame = Frame(main, width=80, height=50, borderwidth=1, relief=SUNKEN)
    rb_frame.place(x=520, y=20)

    Radiobutton(rb_frame, text="DFT", variable=rb_fft, value=0).place(x=0, y=0)
    Radiobutton(rb_frame, text="FFT", variable=rb_fft, value=1).place(x=0, y=20)

    Checkbutton(main, text="Filtrar", variable=cb_filtro, onvalue=1, offvalue=0).place(x=490, y=100)
    
    Checkbutton(main, text="Agregar ruido", variable=cb_ruido, onvalue=1, offvalue=0).place(x=490, y=130)

    Checkbutton(main, text="Notación científica", variable=cb_notacion, onvalue=1, offvalue=0).place(x=490, y=160)

    frec = [0 for i in range(10)]
    amplitud = [0 for i in range(10)]
    dB = [0 for i in range(10)]

    init_frecs = (0, 500, 1000, 1500, 2000, 2500, 3000, 3500, "", "")
    
    rb_frame2 = Frame(main, width=20, height=80, borderwidth=1, relief=SUNKEN)
    rb_frame2.place(x=20,y=100)

    Label(rb_frame2, text="Frecuencia").grid(row=0, column=0)
    Label(rb_frame2, text="Amplitud").grid(row=0, column=1)
    Label(rb_frame2, text="Amplitud (dB)").grid(row=0, column=2)

    for i in range(10):
        frec[i] = Entry(rb_frame2, width=10)
        frec[i].grid(row=i+1, column=0)
        
        frec[i].insert(0, init_frecs[i])

        amplitud[i] = Entry(rb_frame2, width=24)
        amplitud[i].grid(row=i+1, column=1)
        amplitud[i].configure(state="disabled", disabledforeground="black")
    
        dB[i] = Entry(rb_frame2, width=20)
        dB[i].grid(row=i+1, column=2)
        dB[i].configure(state="disabled", disabledforeground="black")

    rb_frame3 = Frame(main, width=12, height=20, borderwidth=1, relief=SUNKEN)
    rb_frame3.place(x=490,y=250)

    Label(rb_frame3, text=" tiempo de cálculo (ms) ").grid(row=0, column=0)
    tx_tiempo = Entry(rb_frame3, width=18)
    tx_tiempo.grid(row=1, column=0)
    tx_tiempo.configure(state="disabled", disabledforeground="black", justify="center")
        
matplotlib.use('TkAgg')  # Para poder usar pyplot desde una GUI

winMain()
