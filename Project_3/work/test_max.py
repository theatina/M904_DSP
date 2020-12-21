#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 21:17:29 2020

@author: theatina
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wf
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import sounddevice as sd

sr_a, a = wf.read('./data/a_audio.wav')
sr_o, o = wf.read('./data/o_audio.wav')
sr_e, e = wf.read('./data/e_audio.wav')
sr_s, s = wf.read('./data/s_audio.wav')
sr_k, k = wf.read('./data/k_audio.wav')

fs = 44100
power_of_2 = 11
time_steps = 2**power_of_2
# time_steps = 4410

fft_a = np.fft.fft( a[:time_steps,0] )
fft_o = np.fft.fft( o[:time_steps,0] )
fft_e = np.fft.fft( e[:time_steps,0] )
fft_s = np.fft.fft( s[:time_steps,0] )
fft_k = np.fft.fft( k[:time_steps,0] )

# print(len(freq_bins))
mag_a = np.sqrt( np.power( fft_a.real , 2 ) + np.power( fft_a.imag , 2 ) )
mag_a = mag_a[:time_steps//2]

mag_o = np.sqrt( np.power( fft_o.real , 2 ) + np.power( fft_o.imag , 2 ) )
mag_o = mag_o[:time_steps//2]

mag_e = np.sqrt( np.power( fft_e.real , 2 ) + np.power( fft_e.imag , 2 ) )
mag_e = mag_e[:time_steps//2]

mag_s = np.sqrt( np.power( fft_s.real , 2 ) + np.power( fft_s.imag , 2 ) )
mag_s = mag_s[:time_steps//2]

mag_k = np.sqrt( np.power( fft_k.real , 2 ) + np.power( fft_k.imag , 2 ) )
mag_k = mag_k[:time_steps//2]

mag_list = [mag_a,mag_o,mag_e,mag_s,mag_k]

freq_bins = np.linspace(0,fs,time_steps)
freq_bins = freq_bins[:time_steps//2]

#testing
plt.plot( freq_bins[:mag_o.size//2], mag_o[:mag_o.size//2])
print(np.argmax(max(mag_o)))

phonemes = ['a', 'o', 'e', 's', 'k']

fig, plots = plt.subplots(3,2,figsize=(27,21))#, sharex=True)
counter=0
for i in range(3):
    for j in range(2):
        if counter==5:
            break
        
        mag_signal = mag_list[counter]
#         print(mag_signal)
        
        plots[i,j].grid(True)
        plots[i,j].title.set_text(r"$\bf{2^{%d} = %d\ bins}$"%(power_of_2,time_steps))
        plots[i,j].plot( freq_bins , mag_signal )        

        #bins of max magnitude (first 4)
        mag_signal_values = np.sort(mag_signal)
#         print(np.where(mag_signal_values in mag_signal_values[:4]))
#         print(f" last 4 {mag_signal_values[:4]}")
        bins_max_mag = [i for i in range(mag_signal.size) if mag_signal[i] in mag_signal_values[-4:]]
        print(bins_max_mag,freq_bins[bins_max_mag])

        x_axis_to_plot = []
        x_axis_to_plot.extend(bins_max_mag)
        x_axis_to_plot.append(bins_max_mag[0]-1)
        x_axis_to_plot.append(bins_max_mag[-1]+1)
#         print(x_axis_to_plot,len(x_axis_to_plot))
        
        x_values_to_plot = []
        x_values_to_plot.extend(freq_bins[x_axis_to_plot].tolist())
        x_values_to_plot.sort()
        print(x_values_to_plot)
    
        plots[i,j].set_xlim(x_values_to_plot[0], x_values_to_plot[-1])#freq_bins_signal[x_axis_to_plot[-1]])
        plots[i,j].set_ylim(0,1.2*max(mag_signal))
        plots[i,j].set_xticks(x_values_to_plot)
        
        #plot points
        f1,f2,f3,f4 = mag_signal_values[-4:]
        print(freq_bins[bins_max_mag[1]])
        p_f1, = plots[i,j].plot(freq_bins[bins_max_mag[0]] , mag_signal[bins_max_mag[0]], "ro")#, label="$f_1\ bin$" )
        p_f2, = plots[i,j].plot(freq_bins[bins_max_mag[1]] , mag_signal[bins_max_mag[1]], "go")#, label="$f_2\ bin$" )
        p_f3, = plots[i,j].plot(freq_bins[bins_max_mag[2]] , mag_signal[bins_max_mag[2]], "bo")#, label="$f_3\ bin$" )
        p_f4, = plots[i,j].plot(freq_bins[bins_max_mag[3]] , mag_signal[bins_max_mag[3]], "o")#, label="$f_4\ bin$" )

        
        x = [freq_bins[bins_max_mag[0]], freq_bins[bins_max_mag[1]], freq_bins[bins_max_mag[2]],freq_bins[bins_max_mag[3]]]
        y = [mag_signal[bins_max_mag[0]], mag_signal[bins_max_mag[1]], mag_signal[bins_max_mag[2]],mag_signal[bins_max_mag[3]]]
        
        points = [freq_bins[bins_max_mag[0]], freq_bins[bins_max_mag[1]], freq_bins[bins_max_mag[2]],freq_bins[bins_max_mag[3]]]
        
        for k, freq in enumerate(points):
            plots[i,j].annotate(int(freq), (x[k], y[k]))
        
        phoneme = phonemes[counter]
        plots[i,j].title.set_text(r"$\bf{Phoneme\ /\ %s\ / }$"%(phoneme))

        
        counter+=1
        
       
        
# fig.legend(handles=[p_f1,p_f2,p_f3])    
       
