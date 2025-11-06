import importlib
import ToolKit.Tollbox_TG_shorten  # Importer le module

# Modifie le code dans Tollbox_TG_shorten.py ici

# Recharge le module
importlib.reload(ToolKit.Tollbox_TG_shorten)

import pyabf
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
import os
import scipy.stats as stat
import seaborn as sns
from scipy.signal import find_peaks, butter, filtfilt
from pyabf.filter import gaussian
from ToolKit.Tollbox_TG_shorten import toto_filter

folder = r"C:\Users\33666\Desktop\Spont CeA 2"
sub_folders = [x for x in glob.glob(rf'{folder}\*') if x.split('\\')[-1] not in ('analysis', 'ISteps', 'Raw')]
if not os.path.exists(rf'{folder}\analysis'): os.makedirs(rf'{folder}\analysis')

n_bin, Raw_bl = 60, []
Output_spike, Output_vm = {}, {}
Files = []

# Définir le seuil de courant en pA
current_threshold = 10e-12  # 10 pA

def highpass_filter(data, cutoff, fs, order=5):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return filtfilt(b, a, data)

def generate_analysis_pdf(raw, envlp, y_ax, binary, matrix_spike, file_name, drug, sub_folder_analysis):
    plt.figure(figsize=(12, 10))
    plt.suptitle(f'Analyse de {file_name} - {drug}', fontsize=16)

    plt.subplot(511)
    plt.plot(raw[::2])
    plt.title('Raw trace')

    plt.subplot(512)
    plt.plot(envlp[::2])
    plt.title('Enveloppe')

    plt.subplot(513)
    plt.plot(y_ax[::2])
    plt.title('Spike trace')

    plt.subplot(514)
    plt.plot(binary)
    plt.title('Binary')

    plt.subplot(515)
    plt.plot(matrix_spike)
    plt.title('Normalized Spike')
    plt.tight_layout()

    # Enregistrer le PDF
    pdf_path = rf'{sub_folder_analysis}\{file_name}_analysis.pdf'
    plt.savefig(pdf_path)
    plt.close()
    print(f'Saved analysis PDF for {file_name} at {pdf_path}')  # Message de confirmation

for sub_folder in sub_folders:
    drug = sub_folder.split('\\')[-1]
    print('\n'*2, '='*30, '\n', drug, '\n', '='*30, '\n'*2)
    files = glob.glob(rf'{sub_folder}\*.abf')
    Files.append([x.split('\\')[-1] for x in files])
    n = len(files)

    sub_folder_analysis = rf'{sub_folder}\analysis'
    if not os.path.exists(sub_folder_analysis):
        os.makedirs(sub_folder_analysis)
        print(f'Created analysis folder: {sub_folder_analysis}')  # Message de création

    matrix_spike, matrix_vm, raw_bl = np.zeros((n, n_bin)), np.zeros((n, n_bin)), []
    mean_frequencies = []  # Liste pour stocker les fréquences moyennes des fichiers
    
    for f, file in enumerate(files):
        file_name = (file.split('\\')[-1]).split('.')[0]
        print(file_name, '\n'*2)

        abf = pyabf.ABF(file)
        raw = list(abf.sweepY)[:600 * abf.sampleRate]

        # Appliquer le filtre passe-haut sur raw
        filtered_raw = highpass_filter(raw, cutoff=10, fs=abf.sampleRate)  # Ajuste le cutoff selon le bruit
        
        y_ax = toto_filter(filtered_raw, sample_rate=abf.sampleRate, freq_low=1, freq_high=500)

        gaussian(abf, 100)
        abf.setSweep(0)
        envlp = abf.sweepY[:600 * abf.sampleRate]

        # Trouver les pics vers le bas (négatifs) avec ajustement des paramètres
        indexes, _ = find_peaks(-y_ax, height=current_threshold, distance=1, prominence=30)
        binary = np.zeros(len(y_ax))
        for ind in indexes:
            binary[ind] = 1
        splt = np.sum(np.split(binary, n_bin), axis=1)
        matrix_spike[f, :] = splt  # Ne pas normaliser pour conserver les valeurs brutes
        raw_bl.append(np.mean(splt[:3]) / 5)

        # Calculer la fréquence des événements détectés pour des bins de 10 secondes
        bin_duration_sec = 10  # Durée du bin en secondes
        points_per_bin = bin_duration_sec * abf.sampleRate  # Nombre de points dans un bin de 10 secondes

        # Diviser le signal binaire en bins de 10 secondes
        n_bins = len(binary) // points_per_bin  # Nombre total de bins de 10 secondes
        binary_binned = np.array_split(binary, n_bins)  # Diviser en bins

        # Calculer la fréquence des événements dans chaque bin en spikes/10 secondes
        event_counts = [np.sum(bin_segment) for bin_segment in binary_binned]  # Nombre de spikes dans chaque bin

        # La fréquence des événements par bin est simplement le nombre de spikes par bin de 10 secondes
        event_frequencies = event_counts  # Fréquence en spikes/10 secondes

        # Calculer la fréquence moyenne pour le fichier actuel
        mean_frequency = np.mean(event_frequencies)  # Fréquence moyenne sur tous les bins
        mean_frequencies.append(mean_frequency)  # Ajouter à la liste
        print(f"{file_name}: Fréquence moyenne des spikes sur 10 secondes = {mean_frequency:.2f} spikes/10s")

        # Analyser le potentiel de membrane
        vm_values = np.nanmean(np.split(envlp, n_bin), axis=1)

        # Toujours générer le PDF
        matrix_vm[f, :] = vm_values
        print(f'{file_name}: Generating analysis PDF regardless of membrane voltage.')  # Message pour indiquer que l'analyse va se faire
        generate_analysis_pdf(raw, envlp, y_ax, binary, matrix_spike[f, :], file_name, drug, sub_folder_analysis)

    Output_spike[drug] = matrix_spike
    Output_vm[drug] = matrix_vm
    Raw_bl.append(raw_bl)

    # Sauvegarder les fréquences moyennes dans un fichier Excel
    df_mean_freq = pd.DataFrame({
        'File_Name': [file.split('\\')[-1] for file in files],
        'Mean_Frequency (spikes/10s)': mean_frequencies
    })
    df_mean_freq_path = f"{folder}/analysis/{drug}_mean_frequencies.xlsx"
    df_mean_freq.to_excel(df_mean_freq_path, index=False)
    print(f"Saved mean frequencies to {df_mean_freq_path}")

    # Création des histogrammes et des fichiers Excel
    plt.figure(figsize=(10, 6))
    plt.suptitle(drug)
    for i, (matrix, label) in enumerate(zip((matrix_spike, matrix_vm), ('AP', 'Vm'))):
        plt.subplot(4, 1, (i * 2) + 1)
        plt.ylabel(f'Time course of {label}')
        mean, sem = np.nanmean(matrix, axis=0), stat.sem(matrix, axis=0)
        plt.fill_between(np.arange(len(mean)), mean - sem, mean + sem, color='lightblue', alpha=0.25, zorder=1)
        plt.plot(mean, c='b', zorder=2)
        plt.xlabel('Time (10s bins)')
        plt.axvline(6, c='gold', lw=2)

        plt.subplot(4, 1, (i * 2) + 2)
        sns.heatmap(matrix, cmap="coolwarm")
    plt.savefig(rf'{sub_folder_analysis}\Heat course and Time map.pdf')
    plt.close()

    # Exporter les résultats dans Excel
    writer = pd.ExcelWriter(f'{folder}/analysis/{drug}_data_spike.xlsx')
    data_spike = [(cell[0:60], cell[0:60], cell[0:60]) for cell in Output_spike[drug]]
    data_spike = np.asarray([[np.nanmean(x) for x in cell] for cell in data_spike])
    
    Cell_id = [(i,) * 3 for i in range(data_spike.shape[0])]
    df_spike = pd.DataFrame({'Cell_ID': [x for y in Cell_id for x in y],
                             'Time': ('Bl', drug, 'Wsh') * data_spike.shape[0],
                             'Score': [x for y in data_spike for x in y]})

    df_spike.to_excel(writer, sheet_name='Spike Data')
    writer.save()

    writer = pd.ExcelWriter(f'{folder}/analysis/{drug}_data_vm.xlsx')
    data_vm = [(cell[:10], cell[10:20], cell[-10:]) for cell in Output_vm[drug]]
    data_vm = np.asarray([[np.nanmean(x) for x in cell] for cell in data_vm])
    
    Cell_id = [(i,) * 3 for i in range(data_vm.shape[0])]
    df_vm = pd.DataFrame({'Cell_ID': [x for y in Cell_id for x in y],
                          'Time': ('Bl', drug, 'Wsh') * data_vm.shape[0],
                          'Score': [x for y in data_vm for x in y]})

    df_vm.to_excel(writer, sheet_name='Vm Data')
    writer.save()

# Fin du script
