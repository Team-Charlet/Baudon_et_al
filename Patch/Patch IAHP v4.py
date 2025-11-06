# -*- coding: utf-8 -*-
"""
Patch Clamp Analyzer for IAHP Analysis - Version Corrigée Détection
Détection automatique de la stimulation et mesure de l'IAHP
CORRECTION: Détection basée sur courant entrant négatif pendant stimulation
Created on Fri Jul 18 14:38:16 2025

@author: 33666
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os
from datetime import datetime
import warnings
import gc  # Pour la gestion de la mémoire
warnings.filterwarnings('ignore')

# Configuration matplotlib pour optimiser la mémoire
plt.rcParams['figure.max_open_warning'] = 50
plt.rcParams['agg.path.chunksize'] = 10000

class PatchClampAnalyzer:
    """
    Analyseur corrigé pour traces patch-clamp
    CORRECTION: Détecte stimulation via courant entrant négatif, puis mesure IAHP sortant positif
    """
    
    def __init__(self, sampling_rate=20000, max_figure_size=(16, 12)):
        """
        Initialise l'analyseur
        
        Parameters:
        -----------
        sampling_rate : int
            Fréquence d'échantillonnage (Hz)
        max_figure_size : tuple
            Taille maximale des figures (largeur, hauteur) en inches
        """
        self.sampling_rate = sampling_rate
        self.max_figure_size = max_figure_size
        self.results = []
        
        # Créer dossier pour les graphiques individuels
        self.output_dir = Path("patch_clamp_analysis_output")
        self.plots_dir = self.output_dir / "individual_cell_plots"
        self.output_dir.mkdir(exist_ok=True)
        self.plots_dir.mkdir(exist_ok=True)
        
    def load_wcp_file(self, file_path):
        """
        Charge un fichier WCP (version simplifiée)
        """
        try:
            with open(file_path, 'rb') as f:
                # Saut de l'en-tête (1024 bytes)
                f.seek(1024)
                
                # Lecture de toutes les données
                raw_data = f.read()
                
                # Conversion en entiers 16-bit
                data_array = np.frombuffer(raw_data, dtype='<i2').astype(float)
                
                # Estimation : 60 secondes à 20kHz = 1.2M points par sweep
                estimated_points_per_sweep = 60 * self.sampling_rate
                n_sweeps = len(data_array) // estimated_points_per_sweep
                
                if n_sweeps > 0:
                    # Reshape en sweeps
                    total_points = n_sweeps * estimated_points_per_sweep
                    data_matrix = data_array[:total_points].reshape(n_sweeps, estimated_points_per_sweep)
                    
                    # Vecteur temps
                    time_array = np.arange(estimated_points_per_sweep) / self.sampling_rate
                    
                    return {
                        'data': data_matrix,
                        'time': time_array,
                        'n_sweeps': n_sweeps,
                        'n_samples': estimated_points_per_sweep
                    }
                else:
                    print(f"Impossible de déterminer la structure des données dans {file_path}")
                    return None
                    
        except Exception as e:
            print(f"Erreur lors du chargement de {file_path}: {e}")
            return None
    
    def detect_stimulation(self, current_trace):
        """
        DÉTECTION CORRIGÉE : Cherche un courant sortant POSITIF (stimulation dépolarisante)
        Stim ~200ms de -20mV à +20mV crée un courant positif
        """
        # Ignorer les 2 premières secondes pour éviter l'artefact initial
        start_search = int(2 * self.sampling_rate)
        end_search = int(15 * self.sampling_rate)  # Chercher jusqu'à 15s max
        search_trace = current_trace[start_search:end_search]
        
        # Calcul baseline locale pour détection
        baseline_level = np.median(search_trace)
        noise_std = np.std(search_trace)
        
        # Seuil de détection : doit être significativement au-dessus de la baseline
        # Sur votre trace, le pic est à ~3000 pA vs baseline ~170 pA
        threshold = baseline_level + 5 * noise_std  # Seuil robuste
        
        # Détecter tous les pics au-dessus du seuil
        above_threshold = search_trace > threshold
        
        if not np.any(above_threshold):
            print(f"    ATTENTION: Aucun pic détecté au-dessus de {threshold:.0f} pA")
            # Fallback : prendre le pic maximum
            peak_idx = np.argmax(search_trace)
            stim_start = start_search + peak_idx
            # Fin = début + 200ms (durée attendue de la stim)
            stim_end = stim_start + int(0.2 * self.sampling_rate)
            return stim_start, min(stim_end, len(current_trace) - 1)
        
        # Trouver le début et fin de la stimulation
        stim_indices = np.where(above_threshold)[0]
        stim_start_idx = stim_indices[0]
        stim_start = start_search + stim_start_idx
        
        # CORRECTION PRINCIPALE : Chercher la vraie fin de la stimulation
        # Ne pas prendre le dernier point au-dessus du seuil (faux)
        # Mais chercher où le signal revient proche de la baseline
        
        # Partir du pic principal et chercher la descente
        peak_in_search = np.argmax(search_trace[stim_start_idx:stim_start_idx + int(1 * self.sampling_rate)])
        peak_absolute = stim_start + peak_in_search
        
        # Chercher après le pic où le signal revient vers la baseline
        post_peak_trace = current_trace[peak_absolute:peak_absolute + int(1 * self.sampling_rate)]
        
        # Critère de fin : signal revient dans les 20% de la baseline + 3*sigma
        end_threshold = baseline_level + 3 * noise_std
        
        for i, value in enumerate(post_peak_trace):
            if value <= end_threshold:
                stim_end = peak_absolute + i
                break
        else:
            # Si pas trouvé, prendre 200ms après le pic (durée attendue)
            stim_end = peak_absolute + int(0.2 * self.sampling_rate)
        
        # Vérification : la stimulation ne doit pas durer plus de 500ms
        max_stim_duration = int(0.5 * self.sampling_rate)
        if (stim_end - stim_start) > max_stim_duration:
            stim_end = stim_start + max_stim_duration
            print(f"    CORRECTION: Durée de stimulation limitée à 500ms")
        
        return stim_start, min(stim_end, len(current_trace) - 1)
    
    def calculate_baseline(self, current_trace, stim_start, duration=2.0):
        """
        Baseline plus courte : 2 secondes avant stimulation
        """
        baseline_points = int(duration * self.sampling_rate)
        baseline_start = max(0, stim_start - baseline_points)
        baseline_end = stim_start
        
        baseline_trace = current_trace[baseline_start:baseline_end]
        return np.mean(baseline_trace)
    
    def find_iahp_window(self, current_trace, stim_end, baseline, max_duration=1.5):
        """
        Trouve la fin de l'IAHP quand le courant retourne à la baseline
        """
        # Section post-stimulation
        post_stim_trace = current_trace[stim_end:]
        
        # Seuil de retour à la baseline : 5% de l'amplitude max de l'IAHP
        iahp_amplitude = np.max(post_stim_trace) - baseline
        threshold = 0.05 * iahp_amplitude  # 5% de l'amplitude
        
        # Chercher le retour à la baseline (avec un minimum de 200ms)
        min_duration = int(0.01 * self.sampling_rate)  # 200ms minimum
        max_points = int(max_duration * self.sampling_rate)
        
        for i in range(min_duration, min(len(post_stim_trace), max_points)):
            if abs(post_stim_trace[i] - baseline) <= threshold:
                return stim_end + i
        
        # Si pas trouvé, prendre la durée maximale
        return min(stim_end + max_points, len(current_trace) - 1)
    
    def calculate_iahp_auc(self, current_trace, baseline, start_idx, end_idx):
        """
        Calcule l'AUC de l'IAHP : courant sortant positif qui diminue
        """
        iahp_trace = current_trace[start_idx:end_idx]
        baseline_corrected = iahp_trace - baseline
        
        # IAHP = courant sortant positif au-dessus de la baseline
        positive_current = np.where(baseline_corrected > 0, baseline_corrected, 0)
        
        # Calcul de l'amplitude maximale
        max_amplitude = np.max(baseline_corrected) if len(baseline_corrected) > 0 else 0
        
        # Intégration
        time_step = 1 / self.sampling_rate
        auc = np.trapz(positive_current, dx=time_step)
        
        return auc, baseline_corrected, positive_current, max_amplitude

    def plot_single_sweep_analysis(self, wcp_data, sweep_result, cell_info, sweep_idx, save_path):
        """
        Crée un graphique détaillé pour UN SEUL sweep (optimisé mémoire)
        """
        # Force la fermeture des figures précédentes
        plt.close('all')
        gc.collect()
        
        # Figure de taille raisonnable
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'{cell_info["condition"].upper()} - {cell_info["cell_id"]} - Sweep {sweep_idx+1}', 
                     fontsize=14, weight='bold')
        
        current_trace = wcp_data['data'][sweep_result['sweep']]
        time_array = wcp_data['time']
        
        # Données de base
        baseline = sweep_result['baseline']
        stim_start = sweep_result['stim_start']
        stim_end = sweep_result['stim_end']
        analysis_end = sweep_result['analysis_end']
        
        # Temps correspondants
        stim_start_time = stim_start / self.sampling_rate
        stim_end_time = stim_end / self.sampling_rate
        analysis_end_time = analysis_end / self.sampling_rate
        
        # **GRAPHIQUE 1 : Vue complète**
        ax1 = axes[0, 0]
        ax1.plot(time_array, current_trace, 'b-', linewidth=1, alpha=0.8)
        ax1.axhline(y=baseline, color='red', linestyle='--', linewidth=2, label=f'Baseline ({baseline:.1f} pA)')
        ax1.axvspan(stim_start_time, stim_end_time, alpha=0.3, color='orange', label='Stimulation')
        ax1.axvspan(stim_end_time, analysis_end_time, alpha=0.2, color='green', label='IAHP')
        ax1.set_title('Vue complète')
        ax1.set_xlabel('Temps (s)')
        ax1.set_ylabel('Courant (pA)')
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)
        
        # **GRAPHIQUE 2 : Zoom sur stimulation + IAHP**
        ax2 = axes[0, 1]
        zoom_start = max(0, stim_start - int(1 * self.sampling_rate))
        zoom_end = min(len(current_trace), analysis_end + int(0.5 * self.sampling_rate))
        zoom_time = time_array[zoom_start:zoom_end]
        zoom_current = current_trace[zoom_start:zoom_end]
        
        ax2.plot(zoom_time, zoom_current, 'b-', linewidth=2)
        ax2.axhline(y=baseline, color='red', linestyle='--', linewidth=2)
        ax2.axvspan(stim_start_time, stim_end_time, alpha=0.4, color='orange')
        
        # Aire IAHP
        iahp_trace = current_trace[stim_end:analysis_end]
        iahp_time = time_array[stim_end:analysis_end]
        iahp_corrected = iahp_trace - baseline
        positive_mask = iahp_corrected > 0
        
        if np.any(positive_mask):
            ax2.fill_between(iahp_time[positive_mask], baseline, 
                           iahp_trace[positive_mask], alpha=0.5, color='green')
        
        ax2.set_title(f'Zoom Stimulation + IAHP\nAUC = {sweep_result["auc"]:.2f} pA.s')
        ax2.set_xlabel('Temps (s)')
        ax2.set_ylabel('Courant (pA)')
        ax2.grid(True, alpha=0.3)
        
        # **GRAPHIQUE 3 : IAHP corrigé**
        ax3 = axes[1, 0]
        ax3.plot(iahp_time, iahp_corrected, 'g-', linewidth=2, label='IAHP corrigé')
        ax3.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        ax3.fill_between(iahp_time, 0, iahp_corrected, where=(iahp_corrected > 0), 
                        alpha=0.4, color='green', interpolate=True)
        
        # Peak
        if len(iahp_corrected) > 0:
            peak_iahp = sweep_result['max_amplitude']
            peak_idx = np.argmax(iahp_corrected)
            if peak_idx < len(iahp_time):
                peak_time = iahp_time[peak_idx]
                ax3.plot(peak_time, peak_iahp, 'ro', markersize=8, 
                        label=f'Peak = {peak_iahp:.1f} pA')
        
        ax3.set_title('IAHP corrigé de la baseline')
        ax3.set_xlabel('Temps (s)')
        ax3.set_ylabel('Courant (pA)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # **GRAPHIQUE 4 : Résumé métriques**
        ax4 = axes[1, 1]
        ax4.axis('off')
        
        # Tableau des métriques
        metrics_text = f"""
MÉTRIQUES D'ANALYSE:

Baseline: {baseline:.2f} pA
Peak IAHP: {sweep_result['max_amplitude']:.2f} pA
AUC IAHP: {sweep_result['auc']:.3f} pA.s

TIMING:
Début stimulation: {stim_start_time:.2f} s
Fin stimulation: {stim_end_time:.2f} s
Durée stimulation: {(stim_end_time-stim_start_time):.2f} s
Durée analyse IAHP: {sweep_result['analysis_duration']:.2f} s

QUALITÉ:
Points analysés: {analysis_end - stim_end:,}
Fréquence: {self.sampling_rate:,} Hz
        """
        
        ax4.text(0.1, 0.9, metrics_text, transform=ax4.transAxes, fontsize=11,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()  # Important: fermer explicitement
        gc.collect()  # Forcer le nettoyage mémoire
        
        return save_path

    def plot_cell_summary(self, wcp_data, result, save_path):
        """
        Crée un graphique résumé pour toute la cellule (tous les sweeps)
        """
        plt.close('all')
        gc.collect()
        
        n_sweeps = result['n_sweeps']
        
        # Adapter la taille selon le nombre de sweeps (mais limiter)
        height = min(4 + 2 * n_sweeps, 16)  # Max 16 inches de hauteur
        fig, axes = plt.subplots(n_sweeps + 1, 1, figsize=(14, height))
        
        if n_sweeps == 1:
            axes = [axes[0], axes[1]]  # S'assurer qu'on a une liste
        
        colors = plt.cm.Set1(np.linspace(0, 1, n_sweeps))
        
        # **Graphique de résumé en haut**
        ax_summary = axes[0]
        
        for i, sweep_result in enumerate(result['sweep_results']):
            current_trace = wcp_data['data'][sweep_result['sweep']]
            time_array = wcp_data['time']
            
            # Sous-échantillonner pour la vue d'ensemble (économie mémoire)
            step = max(1, len(current_trace) // 10000)  # Max 10k points par trace
            ax_summary.plot(time_array[::step], current_trace[::step], 
                          color=colors[i], alpha=0.7, linewidth=1,
                          label=f'Sweep {i+1} (AUC={sweep_result["auc"]:.2f})')
            ax_summary.axhline(y=sweep_result['baseline'], color=colors[i], 
                             linestyle='--', alpha=0.5)
        
        ax_summary.set_title(f'Résumé {result["condition"].upper()} - {result["cell_id"]}\n'
                           f'AUC moyenne: {result["mean_auc"]:.3f} ± {result["std_auc"]:.3f} pA.s')
        ax_summary.set_xlabel('Temps (s)')
        ax_summary.set_ylabel('Courant (pA)')
        ax_summary.legend(fontsize=9, loc='best')
        ax_summary.grid(True, alpha=0.3)
        
        # **Graphiques individuels pour chaque sweep**
        for i, sweep_result in enumerate(result['sweep_results']):
            ax = axes[i + 1]
            current_trace = wcp_data['data'][sweep_result['sweep']]
            time_array = wcp_data['time']
            
            # Zoom sur la région intéressante
            stim_start = sweep_result['stim_start']
            analysis_end = sweep_result['analysis_end']
            zoom_start = max(0, stim_start - int(2 * self.sampling_rate))
            zoom_end = min(len(current_trace), analysis_end + int(1 * self.sampling_rate))
            
            zoom_time = time_array[zoom_start:zoom_end]
            zoom_current = current_trace[zoom_start:zoom_end]
            
            ax.plot(zoom_time, zoom_current, color=colors[i], linewidth=2)
            ax.axhline(y=sweep_result['baseline'], color='red', linestyle='--', alpha=0.8)
            
            # Marquer les zones
            stim_start_time = stim_start / self.sampling_rate
            stim_end_time = sweep_result['stim_end'] / self.sampling_rate
            analysis_end_time = analysis_end / self.sampling_rate
            
            ax.axvspan(stim_start_time, stim_end_time, alpha=0.3, color='orange')
            ax.axvspan(stim_end_time, analysis_end_time, alpha=0.2, color='green')
            
            ax.set_title(f'Sweep {i+1}: AUC = {sweep_result["auc"]:.3f} pA.s, '
                        f'Durée = {sweep_result["analysis_duration"]:.1f}s')
            ax.set_xlabel('Temps (s)')
            ax.set_ylabel('Courant (pA)')
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        gc.collect()
        
        return save_path

    def analyze_single_file(self, file_path, condition, cell_id):
        """
        Analyse un fichier WCP complet avec gestion mémoire optimisée
        """
        wcp = self.load_wcp_file(file_path)
        if wcp is None:
            return None
        
        print(f"\n=== Analyse {condition.upper()} - {cell_id} ===")
        print(f"Sweeps détectés: {wcp['n_sweeps']}")
        
        sweep_results = []
        
        for sweep in range(wcp['n_sweeps']):
            print(f"  Sweep {sweep + 1}/{wcp['n_sweeps']}:")
            current_trace = wcp['data'][sweep]
            
            # 1. Détecter la stimulation
            stim_start, stim_end = self.detect_stimulation(current_trace)
            print(f"    Stimulation: {stim_start/self.sampling_rate:.3f}s à {stim_end/self.sampling_rate:.3f}s")
            
            # 2. Calculer la baseline
            baseline = self.calculate_baseline(current_trace, stim_start)
            print(f"    Baseline: {baseline:.2f} pA")
            
            # 3. Trouver la fin de l'IAHP
            analysis_end = self.find_iahp_window(current_trace, stim_end, baseline)
            print(f"    Fin IAHP: {analysis_end/self.sampling_rate:.3f}s")
            
            # 4. Calculer l'AUC de l'IAHP
            auc, baseline_corrected, positive_current, max_amplitude = self.calculate_iahp_auc(
                current_trace, baseline, stim_end, analysis_end
            )
            print(f"    AUC IAHP: {auc:.2f} pA.s")
            print(f"    Amplitude max: {max_amplitude:.2f} pA")
            
            sweep_results.append({
                'sweep': sweep,
                'baseline': baseline,
                'auc': auc,
                'max_amplitude': max_amplitude,
                'stim_start': stim_start,
                'stim_end': stim_end,
                'analysis_end': analysis_end,
                'analysis_duration': (analysis_end - stim_end) / self.sampling_rate
            })
        
        # Résultats moyens
        aucs = [r['auc'] for r in sweep_results]
        amplitudes = [r['max_amplitude'] for r in sweep_results]
        mean_auc = np.mean(aucs)
        std_auc = np.std(aucs)
        mean_amplitude = np.mean(amplitudes)
        std_amplitude = np.std(amplitudes)
        
        result = {
            'file_path': file_path,
            'condition': condition,
            'cell_id': cell_id,
            'n_sweeps': len(sweep_results),
            'mean_auc': mean_auc,
            'std_auc': std_auc,
            'mean_amplitude': mean_amplitude,
            'std_amplitude': std_amplitude,
            'sweep_results': sweep_results
        }
        
        # Créer le graphique résumé de la cellule
        summary_filename = f"{condition}_{cell_id}_summary.png"
        summary_path = self.plots_dir / summary_filename
        self.plot_cell_summary(wcp, result, summary_path)
        
        # Créer des graphiques détaillés pour chaque sweep (si peu de sweeps)
        if len(sweep_results) <= 5:  # Limiter pour éviter trop de fichiers
            for i, sweep_result in enumerate(sweep_results):
                sweep_filename = f"{condition}_{cell_id}_sweep_{i+1}_detailed.png"
                sweep_path = self.plots_dir / sweep_filename
                cell_info = {'condition': condition, 'cell_id': cell_id}
                self.plot_single_sweep_analysis(wcp, sweep_result, cell_info, i, sweep_path)
        
        print(f"  Résultats: AUC = {mean_auc:.2f} ± {std_auc:.2f} pA.s")
        print(f"  Graphiques sauvegardés dans: {self.plots_dir}")

        return result
    
    def analyze_batch(self, data_directory):
        """
        Analyse par lot selon l'organisation des dossiers
        """
        data_path = Path(data_directory)
        
        if not data_path.exists():
            print(f"Répertoire non trouvé: {data_directory}")
            return
        
        all_results = []
        
        # Parcours des conditions
        for condition_dir in data_path.iterdir():
            if not condition_dir.is_dir():
                continue
                
            condition = condition_dir.name
            print(f"\nAnalyse de la condition: {condition}")
            
            # Parcours des fichiers WCP
            wcp_files = list(condition_dir.glob("*.wcp"))
            for wcp_file in wcp_files:
                cell_id = wcp_file.stem
                
                result = self.analyze_single_file(str(wcp_file), condition, cell_id)
                
                if result:
                    all_results.append(result)
                    self.results.append(result)
                
                # Nettoyage mémoire entre chaque fichier
                gc.collect()
        
        print(f"\nAnalyse terminée. {len(all_results)} cellules analysées.")
        return all_results
    
    def export_results(self):
        """
        Exporte les résultats vers Excel
        """
        if not self.results:
            print("Aucun résultat à exporter")
            return
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = self.output_dir / f"iahp_results_{timestamp}.xlsx"
        
        # DataFrame principal (par cellule)
        main_data = []
        for result in self.results:
            main_data.append({
                'condition': result['condition'],
                'cell_id': result['cell_id'],
                'n_sweeps': result['n_sweeps'],
                'mean_auc_pA_s': result['mean_auc'],
                'std_auc_pA_s': result['std_auc'],
                'sem_auc_pA_s': result['std_auc'] / np.sqrt(result['n_sweeps']),
                'mean_amplitude_pA': result['mean_amplitude'],    
                'std_amplitude_pA': result['std_amplitude'],   
                'sem_amplitude_pA': result['std_amplitude'] / np.sqrt(result['n_sweeps'])  
            })
        df_main = pd.DataFrame(main_data)
        
        # DataFrame détaillé (par sweep)
        detailed_data = []
        for result in self.results:
            for sweep_data in result['sweep_results']:
                detailed_data.append({
                    'condition': result['condition'],
                    'cell_id': result['cell_id'],
                    'sweep': sweep_data['sweep'],
                    'baseline_pA': sweep_data['baseline'],
                    'auc_pA_s': sweep_data['auc'],
                    'max_amplitude_pA': sweep_data['max_amplitude'],
                    'analysis_duration_s': sweep_data['analysis_duration']
                })
        df_detailed = pd.DataFrame(detailed_data)
        
        # Export vers Excel
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df_main.to_excel(writer, sheet_name='Resume_Cellules', index=False)
            df_detailed.to_excel(writer, sheet_name='Donnees_Detaillees', index=False)
            
            # Statistiques par condition (AUC ET Amplitude)
            stats_auc = df_main.groupby('condition')['mean_auc_pA_s'].agg(['count', 'mean', 'std', 'sem']).round(4)
            stats_amplitude = df_main.groupby('condition')['mean_amplitude_pA'].agg(['count', 'mean', 'std', 'sem']).round(4)
            
            stats_auc.to_excel(writer, sheet_name='Stats_AUC')
            stats_amplitude.to_excel(writer, sheet_name='Stats_Amplitude')  # NOUVEAU
        
        print(f"\n=== RÉSULTATS EXPORTÉS ===")
        print(f"Fichier Excel: {output_file}")
        
        # Statistiques résumées
        print(f"\n=== STATISTIQUES AUC ===")
        print(stats_auc)
        print(f"\n=== STATISTIQUES AMPLITUDE ===")  # NOUVEAU
        print(stats_amplitude)
        
        return output_file

# Utilisation
if __name__ == "__main__":
    # Initialisation avec limitation de taille des figures
    analyzer = PatchClampAnalyzer(sampling_rate=20000, max_figure_size=(16, 12))
    
    # Analyse
    data_directory = r"D:\Valentin.GRELOT\Desktop\Revision papier\round 2\Moi\IAHP\IAHP 2"
    
    print("=== ANALYSE PATCH-CLAMP IAHP OPTIMISÉE MÉMOIRE ===\n")
    print("Version corrigée pour éviter les erreurs de mémoire matplotlib\n")
    
    results = analyzer.analyze_batch(data_directory)
    
    if results:
        analyzer.export_results()
        print(f"\n=== ANALYSE TERMINÉE ===")
        print(f"Fichiers générés dans: {analyzer.output_dir}")
    else:
        print("Aucune donnée analysée. Vérifiez le chemin et la structure des fichiers.")
    
    # Nettoyage final
    plt.close('all')
    gc.collect()
    print("Nettoyage mémoire terminé.")
