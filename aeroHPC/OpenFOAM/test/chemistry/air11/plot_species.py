import os
import matplotlib.pyplot as plt
import sys

def plot_species_from_folder(folder_path):
    # Pulisce il path rimuovendo eventuali slash finali doppi o spazi
    folder_path = folder_path.rstrip(os.sep)

    if not os.path.exists(folder_path):
        print(f"Errore: La cartella '{folder_path}' non esiste.")
        print(f"Uso: python3 plot_species.py <nome_cartella>")
        return

    print(f"Analisi della cartella: {folder_path} ...")

    plt.figure(figsize=(10, 7))

    # Filtra file nascosti e ordina
    files = sorted([f for f in os.listdir(folder_path) if not f.startswith('.')])
    files_found = False

    for filename in files:
        filepath = os.path.join(folder_path, filename)

        if not os.path.isfile(filepath):
            continue

        times = []
        concentrations = []

        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            times.append(float(parts[0]))
                            concentrations.append(float(parts[1]))
                        except ValueError:
                            continue

            if times:
                # Usa il nome del file come etichetta nella legenda
                plt.plot(times, concentrations, label=filename, linewidth=2)
                files_found = True
                print(f"  -> Caricato: {filename}")

        except Exception as e:
            print(f"  -> Errore lettura {filename}: {e}")

    if not files_found:
        print("Nessun dato valido trovato.")
        return

    plt.xlabel('Time [s]')
    plt.ylabel('Concentration')
    plt.title(f'Species Evolution (Folder: {folder_path})')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left') # Legenda fuori dal grafico
    plt.grid(True, which="both", linestyle='--', alpha=0.6)
    plt.tight_layout()

    output_name = f'species_plot_{os.path.basename(folder_path)}.png'
    plt.savefig(output_name)
    print(f"Grafico salvato come '{output_name}'")
    plt.show()

if __name__ == "__main__":
    # Default: cerca nella cartella '0' se non viene specificato altro
    target_folder = '0'

    # Se l'utente passa un argomento (es. python plot.py 0.001/) usa quello
    if len(sys.argv) > 1:
        target_folder = sys.argv[1]

    plot_species_from_folder(target_folder)
