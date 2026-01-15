import re
import matplotlib.pyplot as plt

def extract_and_plot(file_path):
    t_tr_values = []
    t_v_values = []
    iterations = []

    # Pattern migliorato: cerca Ttr= seguito da numeri/punti,
    # ma si ferma prima di trovare uno spazio o la lettera successiva
    pattern = re.compile(r"Ttr=([0-9.]+)\s+Tv=([0-9.]+)")

    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Se il file è molto sporco, proviamo a pulire la riga
                match = pattern.search(line)
                if match:
                    try:
                        # Gestione dell'errore specifico se la stringa è malformata
                        raw_tr = match.group(1)
                        raw_v = match.group(2)

                        # Se per qualche motivo ha concatenato i numeri (es. 5501.245501.24)
                        # prendiamo solo la prima metà o solleviamo un alert
                        t_tr = float(raw_tr)
                        t_v = float(raw_v)

                        t_tr_values.append(t_tr)
                        t_v_values.append(t_v)
                        iterations.append(len(iterations) + 1)
                    except ValueError:
                        # Salta la riga se il dato è corrotto
                        continue

        if not t_tr_values:
            print("Nessun dato valido trovato. Controlla il formato del log.")
            return

        plt.figure(figsize=(10, 6))
        plt.plot(iterations, t_tr_values, label='$T_{tr}$ (Translational-Rotational)', color='blue', linewidth=1.5)
        plt.plot(iterations, t_v_values, label='$T_v$ (Vibrational)', color='red', linestyle='--')

        plt.xlabel('Campione (Time step/Iteration)')
        plt.ylabel('Temperatura [K]')
        plt.title('Evoluzione Temperature - log.foamRun')
        plt.legend()
        plt.grid(True, which='both', linestyle=':', alpha=0.5)

        plt.tight_layout()
        plt.savefig('temperature_plot.png')
        print(f"Processati {len(t_tr_values)} punti. Grafico generato.")
        plt.show()

    except FileNotFoundError:
        print(f"Errore: File '{file_path}' non trovato.")

if __name__ == "__main__":
    extract_and_plot('log.foamRun')
