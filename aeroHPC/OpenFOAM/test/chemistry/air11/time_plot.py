import sys
import re
import matplotlib.pyplot as plt

def parse_speedup_file(file_path):
    cores, real_times, user_times = [], [], []
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        sections = re.split(r'(\d+)\s+core', content)
        for i in range(1, len(sections), 2):
            n_core = int(sections[i])
            block = sections[i+1]
            r_m = re.search(r'real\s+(\d+)m([\d.]+)s', block)
            u_m = re.search(r'user\s+(\d+)m([\d.]+)s', block)
            if r_m and u_m:
                cores.append(n_core)
                real_times.append(int(r_m.group(1))*60 + float(r_m.group(2)))
                user_times.append(int(u_m.group(1))*60 + float(u_m.group(2)))
        return cores, real_times, user_times
    except Exception as e:
        print(f"Errore: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <file>")
        sys.exit(1)

    cores, real_t, user_t = parse_speedup_file(sys.argv[1])

    # Calcolo User Time / Cores
    user_per_core = [u / c for u, c in zip(user_t, cores)]

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.set_xlabel('Number of Cores')
    ax1.set_ylabel('Time (seconds)')

    # Plot Real Time
    ax1.plot(cores, real_t, 'o-', label='Real Time (Wall clock)', color='tab:blue')
    # Plot User Time / Cores
    #ax1.plot(cores, user_per_core, 's--', label='User Time / Cores (Avg Workload)', color='tab:green')

    # Secondo asse per lo User Time totale (che ha scala molto diversa)
    # ax2 = ax1.twinx()
    # ax2.set_ylabel('Total User Time (s)', color='tab:red')
    # ax2.plot(cores, user_t, '^-.', label='Total User Time', color='tab:red')
    # ax2.tick_params(axis='y', labelcolor='tab:red')

    plt.title('Speedup analysis for multicell configuration')
    ax1.legend(loc='upper left')
    # ax2.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


    # Salvataggio in PNG
    output_name = "speedup_plot.png"
    plt.savefig(output_name, dpi=300)
    print(f"Grafico salvato correttamente come: {output_name}")

if __name__ == "__main__":
    main()
