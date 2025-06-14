import subprocess
import psutil
import time
from tabulate import tabulate

PEARSON_BINARY = "./pearson"
DATA_PATH = "./data/"
DATA_FILES = ["128.data", "256.data", "512.data", "1024.data"]

def monitor_process(pid):
    p = psutil.Process(pid)
    cpu_percents = []
    memory_usage = []
    max_mem_rss =max_mem_vms = 0
    stats_before = psutil.cpu_stats()
    initial_io = psutil.disk_io_counters()

    # Poll until process ends
    while p.is_running():
        try:
            cpu = p.cpu_percent(interval=0.1)
            mem_rss = p.memory_info().rss
            mem_vms = p.memory_info().vms
            cpu_percents.append(cpu)
            memory_usage.append(round(p.memory_percent(), 2))
            max_mem_rss = max(max_mem_rss, mem_rss)
            max_mem_vms = max(max_mem_vms, mem_vms)

        except psutil.NoSuchProcess:
            break

    stats_after = psutil.cpu_stats()

    final_io = psutil.disk_io_counters()
    read_bytes_total = (final_io.read_bytes - initial_io.read_bytes) / (1024*1024)
    write_bytes_total = (final_io.write_bytes - initial_io.write_bytes) / (1024*1024)

    total_syscalls = abs(stats_after.syscalls - stats_before.syscalls)
    ctx_switches = abs(stats_after.ctx_switches - stats_before.ctx_switches)

    avg_cpu = sum(cpu_percents) / len(cpu_percents) if cpu_percents else 0
    avg_mem = sum(memory_usage) / len(memory_usage) if memory_usage else 0
    max_mem_rss = max_mem_rss / (1024 * 1024)
    max_mem_vms = max_mem_vms / (1024 * 1024)
    return avg_cpu, ctx_switches, avg_mem, max_mem_rss, max_mem_vms, total_syscalls, read_bytes_total, write_bytes_total

def run_pearson_and_monitor(file):
    input_path = f"{DATA_PATH}{file}"
    output_file = f"output_{file.split('.')[0]}_pearson.data"
    output_path = f"{DATA_PATH}{output_file}"

    print(f"▶ Running pearson on {file}...")

    pearson_proc = subprocess.Popen([PEARSON_BINARY, input_path, output_path])

    start_time = time.time()

    # Monitor the process
    avg_cpu, ctx_switches, avg_mem, max_mem_rss, max_mem_vms, total_syscalls, read_bytes_total, write_bytes_total = monitor_process(pearson_proc.pid)

    pearson_proc.wait()

    end_time = time.time()

    elapsed = end_time - start_time

    return {
        'Input Size': file.split('.')[0],
        'Time (s)': f"{elapsed:.3f}",
        'CPU (%)': f"{avg_cpu:.1f}",
        'Context Switches Count': f"{ctx_switches:.1f}",
        'Memory (%)': f"{avg_mem:.3f}",
        'Max Memory RSS (MiB)': f"{max_mem_rss:.2f}",
        'Max Memory VMS (MiB)': f"{max_mem_vms:.2f}",
        'Total Syscalls': total_syscalls,
        'Read Bytes Total (MiB)': f"{read_bytes_total:.2f}",
        'Write Bytes Total (MiB)': f"{write_bytes_total:.2f}"
    }

def main():
    results = []

    for file in DATA_FILES:
        result = run_pearson_and_monitor(file)
        results.append(result)

    headers = [
        'Input Size', 'Time (s)', 'CPU (%)', 'Context Switches Count', 'Memory (%)', 'Max Memory RSS (MiB)', 'Max Memory VMS (MiB)', 'Total Syscalls', 'Read Bytes Total (MiB)', 'Write Bytes Total (MiB)'
    ]

    # Convert dicts to lists to enforce header order
    rows = []
    for d in results:
        rows.append([d[h] for h in headers])

    print("\n" + tabulate(rows, headers=headers, tablefmt="github"))

if __name__ == "__main__":
    main()
