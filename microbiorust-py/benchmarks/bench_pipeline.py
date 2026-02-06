import os
import time
import subprocess
from Bio import SeqIO
from microbiorust import gbk


#warnings.simplefilter('ignore', BioPythonWarning)

class PipelineSuite:
    """
    Benchmarks for microbiorust-py vs BioPython.
    Measures:
      - Parsing Time
      - Peak Memory
      - Latency
      - Energy Consumption (CodeCarbon) using loop iterations
    """

    timeout = 300  # 5 minutes per benchmark

    # Parameterizing creates 4 lines on one graph: Rust/Python × Interactive/Pipeline
    params = [['rust', 'python'], ['interactive', 'pipeline']]
    param_names = ['engine', 'context']

    def setup(self, engine, context):
        bench_dir = os.path.dirname(__file__)
        self.filepath = os.path.join(bench_dir, "Rhiz3841.gbk.gb")

        # CLI scripts for pipeline mode
        self.rust_cli = os.path.join(bench_dir, "rust_via_python_countgbk2faa.py")
        self.python_cli = os.path.join(bench_dir, "bp_equivalent_gbktofaacount.py")

        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"Missing benchmark file: {self.filepath}")

        # Initialize energy dictionary if not already present
        if not hasattr(self, "_energy_joules"):
            self._energy_joules = {}
 
    def _run_once(self, engine, context):
        if context == 'interactive':
            if engine == 'rust':
                return gbk.gbk_to_faa_count(self.filepath)
            else:
                count = 0
                for record in SeqIO.parse(self.filepath, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS":
                            _ = feature.extract(record.seq).translate()
                            count += 1
                return count
        else: # pipeline context
            cli_cmd = self.rust_cli if engine == 'rust' else self.python_cli
            return subprocess.run(
                ["python", cli_cmd, self.filepath],
                capture_output=True, check=True, timeout=self.timeout
            ).returncode

    # --- RUN REPEATEDLY for codecarbon otherwise the script was finishing before it could be measured ---
    def _run_repeatedly(self, engine, context, iterations):
        """Calls the logic many times to make energy measurable."""
        last_result = None
        for _ in range(iterations):
            last_result = self._run_once(engine, context)
        return last_result

    # --- CORE LOGIC WITH CODECARBON ---
    def track_energy(self, engine, context):
        """
        Routes execution based on engine/context and tracks energy with CodeCarbon.
        Stores last measured energy per engine in self._energy_joules.
        """
        is_ci = os.getenv('CI') or os.getenv('GITHUB_ACTIONS')
        start_cpu = time.process_time()
        start_wall = time.perf_counter()
        os.environ["CODECARBON_CARBON_INTENSITY"] = "475"
        if not is_ci:
            try:
              from codecarbon import OfflineEmissionsTracker
              tracker = OfflineEmissionsTracker(measure_power_secs=0.1, log_level="CRITICAL", country_iso_code="USA")
              tracker.start()
              iterations = 500 if engine == 'rust' else 50
              result = None
              try:
                   self._run_repeatedly(engine, context, iterations)
              finally:
                   tracker.stop()

              energy_kwh = getattr(tracker, "total_energy", 0)
            except Exception:
                pass
            # Store energy per engine in Joules
        else:
            self._run_repeatedly(engine, context, iterations)
        cpu_time = time.process_time() - start_cpu
        wall_time = time.perf_counter() - start_wall
        
        # Use CPU-based estimation if CodeCarbon didn't work or returned 0
        if energy_kwh == 0 or energy_kwh is None:
            # CPU utilization factor (how much of wall time was CPU work)
            cpu_utilization = min(cpu_time / wall_time, 1.0) if wall_time > 0 else 1.0
            
            # Power estimates for GitHub Actions runner (2-core VM)
            base_power_w = 10   # System idle power
            cpu_power_w = 30 * cpu_utilization  # CPU power under load
            avg_power_w = base_power_w + cpu_power_w
            
            # Energy in kWh
            energy_kwh = (avg_power_w * wall_time / 3600) / 1000
        
        # Convert to Joules per iteration
        energy_joules_per_iteration = (energy_kwh * 3_600_000) / iterations
        
        return energy_joules_per_iteration
    track_energy.unit = "J"

    # --- 1. PRIMARY TIME BENCHMARK (ASV automatic) ---
    def time_process_all(self, engine, context):
        """Measures parsing speed (seconds)."""
        self._run_once(engine, context)

    # --- 2. PEAK MEMORY BENCHMARK ---
    def peakmem_process_all(self, engine, context):
        """Measures maximum resident set size (RAM)."""
        self._run_once(engine, context)

    # --- 3. PARSING LATENCY (single-pass) ---
    def track_parsing_latency(self, engine, context):
        """Records single-pass parsing latency (seconds)."""
        start = time.perf_counter()
        self._run_once(engine, context)
        end = time.perf_counter()
        return end - start
    track_parsing_latency.unit = "s"

