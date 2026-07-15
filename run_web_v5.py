from pathlib import Path
import argparse,subprocess,sys
ROOT=Path(__file__).resolve().parent
def run(s):subprocess.run([sys.executable,str(ROOT/s),"--root",str(ROOT)],check=True)
def main():
 p=argparse.ArgumentParser();p.add_argument("--skip-observer",action="store_true");p.add_argument("--skip-control",action="store_true");a=p.parse_args()
 run("identification_v4_primary/select_primary_model.py")
 run("digitization/digitize_external_figures.py")
 run("experiments/run_distributed_model_validation.py")
 if not a.skip_observer:
  run("experiments/run_observer_decisive_validation.py");run("experiments/finalize_observer_decisive.py")
 if not a.skip_control:run("experiments/run_control_fair_comparison.py")
if __name__=="__main__":main()
