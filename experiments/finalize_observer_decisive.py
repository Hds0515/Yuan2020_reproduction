"""Finalize observer statistics after the 100x10 off-grid batch is complete."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import sys
import numpy as np,pandas as pd,yaml
PROJECT_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(PROJECT_ROOT))
from experiments.run_observer_decisive_validation import simulate_truth,run_kf,run_bank,metric
from models.model_registry import load_primary_model
from models.three_node_model import load_frozen_parameters,simulate
from models.five_node_model import FiveNodeParameters,interpolate_three_to_five

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=PROJECT_ROOT);ap.add_argument('--output-dir',type=Path);a=ap.parse_args();root=a.root.resolve();outroot=(a.output_dir or root/'outputs_v5').resolve();out=outroot/'observer_decisive'
 primary,ppath,phash=load_primary_model(root);cfg=yaml.safe_load((root/'configs/reproduction_config.yaml').read_text());seed=int(cfg['project']['seed']);heat=str(primary['heat_model']);thermal=load_frozen_parameters(ppath);coolant3=np.array(cfg['model']['coolant_node_temperatures_forward_C'],float);base=FiveNodeParameters.from_three_node(thermal,coolant3)
 d=root/'outputs_v4/digitization/regenerated';f15=pd.read_csv(d/'Fig15_digitized.csv');f16=pd.read_csv(d/'Fig16_digitized.csv');cur1=np.interp(f16.time_s,f15.time_s,f15.paper_simulation_current_A);init3=f16.loc[0,['paper_simulation_T2_C','paper_simulation_middle_C','paper_simulation_T1_C']].to_numpy(float);three=simulate(thermal,cur1,init3,dt_s=1,bidirectional=True,initial_airflow_direction=-1,fan_initially_enabled=True,fan_enable_temperature_C=53,reference_temperature_C=55,smc_deadband_C=.5,smc_update_period_s=10,heat_model=heat,coolant_node_temperatures_forward_C=coolant3)
 dt=float(cfg['observer_v5']['evaluation_dt_s']);stride=int(dt);current=cur1[::stride];duty=three.fan_duty[::stride];direction=three.airflow_direction[::stride];initial=interpolate_three_to_five(init3);sigma=.12;sensors=tuple(int(v)-1 for v in pd.read_csv(out/'sensor_design_only_metrics.csv').sort_values(['design_mean_hotspot_rmse_C','design_worst_hotspot_rmse_C']).iloc[0].sensor_nodes_1_based.split(','));bank=np.array([[1,1,1,1,0],[.85,1,1,1,0],[1.15,1,1,1,0],[1,.75,1,1,0],[1,1,1,.85,0]],float)
 dr=[]
 for idx,name in [(1,'capacity_low'),(2,'capacity_high'),(3,'cooling_low'),(4,'input_low')]:
  sc=bank[idx];truth=simulate_truth(base,sc,current,duty,direction,dt,heat,initial)
  for ns in range(3):
   rng=np.random.default_rng(seed+300000+idx*100+ns);meas=truth[:,sensors]+rng.normal(0,sigma,(len(truth),2));ini=initial+np.linspace(.25,-.25,5);single=run_kf(base,bank[0],1,sensors,meas,current,duty,direction,dt,heat,ini,sigma);exact,prob=run_bank(base,bank,sensors,meas,current,duty,direction,dt,heat,ini,sigma);keep=np.arange(len(bank))!=idx;loo,_=run_bank(base,bank[keep],sensors,meas,current,duty,direction,dt,heat,ini,sigma)
   for typ,est in [('single',single),('exact_in_bank',exact),('leave_one_out',loo)]:dr.append({'scenario':name,'noise_seed':ns,'test':typ,**metric(truth,est),'exact_model_final_probability':prob[-1,idx] if typ=='exact_in_bank' else np.nan})
 pd.DataFrame(dr).to_csv(out/'exact_and_leave_one_out.csv',index=False)
 audit=[];nomtruth=simulate_truth(base,bank[0],current,duty,direction,dt,heat,initial)
 for sig,label in [(.12,'nominal'),(.30,'pure_measurement_noise')]:
  for ns in range(5):
   rng=np.random.default_rng(seed+500000+ns+(100 if sig>.2 else 0));meas=nomtruth[:,sensors]+rng.normal(0,sig,(len(nomtruth),2));ini=initial+np.linspace(.25,-.25,5);sm=metric(nomtruth,run_kf(base,bank[0],1,sensors,meas,current,duty,direction,dt,heat,ini,sig));bm=metric(nomtruth,run_bank(base,bank,sensors,meas,current,duty,direction,dt,heat,ini,sig)[0]);audit.append({'scenario':label,'seed':ns,'single':sm['hotspot_rmse_C'],'bank':bm['hotspot_rmse_C'],'bank_loss_percent':100*(bm['hotspot_rmse_C']-sm['hotspot_rmse_C'])/sm['hotspot_rmse_C']})
 af=pd.DataFrame(audit);af.to_csv(out/'nominal_and_noise_loss.csv',index=False);off=pd.read_csv(out/'continuous_off_grid_100x10.csv');imps=off.hotspot_improvement_percent.to_numpy();rng=np.random.default_rng(seed+88);boot=[]
 for _ in range(5000):boot.append(float(np.mean(rng.choice(imps,size=len(imps),replace=True))))
 ci=[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))];summary={'primary_model_sha256':phash,'selected_sensor_nodes_1_based':[v+1 for v in sensors],'sensor_selection_used_test_set':False,'off_grid_object_count':int(off.object_index.nunique()),'noise_seed_count_per_object':int(off.noise_seed.nunique()),'paired_sample_count':len(off),'median_hotspot_rmse_improvement_percent':float(np.median(imps)),'mean_hotspot_rmse_improvement_percent':float(np.mean(imps)),'worst_hotspot_rmse_improvement_percent':float(np.min(imps)),'q05_hotspot_improvement_percent':float(np.quantile(imps,.05)),'improved_sample_fraction':float(np.mean(imps>0)),'paired_mean_improvement_95pct_CI_percent':ci,'nominal_median_performance_loss_percent':float(af[af.scenario=='nominal'].bank_loss_percent.median()),'pure_noise_median_performance_loss_percent':float(af[af.scenario=='pure_measurement_noise'].bank_loss_percent.median()),'linearization_note':'Heat generation evaluated at 55 C; cooling/conduction dynamics exact in T_i>Tcoolant regime.','exact_and_leave_one_out_noise_seeds':3}
 summary['core_innovation_acceptance']={'off_grid_majority_improves':summary['improved_sample_fraction']>.5,'paired_CI_lower_above_zero':ci[0]>0,'median_improvement_above_10pct':summary['median_hotspot_rmse_improvement_percent']>10,'nominal_and_noise_loss_below_5pct':max(summary['nominal_median_performance_loss_percent'],summary['pure_noise_median_performance_loss_percent'])<5};summary['MM_EKF_core_innovation_supported']=all(summary['core_innovation_acceptance'].values());(out/'observer_decisive_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n');print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
