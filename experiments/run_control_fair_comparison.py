"""Equal-resource PI/MPC comparison and hotspot-energy-switch Pareto audit."""
from __future__ import annotations
import argparse,json,itertools
from dataclasses import replace
from pathlib import Path
import sys
import numpy as np,pandas as pd,yaml
import matplotlib.pyplot as plt
PROJECT_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(PROJECT_ROOT))
from controllers.hotspot_mpc import HotspotMPC
from controllers.pi_smc import PISMCController
from models.five_node_model import FiveNodeParameters,interpolate_three_to_five,step_five_node
from models.model_registry import load_primary_model,write_model_runtime_record
from models.three_node_model import load_frozen_parameters
from observers.ekf import FastExtendedKalmanFilter

def variant(base,capacity=1,cooling=1,conduction=1,offset=0):return replace(base,total_heat_capacity_J_K=base.total_heat_capacity_J_K*capacity,cooling_coefficient_per_node_W_K=base.cooling_coefficient_per_node_W_K*cooling,conduction_edge_W_K=base.conduction_edge_W_K*conduction,coolant_forward_C=base.coolant_forward_C+offset)
def simulate_case(name,true_model,control_model,current,preview,initial,sensors,sigma,seed,dt,heat,reference,limit,mpc_cfg=None):
 rng=np.random.default_rng(seed);noise=rng.normal(0,sigma,(len(current),len(sensors)));T=np.empty((len(current),5));T[0]=initial;estimate=np.full_like(T,np.nan);duty=np.zeros(len(current));direction=np.full(len(current),-1,int);residual=np.zeros(len(current));observer=None
 if name!='AuthorMeasured-PI-SMC':
  observer=FastExtendedKalmanFilter(initial+np.linspace(.25,-.25,5),np.eye(5)*.8,np.eye(5)*.004,sigma**2,sensors,control_model,heat);observer.update(T[0,list(sensors)]+noise[0]);estimate[0]=observer.state_C
 if name in ('AuthorMeasured-PI-SMC','SparseSensor-PI-SMC'):
  controller=PISMCController(control_model.controller_parameters.Kp_1_K,control_model.controller_parameters.Ki_1_Ks,reference_C=reference,feedback_mode='middle',maximum_duty=limit,update_period_s=10)
 elif name in ('MeanTemperature-MPC','Hotspot-MPC'):
  c=mpc_cfg or {};controller=HotspotMPC(control_model,heat,reference_C=reference,maximum_duty=limit,control_dt_s=dt,prediction_dt_s=dt,horizon_steps=6,minimum_direction_dwell_s=c.get('dwell',30),direction_switch_weight=c.get('switch_weight',4),fan_weight=c.get('fan_weight',.35),hotspot_weight=c.get('hotspot_weight',8 if name=='Hotspot-MPC' else 0),mean_temperature_weight=c.get('mean_weight',8 if name=='MeanTemperature-MPC' else 0),gradient_weight=c.get('gradient_weight',2),maximum_duty_step=c.get('maximum_duty_step',.15),duty_move_weight=c.get('duty_move_weight',.8),reversal_dead_time_s=dt)
 else: raise ValueError(name)
 for k in range(len(current)-1):
  if name=='AuthorMeasured-PI-SMC':cmd,flow=controller.command(T[k,[0,2,4]],k*dt,dt)
  elif name=='SparseSensor-PI-SMC':cmd,flow=controller.command(observer.state_C,k*dt,dt)
  else:
   pv=preview[k:k+controller.horizon_steps];cmd,flow,_=controller.command(observer.state_C,pv)
  result=step_five_node(T[k],current[k],min(cmd,limit),int(flow),dt,true_model,heat);T[k+1]=result.next_temperature_C;duty[k]=min(cmd,limit);direction[k]=flow;residual[k]=result.energy_residual_W
  if observer is not None:observer.predict(preview[k],duty[k],int(flow),dt);observer.update(T[k+1,list(sensors)]+noise[k+1]);estimate[k+1]=observer.state_C
 duty[-1]=duty[-2];direction[-1]=direction[-2];hot=T.max(1);grad=T.max(1)-T.min(1);mean=T.mean(1);excess=np.maximum(hot-reference,0);switches=int(np.sum(direction[1:]!=direction[:-1]));metrics={'controller':name,'maximum_hotspot_C':float(hot.max()),'maximum_hotspot_excess_C':float(excess.max()),'hotspot_excess_rmse_C':float(np.sqrt(np.mean(excess**2))),'gradient_rmse_C':float(np.sqrt(np.mean(grad**2))),'fan_energy_proxy_duty2_h':float(np.sum(duty**2)*dt/3600),'direction_switch_count':switches,'mean_temperature_tracking_rmse_C':float(np.sqrt(np.mean((mean-reference)**2))),'maximum_energy_residual_W':float(np.abs(residual).max()),'reference_unattainable_at_saturation_fraction':float(np.mean((hot>reference)&(duty>=limit-1e-12)))}
 return pd.DataFrame({'time_s':np.arange(len(current))*dt,'current_A':current,'hotspot_C':hot,'gradient_C':grad,'mean_C':mean,'duty':duty,'direction':direction}),metrics

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=PROJECT_ROOT);ap.add_argument('--output-dir',type=Path);a=ap.parse_args();root=a.root.resolve();outroot=(a.output_dir or root/'outputs_v5').resolve();out=outroot/'control_fair';out.mkdir(parents=True,exist_ok=True);primary,ppath,phash=load_primary_model(root);write_model_runtime_record(root,out,'control_fair');cfg=yaml.safe_load((root/'configs/reproduction_config.yaml').read_text());thermal=load_frozen_parameters(ppath);heat=str(primary['heat_model']);five=FiveNodeParameters.from_three_node(thermal,np.array(cfg['model']['coolant_node_temperatures_forward_C']));seed=int(cfg['project']['seed']);reference=55.;dt=float(cfg['control_v5']['evaluation_dt_s']);stride=int(dt);d=root/'outputs_v4/digitization/regenerated';f15=pd.read_csv(d/'Fig15_digitized.csv');f16=pd.read_csv(d/'Fig16_digitized.csv');base=np.interp(f16.time_s,f15.time_s,f15.paper_simulation_current_A)[::stride];initial=interpolate_three_to_five(f16.loc[0,['paper_simulation_T2_C','paper_simulation_middle_C','paper_simulation_T1_C']].to_numpy(float));sensors=(0,4)
 scenarios={'nominal':(five,base,base,.12,1,True),'critical':(five,np.minimum(base*1.25,42),np.minimum(base*1.25,42),.12,1,True),'mismatch':(variant(five,.88,.82,1.2),base,base,.12,1,True),'noise':(five,base,base,.30,1,True),'preview_minus10':(five,base,base*.9,.12,1,True),'ambient_plus5':(variant(five,offset=5),base,base,.12,1,True),'actuator_limited':(variant(five,cooling=.65),np.minimum(base*1.75,42),np.minimum(base*1.75,42),.12,.15,False)}
 # Pareto sweep on the critical scenario only.  18 deterministic combinations span all requested weights and dwell times.
 allcomb=list(itertools.product(cfg['control_v5']['minimum_dwell_s'],cfg['control_v5']['direction_switch_weight'],cfg['control_v5']['fan_weight'],cfg['control_v5']['hotspot_weight'],cfg['control_v5']['gradient_weight']));indices=np.linspace(0,len(allcomb)-1,int(cfg['control_v5']['pareto_max_configs'])).astype(int);configs=[]
 for i in indices:
  dwell,sw,fan,hot,grad=allcomb[i];configs.append({'config_id':len(configs),'dwell':dwell,'switch_weight':sw,'fan_weight':fan,'hotspot_weight':hot,'gradient_weight':grad})
 true,current,preview,sigma,limit,_=scenarios['critical'];pareto=[]
 _,baseline=simulate_case('AuthorMeasured-PI-SMC',true,five,current,preview,initial,sensors,sigma,seed,dt,heat,reference,limit)
 for c in configs:
  _,m=simulate_case('Hotspot-MPC',true,five,current,preview,initial,sensors,sigma,seed,dt,heat,reference,limit,c);pareto.append({**c,**m})
 pf=pd.DataFrame(pareto);pf.to_csv(out/'hotspot_mpc_pareto_sweep.csv',index=False)
 # Resource selections.
 equal_energy=pf[(pf.fan_energy_proxy_duty2_h/baseline['fan_energy_proxy_duty2_h']-1).abs()<=.05]
 equal_energy_match_found=not equal_energy.empty
 if equal_energy.empty:equal_energy=pf.iloc[[(pf.fan_energy_proxy_duty2_h-baseline['fan_energy_proxy_duty2_h']).abs().argmin()]]
 selected_energy=equal_energy.sort_values(['hotspot_excess_rmse_C','gradient_rmse_C']).iloc[0]
 selected_switch=pf.iloc[(pf.direction_switch_count-baseline['direction_switch_count']).abs().argmin()]
 selected_original=pf.sort_values(['hotspot_excess_rmse_C','gradient_rmse_C']).iloc[0]
 selections={'original_objective':selected_original.to_dict(),'equal_fan_energy':selected_energy.to_dict(),'equal_switch_budget':selected_switch.to_dict()};(out/'selected_engineering_points.json').write_text(json.dumps(selections,ensure_ascii=False,indent=2)+'\n')
 # Full scenario comparison uses equal-energy hotspot MPC and a mean-temperature MPC with matched resources.
 rows=[]
 trajectories={}
 ecfg={k:selections['equal_fan_energy'][k] for k in ['dwell','switch_weight','fan_weight','hotspot_weight','gradient_weight']};mcfg={**ecfg,'hotspot_weight':0,'mean_weight':8}
 for si,(scenario,(true,current,preview,sigma,limit,feasible)) in enumerate(scenarios.items()):
  for name,c in [('AuthorMeasured-PI-SMC',None),('SparseSensor-PI-SMC',None),('MeanTemperature-MPC',mcfg),('Hotspot-MPC',ecfg)]:
   frame,m=simulate_case(name,true,five,current,preview,initial,sensors,sigma,seed+si,dt,heat,reference,limit,c);rows.append({'scenario':scenario,'feasible':feasible,**m});trajectories[(scenario,name)]=frame
 metrics=pd.DataFrame(rows);metrics.to_csv(out/'controller_comparison_metrics.csv',index=False);agg=metrics[metrics.feasible].groupby('controller',as_index=False).agg(mean_max_hotspot_excess_C=('maximum_hotspot_excess_C','mean'),mean_hotspot_excess_rmse_C=('hotspot_excess_rmse_C','mean'),mean_gradient_rmse_C=('gradient_rmse_C','mean'),mean_fan_energy=('fan_energy_proxy_duty2_h','mean'),total_switches=('direction_switch_count','sum'),mean_tracking_rmse_C=('mean_temperature_tracking_rmse_C','mean'));agg.to_csv(out/'controller_aggregate_feasible.csv',index=False)
 author=agg[agg.controller=='AuthorMeasured-PI-SMC'].iloc[0];hot=agg[agg.controller=='Hotspot-MPC'].iloc[0];summary={'primary_model_sha256':phash,'comparison_resource_point':('equal_fan_energy_within_5pct' if equal_energy_match_found else 'nearest_fan_energy_no_match'),'equal_energy_relative_difference_percent':100*(selected_energy.fan_energy_proxy_duty2_h/baseline['fan_energy_proxy_duty2_h']-1),'equal_switch_count_difference':int(selected_switch.direction_switch_count-baseline['direction_switch_count']),'hotspot_mpc_vs_author':{'maximum_excess_improvement_percent':100*(author.mean_max_hotspot_excess_C-hot.mean_max_hotspot_excess_C)/max(author.mean_max_hotspot_excess_C,1e-12),'excess_rmse_improvement_percent':100*(author.mean_hotspot_excess_rmse_C-hot.mean_hotspot_excess_rmse_C)/author.mean_hotspot_excess_rmse_C,'gradient_rmse_improvement_percent':100*(author.mean_gradient_rmse_C-hot.mean_gradient_rmse_C)/author.mean_gradient_rmse_C,'fan_energy_change_percent':100*(hot.mean_fan_energy/author.mean_fan_energy-1),'switch_count_change':int(hot.total_switches-author.total_switches),'mean_tracking_rmse_change_percent':100*(hot.mean_tracking_rmse_C/author.mean_tracking_rmse_C-1)},'equal_energy_match_found_within_5pct':bool(equal_energy_match_found),'hotspot_MPC_equal_energy_advantage_holds':bool(equal_energy_match_found and hot.mean_hotspot_excess_rmse_C<author.mean_hotspot_excess_rmse_C),'hotspot_MPC_equal_switch_advantage_on_critical':bool(selected_switch.hotspot_excess_rmse_C<baseline['hotspot_excess_rmse_C']),'actuator_limited_excluded_from_feasible_aggregate':True,'safe_temperature_configured':False}
 (out/'control_fair_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
 fig,ax=plt.subplots(figsize=(7,5));sc=ax.scatter(pf.fan_energy_proxy_duty2_h,pf.hotspot_excess_rmse_C,c=pf.direction_switch_count);ax.scatter([baseline['fan_energy_proxy_duty2_h']],[baseline['hotspot_excess_rmse_C']],marker='x',s=100,label='PI-SMC');ax.set(xlabel='Fan energy proxy duty²·h',ylabel='Hotspot excess RMSE (°C)',title='Hotspot MPC resource Pareto sweep');ax.legend();fig.colorbar(sc,ax=ax,label='Direction switches');fig.tight_layout();fig.savefig(out/'pareto_energy_hotspot_switches.png',dpi=190);plt.close(fig);print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
