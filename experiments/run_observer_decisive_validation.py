"""Decisive parameter-model-bank validation on 100 off-grid plants x 10 seeds.

For tractable Monte Carlo execution, the five-node model is represented by its
exact linear form for the operating regime T_i > T_coolant,i.  The only omitted
nonlinearity is the very small quadratic temperature term in the polarization
proxy; heat generation is evaluated at the 55 °C reference.  The same
linearization is used by truth and observers, while plant parameters and input
bias remain off-grid.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import qmc
import yaml
PROJECT_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(PROJECT_ROOT))
from models.model_registry import load_primary_model, write_model_runtime_record
from models.three_node_model import load_frozen_parameters, heat_generation_W, simulate
from models.five_node_model import FiveNodeParameters, interpolate_three_to_five


def lhs(count,seed,ranges):
    return qmc.scale(qmc.LatinHypercube(d=5,seed=seed).random(count),ranges[:,0],ranges[:,1])

def system_matrices(base, scales, current, duty, direction, dt, heat_model):
    c_scale,kc_scale,kn_scale,input_scale,offset=scales
    C=base.total_heat_capacity_J_K*c_scale/5
    kc=base.cooling_coefficient_per_node_W_K*kc_scale
    kn=base.conduction_edge_W_K*kn_scale
    coolant=(base.coolant_forward_C+offset) if direction==1 else (base.coolant_forward_C+offset)[::-1]
    A=np.zeros((5,5))
    for i in range(4):
        A[i,i]-=kn;A[i,i+1]+=kn;A[i+1,i]+=kn;A[i+1,i+1]-=kn
    A[np.arange(5),np.arange(5)]-=kc*duty
    F=np.eye(5)+dt*A/C
    q=heat_generation_W(current*input_scale,55.0,base.controller_parameters,heat_model)/5
    b=dt*(np.full(5,q)+kc*duty*coolant)/C
    return F,b

def simulate_truth(base,scales,current,duty,direction,dt,heat_model,initial):
    x=np.empty((len(current),5));x[0]=initial
    for k in range(len(current)-1):
        F,b=system_matrices(base,scales,current[k],duty[k],int(direction[k]),dt,heat_model);x[k+1]=F@x[k]+b
    return x

def run_kf(base,model_scales,input_scale_override,sensors,measurements,current,duty,direction,dt,heat_model,initial,sigma):
    x=initial.copy();P=np.eye(5)*.8;Q=np.eye(5)*.004;R=np.eye(len(sensors))*sigma**2;H=np.zeros((len(sensors),5));H[np.arange(len(sensors)),sensors]=1
    out=np.empty((len(current),5));out[0]=x
    for k in range(len(current)-1):
        scales=model_scales.copy();scales[3]=input_scale_override
        F,b=system_matrices(base,scales,current[k],duty[k],int(direction[k]),dt,heat_model);x=F@x+b;P=F@P@F.T+Q
        innovation=measurements[k+1]-H@x;S=H@P@H.T+R;K=P@H.T@np.linalg.inv(S);x=x+K@innovation;P=(np.eye(5)-K@H)@P
        out[k+1]=x
    return out

def run_bank(base,bank_scales,sensors,measurements,current,duty,direction,dt,heat_model,initial,sigma):
    m=len(bank_scales);xs=np.repeat(initial[None,:],m,axis=0);Ps=np.repeat((np.eye(5)*.8)[None,:,:],m,axis=0);prob=np.full(m,1/m);Q=np.eye(5)*.004;R=np.eye(len(sensors))*sigma**2;H=np.zeros((len(sensors),5));H[np.arange(len(sensors)),sensors]=1
    estimate=np.empty((len(current),5));history=np.empty((len(current),m));estimate[0]=initial;history[0]=prob
    for k in range(len(current)-1):
        logs=np.empty(m)
        for j,scales in enumerate(bank_scales):
            F,b=system_matrices(base,scales,current[k],duty[k],int(direction[k]),dt,heat_model);xs[j]=F@xs[j]+b;Ps[j]=F@Ps[j]@F.T+Q
            innovation=measurements[k+1]-H@xs[j];S=H@Ps[j]@H.T+R;K=Ps[j]@H.T@np.linalg.inv(S);xs[j]=xs[j]+K@innovation;Ps[j]=(np.eye(5)-K@H)@Ps[j]
            sign,ld=np.linalg.slogdet(S);logs[j]=-.5*(innovation@np.linalg.solve(S,innovation)+ld+len(sensors)*np.log(2*np.pi))
        lw=np.log(np.maximum(prob,1e-6))+logs;lw-=lw.max();prob=np.exp(lw);prob=np.maximum(prob/prob.sum(),1e-4);prob/=prob.sum();estimate[k+1]=prob@xs;history[k+1]=prob
    return estimate,history

def metric(truth,estimate):
    e=estimate-truth;he=estimate.max(1)-truth.max(1)
    return {'hotspot_rmse_C':float(np.sqrt(np.mean(he**2))),'mean_node_rmse_C':float(np.sqrt(np.mean(e**2))),'maximum_absolute_error_C':float(np.abs(e).max()),'q95_absolute_error_C':float(np.quantile(np.abs(e),.95)),'hotspot_position_accuracy':float(np.mean(estimate.argmax(1)==truth.argmax(1)))}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=PROJECT_ROOT);ap.add_argument('--output-dir',type=Path);a=ap.parse_args();root=a.root.resolve();outroot=(a.output_dir or root/'outputs_v5').resolve();out=outroot/'observer_decisive';out.mkdir(parents=True,exist_ok=True)
    primary,ppath,phash=load_primary_model(root);write_model_runtime_record(root,out,'observer_decisive');cfg=yaml.safe_load((root/'configs/reproduction_config.yaml').read_text());settings=cfg['observer_v5'];seed=int(cfg['project']['seed']);heat=str(primary['heat_model']);thermal=load_frozen_parameters(ppath);coolant3=np.array(cfg['model']['coolant_node_temperatures_forward_C'],float);base=FiveNodeParameters.from_three_node(thermal,coolant3)
    d=root/'outputs_v4/digitization/regenerated';f15=pd.read_csv(d/'Fig15_digitized.csv');f16=pd.read_csv(d/'Fig16_digitized.csv');cur1=np.interp(f16.time_s,f15.time_s,f15.paper_simulation_current_A);init3=f16.loc[0,['paper_simulation_T2_C','paper_simulation_middle_C','paper_simulation_T1_C']].to_numpy(float);three=simulate(thermal,cur1,init3,dt_s=1,bidirectional=True,initial_airflow_direction=-1,fan_initially_enabled=True,fan_enable_temperature_C=53,reference_temperature_C=55,smc_deadband_C=.5,smc_update_period_s=10,heat_model=heat,coolant_node_temperatures_forward_C=coolant3)
    stride=int(settings['evaluation_dt_s']);dt=float(stride);current=cur1[::stride];duty=three.fan_duty[::stride];direction=three.airflow_direction[::stride];initial=interpolate_three_to_five(init3);r=settings['parameter_ranges'];ranges=np.array([r['C_th_scale'],r['K_cool_scale'],r['K_node_scale'],r['load_input_scale'],r['coolant_temperature_offset_C']],float);sigma=.12
    design=lhs(int(settings['lhs_design_objects']),seed+10,ranges);layout=[]
    for sensors in [(i,j) for i in range(5) for j in range(i+1,5)]:
        vals=[]
        for z,sc in enumerate(design):
            truth=simulate_truth(base,sc,current,duty,direction,dt,heat,initial);rng=np.random.default_rng(seed+1000+z);meas=truth[:,sensors]+rng.normal(0,sigma,(len(truth),2));est=run_kf(base,np.array([1,1,1,1,0.]),1.0,sensors,meas,current,duty,direction,dt,heat,initial+np.linspace(.25,-.25,5),sigma);vals.append(metric(truth,est)['hotspot_rmse_C'])
        layout.append({'sensor_nodes_1_based':f'{sensors[0]+1},{sensors[1]+1}','design_mean_hotspot_rmse_C':np.mean(vals),'design_worst_hotspot_rmse_C':np.max(vals)})
    lf=pd.DataFrame(layout).sort_values(['design_mean_hotspot_rmse_C','design_worst_hotspot_rmse_C']);lf.to_csv(out/'sensor_design_only_metrics.csv',index=False);sensors=tuple(int(v)-1 for v in lf.iloc[0].sensor_nodes_1_based.split(','))
    bank=np.array([[1,1,1,1,0],[.85,1,1,1,0],[1.15,1,1,1,0],[1,.75,1,1,0],[1,1,1,.85,0]],float);test=lhs(int(settings['lhs_test_objects']),seed+20,ranges);rows=[];prows=[];start=time.perf_counter()
    for oi,sc in enumerate(test):
        truth=simulate_truth(base,sc,current,duty,direction,dt,heat,initial)
        for ns in range(int(settings['noise_seeds_per_test_object'])):
            rng=np.random.default_rng(seed+100000+oi*100+ns);meas=truth[:,sensors]+rng.normal(0,sigma,(len(truth),2));ini=initial+np.linspace(.25,-.25,5);single=run_kf(base,np.array([1,1,1,1,0.]),1,sensors,meas,current,duty,direction,dt,heat,ini,sigma);mm,prob=run_bank(base,bank,sensors,meas,current,duty,direction,dt,heat,ini,sigma);sm=metric(truth,single);bm=metric(truth,mm);imp=100*(sm['hotspot_rmse_C']-bm['hotspot_rmse_C'])/sm['hotspot_rmse_C'];rows.append({'object_index':oi,'noise_seed':ns,'C_th_scale':sc[0],'K_cool_scale':sc[1],'K_node_scale':sc[2],'load_input_scale':sc[3],'coolant_offset_C':sc[4],**{f'single_{k}':v for k,v in sm.items()},**{f'bank_{k}':v for k,v in bm.items()},'hotspot_improvement_percent':imp});prows.append({'object_index':oi,'noise_seed':ns,'final_max_probability':prob[-1].max(),'mean_entropy':np.mean(-np.sum(prob*np.log(np.maximum(prob,1e-12)),1))})
    runtime=time.perf_counter()-start;off=pd.DataFrame(rows);off.to_csv(out/'continuous_off_grid_100x10.csv',index=False);pd.DataFrame(prows).to_csv(out/'off_grid_probability_calibration.csv',index=False)
    # exact-in-bank and leave-one-out
    dr=[]
    for idx,name in [(1,'capacity_low'),(2,'capacity_high'),(3,'cooling_low'),(4,'input_low')]:
        sc=bank[idx];truth=simulate_truth(base,sc,current,duty,direction,dt,heat,initial)
        for ns in range(10):
            rng=np.random.default_rng(seed+300000+idx*100+ns);meas=truth[:,sensors]+rng.normal(0,sigma,(len(truth),2));ini=initial+np.linspace(.25,-.25,5);single=run_kf(base,bank[0],1,sensors,meas,current,duty,direction,dt,heat,ini,sigma);exact,prob=run_bank(base,bank,sensors,meas,current,duty,direction,dt,heat,ini,sigma);keep=np.arange(len(bank))!=idx;loo,_=run_bank(base,bank[keep],sensors,meas,current,duty,direction,dt,heat,ini,sigma)
            for typ,est in [('single',single),('exact_in_bank',exact),('leave_one_out',loo)]:dr.append({'scenario':name,'noise_seed':ns,'test':typ,**metric(truth,est),'exact_model_final_probability':prob[-1,idx] if typ=='exact_in_bank' else np.nan})
    pd.DataFrame(dr).to_csv(out/'exact_and_leave_one_out.csv',index=False)
    audit=[];nomtruth=simulate_truth(base,bank[0],current,duty,direction,dt,heat,initial)
    for sig,label in [(.12,'nominal'),(.30,'pure_measurement_noise')]:
        for ns in range(20):
            rng=np.random.default_rng(seed+500000+ns+(100 if sig>.2 else 0));meas=nomtruth[:,sensors]+rng.normal(0,sig,(len(nomtruth),2));ini=initial+np.linspace(.25,-.25,5);sm=metric(nomtruth,run_kf(base,bank[0],1,sensors,meas,current,duty,direction,dt,heat,ini,sig));bm=metric(nomtruth,run_bank(base,bank,sensors,meas,current,duty,direction,dt,heat,ini,sig)[0]);audit.append({'scenario':label,'seed':ns,'single':sm['hotspot_rmse_C'],'bank':bm['hotspot_rmse_C'],'bank_loss_percent':100*(bm['hotspot_rmse_C']-sm['hotspot_rmse_C'])/sm['hotspot_rmse_C']})
    af=pd.DataFrame(audit);af.to_csv(out/'nominal_and_noise_loss.csv',index=False);imps=off.hotspot_improvement_percent.to_numpy();rng=np.random.default_rng(seed+88);boot=np.mean(imps[rng.integers(0,len(imps),(5000,len(imps)))],1);ci=[np.quantile(boot,.025),np.quantile(boot,.975)];summary={'primary_model_sha256':phash,'selected_sensor_nodes_1_based':[v+1 for v in sensors],'sensor_selection_used_test_set':False,'off_grid_object_count':len(test),'noise_seed_count_per_object':int(settings['noise_seeds_per_test_object']),'paired_sample_count':len(off),'median_hotspot_rmse_improvement_percent':float(np.median(imps)),'mean_hotspot_rmse_improvement_percent':float(np.mean(imps)),'worst_hotspot_rmse_improvement_percent':float(np.min(imps)),'q05_hotspot_improvement_percent':float(np.quantile(imps,.05)),'improved_sample_fraction':float(np.mean(imps>0)),'paired_mean_improvement_95pct_CI_percent':[float(ci[0]),float(ci[1])],'nominal_median_performance_loss_percent':float(af[af.scenario=='nominal'].bank_loss_percent.median()),'pure_noise_median_performance_loss_percent':float(af[af.scenario=='pure_measurement_noise'].bank_loss_percent.median()),'runtime_s':runtime,'linearization_note':'Heat generation evaluated at 55 C; cooling/conduction dynamics exact in T_i>Tcoolant regime.'};summary['core_innovation_acceptance']={'off_grid_majority_improves':summary['improved_sample_fraction']>.5,'paired_CI_lower_above_zero':ci[0]>0,'median_improvement_above_10pct':summary['median_hotspot_rmse_improvement_percent']>10,'nominal_and_noise_loss_below_5pct':max(summary['nominal_median_performance_loss_percent'],summary['pure_noise_median_performance_loss_percent'])<5};summary['MM_EKF_core_innovation_supported']=all(summary['core_innovation_acceptance'].values());(out/'observer_decisive_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n');print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
