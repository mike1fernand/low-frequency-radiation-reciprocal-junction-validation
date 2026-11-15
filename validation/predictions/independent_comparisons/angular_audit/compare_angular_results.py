from pathlib import Path
import json,csv,hashlib
import numpy as np
H=Path(__file__).resolve().parent;B=H.parent
sup=json.loads((B/'input/rev/angular_2400_600_96.json').read_text());run=json.loads((H/'angular_2400_600_96.json').read_text());coarse=json.loads((B/'input/rev/angular_300_150_12.json').read_text());impl=json.loads((H/'angular_implementation_results.json').read_text())
rows=list(csv.DictReader((B/'portable_checks/comparison_inputs/cone_results.csv').open()))
orig={(r['state'],r['direction_kind'],int(float(r['half_angle_deg']))):r for r in rows if r['side']=='forward' and float(r['half_angle_deg']) in [5,10,20,30]}
out=[];absolute=[];relative=[]
for model in sup:
 for key,v in sup[model].items():
  if isinstance(v,dict):
   for cone,x in v.items():
    y=run[model][key][cone];absolute.append(abs(x-y));relative.append(abs(x-y)/abs(x) if x else 0)
  else:
   absolute.append(abs(v-run[model][key]));relative.append(abs(v-run[model][key])/abs(v) if v else 0)
 for angle in [5,10,20,30]:
  old=orig[(model.replace('-','_'),'group',angle)];new=run[model]
  ts=new['qT'][str(angle)];ls=new['qL'][str(angle)]; cs=ts/(ts+ls)
  wn=new['qT_wn'][str(angle)]/(new['qT_wn'][str(angle)]+new['qL_wn'][str(angle)])
  out.append({'model':model,'half_angle_deg':angle,'uploaded_replayed_qT':ts,'uploaded_replayed_qL':ls,'qT_share_pct':100*cs,'wave_normal_qT_share_pct':100*wn,
   'current_qT':float(old['incident_qT']),'current_qL':float(old['incident_qL']),'current_share_pct':100*float(old['qT_share']),
   'qT_relative_difference':ts/float(old['incident_qT'])-1,'qL_relative_difference':ls/float(old['incident_qL'])-1,
   'share_percentage_point_difference':100*(cs-float(old['qT_share'])),
   'coarse_to_fine_qT_relative_change':ts/coarse[model]['qT'][str(angle)]-1,'coarse_to_fine_qL_relative_change':ls/coarse[model]['qL'][str(angle)]-1})
shell={}
for model in run:
 shell[model]={lab:run[model]['tot_'+lab]/impl['new_kernel_states'][model]['T'+lab]-1 for lab in ['qT','qL']}
res={'fresh_replay_exit_code':0,'fresh_replay_reported_runtime_s':109,'directions':288000,'max_abs_difference_from_supplied_fine':max(absolute),'max_relative_difference_from_supplied_fine':max(relative),'full_shell_relative_difference_from_solver':shell,'cone_comparison':out,'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [B/'input/rev/angular_check.py',B/'portable_checks/rev/angular_check.py',H/'angular_2400_600_96.json']}}
(H/'angular_comparison.json').write_text(json.dumps(res,indent=2)+'\n')
print('replay_max_abs',max(absolute),'replay_max_rel',max(relative));print(json.dumps(shell,indent=2))
for r in out:
 print(r['model'],r['half_angle_deg'],'qT/qL relative diffs',r['qT_relative_difference'],r['qL_relative_difference'],'share pp',r['share_percentage_point_difference'])
