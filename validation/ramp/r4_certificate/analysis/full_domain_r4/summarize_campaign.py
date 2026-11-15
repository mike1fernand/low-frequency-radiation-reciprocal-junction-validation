"""Execution diagnostics only; scientific acceptance is in campaign_ledger."""
import json
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent

def main():
 m=json.loads((ROOT/'CAMPAIGN.json').read_text());out={'schema':'paper1-r4-execution-summary-v1','claim':'Execution diagnostics, not an independent mathematical certificate','groups':{},'job_count':len(m['jobs'])}
 for kind in ['Q','g']:
  rows=[]
  for j in m['jobs']:
   if j['kind']!=kind:continue
   p=ROOT/j['state'];s=json.loads(p.read_text()) if p.exists() else None
   rows.append((j,s))
  complete=[(j,s) for j,s in rows if s and s['status']=='rigorous_complete_tile_enclosure']
  out['groups'][kind]={'planned_jobs':len(rows),'complete_jobs':len(complete),'not_started_jobs':sum(s is None for _,s in rows),'partial_jobs':sum(s is not None and s['status']!='rigorous_complete_tile_enclosure' for _,s in rows),'new_complete_physical_parameter_area':str(sum(F(j['area'])*j['multiplicity'] for j,s in complete)),'accepted_leaves':sum(len(s['accepted_ids']) for j,s in rows if s),'stages':sum(len(s['stages']) for j,s in rows if s),'summed_integration_loop_wall_seconds_not_CPU_time':sum(sum(t['seconds'] for t in s['stages']) for j,s in rows if s),'maximum_complete_tile_component_radius_density':max((s['normalized_maximum_radius_diagnostic'] for j,s in complete),default=None)}
 out['coordinator_runs']={p.name:json.loads(p.read_text()) for p in sorted((ROOT/'campaign').glob('RUN_*.json'))}
 (ROOT/'campaign/EXECUTION_SUMMARY.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
