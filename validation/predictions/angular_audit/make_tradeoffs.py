from pathlib import Path
import json,csv,math
P=Path(__file__).resolve().parent
B=P.parent/'baseline/Prediction_Calculations_2026-09-26/angular'
f=json.loads((B/'integral_fine.json').read_text());total=f['total_branch_powers']['full'];tt=sum(total[:2]);tl=total[2]
rows=[]
for a in [10,20,30]:
 fw=f['cone_branch_powers'][f'group_forward_{a}']['full'];bw=f['cone_branch_powers'][f'group_backward_{a}']['full'];ft=sum(fw[:2]);bt=sum(bw[:2]);fl=fw[2];bl=bw[2]
 rows.append({'half_angle_deg':a,'forward_qT_per_incident':ft,'forward_qL_per_incident':fl,'backward_qT_per_incident':bt,'backward_qL_per_incident':bl,'forward_qT_capture_fraction':ft/tt,'forward_qL_capture_fraction':fl/tl,'forward_to_backward_qT':ft/bt,'forward_to_backward_qL':fl/bl,'forward_qT_share':ft/(ft+fl),'backward_qT_share':bt/(bt+bl)})
out={'description':'Derived from archived finite-cone/full-shell powers with the same frozen selected source. Capture fraction divides by each mode cluster total. F/B compares equal cones; neither is physical loaded receiver efficiency.','rows':rows,'exact_opposite_hemisphere_exclusion_half_angle_deg':math.degrees(math.asin(2/(5*math.sqrt(15))))}
(P/'COLLECTION_TRADEOFFS.json').write_text(json.dumps(out,indent=2)+'\n')
with (P/'collection_tradeoffs.csv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
c=json.loads((P/'INDEPENDENT_ANGULAR_CHECK.json').read_text())
print('max independent relative differences',max(abs(r['qT_relative_difference']) for r in c['archive_comparison']),max(abs(r['qL_relative_difference']) for r in c['archive_comparison']))
print(json.dumps(out,indent=2))
