from pathlib import Path
import hashlib,json,platform,shutil,sys,difflib
import numpy as np
ROOT=Path(__file__).resolve().parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while True:
   b=f.read(4*1024*1024)
   if not b:break
   h.update(b)
 return h.hexdigest()
def J(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
cases=[('base80','CANDIDATE_PACKET_BOX80_RESULT.json','candidate_dt004_box80'),('half_step80','CONTROL_DT002_BOX80_RESULT.json','control_dt002_box80'),('larger96','CONTROL_DT004_BOX96_RESULT.json','control_dt004_box96')]
data={};cadence={}
for name,file,folder in cases:
 p=ROOT/file;d=json.loads(p.read_text());B=np.asarray(d['branch_samples']);qT=float(np.mean(d['mean_low_high_total'][:2]));lo,hi=d['longest']
 assert [lo,hi]==[57.4,64.2] and d['frozen_window_all_predicates_pass']
 rows=[x for x in d['rows'] if lo<=x['t']<=hi]
 mean=lambda a:(np.trapezoid(a[:,1:],x=a[:,0],axis=0)/(a[-1,0]-a[0,0])).tolist()
 i6=np.unique(np.r_[np.arange(0,len(B),3),len(B)-1]);i4=np.unique(np.r_[np.arange(0,len(B),2),len(B)-1])
 c6=mean(B[i6]);c4=mean(B[i4]);old6=mean(B[::3]);cmid=lambda x:float(np.mean(x[:2]))
 cadence[name]={'source_file':file,'source_sha256':sha(p),'primary_cadence_0p2_mean':d['mean_low_high_total'],'current_cadence_0p4_same_window_mean':c4,'current_cadence_0p4_relative_qT_difference_percent':100*(cmid(c4)/qT-1),'current_cadence_0p6_plus_final_endpoint_mean':c6,'current_common_grid_times':B[i6,0].tolist(),'legacy_everythird_truncated_mean':old6,'legacy_everythird_actual_window':[float(B[0,0]),float(B[::3][-1,0])],'current_coarse_grid_window':[float(B[0,0]),float(B[-1,0])],'supersedes_old_raw_common_field':file in ['CANDIDATE_PACKET_BOX80_RESULT.json','CONTROL_DT002_BOX80_RESULT.json']}
 evofile=ROOT/folder/('EVOLUTION_CHECK.json' if name!='base80' else 'EVOLUTION_CHECK.json')
 evo=json.loads(evofile.read_text())
 line=np.memmap(ROOT/folder/'line.bin',dtype=float,mode='r').reshape(-1,2,120)
 dt=d['config']['dt'];indices=np.arange(round(lo/dt),round(hi/dt)+1,round(.2/dt));assert len(indices)==len(rows)
 q=line[indices,0];v=line[indices,1];node=.5*v*v;spring=.25*1.69*np.diff(q,axis=1)**2;node[:,:-1]+=spring;node[:,1:]+=spring
 t=np.asarray([x['t'] for x in rows]);E0=d['open_line_E0'];avg=lambda x:float(np.trapezoid(x,x=t)/(hi-lo))
 energy={'line_total':avg(node.sum(1)/E0),'line_left_nodes_0_to_59':avg(node[:,:60].sum(1)/E0),'line_contact_nodes_60_to_66':avg(node[:,60:67].sum(1)/E0),'line_right_nodes_67_to_119':avg(node[:,67:].sum(1)/E0),'open_bulk':avg(np.asarray([x['open_bulk'] for x in rows])),'contact_energy':avg(np.asarray([x['link'] for x in rows]))}
 energy['sum']=energy['line_total']+energy['open_bulk']+energy['contact_energy']
 data[name]={'result_file':file,'result_sha256':sha(p),'configuration_file':folder+'/configuration.json','configuration_sha256':sha(ROOT/folder/'configuration.json'),'window':[lo,hi],'all_five_predicates_pass':d['frozen_window_all_predicates_pass'],'qT_allocation_interval_average':d['mean_low_high_total'][:2],'qT_midpoint_average':qT,'periodic_bulk_average':d['mean_low_high_total'][2],'periodic_qL_midpoint_average':d['mean_low_high_total'][2]-qT,'qT_share_of_periodic_bulk':qT/d['mean_low_high_total'][2],'max_total_six_face_energy_shell_fraction':d['postselection_total_node_energy_shell_max'],'max_axial_energy_shell_fraction':d['postselection_axial_node_energy_shell_max'],'max_transverse_energy_shell_fraction':d['postselection_transverse_node_energy_shell_max'],'current_balance_residual_percent':d['percent_current_balance_residual'],'max_selected_window_energy_drift_at_cadence_0p2':max(x['drift'] for x in rows),'max_whole_run_energy_drift_at_saved_checkpoints_only':evo['sampled_max_relative_energy_drift'],'whole_run_drift_scope':'Monitorcheckpoints,notmaximumovereveryintegratorlevel. Selectedwindowdriftsareindependentlyevaluatedeach0.2cadence.','physical_open_energy_partition_average':energy,'evolution_wall_seconds':evo['wall_seconds'],'evolution_peak_RSS_KiB':evo.get('peak_rss_kib',evo.get('max_rss_kib_process_tree'))}
base=data['base80']['qT_midpoint_average']
for name in data:data[name]['relative_qT_difference_vs_base_percent']=100*(data[name]['qT_midpoint_average']/base-1)
J(ROOT/'CURRENT_CADENCE_CHECKS.json',{'status':'PASS_CORRECTED_SECONDARY_CADENCE_FIELDS','scope':'Primarymeansunchanged. Rawbase/halfstepJSONpreservedasfrozenhistory; theircommon_mean_low_high_totalfieldistruncatedto64.0andsupersededherebytheexplicit64.2endpointmeans.','cases':cadence})
prediction=.7163064796290152
summary={'status':'PASS','date':'2026-11-08','model':'OriginalP1passiveharmonicopenfinitejunction;newrampcontactg=(1,2,3,4,5,6)/7;all48physicalcontactsretained.','prepared_state':'Exactlysame120-siteGaussian/spectralvelocityasreference;bulkexactlyquiet.','normalization':'DirectmodeenergiesdividedbyinitialOPENscalarlineHq=E0; stationaryweightedpredictionusesall-gridperiodicParsevalenergy. Differenceis9.184e-9relativeandmustnotbeconfusedwithaunit-fluxcarriercoefficient.','carrier_stationary_qT_fraction':.7358742883948081,'stationary_full_grid_inband_packet_qT_contribution':prediction,'stationary_tail_only_interval':[prediction,prediction+1-.9996834927581062],'stationary_interval_scope':'Omittedspectralweightonly;notquadrature,finite-time,dynamicsorboxerror.','initial_open_line_E0':json.loads((ROOT/'PREDECLARED_PACKET_CONFIGURATION.json').read_text())['initial_open_line_energy'],'initial_contact_energy':json.loads((ROOT/'PREDECLARED_PACKET_CONFIGURATION.json').read_text())['initial_contact_energy'],'frozen_window_source':'FROZEN_WINDOW.json','frozen_window_sha256':sha(ROOT/'FROZEN_WINDOW.json'),'initial_box56_result':'NO_ACCEPTED_WINDOW;only53.8passesallpredicates;currentringdownandlaterboundaryexposurepreventduration>=2. Preservedunmodified.','cases':data,'cadence_authority':'CURRENT_CADENCE_CHECKS.json supersedes oldsecondarycommon_mean fields;primaryall0.2meansareunchanged.','primary_result_policy':'Base80dt0.04isthepredeclaredmainresult;larger96isafixed-windowboxcontrol,notselectedbecauseitisclosertofrequency-domainprediction.','new_results_are':'Finite-timestoredlower-clusterenergycorroboration,notmonochromatic73.6%measurements,localremoteqTflux,formalexactlatticeerrorenclosureornewglobalcertificate.','no_further_dynamics_selected':True,'runtime':{'python':platform.python_version(),'numpy':np.__version__}}
for name in data:
 q=data[name]['qT_midpoint_average'];data[name]['difference_vs_stationary_packet_contribution_percentage_points']=100*(q-prediction);data[name]['relative_difference_vs_stationary_packet_contribution_percent']=100*(q/prediction-1)
J(ROOT/'FINAL_SUMMARY.json',summary)
# The source archive's own payload hashes authenticate the preserved implementation.
orig=ROOT/'recovered/p1_packet_reproduction';original_hashes=json.loads((orig/'SHA256.json').read_text())
assert all(sha(orig/f)==h for f,h in original_hashes.items())
J(ROOT/'ORIGINAL_PRODUCER_IDENTITY.json',{'status':'PASS','source_archive_sha256':'67b63e05d30857e8fe6a0e206d117d4d814482130f63a2371b0326e606493bd2','version':0,'archive':'inputs/Paper1_Packet_Reproduction_and_Erratum_2026-11-01.zip','archive_sha256':sha(ROOT/'inputs/Paper1_Packet_Reproduction_and_Erratum_2026-11-01.zip'),'all_public_payload_hashes_match':True,'payload_count':len(original_hashes),'original_cpp_sha256':sha(orig/'device/evolve.cpp'),'original_analysis_sha256':sha(orig/'device/analyze.py'),'new_cpp_sha256':sha(ROOT/'evolve_candidate.cpp'),'new_analysis_sha256':sha(ROOT/'analyze_candidate.py'),'public_edition':'The producer code, data and scientific records match the original. Four private links were replaced by their original labels in START_HERE.md, and the public SHA256.json records that change.'})
for old,new,label in [(orig/'device/evolve.cpp',ROOT/'evolve_candidate.cpp','evolution'),(orig/'device/analyze.py',ROOT/'analyze_candidate.py','analysis')]:
 (ROOT/(label+'_changes.patch')).write_text(''.join(difflib.unified_diff(old.read_text().splitlines(keepends=True),new.read_text().splitlines(keepends=True),fromfile=str(old.relative_to(ROOT)),tofile=str(new.relative_to(ROOT)))))
raw=[]
for folder in ['candidate_dt004_box56','candidate_dt004_box80','control_dt002_box80','control_dt004_box96','resume_control_continuous','resume_control_staged']:
 for p in sorted((ROOT/folder).glob('*.bin')):raw.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
J(ROOT/'RAW_DATA_MANIFEST.json',{'status':'RECORDED','raw_files':raw,'core_zip_exclusion_policy':['*/snap.bin','*/snap_times.bin','*/line.bin','*/corner.bin','*/checkpoint_*.bin','evolve_candidate (generatedbinary)','private_helpers/'], 'optional_core_continuation_checkpoint':'candidate_dt004_box80/checkpoint_1650.bin','checkpoint_scope':'Complete576120-coordinatequarterbulk+lineu/vwithtime/params/bondinventory/integratedcurrentstate;equivalenttofull120x80x80C4vstate. Noaxialsymmetryassumed.','raw_data_reproducibility':'Unalteredrecordedinitialstate+exactnewsource/CLI/configurations regenerateallarraysinseconds-to-minutes. Hashesidentifythisexecution;cross-platformbitwiseidentityisnotguaranteed.','raw_total_bytes':sum(x['bytes'] for x in raw)})
print(json.dumps({'status':'PASS','frozen_window':[57.4,64.2],'cases':{k:{'qT':v['qT_midpoint_average'],'relative_difference_vs_base_percent':v['relative_qT_difference_vs_base_percent'],'selected_window_drift':v['max_selected_window_energy_drift_at_cadence_0p2']} for k,v in data.items()},'raw_bytes':sum(x['bytes'] for x in raw)},indent=2))
