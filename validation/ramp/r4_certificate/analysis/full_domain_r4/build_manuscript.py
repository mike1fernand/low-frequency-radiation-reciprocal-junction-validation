"""Compile the revised manuscript; optionally compare the preserved RC4 parent."""
import argparse,concurrent.futures,hashlib,json,os,re,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def tex_environment():
 env=dict(os.environ)
 tree=ROOT/'manuscript/dependencies/revtex-tds'
 if tree.exists():
  # A trailing separator preserves the TeX distribution's default paths.
  env['TEXINPUTS']=str(tree/'tex')+'//:'+env.get('TEXINPUTS','')
  env['BSTINPUTS']=str(tree/'bibtex/bst')+'//:'+env.get('BSTINPUTS','')
 return env

def build(version,stem,out_override=None):
 source=ROOT/'manuscript'/version;out=Path(out_override) if out_override is not None else ROOT/'manuscript/build'/version;out.mkdir(parents=True,exist_ok=True)
 environment=tex_environment()
 class_path=subprocess.run(['kpsewhich','revtex4-2.cls'],env=environment,capture_output=True,text=True,check=True,timeout=10).stdout.strip()
 rec={'version':version,'document':stem,'passes':[],'source_sha256':hashlib.sha256((source/(stem+'.tex')).read_bytes()).hexdigest(),'revtex_class_path':class_path,'revtex_class_sha256':hashlib.sha256(Path(class_path).read_bytes()).hexdigest()}
 for index in range(2):
  start=time.monotonic();cmd=['pdflatex','-interaction=nonstopmode','-halt-on-error','-file-line-error','-output-directory='+str(out),stem+'.tex']
  p=subprocess.run(cmd,cwd=source,env=environment,capture_output=True,text=True,timeout=45)
  (out/f'{stem}.pass{index+1}.stdout.txt').write_text(p.stdout);(out/f'{stem}.pass{index+1}.stderr.txt').write_text(p.stderr)
  rec['passes'].append({'exit':p.returncode,'wall_seconds':time.monotonic()-start})
  if p.returncode:raise RuntimeError(str(source/(stem+'.tex'))+' failed; see build outputs')
 pdf=out/(stem+'.pdf');info=subprocess.run(['pdfinfo',str(pdf)],capture_output=True,text=True,check=True,timeout=15).stdout
 rec['pages']=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1));rec['pdf_sha256']=hashlib.sha256(pdf.read_bytes()).hexdigest()
 log=(out/(stem+'.log')).read_text(errors='replace');rec['overfull_boxes']=re.findall(r'Overfull [^\n]*',log);rec['reference_warnings']=[l for l in log.splitlines() if ('undefined' in l.lower() or 'multiply-defined' in l.lower()) and 'Warning' in l]
 subprocess.run(['pdftotext','-layout',str(pdf),str(out/(stem+'.txt'))],check=True,timeout=15)
 return rec

def main():
 p=argparse.ArgumentParser();p.add_argument('--include-parent',action='store_true');a=p.parse_args()
 jobs=[(v,s) for v in (['parent_rc4','candidate'] if a.include_parent else ['candidate']) for s in ['paper1_main','paper1_supplement']]
 records=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for rec in pool.map(lambda x:build(*x),jobs):records.append(rec)
 (ROOT/'manuscript/BUILD_RECORDS.json').write_text(json.dumps({'status':'compiled','records':records,'visual_review':'separate requirement'},indent=2)+'\n');print(json.dumps(records))
if __name__=='__main__':main()
