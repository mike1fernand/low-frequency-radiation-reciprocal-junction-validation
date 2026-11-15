"""Freeze current atomic evidence with its exact inherited R3 archive."""
import argparse,hashlib,json,os,shutil,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BASE_SHA='d1b3dcac19297522b6e271711a956fce24f681662c7c978c25fb679289039c88'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--destination',type=Path,required=True);p.add_argument('--zip',type=Path,required=True);p.add_argument('--readme',type=Path,required=True);a=p.parse_args()
 dest=a.destination.resolve()
 if dest.exists():raise ValueError('fresh destination required')
 dest.mkdir(parents=True)
 for source in sorted(ROOT.rglob('*')):
  if not source.is_file() or '__pycache__' in source.parts or source.suffix in ('.pyc','.partial'):continue
  rel=source.relative_to(ROOT)
  # Preserve recovery records, but never ship invalid or duplicate transient PDFs.
  if 'candidate_supplement_verified_tmp' in rel.parts:continue
  if 'interrupted_pdf_observation' in rel.parts and source.suffix=='.pdf':continue
  target=dest/'analysis/full_domain_r4'/rel;target.parent.mkdir(parents=True,exist_ok=True)
  raw=source.read_bytes()
  if source.name=='RUN_LOG.jsonl':
   lines=raw.splitlines(keepends=True);valid=[]
   for i,line in enumerate(lines):
    try:json.loads(line)
    except json.JSONDecodeError:
     if i!=len(lines)-1:raise
     break
    valid.append(line if line.endswith(b'\n') else line+b'\n')
   raw=b''.join(valid)
  target.write_bytes(raw)
 baseline=ROOT.parents[1]/'output_r3/PAPER1_COMPLETE_DOMAIN_VALIDATION_R3.zip'
 if sha(baseline)!=BASE_SHA:raise ValueError('inherited baseline archive identity mismatch')
 (dest/'baseline').mkdir();shutil.copyfile(baseline,dest/'baseline'/baseline.name)
 shutil.copyfile(a.readme,dest/'README.md')
 payloads={p.relative_to(dest).as_posix():sha(p) for p in sorted(dest.rglob('*')) if p.is_file()}
 (dest/'MANIFEST.json').write_text(json.dumps({'schema':'paper1-r4-export-manifest-v1','files':payloads},indent=2)+'\n')
 a.zip.parent.mkdir(parents=True,exist_ok=True)
 temporary=a.zip.with_name(a.zip.name+'.partial')
 with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for source in sorted(dest.rglob('*')):
   if source.is_file():z.write(source,source.relative_to(dest).as_posix(),compress_type=zipfile.ZIP_STORED if source.suffix=='.zip' else zipfile.ZIP_DEFLATED)
 with zipfile.ZipFile(temporary) as z:
  if z.testzip() is not None:raise ValueError('new archive CRC check failed')
 with temporary.open('rb') as completed:os.fsync(completed.fileno())
 os.replace(temporary,a.zip)
 print(json.dumps({'directory':str(dest),'zip':str(a.zip.resolve()),'files':len(payloads),'zip_bytes':a.zip.stat().st_size,'zip_sha256':sha(a.zip)}))
if __name__=='__main__':main()
