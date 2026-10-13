#!/usr/bin/env python3
"""Build and record the all-C/Arb CAP. Python standard library ONLY.
All proof inequalities, quadrature, refinement and scalar closure run in C.
This runner checks inputs, process outcomes, coverage records and provenance;
it never imports NumPy/python-flint or accepts a numerical sample as a proof.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal
import argparse, fcntl, hashlib, json, os, platform, re, subprocess, sys, time
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'artifacts'
JOBS={
 'existence':(['build/verify','0.05','compact'],0,'u0/F0 compact coverage: PASS'),
 'T0':(['build/norm'],0,'T0 compact coverage PASS'),
 'supplemental_scalars':(['build/supplemental_scalars'],0,'PASS supplemental scalar comparisons'),
 'scalars':(['build/scalars'],0,'SCALAR CONDITIONS PASS'),
 'evaluators':(['build/test_evaluators'],0,'regression checks: PASS'),
 'derivative_cusp':(['build/derivative_cusp','0.05'],0,'DERIVATIVE CUSP PASS'),
 'derivative_compact':(['build/derivative_bounds','0.05'],0,'DERIVATIVE COMPACT PASS'),
 'local1':(['build/local_t1','0.002'],0,'LOCAL T1 PASS'),
 'local2':(['build/local_t2','0.002'],0,'LOCAL T2 PASS'),
 'operator_regression':(['build/operators','crosscheck'],0,'global kernel crosscheck PASS'),
 'T1':(['build/operators','1','0.002','cover','0.998'],0,'OP1 COMPACT PASS'),
 'T2':(['build/operators','2','0.002','cover','0.998'],0,'OP2 COMPACT PASS'),
 'sign':(['build/sign'],0,'PASS kernel-sign'),
 'refinement_tests':(['build/test_refine'],0,'PASS 300'),
 'kernel_tests':(['build/test_kernel'],0,'PASS ball conversion tests'),
 'closure':(['build/closure'],0,'PASS scalar closure CONDITIONAL'),
 'matrix':(['build/matrix','inputs/grid720.bin','720','artifacts/moments.bin'],0,'PASS: kernel moment'),
 'reject_defect':(['build/verify','0.05','compact','--reject-test'],2,'UNRESOLVED'),
 'reject_derivatives':(['build/derivative_bounds','0.05','--reject-test'],2,'UNRESOLVED'),
 'reject_T1':(['build/operators','1','.1','reject'],2,'UNRESOLVED'),
 'reject_T2':(['build/operators','2','.1','reject'],2,'UNRESOLVED'),
 'reject_local1':(['build/local_t1','0.5'],2,'UNRESOLVED'),
 'reject_local2':(['build/local_t2','0.003'],2,'UNRESOLVED'),
 'reject_sign':(['build/sign','reject'],2,'UNRESOLVED'),
}
LATE={
 'uniqueness':(['build/uniqueness','inputs/grid720.bin','artifacts/moments.bin'],0,'PASS uniqueness full cusp and 719 compact cells'),
 'reject_uniqueness':(['build/uniqueness','inputs/grid720.bin','artifacts/moments.bin','reject'],2,'UNRESOLVED: cusp comparison target')}
def require(ok,message):
 if not ok:raise RuntimeError(message)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atomic(path,obj):
 t=path.with_suffix('.tmp');t.write_text(json.dumps(obj,indent=2)+'\n');t.replace(path)
def now():return datetime.now(timezone.utc).isoformat()
def sources():
 paths=[ROOT/'Makefile',ROOT/'precision.h',ROOT/'PROOF.md',ROOT/'import-manifest.json',Path(__file__)]
 for part in ['core','src','existence','convexity','inputs']:
  paths += [p for p in (ROOT/part).rglob('*') if p.is_file()]
 return {str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}
def check_inputs():
 text=(ROOT/'precision.h').read_text();require(re.search(r'^#define CAP_PREC 128$',text,re.M),'Common precision must be 128')
 candidate=json.loads((ROOT/'inputs/approximation.json').read_text())
 text=(ROOT/'core/basic_functions.c').read_text()
 actual=dict(re.findall(r'arb_set_str\(up->([^,]+),"([^"]+)",prec\);',text))
 require(actual==candidate['assignments'],'Exact approximation manifest differs from core/basic_functions.c')
 closure=json.loads((ROOT/'inputs/closure.json').read_text())
 declared=dict(re.findall(r'(\w+)\s*=\s*decimal\("([^" ]+)"\)',(ROOT/'src/closure.c').read_text()))
 for c_name,key in {'A':'A','r':'r','q0':'q0','d0':'delta0','d1':'delta1','d2':'delta2','a1':'alpha1','a2':'alpha2','af':'alphaf'}.items():
  require(Decimal(declared[c_name])==Decimal(closure[key]),'Scalar input mismatch '+key)
 require(closure['operator_mode']=='raw' and all(Decimal(closure[k])==Decimal('.998') for k in ['q0','q1','q2']) and Decimal(closure['curvature_budget'])==Decimal('.089'),'Unsupported scalar contract')
 manifest=json.loads((ROOT/'inputs/kernel-sign.json').read_text());code=(ROOT/'src/kernel.c').read_text()
 for side in ['minus','plus']:
  found=re.findall(r'arb_set_str\('+side+r'\s*\+\s*(\d+),\s*"([^"]+)",\s*CAP_PREC\)',code)
  require([int(i) for i,_ in found]==list(range(13)) and [v for _,v in found]==manifest['coeffs'][side],'C root coefficients differ')
 require(sha(ROOT/'inputs/grid720.bin')=='b0d0c1783ca0700ff2cf65e520d2e59fe7513718103e00671181113bbe6f26d2','Exact grid input')
 # Bind mathematical inputs to the current paper, not a whole-file frozen hash.
 # In the ancillary package distributed with the paper there is no paper source;
 # the exact targets are then taken from inputs/kernel-sign.json alone.
 paper=ROOT.parent/'paper/appendices/c-uniqueness-kernel-estimates.tex'
 if not paper.exists():
  return {'path':str(paper),'status':'absent (standalone package); exact kernel-sign coefficients and targets taken from inputs/kernel-sign.json'}
 text='\n'.join(l.split('%',1)[0] for l in paper.read_text().splitlines())
 for side,label,symbol in [('minus','eq:gamma-','y^*_-(x)'),('plus','eq:gamma+','y^*_+(x)')]:
  block=text.split('\\label{'+label+'}',1)[1].split('\\end{multline}',1)[0]
  body=block.split(symbol+':=',1)[1].split('\\,',1)[0].replace('\\\\','')
  pattern=r'([+-]?)\s*(\d+\.\d+)\s*x(?:\^\{?(\d+)\}?)?'
  found=re.findall(pattern,body);require(not re.sub(pattern,'',body).strip(),'Unparsed paper polynomial')
  require([int(p or 1) for _,_,p in found]==list(range(1,14)),'Polynomial powers')
  require([(sgn if sgn=='-' else '')+v for sgn,v,_ in found]==manifest['coeffs'][side],'Paper root coefficients differ')
 claim=re.sub(r'\s+','',text.split('\\label{L.signK}',1)[1].split('\\end{proof}',1)[0])
 patterns={'epsilon':r'\\ep:=([0-9.]+)','epsilon_star':r'\\ep\^\*:=([0-9.]+)\\cdot10\^\{(-?\d+)\}','epsilon_K':r'\\ep_\\cK:=([0-9.]+)\\cdot10\^\{(-?\d+)\}','r_minus':r'r_-:=([0-9.]+)','R_minus':r'R_-:=([0-9.]+)','r_plus':r'r_\+:=([0-9.]+)','R_plus':r'R_\+:=([0-9.]+)'}
 for name,pattern in patterns.items():
  m=re.search(pattern,claim);require(m is not None,'Missing target '+name)
  value=Decimal(m[1])*(Decimal(10)**int(m[2]) if len(m.groups())==2 else 1)
  require(value==Decimal(manifest['targets'][name]),'Changed target '+name)
 return {'path':str(paper),'sha256_at_run':sha(paper),'binding':'Parsed coefficient powers, exact decimal strings and seven exact targets. Other manuscript text is not a computational premise.'}
def coverage(name,text):
 patterns={'existence':r'u0/F0 compact coverage: PASS leaves=(\d+) nodes=(\d+)','derivative_compact':r'DERIVATIVE COMPACT PASS leaves=(\d+) nodes=(\d+)','T0':r'T0 compact coverage PASS leaves=(\d+) nodes=(\d+)','T1':r'OP1 COMPACT PASS leaves=(\d+) nodes=(\d+)','T2':r'OP2 COMPACT PASS leaves=(\d+) nodes=(\d+)'}
 if name in patterns:
  m=re.search(patterns[name],text);require(m is not None,'Missing coverage counts '+name);leaves,nodes=map(int,m.groups());require(leaves>0 and nodes==2*leaves-1,'Incomplete recursive cover '+name)
  if name in ['T1','T2']:require(len(re.findall(r'^OP[12] LEAF ',text,re.M))==leaves,'Missing operator leaves')
 if name=='uniqueness':require([int(i) for i in re.findall(r'^COMPARISON_CELL (\d+) ',text,re.M)]==list(range(719)),'Missing comparison interval');require(len(re.findall(r'^BOUND ',text,re.M))==720,'Missing profile bounds')
 if name=='sign':
  m=re.search(r'PASS kernel-sign.*leaves=(\d+) nodes=(\d+)',text);require(m is not None,'Missing sign coverage');leaves,nodes=map(int,m.groups());require(nodes==2*leaves-448 and len(re.findall(r'^SIGN_LEAF ',text,re.M))==leaves,'Incomplete sign forest')
def run(name,spec):
 cmd,expected,marker=spec;start=time.monotonic();binary=sha(ROOT/cmd[0]);log=OUT/(name+'.log');print('START',name,flush=True)
 with log.open('w') as f:r=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=7200)
 text=log.read_text();require(r.returncode==expected and marker in text,'Stage failed: '+name)
 require(binary==sha(ROOT/cmd[0]),'Binary changed during '+name);coverage(name,text)
 result={'status':'PASS','precision_bits':128,'command':cmd,'exit_code':r.returncode,'success_marker':marker,'elapsed_seconds':time.monotonic()-start,'binary_sha256':binary,'log_sha256':sha(log)}
 atomic(OUT/(name+'.json'),result);print('PASS',name,round(result['elapsed_seconds'],2),'s',flush=True);return result
def libraries():
 text=subprocess.check_output(['ldd',str(ROOT/'build/uniqueness')],text=True)
 paths=re.findall(r'(?:=> )?(/\S+)\s+\(',text)
 return {p:sha(Path(p)) for p in paths}
def headers():
 paths=set()
 for dep in (ROOT/'build').glob('*.d'):
  paths.update(re.findall(r'(?<!\S)(/[^\s\\]+\.h)(?=\s|:|$)',dep.read_text()))
 return {p:sha(Path(p)) for p in sorted(paths)}
def existing():
 c=json.loads((OUT/'certificate.json').read_text());require(c['status']=='PASS','No aggregate PASS');require(c['source_sha256']==sources(),'Source/input mismatch');check_inputs()
 for name,d in c['jobs'].items():
  spec=(JOBS|LATE)[name];require(d['command']==spec[0] and d['exit_code']==spec[1],'Changed command/contract')
  require(sha(ROOT/spec[0][0])==d['binary_sha256'],'Changed binary '+name);log=OUT/(name+'.log');require(sha(log)==d['log_sha256'],'Changed log '+name);coverage(name,log.read_text())
 require(set(c['jobs'])==set(JOBS|LATE),'Missing stages')
 for n,h in c['artifacts_sha256'].items():require(sha(ROOT/n)==h,'Changed artifact '+n)
 for n,h in (c['libraries_sha256']|c['headers_sha256']).items():require(sha(Path(n))==h,'Changed arithmetic dependency '+n)
 print('PASS existing all-C certificate: hashes, inputs and coverage contracts')
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--check-existing',action='store_true');args=parser.parse_args()
 if args.check_existing:existing();return
 OUT.mkdir(exist_ok=True);lock=(OUT/'.run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 start=time.monotonic();record={'schema':'whitham-unified-c-v1','status':'RUNNING','started_utc':now(),'working_precision_bits':128,'human_review':'pending','scope':'Existence, derivative/convexity and uniqueness, conditional on the documented analytic arguments and FLINT/Arb arithmetic base.'};atomic(OUT/'certificate.json',record)
 try:
  record['manuscript_inputs']=check_inputs();before=sources();record['source_sha256']=before
  with (OUT/'build.log').open('w') as log:subprocess.run(['make','-B','-j2','all'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
  record['compiler']=subprocess.check_output(['gcc','--version'],text=True).splitlines()[0];record['environment']={'platform':platform.platform(),'python_orchestrator':sys.version}
  libs=libraries();header_hashes=headers();jobs={}
  def work(item):name,spec=item;return name,run(name,spec)
  with ThreadPoolExecutor(max_workers=2) as pool:
   for name,result in pool.map(work,JOBS.items()):jobs[name]=result
  for name,spec in LATE.items():jobs[name]=run(name,spec)
  require(before==sources(),'Source/input change during full run');require(libs==libraries() and header_hashes==headers(),'Arithmetic dependency change during run');check_inputs()
  record.update(status='PASS',finished_utc=now(),elapsed_seconds=time.monotonic()-start,jobs=jobs,libraries_sha256=libs,headers_sha256=header_hashes,artifacts_sha256={'artifacts/moments.bin':sha(OUT/'moments.bin')},arithmetic='All proof arithmetic in C/Arb at 128 working bits and exact signed 128-bit lattice refinement. Python standard library orchestrates only.',remaining=['Independent human proof/code review','Approved manuscript exposition integration'])
  atomic(OUT/'certificate.json',record);print('PASS complete all-C CAP at 128 bits',flush=True)
 except BaseException as error:
  record.update(status='UNRESOLVED',finished_utc=now(),error=repr(error));atomic(OUT/'certificate.json',record);raise
if __name__=='__main__':main()
