#!/usr/bin/env python3
"""New static-algebra check of the fixed-resource concavity corollary.

Uses the actual archived P1 complete-motif blocks with synthetic positive
compliances. It computes no physical bulk kernel, scattering state or cutoff.
The proof, not finite fixtures, establishes the global statement.
"""
from pathlib import Path
import hashlib, json, sys, tempfile
import numpy as np

HERE=Path(__file__).resolve().parent
INPUTS=HERE.parents[1]/'recheck/purity_inputs'

def main():
    metadata=json.loads((INPUTS/'INPUT_PROVENANCE.json').read_text())
    known={x['path']:x for x in metadata['files']}
    payload={n:(INPUTS/n).read_bytes() for n in ['model_api_v1_1.py','reduced_core.py']}
    for n,data in payload.items():assert hashlib.sha256(data).hexdigest()==known[n]['sha256']
    records=[]
    sys.dont_write_bytecode=True
    with tempfile.TemporaryDirectory(prefix='p1-concavity-') as td:
        for n,data in payload.items():(Path(td)/n).write_bytes(data)
        sys.path.insert(0,td)
        from model_api_v1_1 import ModelParameters
        from reduced_core import reduced_static_matrices
        par=ModelParameters()
        for n in [1,2,6]:
            L=np.zeros((n+1,n+1))
            for i in range(n):
                v=np.eye(n+1)[i+1]-np.eye(n+1)[i]
                L+=par.k_spine*np.outer(v,v)
            profiles=[np.ones(n)/n,np.arange(1,n+1,dtype=float)/sum(range(1,n+1))]
            if n>1:
                edge=np.zeros(n);edge[0]=edge[-1]=.5
                profiles.append(edge)
            for scale in [.2,1.,3.]:
                i=np.arange(2*n)
                # Correlated SPD fixture in the same fixed core-coordinate space.
                G=scale*(np.eye(2*n)+.4**np.abs(i[:,None]-i[None,:]))
                Ginv=np.linalg.inv(G)
                def A(g):
                    U,V,K=reduced_static_matrices(g,par)
                    out=L+K-V.conj().T@np.linalg.solve(Ginv+U,V)
                    return (out+out.conj().T)/2
                vertices=[A(np.eye(n)[j]) for j in range(n)]
                vertex_gap=min(float(np.linalg.eigvalsh(v).min()) for v in vertices)
                assert vertex_gap>0
                for p in profiles:
                    actual=A(p)
                    chord=sum(p[j]*vertices[j] for j in range(n))
                    slack=np.linalg.eigvalsh(actual-chord)
                    gap=float(np.linalg.eigvalsh(actual).min())
                    assert slack.min()>-2e-12
                    assert gap>=vertex_gap-2e-12
                    records.append({'motifs':n,'compliance_scale':scale,'profile':p.tolist(),
                        'min_concavity_slack_eigenvalue':float(slack.min()),
                        'max_concavity_slack_eigenvalue':float(slack.max()),
                        'minimum_vertex_static_gap':vertex_gap,'profile_static_gap':gap})
                if n==6 and scale==1.:
                    weak=[{'total_strength':s,'static_gap':float(np.linalg.eigvalsh(A(s*np.ones(n)/n)).min())}
                          for s in [.1,.01,.001]]
                    assert all(weak[j+1]['static_gap']<weak[j]['static_gap'] for j in range(2))
        assert max(r['max_concavity_slack_eigenvalue'] for r in records)>.01
    result={'status':'PASS','scope':'Static concavity with actual motif blocks and synthetic SPD compliance; no physical bulk/scattering calculation',
            'cases':records,'weak_total_coupling_adverse_control':weak,
            'input_hashes':{n:hashlib.sha256(v).hexdigest() for n,v in payload.items()},
            'source_hash':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'certification':'Floating-point diagnostics; no physical static gap or finite-frequency accuracy window inferred'}
    path=HERE/'static_concavity_results.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','cases':len(records),
        'minimum_concavity_slack':min(r['min_concavity_slack_eigenvalue'] for r in records),
        'scope':result['scope']},indent=2))

if __name__=='__main__':main()
