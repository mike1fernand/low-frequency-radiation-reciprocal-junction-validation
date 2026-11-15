"""Resolve all 85 one-factor cases and 2,702 tolerance-screen rows on archived bulk responses."""
import argparse
from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import numpy as np
from verify_packet_spectrum import contact_hessian

REFERENCE = (.1,.2,.2,.2,.2,.1)
COUNTS = dict(frequency=49, coupling_scale=10, length=6, profile=5, source_speed=10, stiffness_ratio=5)
SCREEN_FREQUENCIES = [1.10,1.1125,1.125,1.15,1.175,1.2,1.225,1.2375,1.25,1.275,1.30]
TETRA_FREQUENCIES = [1.1125,1.2375,1.3]
FIELDS = ["Rq","Tq","T_qT_modal","T_qL_modal"]


def require(test, message):
    if not test: raise ValueError(message)


def read_record(path, expected_digest):
    data = path.read_bytes()
    require(sha256(data).hexdigest()==expected_digest, f"Original numerical record differs: {path}")
    return json.loads(data)


@lru_cache(maxsize=240)
def contact(strengths):
    return contact_hessian(strengths)


def solve(omega, strengths, kernel, line_stiffness, eta):
    """Solve simultaneous line/corner displacements. The 85-row ledger uses eta=1e-12."""
    n = len(strengths)+1
    real,lower_0,lower_1,upper = kernel[:,:2*(n-1),:2*(n-1)]
    lower = lower_0+lower_1
    green = real+1j*(lower+upper)
    z = complex(omega,eta)
    ka = 2*np.arcsin(z/(2*np.sqrt(line_stiffness)))
    sigma = -line_stiffness*np.exp(1j*ka)
    gamma = -2*sigma.imag
    require(gamma > 0, "Incident channel is closed")
    hess = contact(tuple(strengths))
    free = line_stiffness*(2*np.eye(n)-np.eye(n,k=1)-np.eye(n,k=-1))
    device = z*z*np.eye(n)-free-hess[:n,:n]
    device[[0,-1],[0,-1]] -= sigma
    coupling, local = hess[n:,:n], hess[n:,n:]
    system = np.block([[device, -coupling.T],[green@coupling,np.eye(green.shape[0])+green@local]])
    incident = np.zeros(len(hess)); incident[0] = np.sqrt(gamma)
    state = np.linalg.solve(system, incident)
    force = -coupling@state[:n] - local@state[n:]
    require(np.linalg.norm(state[n:]-green@force)<2e-12, "Bulk force/displacement residual")
    return np.array([abs(1-1j*np.sqrt(gamma)*state[0])**2,gamma*abs(state[n-1])**2,
                     2*np.vdot(force,lower@force).real,2*np.vdot(force,upper@force).real])


def verify_powers(calculated, recorded):
    diff=float(max(abs(calculated-np.asarray(recorded))))
    require(diff<2e-12, f"Power replay differs by {diff}")
    require(min(calculated)>=0 and abs(sum(calculated)-1)<1.e-9,"Passive power accounting")
    return diff


def one_factor(ledger, kernels):
    families = ledger["families"]
    require({k:len(v) for k,v in families.items()} == COUNTS, "One-factor case counts")
    expected_rows = {"frequency": [((24+i)/40) for i in range(49)],
                     "coupling_scale":[.05,.075,.1,.15,.2,.25,.5,1,1.5,2],
                     "length":[2,3,4,6,8,10],"profile":["reference","uniform","soft_ramp","short_core","edge_heavy"],
                     "source_speed":[.65,.8,.9,1.,1.1,1.2,1.3,1.5,1.7,2.],
                     "stiffness_ratio":[.75,.9,1.,1.1,1.25]}
    keys=dict(frequency="omega",coupling_scale="scale",length="nseg",profile="name",
              source_speed="c_spine_longwave",stiffness_ratio="k_diag_over_k_axis")
    err=0.; allrows=[]
    for family,rows in families.items():
        require([r[keys[family]] for r in rows]==expected_rows[family],f"Changed {family} settings")
        for i,r in enumerate(rows):
            w=r.get("omega",1.2); kq=r.get("c_spine_longwave",1.3)**2
            g=np.array(r.get("g", REFERENCE))*r.get("scale",1.)
            if family=="length":
                g=np.full(r["nseg"],.2);g[[0,-1]]=.1
            kernel=kernels["frequency_24"]
            if family=="frequency": kernel=kernels[f"frequency_{i}"]
            if family=="stiffness_ratio": kernel=kernels[f"stiffness_{i}"]
            values=solve(w,g,kernel,kq,1e-12)
            err=max(err,verify_powers(values,[r[k] for k in FIELDS]))
            require(values[2]>0,"Zero qT in positive-coupling record")
            require(abs(r["Tbulk"]-sum(values[2:]))<2e-12,"Bulk total")
            require(abs(r["qT_fraction_of_bulk"]-values[2]/sum(values[2:]))<2e-12,"Emitted denominator")
            closure=abs((r["Rq"]+r["Tq"]+r["Tbulk"])-1)
            require(abs(closure-r["flux_closure_abs"])<4e-16,"Archived flux residual")
            allrows.append(r)
    # The original producer sums all three branch terms before the modal residual;
    # this may round differently from comparing already rounded lower-group and upper entries.
    maxima=[max(r[k] for r in allrows) for k in ("flux_closure_abs","modal_sum_residual_abs")]
    require([format(x,".2e") for x in maxima]==["3.97e-11","8.33e-17"],"Printed numerical residuals")
    freq = families["frequency"]
    ordering = sorted(freq,key=lambda r:r["T_qT_modal"])
    ends = [[r["omega"],r["T_qT_modal"]] for r in [ordering[0], ordering[-1]]]
    require([round(x[1]*100,1) for x in ends]==[10.5,21.8],"Sampled frequency extrema")
    scales=families["coupling_scale"]
    require(all(p["Tq"]>q["Tq"] and p["Tbulk"]<q["Tbulk"] for p,q in zip(scales,scales[1:])),
            "Sampled trend under taper scaling")
    require(min(r["qT_fraction_of_bulk"] for r in scales)>.91,"Bulk share under taper scaling")
    selected={str(r["scale"]):[r[k] for k in FIELDS] for r in scales if r["scale"] in (1,1.5,2)}
    require(np.allclose(selected["2.0"], [.0445,.3837,.5222,.0496],atol=5e-5,rtol=0),"Doubled taper row")
    dims=ledger["dimensional_unit_covariance"]
    require([r['a'] for r in dims]==[.25,.5,1,2,4] and
            all(all(r[k]==dims[0][k] for k in FIELDS) for r in dims),"Lattice spacing re-expression")
    zero=solve(1.2,[0.]*6,kernels["frequency_24"],1.69,1e-12)
    verify_powers(zero,[ledger["coupling_including_zero"][0][k] for k in FIELDS])
    require(np.array_equal(zero[2:],[0.,0.]),"Zero-contact control")
    return dict(solves=len(allrows),counts=COUNTS,zero_contact_additional_control=True,
                max_abs_new_minus_archived_channel_power=err,archive_residual_maxima=maxima,
                sampled_frequency_extrema_omega_and_incident_qT=ends,
                three_taper_scales_Rq_Tq_TqT_TqL=selected,
                minimum_positive_case_incident_qT=min(r["T_qT_modal"] for r in allrows),
                min_emitted_qT_share_in_positive_taper_scale_scan=min(r["qT_fraction_of_bulk"] for r in scales))


def tolerance_screen(screen,kernels):
    result={}
    profile_map={}
    for method,rows,frequencies in [("surface",screen["surface_rows"],SCREEN_FREQUENCIES),
                                    ("tetra",screen["independent_rows"],TETRA_FREQUENCIES)]:
        groups=defaultdict(dict);errors=[];minimum=None
        for r in rows:
            w=r["omega"];g=tuple(r["g"]);label=r["setting"]
            require(label not in groups[w],"Duplicate (frequency,setting)")
            groups[w][label]=g
            keys=["Rq","Tline","TqT","TqL"] if method=="surface" else ["Rq_amplitude","Tq","TT_modal","TL_modal"]
            if method=="surface":
                require(r["shell_grid"]==(128 if w in TETRA_FREQUENCIES else 80) and r["PV_order"]==256,
                        "Screen resolution")
            if method=="tetra":
                require(abs(sum(r["T_branches"][:2])-r["TT_modal"])<1e-15 and
                        abs(r["T_branches"][2]-r["TL_modal"])<1e-15,"Branch indexing")
            values=solve(w,g,kernels[f"{method}_{w:.4f}"],1.69,0.)
            errors.append(verify_powers(values,[r[k] for k in keys]))
            require(values[2]>.10,"Declares a sampled threshold violation")
            if minimum is None or values[2]<minimum[0]: minimum=(float(values[2]),w,label)
        require(sorted(groups)==frequencies,"Frequency grid")
        for group in groups.values():
            require(len(group)==len(set(group.values()))==193 and group["nominal"]==REFERENCE,"Profile count")
            corners=set(product(*[(.95*g,1.05*g) for g in REFERENCE]))
            require({g for tag,g in group.items() if tag.startswith("corner_")}==corners,"Reference corners")
            interiors=[g for tag,g in group.items() if tag.startswith("interior_")]
            require(len(interiors)==128 and all(all(.95*a<b<1.05*a for a,b in zip(REFERENCE,g)) for g in interiors),
                    "Strictly interior profiles")
            if not profile_map: profile_map=group
            require(group==profile_map,"Different frequency or method used different trial profiles")
        require(minimum[1:]==(1.3,"corner_0"),"Common minimum setting")
        result[method]=dict(solves=len(rows),frequency_points=len(groups),profiles_per_frequency=193,
                            minimum_qT=minimum[0],minimum_omega=minimum[1],minimum_setting=minimum[2],
                            max_abs_new_minus_archived_channel_power=max(errors))
    return result


def verify(root):
    path=root/"evidence/robustness"
    active=read_record(root/"evidence/vi_results/ACTIVE85.json","93ae058159343720a54fb272b2d9f9dd5b0cce839d652bf1e0e66f69b8174ef5")
    screen=read_record(path/"REFERENCE_FREQUENCY_SCREEN.json","4b5da237ca9fe53908e41012d69c99811575c40e7b95218cab49f0379ae3398a")
    with np.load(path/"KERNELS.npz",allow_pickle=False) as source: kernels=dict(source)
    return dict(status="PASS",one_factor=one_factor(active,kernels),tolerance_screen=tolerance_screen(screen,kernels),
                scope="New coupled solves on retained numerical bulk response matrices. Common-medium surface principal values reconstructed from archived spectral histograms. No new integration of the underlying bulk shells and no new continuous-domain certificate.")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    print(json.dumps(verify(parser.parse_args().root),indent=2,allow_nan=False))
