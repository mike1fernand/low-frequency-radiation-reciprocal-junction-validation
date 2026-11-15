"""Verify the Section VI.A and VI.D results from the archived numerical inputs."""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import numpy as np

INPUTS = Path(__file__).resolve().parent.parent / "evidence" / "vi_results"


def read(name):
    return json.loads((INPUTS / name).read_text())


def one(rows, **fields):
    matches = [r for r in rows if all(r.get(k) == v for k, v in fields.items())]
    if len(matches) != 1:
        raise ValueError(f"Expected one row matching {fields}; got {len(matches)}")
    return matches[0]


def close(actual, expected, atol=2e-10):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=atol)


def loaded_line_response(strengths, compliance, qT, qL):
    """Solve the simultaneous displacement equations in the exact retarded lead limit.

    The real-axis projected bulk Green kernel and shell matrices come from the archive.
    From the two squared extensions on a segment we form the stiffness of the complete
    reciprocal contact. Both the line and corner coordinates participate in the solve.
    """
    nseg = len(strengths)
    nline, ncore = nseg + 1, 2*nseg
    n = nline + ncore
    stiffness = np.zeros((n, n))
    for s, g in enumerate(strengths):
        axial, radial = nline + 2*s, nline + 2*s + 1
        x = np.zeros(n)
        x[[s, s + 1, axial]] = [-1, -1, 1]
        r = np.zeros(n)
        r[[s, s + 1, radial]] = [-1, 1, np.sqrt(2)]
        stiffness += (2*g/3) * (np.outer(x, x) + np.outer(r, r))
    Ks = stiffness[:nline, :nline]
    V = stiffness[nline:, :nline]
    U = stiffness[nline:, nline:]
    omega, Kq = 1.2, 1.69
    ka = 2*np.arcsin(omega/(2*np.sqrt(Kq)))
    gamma = 2*Kq*np.sin(ka)
    line = Kq*(2*np.eye(nline)-np.eye(nline, k=1)-np.eye(nline,k=-1))
    lhs = omega**2*np.eye(nline)-line-Ks+0j
    lhs[0,0] += Kq*np.exp(1j*ka)
    lhs[-1,-1] += Kq*np.exp(1j*ka)
    system = np.block([[lhs, -V.T], [compliance @ V, np.eye(ncore)+compliance @ U]])
    source = np.zeros(n)
    source[0] = np.sqrt(gamma)
    z = np.linalg.solve(system, source)
    line_motion, core_motion = z[:nline], z[nline:]
    force = -V @ line_motion - U @ core_motion
    close(core_motion, compliance @ force, 3e-13)
    R = abs(1 - 1j*np.sqrt(gamma)*line_motion[0])**2
    Tq = gamma*abs(line_motion[-1])**2
    powers = [float(2*np.vdot(force, shell@force).real) for shell in (qT,qL)]
    close(R+Tq+sum(powers), 1, 2e-13)
    return np.array([R, Tq, *powers])


def reference_and_profiles():
    reference = np.array([.1,.2,.2,.2,.2,.1])
    cases = {"reference":reference, "endpoint":np.array([.5,0,0,0,0,.5]),
             "doubled":2*reference, "ramp":np.arange(1.,7.)/7}
    equal_total = read("EQUAL_TOTAL.json")["states"]
    scales = read("SCALE_PROBE.json")["states"]
    profiles = read("PROFILE_COMPARISON.json")["states"]
    archives = {"reference":equal_total[0], "endpoint":equal_total[1],
                "doubled":one(scales, cache="N192_M256", scale=2),
                "ramp":next(r for r in profiles if r["cache"]=="N192_M256"
                       and r["shape"]=="increasing_linear"
                       and np.allclose(r["g"],cases["ramp"],rtol=0,atol=1e-15))}
    result = {}
    with np.load(INPUTS/"REFERENCE_KERNEL.npz", allow_pickle=False) as data:
        matrices = [data[key][:12,:12] for key in ("G0","Im_qT","Im_qL")]
        close(data["Im_branch_0"]+data["Im_branch_1"], data["Im_qT"])
        close(data["Im_branch_2"], data["Im_qL"])
    for name, g in cases.items():
        row = archives[name]
        close(g, row["g"], 1e-15)
        close(8*sum(g), row["physical_total_spring_stiffness"])
        if row["active_springs"] != 8*np.count_nonzero(g):
            raise ValueError(f"Active spring count mismatch for {name}")
        powers = loaded_line_response(g, *matrices)
        close(powers, [row[k] for k in ("Rq","Tline","TqT","TqL")])
        share = powers[2]/sum(powers[2:])
        close(share, row["bulk_qT_share"])
        result[name] = {"Rq_Tline_TqT_TqL": powers.tolist(), "emitted_qT_share":share}
    baseline = read("REFERENCE_KERNEL.json")["scattering"]
    recorded_residual = abs(sum(baseline[k] for k in ("Rq","Tline","qT","qL")) - 1)
    close(recorded_residual, baseline["amplitude_closure"], 3e-16)
    if f"{recorded_residual:.2e}" != "1.30e-11":
        raise ValueError("Archived flux residual fails to round to the reported value")
    methods = read("REFERENCE_METHODS.json")["comparisons"]
    surface,tetra,gb = [row["value"] for row in methods]
    close(surface, result["reference"]["Rq_Tline_TqT_TqL"][2])
    close(100*(surface-tetra)/surface,.0408,5e-5)
    if not gb < tetra < surface:
        raise ValueError("Unexpected central-estimate ordering")
    return {"profiles":result, "archived_flux_residual":recorded_residual,
            "central_estimates_surface_tetra_retained_GB":[surface,tetra,gb],
            "central_estimate_spread":surface-gb}


def check_frequency_scan():
    rows=read("ACTIVE85.json")["families"]["frequency"]
    frequencies=np.array([r["omega"] for r in rows])
    close(frequencies, np.linspace(.6,1.8,49), 1e-14)
    for row in rows:
        if row["method"]["surface_grid"]!=80 or row["method"]["PV_order"]!=256:
            raise ValueError("Unexpected sweep resolution")
        close(row["Rq"]+row["Tq"]+row["Tbulk"],1.,5e-9)
        close(row["T_qT_modal"]+row["T_qL_modal"],row["Tbulk"],5e-14)
        close(row["qT_fraction_of_bulk"],row["T_qT_modal"]/row["Tbulk"],5e-14)
    low,high = [fn(rows,key=lambda r:r["T_qT_modal"]) for fn in (min,max)]
    close([low["omega"],high["omega"]],[1.575,.825])
    close([low["T_qT_modal"],high["T_qT_modal"]],[.104864,.217743],5e-7)
    return {"rows":len(rows),"min_frequency_TqT":[low["omega"],low["T_qT_modal"]],
            "max_frequency_TqT":[high["omega"],high["T_qT_modal"]],
            "samples_above_10_percent":sum(r["T_qT_modal"]>.1 for r in rows),
            "reference_frequency_sweep_TqT":one(rows,omega=1.2)["T_qT_modal"]}


def check_certificate_reporting():
    ref = read("REFERENCE_EXACT_LATTICE_BRIDGE.json")
    cert = read("RAMP_CERTIFICATE_CHAIN.json")
    if not ref["exact_lattice_floor_certified"] or Fraction(ref["rational_floor"])!=Fraction(1037,10000):
        raise ValueError("Reference certificate differs")
    if cert["target_model"]["g"]!=[f"{i}/7" for i in range(1,7)]:
        raise ValueError("Wrong nominal ramp")
    if Fraction(cert["target_model"]["omega"])!=Fraction(6,5):
        raise ValueError("Wrong certified frequency")
    interval = cert["nominal_infinite_lattice_interval"]
    lower, upper = [Fraction(interval[k]) for k in ["lower_rational","upper_rational"]]
    floor = Fraction(cert["whole_six_motif_plus_minus5_percent_infinite_lattice_floor"]["lower_rational"])
    if not (Fraction(".732166")<lower<upper<Fraction(".739652") and floor>Fraction(".69")):
        raise ValueError("The stated bounds do not follow from the saved rationals")
    return {"nominal_interval":[float(lower),float(upper)],"box_floor":float(floor),
            "scope":"Arithmetic check on the endpoints of the existing proofs."}


def check_upstream():
    freq = read("BACKWARD_FREQUENCY_REFINEMENT.json")
    pv = read("BACKWARD_PV_REFINEMENT.json")
    at26 = one(freq,omega=2.6,cq=1.5,Nc=24,s=1,Nperp=2048)
    at29 = one(pv,omega=2.9,cq=1.5,Nc=24,s=1,PV="exact_kx_midpoint_PV_N1536")
    output = {}
    for row,key,printed in [(at26,"powers_by_sorted_branch_forward_backward",[33.2,17.1,51.5]),
                            (at29,"powers",[72.1,71.,98.5])]:
        directional = np.array(row[key])
        bulk = directional.sum()
        upstream = directional[:,1].sum()
        close([bulk,upstream,upstream/bulk], [row["Tbulk"],row["backward_incident"],
                                                row["backward_share"]],2e-12)
        close(1, row["R"]+row["Tline"]+bulk)
        values = 100*np.array([bulk,upstream,upstream/bulk])
        close(values,printed,.05)
        output[str(row["omega"])] = values.tolist()
    close(np.array(at29["powers"])[:2], 0, 1e-15)
    return {"percent_bulk_upstream_upstream_share":output,
            "scope":"Saved directional integrations; forward and backward branch powers summed here."}


def main():
    for filename, record in read("PROVENANCE.json").items():
        if sha256((INPUTS/filename).read_bytes()).hexdigest()!=record["sha256"]:
            raise ValueError("Input digest mismatch for "+filename)
    result={"status":"PASS","stationary_reference":reference_and_profiles(),
            "frequency_scan":check_frequency_scan(),"certificate_reporting":check_certificate_reporting(),
            "directional_example":check_upstream(),
            "scope":"New coupled solves on saved numerical kernel. Scan and integration records rechecked; no fresh bulk quadrature or new error certificate."}
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
