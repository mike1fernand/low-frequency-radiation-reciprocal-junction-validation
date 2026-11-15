#!/usr/bin/env python3
"""Render only four affected Paper1 plots from a selected, hash-bound JSON contract.

No scientific input is embedded. The binding's source paths are relative to the
binding file. A reference is {"input": source_key, "pointer": "/JSON/pointer"}.
All plotted arrays, labels, numerical derivations and source identities are
written beside the figures. Unaffected main figures are checked, never rewritten.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NAMES = {
    "reference": "paper1_sectionVI_fig4_reference_partition_v2.pdf",
    "packet": "paper1_sectionVII_fig6a_qT_window_v2.pdf",
    "alpha": "Halpha_transmission_scan_p06.png",
    "continuum": "continuum_lattice_share_comparison_p06.png",
}
FIXED_DATE = datetime(2026, 9, 20, tzinfo=timezone.utc)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pointer(value, path):
    if path == "":
        return value
    if not path.startswith("/"):
        raise ValueError("Expected RFC6901 JSON pointer")
    for token in path[1:].split("/"):
        token = token.replace("~1", "/").replace("~0", "~")
        value = value[int(token)] if isinstance(value, list) else value[token]
    return value


def finite(x):
    a = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError("Nonfinite scientific figure input")
    return a


def collect(binding_path):
    b = json.loads(binding_path.read_text())
    if b["schema"] != "paper1.selected-figure-inputs.v1" or b["status"] != "selected-current":
        raise ValueError("A selected-current scientific binding is required")
    inputs = {}
    for key, item in b["inputs"].items():
        path = (binding_path.parent / item["path"]).resolve()
        if sha256(path) != item["sha256"]:
            raise ValueError(f"Changed figure input: {key}")
        canonical = Path(__file__).resolve().parents[3] / item["source_path"]
        if sha256(canonical) != item["source_sha256"] or item["source_sha256"] != item["sha256"]:
            raise ValueError(f"Canonical figure input changed: {key}")
        inputs[key] = json.loads(path.read_text())
    for item in b["preserve_unchanged"].values():
        path = (binding_path.parent / item["path"]).resolve()
        if sha256(path) != item["sha256"]:
            raise ValueError(f"Changed unaffected figure: {path.name}")
    if len(b["preserve_unchanged"]) != 2:
        raise ValueError("Bind both unaffected main figures")

    def ref(spec):
        return pointer(inputs[spec["input"]], spec["pointer"])

    r = ref(b["reference"])
    fields = ["Rq", "Tq", "TL_modal", "TT_modal"]
    vals = finite([r[k] for k in fields])
    if min(vals) < -1e-10 or abs(sum(vals)-1) > 1e-7:
        raise ValueError("Reference partition fails its independently supplied normalization")
    if abs(r["TT_modal"]+r["TL_modal"]-r["Tbulk"]) > 1e-9:
        raise ValueError("Reference branch sum disagrees with bulk power")
    comparisons = [{"label": x["label"], "value": float(ref(x)), "status": x["status"]}
                   for x in b["comparisons"]]
    finite([x["value"] for x in comparisons])
    if not comparisons or comparisons[0]["value"] != float(r["TT_modal"]):
        raise ValueError("First comparison must be the controlling reference")
    tetra_idx = b["tetrahedron_comparison_index"]
    if not 0 <= tetra_idx < len(comparisons):
        raise ValueError("Invalid tetrahedron comparison index")
    d4 = {
        "partition": dict(zip(fields, map(float, vals))),
        "partition_sum": float(sum(vals)), "raw_closure": float(r["closure"]),
        "comparisons": comparisons,
        "surface_tetra_relative_percent": 100*abs(comparisons[tetra_idx]["value"]-r["TT_modal"])/r["TT_modal"],
        "relative_difference_denominator": "controlling surface TT_modal",
        "reference_label": b["reference_label"],
    }

    selector = ref(b["selector"])
    times = finite([x["time"] for x in selector["rows"]])
    trace = finite([x["qT_production"] for x in selector["rows"]])
    if np.any(np.diff(times) <= 0):
        raise ValueError("Selector times are not strictly increasing")
    central = float(ref(b["packet"]["central"]))
    variants = [{"label": x["label"], "value": float(ref(x))} for x in b["packet"]["variants"]]
    pv = finite([central]+[x["value"] for x in variants])
    if np.any(pv<0) or np.any(pv>1):
        raise ValueError("Packet predictions outside passive coefficient range")
    if b["packet"]["normalization"] != "full_prepared_packet_energy":
        raise ValueError("Figure5(a) requires the declared full-energy packet normalization")
    d5 = {"time": times.tolist(), "qT": trace.tolist(),
          "accepted_interval": list(map(float, selector["selected_interval"])),
          "prediction_central": central, "prediction_variants": variants,
          "prediction_range": [float(pv.min()),float(pv.max())],
          "central_label": b["packet"]["central_label"],
          "band_label": b["packet"]["band_label"],
          "band_meaning": "Observed selected numerical-method/mesh sensitivity; not a total error bound or omitted-tail interval",
          "normalization": b["packet"]["normalization"]}
    alpha = ref(b["alpha"])
    da = {k: finite([x[k] for x in alpha]).tolist() for k in ["alpha", "TT_modal", "TL_modal", "Tbulk"]}
    if len(alpha) != 21 or not np.allclose(da["alpha"],np.linspace(0,1,21),atol=1e-14,rtol=0):
        raise ValueError("Expected the declared 21-point alpha family")
    if np.max(np.abs(np.array(da["TT_modal"])+da["TL_modal"]-da["Tbulk"]))>1e-9:
        raise ValueError("Alpha branch sum disagrees with bulk power")
    da["method_label"] = b["alpha_label"]
    continuum = ref(b["continuum"])
    dc = {"labels": ["Point force S","Finite X S","Finite R S","Finite X+R S","Lattice modal qT"],
          "values": [float(continuum["point_force_S_fraction"])] +
                    [float(continuum["finite_source"][k]["S_fraction"]) for k in ["X","R","full"]] +
                    [float(r["TT_modal"]/r["Tbulk"])],
          "lattice_reference": b["reference_label"]}
    finite(dc["values"])
    return b, {"reference":d4, "packet":d5, "alpha":da, "continuum":dc}


def configure():
    mpl.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"DejaVu Sans",
        "mathtext.fontset":"dejavusans","axes.unicode_minus":True,"axes.linewidth":.9,
        "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"xtick.labelsize":9,
        "ytick.labelsize":9,"xtick.major.width":.8,"ytick.major.width":.8,
        "savefig.transparent":False,"figure.facecolor":"white","axes.facecolor":"white"})


def save(fig,path):
    if path.suffix == ".pdf":
        fig.savefig(path,metadata={"Title":"Paper1 corrected numerical evidence",
            "Author":"M. Fernandes","Subject":"Selected, source-bound numerical inputs",
            "Creator":"Paper1 current figure producer","Producer":f"Matplotlib {mpl.__version__}",
            "CreationDate":FIXED_DATE,"ModDate":FIXED_DATE},bbox_inches="tight",pad_inches=.02)
    else:
        fig.savefig(path,dpi=220,metadata={"Software":"Paper1 current figure producer"})
    plt.close(fig)


def reference_plot(d,path):
    fig,(ax,bx)=plt.subplots(1,2,figsize=(7,3),gridspec_kw={"wspace":.44})
    labels=[r"$R_q$",r"$T_q$",r"$T_{\rm qL}^{\rm modal}$",r"$T_{\rm qT}^{\rm modal}$"]
    left=0
    for value,label in zip(d["partition"].values(),labels):
        ax.barh([0],[value],left=left,height=.26,label=label);left+=value
    ax.set(xlim=(0,1),ylim=(-.32,.32),yticks=[],xlabel=r"incident energy-flux fraction at $\omega=1.2$")
    ax.set_title("(a) reference flux partition",loc="left",pad=47)
    ax.legend(ncol=2,loc="lower center",bbox_to_anchor=(.5,1.02),frameon=False,fontsize=9)
    for i,(label,v) in enumerate(zip(labels,d["partition"].values())):
        ax.text(.02+.52*(i%2),.17-.14*(i//2),f"{label} = {v:.4f}",transform=ax.transAxes,fontsize=9)
    ax.text(.5,.85,f"raw closure = {d['raw_closure']:.2e}",transform=ax.transAxes,ha="center",fontsize=9)
    ax.grid(axis="x",alpha=.18)
    values=np.array([x["value"] for x in d["comparisons"]]);x=np.arange(len(values))
    bx.axhspan(values.min(),values.max(),alpha=.15,label="central-estimate spread")
    bx.scatter(x,values,s=48,zorder=3)
    for xx,v in zip(x,values):
        ha="left" if xx==0 else "right" if xx==len(values)-1 else "center"
        dx=4 if xx==0 else -4 if xx==len(values)-1 else 0
        dy=24 if xx==1 else 10
        bx.annotate(f"{v:.8f}",(xx,v),xytext=(dx,dy),textcoords="offset points",ha=ha,fontsize=9)
    bx.set_xticks(x);bx.set_xticklabels([z["label"] for z in d["comparisons"]])
    bx.set_ylabel(r"$T_{\rm qT}^{\rm modal}$")
    bx.set_title("(b) cross-formulation\ncentral estimates",loc="left",pad=33)
    pad=max(float(np.ptp(values))*1.2,1.25e-4)
    bx.set_ylim(values.min()-pad,values.max()+pad)
    bx.text(.5,.20,f"tetrahedron difference = {d['surface_tetra_relative_percent']:.4f}%",
            transform=bx.transAxes,ha="center",fontsize=9)
    bx.legend(loc="lower center",frameon=False,fontsize=9);bx.grid(axis="y",alpha=.22)
    save(fig,path)


def packet_plot(d,path):
    fig,ax=plt.subplots(figsize=(3.6,2.65))
    line,=ax.plot(d["time"],d["qT"],marker="o",markersize=2.8,linewidth=1.3,zorder=3,label=r"$\mathcal{E}_{\rm qT}^{\rm TD}(\tau)$")
    accepted=ax.axvspan(*d["accepted_interval"],alpha=.16,zorder=0,label="accepted interval")
    band=ax.axhspan(*d["prediction_range"],alpha=.11,hatch="///",zorder=0,label=d["band_label"])
    central=ax.axhline(d["prediction_central"],linestyle="--",linewidth=1.2,zorder=2,label=d["central_label"])
    lo=min(.122,min(d["qT"]),d["prediction_range"][0]);hi=max(.129,max(d["qT"]),d["prediction_range"][1])
    ax.set(xlim=(d["time"][0],d["time"][-1]),ylim=(lo,hi),xlabel=r"Time $\tau$",ylabel="Initial-energy fraction")
    ax.grid(alpha=.22);ax.legend(handles=[line,accepted,band,central],loc="lower right",fontsize=8.1,frameon=True)
    save(fig,path)


def alpha_plot(d,path):
    fig,ax=plt.subplots(figsize=(7.2,4.6))
    for key,marker,label in [("TT_modal","o",r"$T_{\rm qT}^{\rm modal}$"),("TL_modal","s",r"$T_{\rm qL}^{\rm modal}$"),("Tbulk","^",r"$T_{\rm bulk}$")]:
        ax.plot(d["alpha"],d[key],marker=marker,markersize=3,label=label)
    ax.set(xlabel=r"full radial-sector scale $\alpha$ in $H_X+\alpha H_R$",ylabel="Incident-power fraction",title="Separately dressed surrogate sensitivity family")
    ax.legend();fig.tight_layout();save(fig,path)


def continuum_plot(d,path):
    fig,ax=plt.subplots(figsize=(7.8,4.7));ax.bar(d["labels"],d["values"])
    ax.set(ylim=(min(.84,min(d["values"])-.01),max(.94,max(d["values"])+.01)),ylabel="Generated-bulk share",
           title="Continuum S fractions and lattice modal qT share\n(analogous but nonidentical observables)")
    ax.tick_params(axis="x",rotation=18);fig.tight_layout();save(fig,path)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--binding",type=Path,required=True);ap.add_argument("--output-dir",type=Path,required=True);a=ap.parse_args()
    binding=a.binding.resolve();b,data=collect(binding)
    if a.output_dir.exists() and any(a.output_dir.iterdir()):
        ap.error("Output directory must be new or empty")
    a.output_dir.mkdir(parents=True,exist_ok=True);dd=a.output_dir/"data";dd.mkdir()
    configure()
    for key,producer in [("reference",reference_plot),("packet",packet_plot),("alpha",alpha_plot),("continuum",continuum_plot)]:
        (dd/f"{key}.json").write_text(json.dumps(data[key],indent=2,sort_keys=True)+"\n")
        producer(data[key],a.output_dir/NAMES[key])
    manifest={"schema":"paper1.current-scientific-figures.v1","status":"generated",
        "binding_sha256":sha256(binding),"binding":b,"producer_sha256":sha256(__file__),
        "environment":{"matplotlib":mpl.__version__,"numpy":np.__version__},
        "outputs":{k:{"filename":n,"sha256":sha256(a.output_dir/n),"bytes":(a.output_dir/n).stat().st_size} for k,n in NAMES.items()},
        "data":{k:{"filename":f"data/{k}.json","sha256":sha256(dd/f"{k}.json")} for k in data},
        "unaffected_figures":"Verified by hash; not rewritten"}
    (a.output_dir/"FIGURE_REGENERATION_MANIFEST.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps(manifest,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
