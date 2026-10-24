"""用全谱独立复核小簇能量矩，拒绝大矩阵。"""
import argparse
import json
from pathlib import Path
import numpy as np
import ed_cls_moments as probe


def audit(path, graph):
    data=json.loads(path.read_text(encoding="utf-8"))
    assert data["nsite"] <= 10
    assert data["source_sha256"] == probe.base.sha256_file(Path(probe.__file__))
    assert data["base_source_sha256"] == probe.base.sha256_file(Path(probe.base.__file__))
    assert data["graph_sha256"] == probe.base.sha256_file(graph)
    n,bonds,_,_,_=probe.base.make_lattice(graph,data["lx"],data["ly"])
    nup=int(n/2+data["ground_spin"])
    h,states=probe.sector(n,bonds,nup)
    energies,vectors=np.linalg.eigh(h.toarray())
    ground=vectors[:,0]
    assert abs(energies[0]-data["ground_energy"]) < 1e-10
    results={}
    for channel,delta in [("minus",-1),("z",0),("plus",1)]:
        ht,target=probe.sector(n,bonds,nup+delta)
        e,u=np.linalg.eigh(ht.toarray())
        v=probe.apply_spin(states,ground,target,np.array(data["cls_weights"]),channel)
        strength=np.abs(u.T@v)**2
        omega=e-energies[0]
        excited=omega>1e-9
        norm=float(strength[excited].sum())
        mean=float(np.dot(strength[excited],omega[excited])/norm)
        variance=float(np.dot(strength[excited],(omega[excited]-mean)**2)/norm)
        recorded=data["channels"][channel]["local_cls"]
        errors=[abs(norm-recorded["inelastic_weight"]),
                abs(mean-recorded["mean_energy"]),abs(variance-recorded["energy_variance"])]
        assert max(errors)<1e-10
        poles=[]
        for energy,weight in zip(omega[excited],strength[excited]/norm):
            if weight < 1e-12:
                continue
            if poles and abs(poles[-1]["energy"]-energy)<1e-9:
                poles[-1]["normalized_weight"]+=float(weight)
            else:
                poles.append(dict(energy=float(energy),normalized_weight=float(weight)))
        results[channel]=dict(max_moment_error=max(errors),poles=poles,
                              retained_normalized_weight=sum(p["normalized_weight"] for p in poles))
    return dict(status="passed",method="all eigenpairs; dense eigh; independently formed Lehmann sums",
                result_sha256=probe.base.sha256_file(path),
                audit_source_sha256=probe.base.sha256_file(Path(__file__)),channels=results)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result",type=Path,required=True)
    parser.add_argument("--graph",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    probe.base.require_new_output(args.output)
    result=audit(args.result,args.graph)
    with args.output.open("x",encoding="utf-8") as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write("\n")
    print(json.dumps(result,indent=2))
