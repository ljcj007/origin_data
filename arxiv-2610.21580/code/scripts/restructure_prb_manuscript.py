#!/usr/bin/env python3
"""将现稿重排为方法与结果相邻的 PRB 版本。"""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "paper" / "manuscript_v2.tex"
TARGET = ROOT / "paper" / "manuscript_prb.tex"


def section(text: str, start: str, end: str) -> str:
    left = text.index(start)
    right = text.index(end, left)
    return text[left:right].strip()


def body_without_heading(block: str) -> str:
    return block.split("\n", 1)[1].strip()


def main() -> None:
    source = SOURCE.read_text(encoding="utf-8")

    preamble = section(source, "\\documentclass", "\\begin{abstract}")
    preamble = "\n".join(
        line
        for line in preamble.splitlines()
        if not line.startswith(("\\newcommand{\\todo}", "\\author{", "\\affiliation{"))
    )

    introduction = section(source, "\\section{Introduction}", "\\section{Lattice and Hamiltonians}")
    old_scope = (
        "We analyze the uniform nearest-neighbor model using exact graph and band arguments, "
        "a decoration mapping in the classical Ising limit, linear spin-wave theory (LSWT), "
        "and stochastic-series-expansion (SSE) quantum Monte Carlo. This combination separates "
        "three statements that are easily conflated: rank-nullity fixes the single-particle flat "
        "band, the Lieb--Mattis theorem fixes the ground-state spin, and LSWT predicts the "
        "finite-frequency magnetic flat band. Static SSE estimates low-temperature observables "
        "along a fixed-aspect-ratio sequence, while imaginary-time correlations and stochastic "
        "analytic continuation test the interacting excitation spectrum."
    )
    new_scope = (
        "We organize the analysis method by method, placing each numerical or analytic result "
        "next to the construction that supports it. Exact graph arguments and LSWT identify the "
        "noninteracting and harmonic constraints; an exact Ising mapping provides a classical "
        "benchmark; SSE establishes the interacting ferrimagnetic ground state; and finite-cluster "
        "spectral moments together with frequency-resolved SSE and stochastic analytic continuation "
        "test whether the harmonic flat excitation survives interactions. This separation is "
        "essential: rank--nullity protects a single-particle zero mode, the Lieb--Mattis theorem "
        "fixes the ground-state spin, but neither theorem protects a sharp finite-frequency magnon."
    )
    introduction = introduction.replace(old_scope, new_scope)

    lattice = section(source, "\\subsection{Periodic graph from the POSCAR}", "\\subsection{Single-particle model}")
    spins = section(source, "\\subsection{Spin models}", "\\section{Methods}")
    model = (
        "\\section{Lattice, Hamiltonians, and exact constraints}\n\n"
        "\\subsection{Periodic lattice}\n\n"
        + body_without_heading(lattice)
        + "\n\n\\subsection{Spin Hamiltonians and ground-state spin}\n\n"
        + body_without_heading(spins)
    )
    model = model.replace(
        "The Cartesian first-neighbor shell of the structural coordinate file contains six undirected bonds per cell at distance $1.4225911$ in the file's arbitrary length units.",
        "The Cartesian first-neighbor shell of the crystallographic coordinates contains six undirected bonds per cell at distance $1.4225911$ in the coordinate units.",
    ).replace(
        "the edge-site cell offsets in the exported graph are",
        "the edge-site cell offsets in the periodic embedding are",
    ).replace(
        "follow the site order in the graph",
        "label the sites within the primitive cell",
    )

    band_model = section(source, "\\subsection{Single-particle model}", "\\subsection{Spin models}")
    graph_method = section(source, "\\subsection{Graph, band, and spin-wave checks}", "\\subsection{Exact Ising mapping")
    flat_results = section(source, "\\subsection{One structural null space, two flat spectral sectors}", "\\subsection{Ferrimagnetic Goldstone mode}")
    goldstone_results = section(source, "\\subsection{Ferrimagnetic Goldstone mode}", "\\subsection{Exact and simulated classical transition}")
    lswt_appendix = section(source, "\\section{Analytic origin of the LSWT branches}", "\\section{Numerical validation and provenance boundary}")
    band_lswt = (
        "\\subsection{Bloch analysis and linear spin-wave theory}\n\n"
        "\\paragraph{Construction.}\n"
        + body_without_heading(band_model)
        + "\n\n"
        + body_without_heading(graph_method)
        + "\n\n\\paragraph{Method-resolved results.}\n"
        + body_without_heading(flat_results)
        + "\n\n"
        + body_without_heading(goldstone_results)
        + "\n\n\\paragraph{Analytic origin of the harmonic branches.}\n"
        + body_without_heading(lswt_appendix)
    )
    band_lswt = band_lswt.replace(
        "Nearest-neighbor bonds were selected using the full Cartesian metric of the structural coordinate file under periodic boundary conditions in the two in-plane directions. The exported graph was checked for self-links, duplicate bonds, coordination numbers, bipartiteness, and consistency after supercell replication.",
        "Nearest-neighbor bonds were identified with the full Cartesian metric under periodic boundary conditions in the two in-plane directions. The resulting graph was checked for self-links, duplicate bonds, coordination numbers, bipartiteness, and consistency under supercell replication.",
    ).replace(
        "The POSCAR-derived graph reproduces the edge-centered honeycomb connectivity.",
        "The periodic graph reproduces the edge-centered honeycomb connectivity.",
    )

    ising_method = section(source, "\\subsection{Exact Ising mapping and cluster simulation}", "\\subsection{Stochastic-series-expansion quantum Monte Carlo}")
    ising_figure = section(source, "\\begin{figure}[t]", "\\subsection{One structural null space")
    ising_results = section(source, "\\subsection{Exact and simulated classical transition}", "\\subsection{SSE fixed-aspect-ratio scaling}")
    ising = (
        "\\subsection{Decoration mapping and Wolff-cluster simulation}\n\n"
        "\\paragraph{Method.}\n"
        + body_without_heading(ising_method)
        + "\n\n"
        + ising_figure
        + "\n\n\\paragraph{Results.}\n"
        + body_without_heading(ising_results)
    )

    sse_method = section(source, "\\subsection{Stochastic-series-expansion quantum Monte Carlo}", "\\subsection{Finite-cluster spectral-moment test}")
    sse_results = section(source, "\\subsection{SSE fixed-aspect-ratio scaling}", "\\subsection{Quantum response of a compact spin probe}")
    sse = (
        "\\subsection{Stochastic-series-expansion simulation of quantum ferrimagnetism}\n\n"
        "\\paragraph{Method.}\n"
        + body_without_heading(sse_method)
        + "\n\n\\paragraph{Results.}\n"
        + body_without_heading(sse_results)
    )
    sse = sse.replace(
        "All 14 production cells passed the completion, provenance, and artifact-hash checks. Each size therefore contributes 128 pooled bins from two independent replicas.",
        "Each size contributes 128 pooled bins from two independent replicas.",
    ).replace(
        "An independent split-half audit returned \\texttt{provisional=false} for every cell and found a largest auxiliary blocking ratio of $1.47$.",
        "An independent split-half audit found a largest auxiliary blocking ratio of $1.47$.",
    ).replace(
        "no diagnostic survived the registered family-wise correction",
        "no diagnostic remained significant after the family-wise correction",
    ).replace(
        "The primary finite-size forms were registered after $L\\leq10$ had completed but before the three largest sizes were known:",
        "The primary finite-size forms were fixed after the $L\\leq10$ data were available but before the three largest sizes were known:",
    ).replace(
        "The preregistered $L=4$--16 fits gave",
        "The primary $L=4$--16 fits gave",
    ).replace(
        "the preregistered primary fits",
        "the primary fits",
    ).replace(
        "Every input graph was checked at runtime to ensure that each bond joined opposite colors and that all $J$ values were positive.",
        "Each periodic graph was verified to contain only positive couplings between opposite sublattices.",
    )

    moment_method = section(source, "\\subsection{Finite-cluster spectral-moment test}", "\\subsection{Frequency-resolved SSE and analytic continuation}")
    moment_results = section(source, "\\subsection{Quantum response of a compact spin probe}", "\\subsection{Frequency-resolved interacting spectrum}")
    moments = (
        "\\subsection{Finite-cluster exact diagonalization and spectral moments}\n\n"
        "\\paragraph{Method.}\n"
        + body_without_heading(moment_method)
        + "\n\n\\paragraph{Results.}\n"
        + body_without_heading(moment_results)
    )

    sac_method = section(source, "\\subsection{Frequency-resolved SSE and analytic continuation}", "\\subsection{Two-exchange thermal pilot}")
    sac_results = section(source, "\\subsection{Frequency-resolved interacting spectrum}", "\\subsection{Exploratory two-exchange thermodynamics}")
    sac = (
        "\\subsection{Frequency-resolved SSE and stochastic analytic continuation}\n\n"
        "\\paragraph{Method.}\n"
        + body_without_heading(sac_method)
        + "\n\n\\paragraph{Results.}\n"
        + body_without_heading(sac_results)
    )

    two_exchange_method = section(source, "\\subsection{Two-exchange thermal pilot}", "\\section{Results}")
    two_exchange_results = section(source, "\\subsection{Exploratory two-exchange thermodynamics}", "\\section{Discussion}")
    star = section(source, "\\section{Exact decoupled-star thermodynamics}", "\\bibliographystyle")
    unequal = (
        "\\subsection{Unequal-exchange extension and its decoupled-star limit}\n\n"
        "\\paragraph{Method.}\n"
        + body_without_heading(two_exchange_method)
        + "\n\n\\paragraph{Exact reference limit.}\n"
        + body_without_heading(star)
        + "\n\n\\paragraph{Exploratory results.}\n"
        + body_without_heading(two_exchange_results)
    )
    unequal = unequal.replace(
        "Downloaded inputs, source and executable hashes, return codes, and raw-bin identities were checked before recomputing observables. For energy, uniform susceptibility, and expansion order, diagnostics flagged",
        "For energy, uniform susceptibility, and expansion order, statistical diagnostics flagged",
    ).replace(
        "under the registered rule $|\\Delta|\\leq5\\,{\\rm SE}+10^{-4}$. Only after this gate passed did we run the lattice pilot:",
        "under the criterion $|\\Delta|\\leq5\\,{\\rm SE}+10^{-4}$. Only validated parameter sets were used for the lattice pilot:",
    ).replace(
        "Passing this gate did not resolve every low-temperature response:",
        "This validation did not resolve every low-temperature response:",
    ).replace(
        "This post hoc diagnostic did not replace the registered gate.",
        "This additional diagnostic indicates the remaining low-temperature sampling limitation.",
    )
    band_lswt = band_lswt.replace(
        "the four dispersive positive frequencies in Sec.~IV~A",
        "the four dispersive positive frequencies above",
    )

    methods = (
        "\\section{Methods and method-resolved results}\n\n"
        + band_lswt
        + "\n\n"
        + ising
        + "\n\n"
        + sse
        + "\n\n"
        + moments
        + "\n\n"
        + sac
        + "\n\n"
        + unequal
    )
    methods = methods.replace(
        "The displayed free-spectrum continuation is accompanied by independent-seed, dominant-delta, and rebin/time-grid sensitivity analyses in the archived results.",
        "Independent-seed, dominant-delta, and rebin/time-grid variants quantify continuation sensitivity.",
    )

    discussion = r"""\section{Discussion}

The calculations expose a common geometric origin for phenomena that otherwise look unrelated. The $2:3$ bipartite imbalance makes the intersublattice block rectangular. In the hopping problem, rank--nullity then produces an edge-supported zero mode at every momentum. In the spin problem, the same imbalance fixes an extensive ground-state spin through the Lieb--Mattis theorem. LSWT inherits the null vector as an edge-supported mode, but the local exchange field shifts it to the finite frequency $2JS$. Thus the single-particle flat band, ferrimagnetic moment, and harmonic flat magnon are linked by the same lattice constraint while remaining logically distinct statements.

The interacting results locate the boundary of that connection. Static SSE follows the exact total-spin constraint and yields nonzero staggered order and stiffness, establishing a quantum ferrimagnet with a quadratic type-B Goldstone mode. The dynamical probes, however, show that the harmonic flat line is not protected as a sharp many-body excitation. A compact spin operator has nonzero spectral variance already on the 20-site torus, and the larger-cluster SSE--SAC response carries broad finite-frequency weight. Its centroid disperses only weakly, preserving a memory of destructive interference, whereas the continuation-sensitive maximum and $O(J)$ spectral width preclude a quasiparticle interpretation. Sublattice imbalance therefore organizes spectral weight without forbidding interaction-induced hybridization and broadening.

This distinction separates the present mechanism from kagome antiferromagnetism. Kagome flat modes arise from frustrated triangular constraints and accompany a macroscopically degenerate classical manifold~\cite{harris1992,maksymenko2017}. Here every nearest-neighbor bond can be satisfied classically, the lattice admits a Marshall sign convention, and the imbalance enforces ferrimagnetism. The Goldstone mode is quadratic rather than the linear mode of a compensated collinear antiferromagnet, and the localized harmonic branch occurs at finite frequency rather than at the bottom of a frustrated zero-mode manifold. The shared hexagonal compact pattern is therefore geometric, not evidence of the same magnetic phase.

The classical Ising calculation provides an independent benchmark of the graph rather than a finite-temperature prediction for the quantum model. The exact decoration transformation and Wolff crossings recover the two-dimensional Ising transition and its critical exponents. By contrast, the isotropic short-range Heisenberg model cannot spontaneously break continuous spin symmetry at nonzero temperature~\cite{mermin1966}; its transition is at $T=0$. The quadratic mode nevertheless controls a broad low-temperature regime. For $\omega=Dq^2$, with $D=JSa^2/4$, LSWT gives $u_{\rm th}/(NJ)=[\pi\sqrt3/(60S)]t^2+o(t^2)$ and $C/(Nk_B)=[\pi\sqrt3/(30S)]t+o(t)$, where $t=k_BT/J$. These crossover forms do not imply finite-temperature order, and the infrared magnon-number divergence makes interaction corrections unavoidable.

The unequal-exchange extension clarifies how the uniform ferrimagnet can be approached from local building blocks. In the $J_B=0$ limit, a four-spin star forms a spin-1 triplet while the second hub remains free; weak coupling generates an alternating spin-1--spin-$1/2$ honeycomb model with $J_{\rm eff}=5J_B/12$. The exact star thermodynamics predicts the intermediate-window entropy $\ln6/5$ and $T\chi=11/60$, as well as the leading correction $-J_B/(12T)$. The pilot correlations are consistent with strong-bond cluster formation preceding weaker interstar correlations, but the heat capacity, entropy, and thermodynamic crossover scales are not resolved. For every connected cluster with $J_A,J_B>0$, the Lieb--Mattis spin remains $S_0=L^2/2$; tuning the bond ratio alone does not create a singlet ground state.

There is also a distinct exact flat mode in the fully polarized phase. The one-magnon graph Laplacian has a flat eigenvalue $2$ and maximum eigenvalue $5$, giving $\epsilon_{\rm flat}/J=h/J-2S$ and $h_{\rm sat}/J=5S$. This result is protected by the one-magnon reduction and should not be conflated with the zero-field continuum found by SSE--SAC.

Photonic super-honeycomb lattices already realize the five-site connectivity and its compact single-particle states~\cite{zhong2017,yan2020}, and an optical-lattice proposal provides the same geometry~\cite{lan2012}. Realizing the spin model requires an interacting Mott regime with dominant, nearly uniform antiferromagnetic nearest-neighbor exchange. Existing magnetic candidates instead involve kagome frustration, unequal spin lengths, anisotropy, interlayer coupling, or multi-orbital electronic structure~\cite{misumi2020,pratt2025,yoruk2025,lobocheca2024,wang2023,frezza2024}. The present model is therefore best viewed as a controlled design principle: sublattice imbalance can combine robust ferrimagnetism with localized spectral weight, while interactions determine whether that weight forms a sharp mode or a continuum.
"""

    conclusion = r"""\section{Conclusion}

Sublattice imbalance provides a single organizing principle for the edge-centered honeycomb antiferromagnet. It enforces a single-particle flat band, fixes an extensive ferrimagnetic ground-state spin, and produces an edge-localized flat branch in harmonic spin-wave theory. The many-body dynamics reveals the essential qualification: the localized geometry survives as a weakly dispersive finite-frequency spectral sector, but not as a protected sharp magnon. This conversion of a harmonic flat band into a structured interacting continuum distinguishes imbalance-driven localization from frustration-induced kagome flat modes and identifies a route to ferrimagnets in which robust static order coexists with strongly nontrivial localized dynamics.
"""

    abstract = r"""\begin{abstract}
Flat bands and ferrimagnetism are commonly attributed to different mechanisms---destructive interference for the former and unequal magnetic moments for the latter. We show that both can instead descend from a single lattice constraint. In the spin-$1/2$ edge-centered honeycomb antiferromagnet, a $2:3$ bipartite imbalance simultaneously enforces a single-particle flat band and an extensive Lieb--Mattis ground-state spin. Linear spin-wave theory converts the same compact null vector into a finite-frequency flat magnon above a ferrimagnetic background, while unbiased quantum Monte Carlo establishes the corresponding interacting order and quadratic Goldstone dynamics. Crucially, frequency-resolved simulations and exact spectral moments show that interactions transform the harmonic flat line into a broad, only weakly dispersive spectral sector rather than a protected quasiparticle. Sublattice imbalance therefore unifies flat-band localization and ferrimagnetism without making their dynamical signatures equivalent, providing a design principle for unfrustrated magnets with robust order and localized many-body spectral weight.
\end{abstract}"""

    availability = r"""\section*{Data availability}

The data supporting the conclusions of this study and the numerical implementation are available from the corresponding author upon reasonable request.
"""

    output = "\n\n".join(
        [
            preamble.rstrip(),
            abstract,
            "\\maketitle",
            introduction,
            model,
            methods,
            discussion,
            conclusion,
            availability,
            "\\bibliographystyle{apsrev4-2}\n\\bibliography{references_prb_verified}\n\n\\end{document}",
        ]
    )
    TARGET.write_text(output + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
