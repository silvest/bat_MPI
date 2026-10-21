#!/usr/bin/env python3
"""Production-only MCMC diagnostics for BAT ROOT output, including PT transport.

Uses ArviZ's rank-normalized split/folded R-hat, bulk/tail ESS and raw-mean
ESS/MCSE. No reimplementation of the statistical estimators is used.
Region conditions are data (no Python expression evaluation).
"""
import argparse
import importlib.util
import csv
import json
import math
import operator
from pathlib import Path
import sys
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message=r"\s*ArviZ is undergoing", category=FutureWarning)
    import arviz as az
import numpy as np
import uproot

# Load beside this script, including when imported by an external test runner.
_pt_spec = importlib.util.spec_from_file_location('bat_pt_diagnostics', Path(__file__).with_name('pt_diagnostics.py'))
_pt_module = importlib.util.module_from_spec(_pt_spec)
_pt_spec.loader.exec_module(_pt_module)
analyze_streams = _pt_module.analyze_streams


def analyze_pt_only(args):
    prefixes = args.transport_prefix or [str(p)[:-len('cold_chains.root')] for p in args.inputs if str(p).endswith('cold_chains.root')]
    if not prefixes:
        raise ValueError('--pt-only requires --transport-prefix (or a conventional cold_chains.root filename)')
    stream = analyze_streams(prefixes, args.output, args.no_plots, args.trace_points)
    status = 'attention' if stream['issues'] else 'transport_checked' if stream['available'] else 'incomplete'
    report = dict(status=status, streaming_pt=stream,
                  note='PT transport only: cold-chain R-hat, ESS and posterior region weights were not checked.')
    output = Path(args.output)
    (output/'summary.json').write_text(json.dumps(finite_json(report), indent=2, allow_nan=False)+'\n')
    lines = [f'PT transport: {status}', *stream['issues'], report['note']]
    (output/'summary.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
    return report

RESERVED = {"Chain", "Iteration", "Phase", "Replica", "Temperature"}
OPS = {"<": operator.lt, "<=": operator.le, ">": operator.gt,
       ">=": operator.ge, "==": operator.eq}
TRANSFORMS = {"identity": lambda x: x, "cos": np.cos, "sin": np.sin}


def text(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def load_cold(paths, tree_name=None):
    """Keep chain boundaries and ordering; reject mismatched or irregular draws."""
    groups, fixed_sets, columns = [], [], None
    for file_index, path in enumerate(paths):
        with uproot.open(path) as root:
            names = [k.split(";")[0] for k in root.keys()
                     if k.split(";")[0].endswith("_mcmc")]
            if tree_name is None and len(set(names)) != 1:
                raise ValueError(f"{path}: specify --tree (found {names})")
            selected = tree_name or names[0]
            tree = root[selected]
            arrays = tree.arrays(library="np")
            if not isinstance(arrays, dict):
                arrays = {name: arrays[name] for name in arrays.dtype.names}
            if not {"Chain", "Iteration", "Phase"}.issubset(arrays):
                raise ValueError(f"{path}: Chain, Iteration and Phase branches are required")
            numeric = sorted(name for name, value in arrays.items()
                             if name not in RESERVED and value.ndim == 1
                             and np.issubdtype(value.dtype, np.number))
            if columns is not None and numeric != columns:
                raise ValueError("Input files have different variable columns")
            columns = numeric
            fixed = set()
            parameter_tree = selected.removesuffix("_mcmc") + "_parameters"
            if parameter_tree in root:
                metadata = root[parameter_tree]
                required = {"parameter", "safe_name", "fixed"}
                # BAT also writes unrelated empty C-string branches. Read only the
                # metadata needed here; older ROOT empty-string baskets can fail
                # decoding in uproot even though the required branches are valid.
                meta = metadata.arrays(sorted(required), library="np") if required.issubset(metadata.keys()) else {}
                if required.issubset(meta):
                    fixed = {text(name) for name, par, flag in
                             zip(meta["safe_name"], meta["parameter"], meta["fixed"])
                             if par and flag}
            fixed_sets.append(fixed)
            production = arrays["Phase"] == 1
            if not np.any(production):
                raise ValueError(f"{path}: no production samples")
            chains = sorted(np.unique(arrays["Chain"][production]))
            for chain in chains:
                mask = production & (arrays["Chain"] == chain)
                iterations = arrays["Iteration"][mask]
                steps = np.diff(iterations)
                if len(steps) and (np.any(steps <= 0) or np.any(steps != steps[0])):
                    raise ValueError(f"{path}: chain {chain} has duplicate, unsorted or irregular iterations")
                groups.append((f"{file_index}:{int(chain)}", iterations,
                               {name: arrays[name][mask] for name in numeric}))
    sizes = {len(group[1]) for group in groups}
    if len(sizes) != 1:
        raise ValueError("Unequal chain lengths: analyze consistent production windows; no automatic truncation")
    if not groups or next(iter(sizes)) == 0:
        raise ValueError("No cold samples")
    # Every chain in a BAT cold-only tree is an independent ladder. Catch accidentally
    # providing a combined hot/cold tree if it has an explicit temperature branch.
    for path in paths:
        with uproot.open(path) as root:
            name = tree_name or next(k.split(";")[0] for k in root.keys() if k.split(";")[0].endswith("_mcmc"))
            if "Temperature" in root[name]:
                a = root[name].arrays(["Phase", "Temperature"], library="np")
                if np.any(a["Temperature"][a["Phase"] == 1] != 1):
                    raise ValueError("Input includes heated samples; provide cold-only output")
    fixed = set.intersection(*fixed_sets) if fixed_sets else set()
    data = {name: np.stack([group[2][name] for group in groups]) for name in columns}
    # Fixed metadata must agree with observed values; don't mask a broken output file.
    for name in fixed:
        if name in data and not np.all(data[name] == data[name].flat[0]):
            raise ValueError(f"Declared fixed parameter {name} varies in the output")
    return data, [g[0] for g in groups], np.stack([g[1] for g in groups]), fixed


def diagnostic(values, fixed=False, rhat_limit=1.01, ess_per_chain=100):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or not x.size:
        raise ValueError("Expected nonempty chain x draw array")
    m, n = x.shape
    out = dict(chains=m, draws_per_chain=n, mean=float(np.mean(x)),
               sd=float(np.std(x, ddof=1)) if x.size > 1 else math.nan,
               rhat=math.nan, ess_bulk=math.nan, ess_tail=math.nan,
               ess_mean=math.nan, mcse_mean=math.nan, status="ok")
    if not np.isfinite(x).all():
        out["status"] = "nonfinite"
        return out
    if fixed:
        out.update(status="fixed", mcse_mean=0.)
        return out
    if n < 4:
        out["status"] = "insufficient_draws"
        return out
    if np.all(x == x.flat[0]):
        out["status"] = "constant_unverified"
        return out
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if m >= 2:
            out["rhat"] = float(az.rhat(x, method="rank"))
        for method in ("bulk", "tail", "mean"):
            out["ess_" + method] = float(az.ess(x, method=method))
        out["mcse_mean"] = float(az.mcse(x, method="mean"))
    flags = []
    if m < 2:
        flags.append("insufficient_chains")
    elif not math.isfinite(out["rhat"]) or out["rhat"] > rhat_limit:
        flags.append("rhat_high_or_undefined")
    if np.any(np.all(x == x[:, :1], axis=1)):
        flags.append("stuck_chain")
        # An unchanging chain is not evidence for an exactly known expectation.
        out["mcse_mean"] = math.nan
    for key in ("ess_bulk", "ess_tail"):
        if not math.isfinite(out[key]) or out[key] < ess_per_chain * m:
            flags.append("low_" + key)
    if not math.isfinite(out["mcse_mean"]):
        flags.append("mcse_unavailable")
    out["status"] = ";".join(flags) or "ok"
    return out


def condition_mask(condition, data):
    if set(condition) - {"parameter", "transform", "op", "value"}:
        raise ValueError(f"Unknown condition fields: {condition}")
    name = condition["parameter"]
    if name not in data:
        raise KeyError(name)
    transform = condition.get("transform", "identity")
    op = condition.get("op", ">=")
    if transform not in TRANSFORMS or op not in OPS:
        raise ValueError("Unsupported transform/operator in region condition")
    value = float(condition["value"])
    if not math.isfinite(value) or not np.isfinite(data[name]).all():
        raise ValueError("Region conditions require finite values and samples")
    return OPS[op](TRANSFORMS[transform](data[name]), value)


def classify_regions(config, data):
    """Return (partition, region, mask). Named regions must not overlap.

    Binary partitions are independent views, not mutually exclusive modes
    across partitions. An explicit disjoint 'regions' list defines joint modes.
    """
    if not isinstance(config, dict):
        raise ValueError("Region configuration must be a JSON object")
    for key in ("partitions", "regions"):
        if key in config and (not isinstance(config[key], list) or any(
            not isinstance(item, dict) or not isinstance(item.get("name"), str) or
            not isinstance(item.get("conditions"), list) or
            any(not isinstance(c, dict) for c in item["conditions"]) for item in config[key])):
            raise ValueError(f"{key} must be a list of named objects with condition lists")
    result, skipped = [], []
    shape = next(iter(data.values())).shape
    seen = set()
    for partition in config.get("partitions", []):
        name = partition["name"]
        if name in seen or name == "regions":
            raise ValueError("Duplicate/reserved partition name")
        seen.add(name)
        conditions = partition["conditions"]
        if not conditions:
            raise ValueError("Empty partition conditions")
        try:
            mask = np.logical_and.reduce([condition_mask(c, data) for c in conditions])
        except KeyError as error:
            if not partition.get("optional", False):
                raise ValueError(f"Partition {name}: missing parameter {error}") from error
            skipped.append(name)
            continue
        result.extend([(name, partition.get("true_label", "true"), mask),
                       (name, partition.get("false_label", "false"), ~mask)])
    coverage = np.zeros(shape, dtype=bool)
    names = set()
    for region in config.get("regions", []):
        name = region["name"]
        if name in names or name == "other":
            raise ValueError("Duplicate/reserved region name")
        names.add(name)
        if not region["conditions"]:
            raise ValueError("Empty region conditions")
        try:
            mask = np.logical_and.reduce([condition_mask(c, data) for c in region["conditions"]])
        except KeyError as error:
            raise ValueError(f"Region {name}: missing parameter {error}") from error
        if np.any(mask & coverage):
            raise ValueError("Named regions overlap; define disjoint conditions")
        coverage |= mask
        result.append(("regions", name, mask))
    if names and np.any(~coverage):
        result.append(("regions", "other", ~coverage))
    return result, skipped


def region_diagnostic(mask, **options):
    out = diagnostic(mask.astype(float), **options)
    transitions = np.count_nonzero(np.diff(mask.astype(int), axis=1), axis=1)
    out.update(weight=out.pop("mean"), ess_indicator=out.pop("ess_mean"),
               mcse_weight=out.pop("mcse_mean"), transitions=int(transitions.sum()),
               min_chain_transitions=int(transitions.min()))
    if not mask.any():
        out["status"] = "not_visited"
    elif mask.all():
        out["status"] = "no_transitions"
    elif transitions.min() == 0:
        out["status"] = "chain_without_region_transitions"
        out["mcse_weight"] = math.nan
    return out


def csv_write(path, rows, fields=None):
    rows = list(rows)
    if not rows and fields is None:
        Path(path).write_text("")
        return
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def finite_json(value):
    if isinstance(value, dict):
        return {str(k): finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(v) for v in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_transport(prefixes):
    """Preserve source identity when reporting several independent runs."""
    roundtrips, swaps, missing = [], [], []
    for run, prefix in enumerate(prefixes):
        for suffix, destination in [("pt_roundtrips.csv", roundtrips),
                                    ("pt_swaps_by_ladder.csv", swaps)]:
            path = Path(str(prefix) + suffix)
            if not path.exists():
                missing.append(str(path))
                continue
            with path.open() as stream:
                for row in csv.DictReader(stream):
                    destination.append(dict(run=run, **row))
        # Older runs only contain aggregate swap counts.
        aggregate = Path(str(prefix) + "pt_swaps.csv")
        if not Path(str(prefix) + "pt_swaps_by_ladder.csv").exists() and aggregate.exists():
            with aggregate.open() as stream:
                swaps.extend(dict(run=run, ladder="aggregate", **row) for row in csv.DictReader(stream))
    return roundtrips, swaps, missing


def make_plots(outdir, data, iterations, chain_names, regions, cumulative, transport):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    chosen = [name for name in data if not np.all(data[name] == data[name].flat[0])][:8] or list(data)[:8]
    fig, axes = plt.subplots(len(chosen), 1, figsize=(10, max(3, 2*len(chosen))), squeeze=False)
    stride = max(1, iterations.shape[1] // 2000)
    for ax, name in zip(axes[:, 0], chosen):
        for c, chain in enumerate(chain_names):
            ax.plot(iterations[c, ::stride], data[name][c, ::stride], lw=.6, label=chain)
        ax.set_ylabel(name, fontsize=8)
    axes[0, 0].legend(title="run:chain", fontsize=7, ncol=min(4, len(chain_names)))
    axes[-1, 0].set_xlabel("Recorded chain iteration")
    fig.tight_layout(); fig.savefig(outdir / "traces.png", dpi=140); plt.close(fig)
    rank_variables = [name for name in chosen if not np.all(data[name] == data[name].flat[0])]
    if rank_variables:
        axes = az.plot_rank(az.from_dict(posterior={name: data[name] for name in rank_variables}),
                            var_names=rank_variables, kind="bars", show=False)
        figure = np.asarray(axes).flat[0].get_figure()
        figure.tight_layout(); figure.savefig(outdir / "rank_plots.png", dpi=140); plt.close(figure)
    if regions:
        chosen_regions = regions[:8]
        fig, axes = plt.subplots(len(chosen_regions), 2, figsize=(12, max(3, 2*len(chosen_regions))), squeeze=False)
        for (left, right), (partition, name, mask) in zip(axes, chosen_regions):
            for c, chain in enumerate(chain_names):
                left.step(iterations[c, ::stride], mask[c, ::stride].astype(int)+c*1.3, where="post", lw=.6)
                rows = [row for row in cumulative if row["partition"] == partition and row["region"] == name and row["chain"] == chain]
                right.plot([r["iteration"] for r in rows], [r["weight"] for r in rows], label=chain)
            left.set_ylabel(f"{partition}\n{name}", fontsize=7)
            right.set_ylim(-.02, 1.02)
        axes[0, 0].set_title("Region indicator (chains offset)")
        axes[0, 1].set_title("Cumulative region weight")
        axes[0, 1].legend(title="run:chain", fontsize=7)
        axes[-1, 0].set_xlabel("Recorded chain iteration")
        axes[-1, 1].set_xlabel("Recorded chain iteration")
        fig.tight_layout(); fig.savefig(outdir / "region_traces.png", dpi=140); plt.close(fig)
    trips, swaps = transport
    if swaps:
        fig, ax = plt.subplots(figsize=(10, 4))
        keys = sorted({(r["run"], r["ladder"]) for r in swaps}, key=str)
        for run, ladder in keys:
            rows = [r for r in swaps if (r["run"], r["ladder"]) == (run, ladder)]
            ax.plot([float(r["temperature_low"]) for r in rows],
                    [float(r["acceptance_rate"]) for r in rows], "o-", label=f"{run}:{ladder}")
        ax.set(xscale="log", xlabel="Lower temperature of adjacent pair", ylabel="Swap acceptance", ylim=(0, 1))
        ax.legend(title="run:ladder", fontsize=7); fig.tight_layout()
        fig.savefig(outdir / "swap_acceptance.png", dpi=140); plt.close(fig)
    if trips:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.bar(np.arange(len(trips)), [int(r["round_trips"]) for r in trips])
        ax.set(xlabel="Configuration ID (runs concatenated)", ylabel="Completed cold–hot–cold trips")
        fig.tight_layout(); fig.savefig(outdir / "round_trips.png", dpi=140); plt.close(fig)


def analyze(args):
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    if getattr(args, 'pt_only', False):
        return analyze_pt_only(args)
    data, chain_names, iterations, fixed = load_cold(args.inputs, args.tree)
    if not data:
        raise ValueError("No numerical variables")
    config = json.loads(Path(args.regions).read_text()) if args.regions else {}
    regions, skipped = classify_regions(config, data)
    variables = args.var or list(data)
    missing = set(variables) - set(data)
    if missing:
        raise ValueError(f"Missing requested variables: {sorted(missing)}")
    options = dict(rhat_limit=args.rhat_limit, ess_per_chain=args.ess_per_chain)
    monitored, summaries = {}, []
    for name in variables:
        monitored[name] = data[name]
        summaries.append(dict(variable=name, **diagnostic(data[name], name in fixed, **options)))
        if name.endswith("_arg") or name in args.angular:
            for transform in ("sin", "cos"):
                derived = f"{transform}({name})"
                monitored[derived] = TRANSFORMS[transform](data[name])
                summaries.append(dict(variable=derived, **diagnostic(monitored[derived], name in fixed, **options)))
    region_rows, ladder_rows, blocks, cumulative = [], [], [], []
    m, n = iterations.shape
    for partition, name, mask in regions:
        diag = region_diagnostic(mask, **options)
        if math.isfinite(diag["mcse_weight"]) and diag["mcse_weight"] > args.region_mcse:
            diag["status"] = ("" if diag["status"] == "ok" else diag["status"] + ";") + "region_mcse_high"
        region_rows.append(dict(partition=partition, region=name, **diag))
        for c, chain in enumerate(chain_names):
            ladder_rows.append(dict(partition=partition, region=name, chain=chain,
                                   **region_diagnostic(mask[c:c+1], **options)))
            for k, indices in enumerate(np.array_split(np.arange(n), min(args.blocks, n))):
                blocks.append(dict(partition=partition, region=name, chain=chain, block=k,
                                   first_iteration=int(iterations[c, indices[0]]),
                                   last_iteration=int(iterations[c, indices[-1]]),
                                   draws=len(indices), weight=float(mask[c, indices].mean())))
            indices = np.unique(np.linspace(0, n-1, min(n, args.trace_points), dtype=int))
            weights = np.cumsum(mask[c]) / np.arange(1, n+1)
            for j in indices:
                cumulative.append(dict(partition=partition, region=name, chain=chain,
                                       iteration=int(iterations[c, j]), weight=float(weights[j])))
    prefixes = args.transport_prefix
    if not prefixes:
        prefixes = [str(p)[:-len("cold_chains.root")] for p in args.inputs if str(p).endswith("cold_chains.root")]
    streaming_pt = analyze_streams(prefixes, outdir, args.no_plots, args.trace_points)
    trips, swaps, missing_transport = read_transport(prefixes)
    transport_flags = []
    for run, ladder in sorted({(r["run"], r["ladder"]) for r in trips}):
        count = sum(int(r["round_trips"]) for r in trips if (r["run"], r["ladder"]) == (run, ladder))
        if not count:
            transport_flags.append(f"run {run} ladder {ladder}: no completed round trips")
    for run, ladder in sorted({(r["run"], r["ladder"]) for r in trips}):
        rows = [r for r in trips if (r["run"], r["ladder"]) == (run, ladder)]
        missing = sum(int(r["round_trips"]) == 0 for r in rows)
        if missing and missing < len(rows):
            transport_flags.append(f"run {run} ladder {ladder}: {missing} configurations without completed round trips")
    for row in swaps:
        attempts, accepts = int(row["attempts"]), int(row["accepted"])
        if attempts == 0 or accepts == 0:
            transport_flags.append(f"run {row['run']} ladder {row['ladder']}: no swaps for {row['temperature_low']}–{row['temperature_high']}")
    transport_flags.extend(streaming_pt["issues"])
    csv_write(outdir / "diagnostics.csv", summaries)
    csv_write(outdir / "region_weights.csv", region_rows)
    csv_write(outdir / "region_weights_by_chain.csv", ladder_rows)
    csv_write(outdir / "region_weights_by_block.csv", blocks)
    csv_write(outdir / "region_weights_cumulative.csv", cumulative)
    csv_write(outdir / "round_trips.csv", trips)
    csv_write(outdir / "swaps.csv", swaps)
    issues = [dict(variable=r["variable"], status=r["status"]) for r in summaries if r["status"] not in ("ok", "fixed")]
    issues += [dict(region=f"{r['partition']}/{r['region']}", status=r["status"]) for r in region_rows if r["status"] != "ok"]
    report = dict(status="attention" if issues or transport_flags else ("incomplete" if not regions or missing_transport else "screening_passed"),
                  inputs=[str(p) for p in args.inputs], chains=chain_names, draws_per_chain=n,
                  thresholds=dict(rhat=args.rhat_limit, ess_per_chain=args.ess_per_chain, region_mcse=args.region_mcse),
                  arviz_version=az.__version__, issues=issues, transport_issues=transport_flags,
                  streaming_pt=streaming_pt, missing_transport=missing_transport, skipped_optional_partitions=skipped,
                  regions_config=config, note="These checks cannot establish that every relevant mode was discovered. Missing transport or regions means those checks were not performed.")
    (outdir / "summary.json").write_text(json.dumps(finite_json(report), indent=2, allow_nan=False) + "\n")
    lines = [f"Diagnostic screening: {report['status']}", f"Cold chains: {m}; production draws per chain: {n}",
             f"Flagged variables/regions: {len(issues)}", f"Transport flags: {len(transport_flags)}"]
    if not regions:
        lines.append("Region occupancy NOT checked: provide --regions with scientifically relevant definitions.")
    if missing_transport:
        lines.append("Some transport files unavailable: old ROOT output cannot reconstruct round trips.")
    lines.extend(f"  {issue}" for issue in transport_flags)
    lines.append(report["note"])
    (outdir / "summary.txt").write_text("\n".join(lines) + "\n")
    # Never leave plots from a different, previous configuration in this report.
    for name in ("traces.png", "rank_plots.png", "region_traces.png", "swap_acceptance.png", "round_trips.png"):
        (outdir / name).unlink(missing_ok=True)
    if not args.no_plots:
        flagged_names = [r["variable"] for r in summaries if r["status"] not in ("ok", "fixed")]
        plot_order = list(dict.fromkeys(flagged_names + list(monitored)))
        make_plots(outdir, {name: monitored[name] for name in plot_order}, iterations, chain_names, regions, cumulative, (trips, swaps))
    print("\n".join(lines))
    return report


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", nargs="*", type=Path, help="Independent cold-chain ROOT files")
    p.add_argument("--pt-only", action="store_true", help="Analyze committed streaming PT files without opening cold ROOT output")
    p.add_argument("--output", required=True, help="Directory for diagnostic tables/plots")
    p.add_argument("--tree", help="MCMC TTree name (auto-detected if unique)")
    p.add_argument("--regions", help="JSON with named disjoint regions and/or binary partitions")
    p.add_argument("--var", action="append", help="Monitor only this variable (repeatable)")
    p.add_argument("--angular", action="append", default=[], help="Additional radian variable for sin/cos diagnostics (_arg is automatic)")
    p.add_argument("--transport-prefix", action="append", default=[], help="Run prefix before pt_roundtrips.csv etc.; repeat for independent runs")
    p.add_argument("--blocks", type=int, default=10)
    p.add_argument("--trace-points", type=int, default=2000)
    p.add_argument("--rhat-limit", type=float, default=1.01)
    p.add_argument("--ess-per-chain", type=float, default=100.)
    p.add_argument("--region-mcse", type=float, default=.01)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--fail-on-warning", action="store_true", help="Exit 2 for flagged diagnostics (ordinary mode exits 0 after producing a report)")
    return p


def main():
    p = parser()
    args = p.parse_args()
    if not args.inputs and not args.pt_only:
        p.error("Provide cold-chain ROOT inputs or use --pt-only")
    if args.blocks < 1 or args.trace_points < 2 or not math.isfinite(args.rhat_limit) or args.rhat_limit < 1 or not math.isfinite(args.ess_per_chain) or args.ess_per_chain <= 0 or not math.isfinite(args.region_mcse) or args.region_mcse <= 0:
        p.error("Invalid diagnostic limits/block/trace settings")
    try:
        report = analyze(args)
    except (ValueError, KeyError, OSError) as error:
        p.exit(1, f"Diagnostic input error: {error}\n")
    if args.fail_on_warning and report["status"] != "screening_passed":
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
