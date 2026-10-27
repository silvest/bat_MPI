"""Analyze committed streaming PT diagnostics without retaining full run histories.

Counts are event-level sampler counters, not inferred from decimated traces.
Hot occupancies describe tempered targets and are never pooled into the posterior.
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

FILES = {
    "pt_local_moves.csv": "pt_local_move_blocks.csv",
    "pt_region_occupancy.csv": "pt_occupancy_blocks.csv",
    "pt_region_transitions.csv": "pt_transition_blocks.csv",
    "pt_swaps_by_block.csv": "pt_swap_blocks.csv",
    "pt_roundtrip_regions.csv": "pt_region_return_blocks.csv",
}


def rows(path):
    """Ignore an unterminated final physical line left by an interrupted writer."""
    with Path(path).open(newline="") as stream:
        def complete_lines():
            for line in stream:
                if not line.endswith('\n'):
                    return
                yield line
        yield from csv.DictReader(complete_lines())


def number(row, key, integer=False):
    try:
        return int(row[key]) if integer else float(row[key])
    except (KeyError, ValueError, TypeError) as exc:
        raise ValueError(f"Invalid PT field {key}: {row}") from exc


def write_table(path, records):
    records = list(records)
    if not records:
        Path(path).write_text("")
        return
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def checkpoints(prefix):
    phases = {}
    for row in rows(str(prefix) + 'pt_checkpoints.csv'):
        if any(value is None for value in row.values()) or None in row:
            raise ValueError(f'Invalid PT checkpoint ledger row: {prefix}')
        phase = row['phase']
        block, first, last = (number(row, k, True) for k in ('block', 'first_sweep', 'last_sweep'))
        prior = phases.get(phase)
        if block != (prior['last_block'] + 1 if prior else 0) or first != (prior['last_sweep'] + 1 if prior else 1) or last < first:
            raise ValueError(f"Noncontiguous PT checkpoint ledger: {prefix}")
        phases[phase] = dict(last_block=block, last_sweep=last,
                            interval=prior['interval'] if prior else last, completed=False)
    for phase in ('prerun', 'warmup', 'production'):
        path = Path(str(prefix) + f'pt_progress_{phase}.csv')
        if not path.exists():
            continue
        progress = list(rows(path))
        if len(progress) != 1:
            raise ValueError(f"Invalid PT phase progress: {path}")
        row = progress[0]
        boundary = phases.setdefault(phase, dict(last_block=-1, last_sweep=0, interval=1, completed=False))
        boundary['interval'] = number(row, 'checkpoint_sweeps', True)
        boundary['completed'] = row['completed'] == '1' and number(row, 'last_sweep', True) == boundary['last_sweep']
        if boundary['interval'] <= 0:
            raise ValueError('Invalid PT checkpoint interval')
    return phases


def analyze_streams(prefixes, output, no_plots=False, trace_points=2000):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if trace_points < 2:
        raise ValueError('trace_points must be at least 2')
    for name in list(FILES.values()) + ['pt_local_acceptance.csv', 'pt_occupancy.csv', 'pt_region_crossings.csv', 'pt_region_returns.csv', 'pt_swap_acceptance.csv', 'pt_local_acceptance.png', 'pt_swap_history.png']:
        (output/name).unlink(missing_ok=True)
    for pattern in ('pt_regions_run*.png', 'pt_paths_run*.png'):
        for path in output.glob(pattern): path.unlink()
    writers, streams = {}, []
    local, occupancy, transitions, excursions = {}, {}, defaultdict(int), defaultdict(int)
    block_swaps, trace, swaps = {}, defaultdict(list), {}
    schema_by_run, reports, issues = {}, [], []
    try:
        for run, prefix in enumerate(prefixes):
            if not Path(str(prefix) + 'pt_checkpoints.csv').exists():
                continue  # Backward compatibility: old cold-only diagnostics still work.
            phases = checkpoints(prefix)
            preferred = 'production' if 'production' in phases else ('warmup' if 'warmup' in phases else 'prerun')
            defs = list(rows(str(prefix) + 'pt_partitions.csv'))
            schema = defaultdict(dict)
            for row in defs:
                schema[row['partition']][row['region']] = number(row, 'region_index', True)
            schema_by_run[run] = dict(schema)
            run_info = dict(run=run, prefix=str(prefix), phases=phases, preferred_phase=preferred,
                            partitions=list(schema), ignored_uncommitted_rows=0, observations=[])
            reports.append(run_info)
            status_path = Path(str(prefix) + 'warmup_status.csv')
            if status_path.exists():
                status_rows = list(rows(status_path))
                if len(status_rows) != 1 or any(value is None for value in status_rows[0].values()):
                    raise ValueError(f'Invalid warmup status: {status_path}')
                status = dict(status_rows[0])
                for key in ('adaptation_sweeps', 'settling_sweeps', 'adaptation_efficiency_outliers',
                            'adaptation_zero_acceptance', 'settling_efficiency_outliers', 'settling_zero_acceptance',
                            'factorized_scale_updates', 'joint_scale_updates'):
                    if key in status:
                        status[key] = number(status, key, True)
                        if status[key] < 0: raise ValueError(f'Negative warmup status counter: {key}')
                for key in ('scale_target', 'scale_gain', 'scale_decay', 'scale_max_log_step'):
                    if key in status:
                        status[key] = number(status, key)
                        if not math.isfinite(status[key]): raise ValueError(f'Nonfinite warmup scale setting: {key}')
                for key in ('adaptation_complete', 'settling_complete', 'proposal_loaded',
                            'legacy_rvalue_checked', 'legacy_rvalue_passed', 'posterior_convergence_assessed'):
                    if key in status:
                        if status[key] not in ('0', '1'): raise ValueError(f'Invalid warmup status flag: {key}')
                        status[key] = status[key] == '1'
                run_info['warmup_status'] = status
            if not phases.get(preferred, {}).get('completed', False):
                issues.append(f"run {run}: {preferred} interrupted or still running; using committed sweeps through {phases.get(preferred, {}).get('last_sweep', 0)}")
            # Copy full committed block tables in a streaming pass. Only compact
            # cumulative counters and bounded plot samples remain in memory.
            for suffix, destination in FILES.items():
                path = Path(str(prefix) + suffix)
                if not path.exists():
                    raise ValueError(f"Missing streaming PT file: {path}")
                for row in rows(path):
                    phase = row.get('phase')
                    try:
                        block = int(row['block'])
                    except (KeyError, ValueError, TypeError):
                        run_info['ignored_uncommitted_rows'] += 1
                        continue
                    limit = phases.get(phase)
                    if not limit or block > limit['last_block']:
                        run_info['ignored_uncommitted_rows'] += 1
                        continue
                    if block < 0 or any(v is None for v in row.values()) or None in row:
                        raise ValueError(f"Incomplete committed PT row: {path}")
                    first, last = (number(row, k, True) for k in ('first_sweep', 'last_sweep'))
                    if first != block*limit['interval']+1 or last != min((block+1)*limit['interval'], limit['last_sweep']):
                        raise ValueError(f"PT row does not match committed window: {path}")
                    row = dict(run=run, **row)
                    if destination not in writers:
                        stream = (output / destination).open('w', newline='')
                        streams.append(stream)
                        writer = csv.DictWriter(stream, fieldnames=list(row))
                        writer.writeheader()
                        writers[destination] = writer
                    writers[destination].writerow(row)
                    if suffix == 'pt_local_moves.csv':
                        key = (run, phase, row['ladder'], row['slot'], row['temperature'], row['kernel'])
                        total = local.setdefault(key, dict(attempts=0, accepted=0, invalid_proposals=0, nonfinite_evaluations=0))
                        for k in total:
                            total[k] += number(row, k, True)
                    elif suffix == 'pt_region_occupancy.csv':
                        if row['partition'] not in schema or row['region'] not in schema[row['partition']]:
                            raise ValueError('PT occupancy references an unknown region')
                        if number(row, 'draws', True) != last-first+1:
                            raise ValueError('PT occupancy draw count differs from its window')
                        key = (run, phase, row['ladder'], row['slot'], row['temperature'], row['partition'], row['region'])
                        total = occupancy.setdefault(key, [0, 0])
                        total[0] += number(row, 'samples', True)
                        total[1] += number(row, 'draws', True)
                    elif suffix == 'pt_region_transitions.csv':
                        key = (run, phase, row['ladder'], row['slot'], row['temperature'], row['partition'], row['from_region'], row['to_region'], row['cause'])
                        if row['cause'] not in ('local', 'swap') or row['from_region'] == row['to_region']:
                            raise ValueError('Invalid PT transition cause or unchanged label')
                        transitions[key] += number(row, 'transitions', True)
                    elif suffix == 'pt_roundtrip_regions.csv':
                        key = (run, phase, row['ladder'], row['replica'], row['partition'], row['departure_region'], row['return_region'])
                        excursions[key] += number(row, 'round_trips', True)
                    else:
                        attempts, accepted = (number(row, k, True) for k in ('attempts', 'accepted'))
                        if not 0 <= accepted <= attempts:
                            raise ValueError('Invalid PT swap counters')
                        edge = (run, phase, row['ladder'], row['temperature_low'], row['temperature_high'])
                        total_edge = swaps.setdefault(edge, [0, 0])
                        total_edge[0] += attempts
                        total_edge[1] += accepted
                        stride = max(1, (limit['last_block']+1)//trace_points)
                        if block % stride == 0 or block == limit['last_block']:
                            key = (run, phase, block, last, row['temperature_low'], row['temperature_high'])
                            total = block_swaps.setdefault(key, [0, 0])
                            total[0] += attempts
                            total[1] += accepted
                        if attempts and not accepted and phase == preferred:
                            # A zero in one short block is not itself a diagnosis.
                            pass
            state_path = Path(str(prefix) + 'pt_states.csv')
            if not state_path.exists():
                raise ValueError(f'Missing PT state trace: {state_path}')
            stride = max(1, phases.get(preferred, {}).get('last_sweep', 0)//trace_points)
            for row in rows(state_path):
                if row.get('phase') != preferred or any(v is None for v in row.values()):
                    continue
                sweep = number(row, 'sweep', True)
                if sweep > phases.get(preferred, {}).get('last_sweep', 0):
                    continue
                key = (run, row['ladder'], row['replica'])
                history = trace[key]
                # A bounded deterministic subsample, retaining the last point.
                if not history or sweep-history[-1][0] >= stride:
                    history.append((sweep, number(row, 'temperature'), number(row, 'log_likelihood')))
            # Missing region rows must never quietly become a zero occupancy.
            groups = defaultdict(dict)
            for key, (count, draws) in occupancy.items():
                if key[0] != run:
                    continue
                _, phase, ladder, slot, temp, part, region = key
                if draws != phases[phase]['last_sweep'] or not 0 <= count <= draws:
                    raise ValueError('Missing/invalid committed PT occupancy counts')
                groups[(phase, ladder, slot, temp, part)][region] = count
            for key in local:
                if key[0] == run and phases[key[1]]['last_sweep']:
                    for part in schema:
                        if (key[1], key[2], key[3], key[4], part) not in groups:
                            raise ValueError('Missing committed PT occupancy group')
            for key, counts in groups.items():
                if set(counts) != set(schema[key[-1]]) or sum(counts.values()) != phases[key[0]]['last_sweep']:
                    raise ValueError('PT region occupancies do not partition all committed sweeps')
            for part, labels in schema.items():
                observed = [(key, value) for key, value in occupancy.items()
                            if key[0] == run and key[1] == preferred and key[5] == part]
                if not observed:
                    continue
                hottest = max(float(key[4]) for key, _ in observed)
                local_hot = sum(count for key, count in transitions.items() if key[0] == run and key[1] == preferred and
                                float(key[4]) == hottest and key[5] == part and key[-1] == 'local')
                cold_changes = sum(count for key, count in transitions.items() if key[0] == run and key[1] == preferred and
                                   float(key[4]) == 1 and key[5] == part)
                run_info['observations'].append(dict(partition=part, hottest_temperature=hottest,
                                                      hottest_local_transitions=local_hot, cold_transitions=cold_changes))
                def disagreement(temp):
                    for label in labels:
                        weights = [count/draws for key, (count, draws) in observed if float(key[4]) == temp and key[6] == label and draws]
                        if weights and max(weights)-min(weights) > .5:
                            return True
                    return False
                if not local_hot and disagreement(hottest):
                    issues.append(f'run {run} {part}: hottest ladders disagree and no hottest local region crossings were observed')
                elif local_hot and not cold_changes and disagreement(1):
                    issues.append(f'run {run} {part}: hot local crossings occur, but cold ladders remain separated without crossings')
        # Validate and export cumulative counters after reading all runs.
        move_rows = []
        for key, values in local.items():
            run, phase, ladder, slot, temp, kernel = key
            if min(values.values()) < 0 or values['accepted'] > values['attempts'] or values['invalid_proposals']+values['nonfinite_evaluations'] > values['attempts']-values['accepted']:
                raise ValueError('Inconsistent PT proposal counters')
            move_rows.append(dict(run=run, phase=phase, ladder=ladder, slot=slot, temperature=temp, kernel=kernel,
                                  **values, acceptance_rate=values['accepted']/values['attempts'] if values['attempts'] else math.nan))
            if values['attempts'] and not values['accepted'] and phase == next(r['preferred_phase'] for r in reports if r['run'] == run):
                issues.append(f'run {run} ladder {ladder} T={temp}: no accepted {kernel} proposals')
        swap_rows = []
        for key, (attempts, accepted) in swaps.items():
            run, phase, ladder, low, high = key
            swap_rows.append(dict(run=run, phase=phase, ladder=ladder, temperature_low=low, temperature_high=high,
                                  attempts=attempts, accepted=accepted, acceptance_rate=accepted/attempts if attempts else math.nan))
            preferred = next(r['preferred_phase'] for r in reports if r['run'] == run)
            if phase == preferred and not accepted:
                issues.append(f'run {run} ladder {ladder}: no accepted swaps for T={low}–{high} ({attempts} attempts)')
        write_table(output/'pt_swap_acceptance.csv', swap_rows)
        write_table(output/'pt_local_acceptance.csv', move_rows)
        write_table(output/'pt_occupancy.csv', [dict(run=k[0], phase=k[1], ladder=k[2], slot=k[3], temperature=k[4],
                    partition=k[5], region=k[6], samples=v[0], draws=v[1], weight=v[0]/v[1] if v[1] else math.nan) for k,v in occupancy.items()])
        write_table(output/'pt_region_crossings.csv', [dict(run=k[0], phase=k[1], ladder=k[2], slot=k[3], temperature=k[4],
                    partition=k[5], from_region=k[6], to_region=k[7], cause=k[8], transitions=v) for k,v in transitions.items()])
        write_table(output/'pt_region_returns.csv', [dict(run=k[0], phase=k[1], ladder=k[2], replica=k[3], partition=k[4],
                    departure_region=k[5], return_region=k[6], round_trips=v, changed_region=k[5]!=k[6]) for k,v in excursions.items()])
        report = dict(available=bool(reports), runs=reports, issues=list(dict.fromkeys(issues)),
                      note='Occupancy and transitions use exact committed counters. State/temperature plots are decimated. Hot samples are not posterior samples; changed round-trip endpoints do not identify where the crossing happened.')
        if not no_plots and reports:
            plots(output, reports, move_rows, occupancy, transitions, excursions, block_swaps, trace, schema_by_run)
        return report
    finally:
        for stream in streams:
            stream.close()


def plots(output, reports, move_rows, occupancy, transitions, excursions, block_swaps, trace, schemas):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # Acceptance curves distinguish factorized coordinate proposals and joint proposals.
    grouped = defaultdict(lambda: [0, 0])
    for row in move_rows:
        key = (row['run'], row['phase'], row['kernel'], float(row['temperature']))
        grouped[key][0] += row['attempts']; grouped[key][1] += row['accepted']
    fig, ax = plt.subplots(figsize=(9, 4))
    for series in sorted({key[:3] for key in grouped}):
        points = sorted((key[3], v[1]/v[0]) for key,v in grouped.items() if key[:3] == series and v[0])
        if points:
            ax.plot(*zip(*points), marker='o', label=' / '.join(map(str,series)))
    ax.set(xscale='log', xlabel='Temperature', ylabel='Local proposal acceptance', ylim=(0, 1))
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(output/'pt_local_acceptance.png'); plt.close(fig)
    # Independent ladder curves are kept separate for regional mixing.
    for info in reports:
        run, phase = info['run'], info['preferred_phase']
        parts = info['partitions'][:8]
        if parts:
            fig, axes = plt.subplots(len(parts), 3, figsize=(14, 2.7*len(parts)), squeeze=False)
            for p, part in enumerate(parts):
                label0 = min(schemas[run][part], key=schemas[run][part].get)
                for ladder in sorted({k[2] for k in occupancy if k[0]==run and k[1]==phase and k[5]==part}, key=int):
                    points = sorted((float(k[4]), v[0]/v[1]) for k,v in occupancy.items() if k[0]==run and k[1]==phase and k[2]==ladder and k[5]==part and k[6]==label0 and v[1])
                    if points: axes[p,0].plot(*zip(*points), marker='o', label=f'ladder {ladder}')
                    temps = sorted({float(k[4]) for k in occupancy if k[0]==run and k[1]==phase and k[2]==ladder})
                    for column, cause in ((1,'local'),(2,'swap')):
                        counts = [sum(v for k,v in transitions.items() if k[0]==run and k[1]==phase and k[2]==ladder and float(k[4])==temp and k[5]==part and k[-1]==cause) for temp in temps]
                        axes[p,column].plot(temps, counts, marker='o')
                axes[p,0].set(ylabel=part, ylim=(-.03,1.03), title=f'Occupancy: {label0}')
                axes[p,1].set(title='Local region transitions', yscale='symlog')
                axes[p,2].set(title='Swap-induced region transitions', yscale='symlog')
                for ax in axes[p]: ax.set(xscale='log', xlabel='Temperature')
            axes[0,0].legend(fontsize=6); fig.suptitle(f'Run {run}: {phase} (exact committed counts)')
            fig.tight_layout(); fig.savefig(output/f'pt_regions_run{run}.png'); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4))
    for series in sorted({(k[0],k[1],k[4],k[5]) for k in block_swaps}):
        points = sorted((k[3], v[1]/v[0]) for k,v in block_swaps.items() if (k[0],k[1],k[4],k[5])==series and v[0])
        if points: ax.plot(*zip(*points), label=f'{series[0]} {series[1]} T={series[2]}–{series[3]}')
    ax.set(xlabel='Sweep within phase', ylabel='Swap acceptance in sampled blocks', ylim=(0,1))
    if block_swaps: ax.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(output/'pt_swap_history.png'); plt.close(fig)
    for info in reports:
        run=info['run']
        keys=sorted([k for k in trace if k[0]==run],key=lambda k:(int(k[1]),int(k[2])))
        if not keys: continue
        ladder=keys[0][1]
        fig, axes=plt.subplots(2,1,figsize=(10,7))
        for key in [k for k in keys if k[1]==ladder][:8]:
            history=trace[key]
            axes[0].plot([p[0] for p in history],[p[1] for p in history],lw=.6,label=f'replica {key[2]}')
            finite=[p for p in history if math.isfinite(p[2])]
            axes[1].plot([p[0] for p in finite],[p[2] for p in finite],lw=.6)
        axes[0].set(yscale='log',ylabel='Temperature',title=f'Run {run}, ladder {ladder}: decimated replica trajectories')
        axes[0].legend(fontsize=7)
        axes[1].set(xlabel='Sweep',ylabel='Untempered log likelihood')
        fig.tight_layout();fig.savefig(output/f'pt_paths_run{run}.png');plt.close(fig)
