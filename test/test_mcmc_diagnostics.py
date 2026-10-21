"""Regression tests for scientific diagnostics, ROOT ordering and region checks."""
import csv
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
import uproot

SCRIPT = Path(__file__).resolve().parents[1] / "diagnostics" / "analyze_mcmc.py"
spec = importlib.util.spec_from_file_location("diagnostics", SCRIPT)
diag = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diag)


def write_root(path, x, warmup=0, bad_iterations=False, temperatures=False, fixed=False):
    m, n = x.shape
    arrays = {"Chain": np.tile(np.arange(m, dtype="uint32"), n),
              "Iteration": np.repeat(np.arange(1, n+1, dtype="uint32"), m),
              "Phase": np.ones(m*n, dtype="int32"), "x": x.T.reshape(-1)}
    if fixed:
        arrays["fixed"] = np.full(m*n, 0.4)
    if temperatures:
        arrays["Temperature"] = np.tile(np.arange(1, m+1, dtype=float), n)
    if bad_iterations:
        arrays["Iteration"][m] = 1
    if warmup:
        for k, value in list(arrays.items()):
            if k == "Chain": head = np.tile(np.arange(m, dtype="uint32"), warmup)
            elif k == "Iteration": head = np.repeat(np.arange(1, warmup+1, dtype="uint32"), m)
            elif k == "Phase": head = np.zeros(warmup*m, dtype="int32")
            else: head = np.full(warmup*m, 99., dtype=value.dtype)
            arrays[k] = np.concatenate([head, value])
    with uproot.recreate(path) as root:
        tree = root.mktree("model_mcmc", {k: v.dtype for k, v in arrays.items()})
        tree.extend(arrays)
        if fixed:
            meta = root.mktree("model_parameters", {"safe_name": "string", "parameter": "bool", "fixed": "bool"})
            meta.extend({"safe_name": ["x", "fixed"], "parameter": [True, True], "fixed": [False, True]})


class StatisticalDiagnostics(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(82903)
        self.x = self.rng.normal(size=(4, 4000))

    def test_iid_chains(self):
        result = diag.diagnostic(self.x)
        self.assertLess(result["rhat"], 1.01)
        self.assertGreater(result["ess_bulk"], self.x.size*.7)
        self.assertLess(result["mcse_mean"], .012)
        self.assertEqual(result["status"], "ok")

    def test_shifted_chains_flagged(self):
        self.x[0] += 4
        result = diag.diagnostic(self.x)
        self.assertGreater(result["rhat"], 1.2)
        self.assertIn("rhat_high", result["status"])

    def test_scale_mismatch_detected_by_folded_rhat(self):
        self.x[0] *= 8
        self.assertGreater(diag.diagnostic(self.x)["rhat"], 1.1)

    def test_autocorrelation_reduces_ess(self):
        for j in range(1, self.x.shape[1]):
            self.x[:, j] = .95*self.x[:, j-1] + self.x[:, j]
        result = diag.diagnostic(self.x)
        self.assertLess(result["ess_bulk"], self.x.size*.1)
        self.assertGreater(result["mcse_mean"], .05)

    def test_fixed_and_stuck_are_not_equivalent(self):
        constant = np.ones((4, 100))
        self.assertEqual(diag.diagnostic(constant, fixed=True)["status"], "fixed")
        result = diag.diagnostic(constant)
        self.assertEqual(result["status"], "constant_unverified")
        self.assertTrue(math.isnan(result["mcse_mean"]))

    def test_stuck_modes_do_not_report_zero_weight_error(self):
        mask = np.vstack([np.zeros((2, 100)), np.ones((2, 100))]).astype(bool)
        result = diag.region_diagnostic(mask)
        self.assertEqual(result["weight"], .5)
        self.assertEqual(result["transitions"], 0)
        self.assertTrue(math.isnan(result["mcse_weight"]))
        self.assertIn("without_region_transitions", result["status"])

    def test_unvisited_region_not_precisely_zero(self):
        result = diag.region_diagnostic(np.zeros((4, 100), dtype=bool))
        self.assertEqual(result["status"], "not_visited")
        self.assertTrue(math.isnan(result["mcse_weight"]))
        result = diag.region_diagnostic(np.ones((4, 100), dtype=bool))
        self.assertEqual(result["status"], "no_transitions")

    def test_single_and_short_chains_not_certified(self):
        self.assertIn("insufficient_chains", diag.diagnostic(self.x[:1])["status"])
        self.assertEqual(diag.diagnostic(self.x[:, :2])["status"], "insufficient_draws")

    def test_nonfinite_input(self):
        self.x[0, 0] = np.nan
        self.assertEqual(diag.diagnostic(self.x)["status"], "nonfinite")


class RootAndRegions(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.directory = Path(self.tmp.name)
        self.rng = np.random.default_rng(34031)
        self.x = self.rng.normal(size=(4, 4000))
        self.path = self.directory / "samples.root"

    def tearDown(self):
        self.tmp.cleanup()

    def config(self):
        return {"partitions": [{"name": "sign", "conditions": [{"parameter": "x", "op": ">=", "value": 0}]}]}

    def test_root_order_warmup_filter_and_metadata(self):
        write_root(self.path, self.x, warmup=10, fixed=True)
        data, chains, iterations, fixed = diag.load_cold([self.path])
        np.testing.assert_array_equal(data["x"], self.x)
        self.assertEqual(len(chains), 4)
        self.assertEqual(iterations[0, 0], 1)
        self.assertEqual(fixed, {"fixed"})

    def test_independent_files_preserve_chain_identity(self):
        write_root(self.path, self.x)
        other = self.directory / "other.root"
        write_root(other, self.x + .1)
        data, chains, _, _ = diag.load_cold([self.path, other])
        self.assertEqual(data["x"].shape, (8, 4000))
        self.assertEqual(chains[-1], "1:3")

    def test_bad_iterations_rejected(self):
        write_root(self.path, self.x, bad_iterations=True)
        with self.assertRaisesRegex(ValueError, "irregular"):
            diag.load_cold([self.path])

    def test_unequal_lengths_rejected(self):
        write_root(self.path, self.x)
        other = self.directory / "short.root"
        write_root(other, self.x[:, :100])
        with self.assertRaisesRegex(ValueError, "Unequal"):
            diag.load_cold([self.path, other])

    def test_hot_output_rejected(self):
        write_root(self.path, self.x, temperatures=True)
        with self.assertRaisesRegex(ValueError, "heated"):
            diag.load_cold([self.path])

    def test_partitions_and_boundaries(self):
        data = {"x": np.array([[.5, .75, 1.]])}
        config = {"partitions": [{"name": "cut", "conditions": [{"parameter": "x", "op": "<", "value": .75}]}]}
        regions, _ = diag.classify_regions(config, data)
        np.testing.assert_array_equal(regions[0][2], [[True, False, False]])
        np.testing.assert_array_equal(regions[1][2], [[False, True, True]])

    def test_overlapping_joint_modes_rejected(self):
        config = {"regions": [{"name": name, "conditions": [{"parameter": "x", "op": ">=", "value": 0}]} for name in ["a", "b"]]}
        with self.assertRaisesRegex(ValueError, "overlap"):
            diag.classify_regions(config, {"x": self.x})

    def test_missing_optional_parameters_and_safe_conditions(self):
        config = {"partitions": [{"name": "missing", "optional": True, "conditions": [{"parameter": "absent", "value": 0}]}]}
        regions, skipped = diag.classify_regions(config, {"x": self.x})
        self.assertFalse(regions)
        self.assertEqual(skipped, ["missing"])
        with self.assertRaises(ValueError):
            diag.condition_mask({"parameter": "x", "transform": "__import__", "value": 0}, {"x": self.x})

    def test_other_joint_mode_and_unvisited_declared_mode(self):
        config = {"regions": [{"name": "impossible", "conditions": [{"parameter": "x", "op": ">", "value": 100}]}]}
        regions, _ = diag.classify_regions(config, {"x": self.x})
        self.assertEqual([r[1] for r in regions], ["impossible", "other"])
        self.assertFalse(regions[0][2].any())
        self.assertTrue(regions[1][2].all())

    def test_complete_report_and_plots(self):
        write_root(self.path, self.x, warmup=10, fixed=True)
        config_path = self.directory / "regions.json"
        config_path.write_text(json.dumps(self.config()))
        out = self.directory / "report"
        args = diag.parser().parse_args([str(self.path), "--output", str(out), "--regions", str(config_path), "--angular", "x"])
        report = diag.analyze(args)
        self.assertEqual(report["status"], "screening_passed")
        with (out / "diagnostics.csv").open() as stream:
            summaries = {row["variable"]: row for row in csv.DictReader(stream)}
        self.assertEqual(summaries["fixed"]["status"], "fixed")
        self.assertIn("sin(x)", summaries)
        self.assertIn("cos(x)", summaries)
        self.assertEqual(int(summaries["x"]["draws_per_chain"]), 4000)
        with (out / "region_weights_cumulative.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        last = [row for row in rows if row["chain"] == "0:0" and row["region"] == "true"][-1]
        self.assertAlmostEqual(float(last["weight"]), float((self.x[0] >= 0).mean()))
        self.assertGreater((out / "traces.png").stat().st_size, 1000)
        self.assertGreater((out / "rank_plots.png").stat().st_size, 1000)
        self.assertGreater((out / "region_traces.png").stat().st_size, 1000)
        json.loads((out / "summary.json").read_text())

    def test_invalid_region_schema(self):
        for config in [[], {"regions": "x"}, {"partitions": [{"name": "x"}]}]:
            with self.assertRaises(ValueError):
                diag.classify_regions(config, {"x": self.x})

    def test_transport_and_legacy_missing_data(self):
        prefix = self.directory / "run_"
        (self.directory / "run_pt_swaps.csv").write_text("temperature_low,temperature_high,attempts,accepted,acceptance_rate\n1,4,100,20,0.2\n")
        trips, swaps, missing = diag.read_transport([prefix])
        self.assertEqual(len(swaps), 1)
        self.assertEqual(swaps[0]["ladder"], "aggregate")
        self.assertFalse(trips)
        self.assertEqual(len(missing), 2)

    def test_transport_run_identity(self):
        prefixes = [self.directory / "a_", self.directory / "b_"]
        for prefix in prefixes:
            Path(str(prefix)+"pt_roundtrips.csv").write_text("replica,ladder,round_trips,total_sweeps,mean_sweeps,min_sweeps,max_sweeps,journey_stage\n0,0,2,100,50,40,60,1\n")
            Path(str(prefix)+"pt_swaps_by_ladder.csv").write_text("ladder,temperature_low,temperature_high,attempts,accepted,acceptance_rate\n0,1,4,100,20,0.2\n")
        trips, swaps, missing = diag.read_transport(prefixes)
        self.assertEqual([r["run"] for r in trips], [0, 1])
        self.assertEqual([r["run"] for r in swaps], [0, 1])
        self.assertFalse(missing)

    def test_missing_regions_are_explicit(self):
        write_root(self.path, self.x)
        args = diag.parser().parse_args([str(self.path), "--output", str(self.directory / "report"), "--no-plots"])
        report = diag.analyze(args)
        self.assertEqual(report["status"], "incomplete")


class StreamingPTDiagnostics(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.prefix = self.directory/'run_'
        self.output = self.directory/'report'
        self.pt = diag._pt_module
        self.table('pt_partitions.csv', ['partition_index','partition','region_index','region'],
                   [[0,'sign',0,'negative'],[0,'sign',1,'positive']])
        common = ['phase','block','first_sweep','last_sweep']
        ledger, moves, occ, changes, swaps, trips = [], [], [], [], [], []
        for phase in ('prerun','warmup','production'):
            self.table(f'pt_progress_{phase}.csv', ['phase','last_sweep','completed','checkpoint_sweeps','trace_sweeps'], [[phase,3,1,2,0]])
            for block, first, last in ((0,1,2),(1,3,3)):
                base = [phase,block,first,last]; ledger.append(base)
                for ladder in range(2):
                    swaps.append(base+[ladder,1,4,2,1,.5])
                    for temp_index, temp in enumerate((1,4)):
                        slot = temp_index*2+ladder
                        moves.append(base+[ladder,slot,temp,'factorized',4,2,1,0])
                        for label in ('negative','positive'):
                            occ.append(base+[ladder,slot,temp,'sign',label,last-first+1 if (label=='negative')==(ladder==0) else 0,last-first+1,.5])
                        if temp_index:
                            changes.append(base+[ladder,slot,temp,'sign','negative','positive','local',3])
                        else:
                            changes.append(base+[ladder,slot,temp,'sign','negative','positive','swap',1])
                    trips.append(base+[ladder,ladder,'sign','negative','positive',1])
        self.table('pt_checkpoints.csv',common,ledger)
        self.table('pt_local_moves.csv', common+['ladder','slot','temperature','kernel','attempts','accepted','invalid_proposals','nonfinite_evaluations'], moves)
        self.table('pt_region_occupancy.csv', common+['ladder','slot','temperature','partition','region','samples','draws','weight'],occ)
        self.table('pt_region_transitions.csv',common+['ladder','slot','temperature','partition','from_region','to_region','cause','transitions'],changes)
        self.table('pt_swaps_by_block.csv',common+['ladder','temperature_low','temperature_high','attempts','accepted','acceptance_rate'],swaps)
        self.table('pt_roundtrip_regions.csv',common+['ladder','replica','partition','departure_region','return_region','round_trips'],trips)
        self.table('pt_states.csv',['phase','sweep','ladder','slot','temperature','replica','log_likelihood','log_prior','region_0'],[])

    def table(self, suffix, columns, data):
        with Path(str(self.prefix)+suffix).open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(columns);writer.writerows(data)

    def run_report(self, plots=False):
        return self.pt.analyze_streams([self.prefix],self.output,no_plots=not plots)

    def test_phases_exact_counts_and_event_attribution(self):
        result=self.run_report()
        self.assertTrue(result['available']);self.assertFalse(result['issues'])
        observation=result['runs'][0]['observations'][0]
        self.assertEqual(observation['hottest_local_transitions'],12)
        self.assertEqual(observation['cold_transitions'],4)
        self.assertEqual(len(list(self.pt.rows(self.output/'pt_local_acceptance.csv'))),12)
        returns=list(self.pt.rows(self.output/'pt_region_returns.csv'))
        self.assertEqual(len(returns),6)
        self.assertTrue(all(r['changed_region']=='True' and int(r['round_trips'])==2 for r in returns))
        weights=list(self.pt.rows(self.output/'pt_occupancy.csv'))
        self.assertTrue(all(int(r['draws'])==3 for r in weights))

    def test_interruption_ignores_uncommitted_tail(self):
        path=Path(str(self.prefix)+'pt_checkpoints.csv')
        path.write_text(path.read_text().replace('production,1,3,3\n',''))
        progress=Path(str(self.prefix)+'pt_progress_production.csv')
        progress.write_text('phase,last_sweep,completed,checkpoint_sweeps,trace_sweeps\nproduction,2,0,2,0\n')
        with Path(str(self.prefix)+'pt_local_moves.csv').open('a') as stream:stream.write('production,2,5')
        result=self.run_report()
        self.assertGreater(result['runs'][0]['ignored_uncommitted_rows'],0)
        self.assertTrue(any('interrupted' in issue for issue in result['issues']))
        self.assertEqual(result['runs'][0]['phases']['production']['last_sweep'],2)

    def test_missing_committed_region_rejected(self):
        path=Path(str(self.prefix)+'pt_region_occupancy.csv')
        lines=path.read_text().splitlines(True);path.write_text(''.join(lines[:-1]))
        with self.assertRaisesRegex(ValueError,'occupancy'):self.run_report()

    def test_missing_entire_committed_slot_rejected(self):
        path=Path(str(self.prefix)+'pt_region_occupancy.csv')
        lines=path.read_text().splitlines(True)
        path.write_text(''.join(line for line in lines if not line.startswith('production,') or ',0,0,1,sign,' not in line))
        with self.assertRaisesRegex(ValueError,'occupancy'):self.run_report()

    def test_corrupt_swap_count_rejected(self):
        path=Path(str(self.prefix)+'pt_swaps_by_block.csv')
        path.write_text(path.read_text().replace(',2,1,0.5',',2,3,0.5'))
        with self.assertRaisesRegex(ValueError,'swap'):self.run_report()

    def test_run_ids_remain_separate(self):
        result=self.pt.analyze_streams([self.prefix,self.prefix],self.output,no_plots=True)
        self.assertEqual([r['run'] for r in result['runs']],[0,1])
        self.assertEqual({r['run'] for r in self.pt.rows(self.output/'pt_region_returns.csv')},{'0','1'})

    def test_legacy_files_remain_supported(self):
        result=self.pt.analyze_streams([self.directory/'old_'],self.output,no_plots=True)
        self.assertFalse(result['available']);self.assertFalse(result['issues'])

    def test_pt_only_report_and_plots(self):
        args=diag.parser().parse_args(['--pt-only','--transport-prefix',str(self.prefix),'--output',str(self.output)])
        report=diag.analyze(args)
        self.assertEqual(report['status'],'transport_checked')
        self.assertIn('not checked',report['note'])
        self.assertGreater((self.output/'pt_regions_run0.png').stat().st_size,1000)
        self.run_report()
        self.assertFalse((self.output/'pt_regions_run0.png').exists())

    def test_hot_trapping_warning_requires_disagreement(self):
        path=Path(str(self.prefix)+'pt_region_transitions.csv')
        lines=path.read_text().splitlines(True)
        path.write_text(''.join(line for line in lines if ',local,' not in line))
        result=self.run_report()
        self.assertTrue(any('hottest ladders disagree' in x for x in result['issues']))


if __name__ == "__main__":
    unittest.main()
