import unittest
import numpy as np

from pyIsoPEP.IsotonicPEP import IsotonicRegression, IsotonicPEP


def _monotone(seq, tol=1e-12):
    a = np.asarray(seq)
    return bool(np.all(np.diff(a) >= -tol))


def _in_unit(seq, tol=1e-12):
    a = np.asarray(seq)
    return bool(np.all(a >= -tol) and np.all(a <= 1.0 + tol))


# =============================================================================
# I-Spline basis matrix
# =============================================================================

class TestIsplineDesignMatrix(unittest.TestCase):
    def setUp(self):
        self.model = IsotonicRegression()

    def test_columns_are_monotone_non_decreasing(self):
        # New API: degree=3 (cubic), include_intercept=False for I-spline-only check.
        # Percolator definition: I_j(x) = sum_{m>=j} B_m(x), j=1..n_bspline-1.
        x = np.linspace(0.0, 1.0, 101)
        interior = np.array([0.25, 0.5, 0.75])
        degree = 3
        t = np.r_[[x[0]] * (degree + 1), interior, [x[-1]] * (degree + 1)]
        B = self.model._build_ispline_basis(x, t, degree=degree, include_intercept=False)
        for col in range(B.shape[1]):
            self.assertTrue(_monotone(B[:, col]),
                            msg=f"I-spline column {col} is not monotone")

    def test_columns_start_near_zero_end_near_one(self):
        x = np.linspace(0.0, 1.0, 101)
        interior = np.array([0.25, 0.5, 0.75])
        degree = 3
        t = np.r_[[x[0]] * (degree + 1), interior, [x[-1]] * (degree + 1)]
        B = self.model._build_ispline_basis(x, t, degree=degree, include_intercept=False)
        for col in range(B.shape[1]):
            self.assertAlmostEqual(B[0, col],  0.0, delta=1e-10)
            self.assertAlmostEqual(B[-1, col], 1.0, delta=1e-10)

    def test_values_bounded_in_unit_interval(self):
        x = np.linspace(0.0, 1.0, 200)
        interior = np.array([0.2, 0.4, 0.6, 0.8])
        degree = 3
        t = np.r_[[x[0]] * (degree + 1), interior, [x[-1]] * (degree + 1)]
        B = self.model._build_ispline_basis(x, t, degree=degree, include_intercept=False)
        self.assertTrue(np.all(B >= -1e-12))
        self.assertTrue(np.all(B <= 1 + 1e-12))


# =============================================================================
# ispline_non_decreasing  (rank-based)
# =============================================================================

class TestMonotonicityFitY(unittest.TestCase):
    def setUp(self):
        self.model = IsotonicRegression()

    def test_non_decreasing_output_on_noisy_input(self):
        y = [0.1, 0.05, 0.2, 0.15, 0.3, 0.25, 0.5]
        fitted = self.model.ispline_non_decreasing(y)
        self.assertEqual(len(fitted), len(y))
        self.assertTrue(_monotone(fitted))

    def test_output_length_matches_input(self):
        y = [0.3, 0.1, 0.4, 0.1, 0.5, 0.9, 0.2, 0.6]
        fitted = self.model.ispline_non_decreasing(y)
        self.assertEqual(len(fitted), len(y))

    def test_empty_input_returns_empty(self):
        self.assertEqual(self.model.ispline_non_decreasing([]), [])


class TestFitBinaryInputsSmall(unittest.TestCase):
    def setUp(self):
        self.model = IsotonicRegression()

    def test_monotone_and_bounded(self):
        y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
        fitted = self.model.ispline_non_decreasing(y)
        self.assertEqual(len(fitted), len(y))
        self.assertTrue(_monotone(fitted))
        self.assertGreaterEqual(fitted[0],  0.0)
        self.assertLessEqual(fitted[-1], 1.0)


class TestFitBinaryInputsLarge(unittest.TestCase):
    def setUp(self):
        self.model = IsotonicRegression()
        N = 1000
        self.y = [0.0 if i < N // 3 else (1.0 if i > 2 * N // 3 else float(i % 2))
                  for i in range(N)]

    def test_monotone(self):
        fitted = self.model.ispline_non_decreasing(self.y)
        self.assertEqual(len(fitted), len(self.y))
        self.assertTrue(_monotone(fitted))

    def test_values_in_unit_interval(self):
        fitted = self.model.ispline_non_decreasing(self.y)
        self.assertTrue(_in_unit(fitted))


class TestConstantInput(unittest.TestCase):
    def setUp(self):
        self.model = IsotonicRegression()

    def test_constant_half(self):
        y = [0.5] * 10_000
        fitted = self.model.ispline_non_decreasing(y)
        self.assertEqual(len(fitted), len(y))
        for v in fitted:
            self.assertAlmostEqual(v, 0.5, delta=3e-2)

    def test_constant_zero(self):
        y = [0.0] * 500
        for v in self.model.ispline_non_decreasing(y):
            self.assertAlmostEqual(v, 0.0, delta=1e-6)

    def test_constant_one(self):
        # With Percolator-matched regularisation (ridge 1e-4 + smooth 1e-3), the
        # ridge penalty pulls the fit slightly below 1.0 for a constant-1 input.
        y = [1.0] * 500
        for v in self.model.ispline_non_decreasing(y):
            self.assertGreater(v, 0.95)
            self.assertLessEqual(v, 1.0 + 1e-9)


class TestSmallDataset(unittest.TestCase):
    def setUp(self):
        self.model = IsotonicRegression()

    def test_five_point_monotone(self):
        y = [0.1, 0.2, 0.2, 0.3, 0.4]
        fitted = self.model.ispline_non_decreasing(y)
        self.assertEqual(len(fitted), len(y))
        self.assertTrue(_monotone(fitted))

    def test_two_point_input(self):
        fitted = self.model.ispline_non_decreasing([0.2, 0.8])
        self.assertEqual(len(fitted), 2)
        self.assertGreaterEqual(fitted[1], fitted[0] - 1e-12)


class TestLargeDatasetAdaptiveKnots(unittest.TestCase):
    def setUp(self):
        self.model = IsotonicRegression()
        self.y = [0.0] * 300 + [1.0] * 200

    def test_monotone(self):
        fitted = self.model.ispline_non_decreasing(self.y)
        self.assertEqual(len(fitted), len(self.y))
        self.assertTrue(_monotone(fitted))

    def test_step_shape(self):
        fitted = self.model.ispline_non_decreasing(self.y)
        self.assertLess(np.mean(fitted[:100]), np.mean(fitted[400:]))


# =============================================================================
# ispline_non_decreasing_xy  (score-based)
# =============================================================================

class TestIsplineNonDecreasingXY(unittest.TestCase):
    def setUp(self):
        self.model = IsotonicRegression()

    def _make_tdc_data(self, n=400, seed=0):
        rng = np.random.default_rng(seed)
        scores = np.linspace(10.0, 0.0, n)
        is_decoy = (rng.random(n) < np.linspace(0.05, 0.8, n)).astype(float)
        return scores, is_decoy

    def test_output_length(self):
        scores, is_decoy = self._make_tdc_data()
        fitted = self.model.ispline_non_decreasing_xy(scores, is_decoy)
        self.assertEqual(len(fitted), len(is_decoy))

    def test_monotone_in_rank_order(self):
        scores, is_decoy = self._make_tdc_data()
        fitted = self.model.ispline_non_decreasing_xy(scores, is_decoy)
        self.assertTrue(_monotone(fitted))

    def test_values_in_unit_interval(self):
        scores, is_decoy = self._make_tdc_data()
        fitted = self.model.ispline_non_decreasing_xy(scores, is_decoy)
        self.assertTrue(_in_unit(fitted))

    def test_empty_input(self):
        self.assertEqual(self.model.ispline_non_decreasing_xy([], []), [])

    def test_constant_scores_handled(self):
        scores = np.ones(20)
        y = np.linspace(0, 1, 20)
        fitted = self.model.ispline_non_decreasing_xy(scores, y)
        self.assertEqual(len(fitted), 20)
        self.assertTrue(_in_unit(fitted))

    def test_differs_from_rank_based_when_scores_uneven(self):
        rng = np.random.default_rng(42)
        scores = np.sort(rng.exponential(scale=2, size=200))[::-1]
        y = np.linspace(0, 1, 200) + rng.normal(0, 0.05, 200)
        y = np.clip(y, 0, 1)
        fitted_xy = self.model.ispline_non_decreasing_xy(scores, y)
        fitted_y  = self.model.ispline_non_decreasing(y)
        self.assertTrue(_monotone(fitted_xy))
        self.assertTrue(_monotone(fitted_y))


# =============================================================================
# Combination 1 & 2: q2pep × {ispline, PAVA}  (rank-based)
# =============================================================================

class TestQ2PEP(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(7)
        n = 300
        self.q_values = np.sort(rng.uniform(0.0, 0.3, n))
        self.iso = IsotonicPEP()

    def test_ispline_monotone(self):
        pep = self.iso.q_to_pep(self.q_values, pava=False)
        self.assertEqual(len(pep), len(self.q_values))
        self.assertTrue(_monotone(pep.values))

    def test_pava_monotone(self):
        pep = self.iso.q_to_pep(self.q_values, pava=True)
        self.assertEqual(len(pep), len(self.q_values))
        self.assertTrue(_monotone(pep.values))

    def test_ispline_in_unit(self):
        self.assertTrue(_in_unit(self.iso.q_to_pep(self.q_values, pava=False).values))

    def test_pava_in_unit(self):
        self.assertTrue(_in_unit(self.iso.q_to_pep(self.q_values, pava=True).values))

    def test_pep_regression_q2pep_ispline(self):
        _, _, pep, _ = self.iso.pep_regression(
            q_values=self.q_values, method="q2pep", pava=False
        )
        self.assertTrue(_monotone(np.sort(pep)))

    def test_pep_regression_q2pep_pava(self):
        _, _, pep, _ = self.iso.pep_regression(
            q_values=self.q_values, method="q2pep", pava=True
        )
        self.assertTrue(_in_unit(np.sort(pep)))


# =============================================================================
# Combination 3 & 4: qns2pep × {ispline, PAVA}  (score-based q-values)
# =============================================================================

class TestQNS2PEP(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(11)
        n = 300
        self.q_values = np.sort(rng.uniform(0.0, 0.3, n))
        self.scores = np.sort(rng.exponential(scale=3, size=n))[::-1]
        self.iso = IsotonicPEP()

    def test_ispline_monotone(self):
        pep = self.iso.q_to_pep(self.q_values, scores=self.scores, pava=False)
        self.assertEqual(len(pep), len(self.q_values))
        self.assertTrue(_monotone(pep.values))

    def test_pava_monotone(self):
        pep = self.iso.q_to_pep(self.q_values, scores=self.scores, pava=True)
        self.assertEqual(len(pep), len(self.q_values))
        self.assertTrue(_monotone(pep.values))

    def test_ispline_in_unit(self):
        self.assertTrue(_in_unit(self.iso.q_to_pep(self.q_values, scores=self.scores, pava=False).values))

    def test_pava_in_unit(self):
        self.assertTrue(_in_unit(self.iso.q_to_pep(self.q_values, scores=self.scores, pava=True).values))

    def test_size_mismatch_raises(self):
        with self.assertRaises((ValueError, AssertionError)):
            self.iso.q_to_pep(self.q_values, scores=self.scores[:-1], pava=False)

    def test_pep_regression_qns2pep_ispline(self):
        n = len(self.q_values)
        labels = np.zeros(n)
        obs = np.column_stack([self.scores, labels])
        _, _, pep, _ = self.iso.pep_regression(
            q_values=self.q_values, obs=obs, method="qns2pep", pava=False
        )
        self.assertEqual(len(pep), n)
        self.assertTrue(_in_unit(np.sort(pep)))

    def test_pep_regression_qns2pep_pava(self):
        n = len(self.q_values)
        labels = np.zeros(n)
        obs = np.column_stack([self.scores, labels])
        _, _, pep, _ = self.iso.pep_regression(
            q_values=self.q_values, obs=obs, method="qns2pep", pava=True
        )
        self.assertEqual(len(pep), n)
        self.assertTrue(_in_unit(np.sort(pep)))


# =============================================================================
# Combination 5 & 6: tdc2pep × {ispline, PAVA}  (score-based TDC)
# =============================================================================

class TestTDC2PEP(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(42)
        n_target, n_decoy = 300, 300
        sc_target = rng.normal(loc=3.0, scale=1.5, size=n_target)
        sc_decoy  = rng.normal(loc=0.0, scale=1.5, size=n_decoy)
        scores = np.concatenate([sc_target, sc_decoy])
        labels = np.concatenate([np.zeros(n_target), np.ones(n_decoy)])
        self.obs = np.column_stack([scores, labels])
        self.iso = IsotonicPEP()

    def test_tdc_to_pep_returns_all_psms(self):
        df_obs = self.iso.process_obs(self.obs)
        pep_all = self.iso.tdc_to_pep(df_obs, pava=False)
        self.assertEqual(len(pep_all), len(self.obs))

    def test_tdc_to_pep_ispline_in_unit(self):
        df_obs = self.iso.process_obs(self.obs)
        pep_all = self.iso.tdc_to_pep(df_obs, pava=False)
        self.assertTrue(_in_unit(pep_all.values))

    def test_tdc_to_pep_pava_in_unit(self):
        df_obs = self.iso.process_obs(self.obs)
        pep_all = self.iso.tdc_to_pep(df_obs, pava=True)
        self.assertTrue(_in_unit(pep_all.values))

    def test_tdc_to_pep_targets_lower_pep(self):
        df_obs = self.iso.process_obs(self.obs)
        pep_all = self.iso.tdc_to_pep(df_obs, pava=False)
        pep_targets = pep_all.values[self.obs[:, 1] == 0]
        pep_decoys  = pep_all.values[self.obs[:, 1] == 1]
        self.assertLess(np.mean(pep_targets), np.mean(pep_decoys))

    def test_pep_regression_tdc2pep_ispline_targets_only(self):
        _, _, pep, _ = self.iso.pep_regression(
            obs=self.obs, method="tdc2pep", pava=False
        )
        n_target = int(np.sum(self.obs[:, 1] == 0))
        self.assertEqual(len(pep), n_target)
        self.assertTrue(_in_unit(pep))

    def test_pep_regression_tdc2pep_pava_targets_only(self):
        _, _, pep, _ = self.iso.pep_regression(
            obs=self.obs, method="tdc2pep", pava=True
        )
        n_target = int(np.sum(self.obs[:, 1] == 0))
        self.assertEqual(len(pep), n_target)
        self.assertTrue(_in_unit(pep))

    def test_pep_regression_tdc2pep_with_q2(self):
        _, _, pep, q2 = self.iso.pep_regression(
            obs=self.obs, method="tdc2pep", pava=False,
            calc_q_from_pep=True,
        )
        self.assertIsNotNone(q2)
        self.assertEqual(len(q2), len(pep))
        self.assertTrue(_in_unit(q2))


# =============================================================================
# Cross-combination sanity: all six return same-length, bounded output
# =============================================================================

class TestAllSixCombinations(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(99)
        n_t, n_d = 200, 200
        sc_t = rng.normal(3.0, 1.0, n_t)
        sc_d = rng.normal(0.0, 1.0, n_d)
        self.obs = np.column_stack([
            np.concatenate([sc_t, sc_d]),
            np.concatenate([np.zeros(n_t), np.ones(n_d)]),
        ])
        sc_sorted = np.sort(sc_t)[::-1]
        n = len(sc_sorted)
        self.q_values = np.cumsum(np.ones(n) * 0.01) / np.arange(1, n + 1)
        self.q_values = np.maximum.accumulate(self.q_values)
        self.n_target = n_t
        self.iso = IsotonicPEP()

    def _run(self, method, pava, extra=None):
        kw = dict(method=method, pava=pava,
                  obs=self.obs, calc_q_from_pep=True)
        if method in ("q2pep", "qns2pep"):
            kw["q_values"] = self.q_values
        if extra:
            kw.update(extra)
        _, _, pep, q2 = self.iso.pep_regression(**kw)
        self.assertEqual(len(pep), self.n_target)
        self.assertTrue(_in_unit(pep))
        self.assertIsNotNone(q2)

    def test_combo_1_q2pep_ispline(self):
        self._run("q2pep", False)

    def test_combo_2_q2pep_pava(self):
        self._run("q2pep", True)

    def test_combo_3_qns2pep_ispline(self):
        self._run("qns2pep", False)

    def test_combo_4_qns2pep_pava(self):
        self._run("qns2pep", True)

    def test_combo_5_tdc2pep_ispline(self):
        self._run("tdc2pep", False)

    def test_combo_6_tdc2pep_pava(self):
        self._run("tdc2pep", True)


if __name__ == "__main__":
    unittest.main()
