"""
Unit tests for the I-spline isotonic regression implementation in IsotonicPEP.

Inspired by:
  UnitTest_Percolator_IsplineRegression.cpp
  https://github.com/statisticalbiotechnology/percolator/blob/master/tests/unit_tests/
  percolator/UnitTest_Percolator_IsplineRegression.cpp
"""

import unittest
import numpy as np

from pyIsoPEP.IsotonicPEP import IsotonicRegression, IsotonicPEP


class TestIsplineDesignMatrix(unittest.TestCase):
    """Tests for the _build_ispline_basis helper (mirrors IsplineDesignMatrixMonotonicity)."""

    def setUp(self):
        self.model = IsotonicRegression()

    def test_columns_are_monotone_non_decreasing(self):
        """Each I-spline basis column must be non-decreasing across evaluation points."""
        x = np.linspace(0.0, 1.0, 101)
        interior = np.array([0.25, 0.5, 0.75])
        order = 4
        t = np.r_[[x[0]] * order, interior, [x[-1]] * order]

        B = self.model._build_ispline_basis(x, t, order)

        for col in range(B.shape[1]):
            diffs = np.diff(B[:, col])
            self.assertTrue(
                np.all(diffs >= -1e-12),
                msg=f"I-spline column {col} is not monotone non-decreasing"
            )

    def test_columns_start_near_zero_end_near_one(self):
        """Each I-spline basis column should start near 0 and end near 1."""
        x = np.linspace(0.0, 1.0, 101)
        interior = np.array([0.25, 0.5, 0.75])
        order = 4
        t = np.r_[[x[0]] * order, interior, [x[-1]] * order]

        B = self.model._build_ispline_basis(x, t, order)

        for col in range(B.shape[1]):
            self.assertAlmostEqual(B[0, col],  0.0, delta=1e-10,
                                   msg=f"Column {col} does not start near 0")
            self.assertAlmostEqual(B[-1, col], 1.0, delta=1e-10,
                                   msg=f"Column {col} does not end near 1")

    def test_values_bounded_in_unit_interval(self):
        """All basis values must lie in [0, 1]."""
        x = np.linspace(0.0, 1.0, 200)
        interior = np.array([0.2, 0.4, 0.6, 0.8])
        order = 4
        t = np.r_[[x[0]] * order, interior, [x[-1]] * order]

        B = self.model._build_ispline_basis(x, t, order)

        self.assertTrue(np.all(B >= -1e-12), "Basis values below 0")
        self.assertTrue(np.all(B <=  1 + 1e-12), "Basis values above 1")


class TestMonotonicityFitY(unittest.TestCase):
    """Tests that ispline_non_decreasing produces a non-decreasing sequence."""

    def setUp(self):
        self.model = IsotonicRegression()

    def test_non_decreasing_output_on_noisy_input(self):
        """Mirrors MonotonicityFitY: unordered input must yield monotone output."""
        y = [0.1, 0.05, 0.2, 0.15, 0.3, 0.25, 0.5]
        fitted = self.model.ispline_non_decreasing(y)

        self.assertEqual(len(fitted), len(y))
        for i in range(1, len(fitted)):
            self.assertGreaterEqual(
                fitted[i], fitted[i - 1] - 1e-12,
                msg=f"Monotonicity violated at position {i}"
            )

    def test_output_length_matches_input(self):
        y = [0.3, 0.1, 0.4, 0.1, 0.5, 0.9, 0.2, 0.6]
        fitted = self.model.ispline_non_decreasing(y)
        self.assertEqual(len(fitted), len(y))

    def test_empty_input_returns_empty(self):
        self.assertEqual(self.model.ispline_non_decreasing([]), [])


class TestFitBinaryInputsSmall(unittest.TestCase):
    """Mirrors FitXYBinaryInputs: small binary dataset."""

    def setUp(self):
        self.model = IsotonicRegression()

    def test_monotone_and_bounded(self):
        y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
        fitted = self.model.ispline_non_decreasing(y)

        self.assertEqual(len(fitted), len(y))

        for i in range(1, len(fitted)):
            self.assertGreaterEqual(
                fitted[i], fitted[i - 1] - 1e-12,
                msg=f"Monotonicity violated at position {i}"
            )

        self.assertGreaterEqual(fitted[0],  0.0)
        self.assertLessEqual(fitted[-1], 1.0)


class TestFitBinaryInputsLarge(unittest.TestCase):
    """Mirrors FitXYBinaryInputsLarge: 1000-point dataset with three regions."""

    def setUp(self):
        self.model = IsotonicRegression()
        N = 1000
        y = []
        for i in range(N):
            if i < N // 3:
                y.append(0.0)
            elif i > 2 * N // 3:
                y.append(1.0)
            else:
                y.append(float(i % 2))   # alternating 0/1 in the middle
        self.y = y

    def test_monotone(self):
        fitted = self.model.ispline_non_decreasing(self.y)

        self.assertEqual(len(fitted), len(self.y))
        for i in range(1, len(fitted)):
            self.assertGreaterEqual(
                fitted[i], fitted[i - 1] - 1e-12,
                msg=f"Monotonicity violated at position {i}"
            )

    def test_values_in_unit_interval(self):
        fitted = self.model.ispline_non_decreasing(self.y)

        for i, v in enumerate(fitted):
            self.assertGreaterEqual(v, 0.0, msg=f"Negative value at position {i}")
            self.assertLessEqual(v,   1.0, msg=f"Value > 1 at position {i}")


class TestConstantInput(unittest.TestCase):
    """Mirrors HandlesConstantInput: all-identical input should yield near-constant output."""

    def setUp(self):
        self.model = IsotonicRegression()

    def test_constant_half(self):
        y = [0.5] * 10_000
        fitted = self.model.ispline_non_decreasing(y)

        self.assertEqual(len(fitted), len(y))
        for i, v in enumerate(fitted):
            self.assertAlmostEqual(
                v, 0.5, delta=3e-2,
                msg=f"Fitted value {v} too far from 0.5 at position {i}"
            )

    def test_constant_zero(self):
        y = [0.0] * 500
        fitted = self.model.ispline_non_decreasing(y)
        for v in fitted:
            self.assertAlmostEqual(v, 0.0, delta=1e-6)

    def test_constant_one(self):
        y = [1.0] * 500
        fitted = self.model.ispline_non_decreasing(y)
        for v in fitted:
            self.assertAlmostEqual(v, 1.0, delta=1e-4)


class TestSmallDataset(unittest.TestCase):
    """Mirrors SmallDatasetUniformKnots: few points should still yield monotone output."""

    def setUp(self):
        self.model = IsotonicRegression()

    def test_five_point_monotone(self):
        y = [0.1, 0.2, 0.2, 0.3, 0.4]
        fitted = self.model.ispline_non_decreasing(y)

        self.assertEqual(len(fitted), len(y))
        for i in range(1, len(fitted)):
            self.assertGreaterEqual(
                fitted[i], fitted[i - 1] - 1e-12,
                msg=f"Monotonicity violated at position {i}"
            )

    def test_two_point_input(self):
        """Two-point input is the smallest meaningful monotone problem."""
        fitted = self.model.ispline_non_decreasing([0.2, 0.8])
        self.assertEqual(len(fitted), 2)
        self.assertGreaterEqual(fitted[1], fitted[0] - 1e-12)


class TestLargeDatasetAdaptiveKnots(unittest.TestCase):
    """Mirrors LargeDatasetAdaptiveKnots: 500-point step function."""

    def setUp(self):
        self.model = IsotonicRegression()
        self.y = [0.0] * 300 + [1.0] * 200

    def test_monotone(self):
        fitted = self.model.ispline_non_decreasing(self.y)

        self.assertEqual(len(fitted), len(self.y))
        for i in range(1, len(fitted)):
            self.assertGreaterEqual(
                fitted[i], fitted[i - 1] - 1e-12,
                msg=f"Monotonicity violated at position {i}"
            )

    def test_step_shape(self):
        """Early fitted values should be less than late fitted values."""
        fitted = self.model.ispline_non_decreasing(self.y)
        mean_early = np.mean(fitted[:100])
        mean_late  = np.mean(fitted[400:])
        self.assertLess(mean_early, mean_late)


if __name__ == "__main__":
    unittest.main()
