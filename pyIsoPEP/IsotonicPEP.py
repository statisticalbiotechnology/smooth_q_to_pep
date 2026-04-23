import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy.interpolate import BSpline

DEFAULT_MAX_ITER = 1000
DEFAULT_LAMBDA = 1e-4
DEFAULT_SMOOTH_LAMBDA = 1e-3
PEP_CLIP_LO = 1e-10


class PreProcessing:
    def process_obs(self, obs):
        if not isinstance(obs, np.ndarray) or obs.ndim != 2 or obs.shape[1] != 2:
            raise ValueError("obs must be a numpy array with shape (n, 2).")
        scores = obs[:, 0].astype(float)
        labels = obs[:, 1].astype(int)
        if not set(np.unique(labels)).issubset({0, 1}):
            raise ValueError("Labels must be 0 (target) or 1 (decoy).")
        df = pd.DataFrame({"score": scores, "label": labels})
        df["orig_idx"] = np.arange(len(df))
        return df

    def calc_q_from_fdr(self, obs):
        df = self.process_obs(obs)
        df_sorted = df.sort_values(by="score", ascending=False, kind="mergesort").reset_index(drop=True)
        df_sorted["cumulative_target"] = (df_sorted["label"] == 0).cumsum()
        df_sorted["cumulative_decoy"] = (df_sorted["label"] == 1).cumsum()
        df_sorted["FDR"] = (df_sorted["cumulative_decoy"] + 0.5) / (df_sorted["cumulative_target"] + 0.5)
        q = df_sorted["FDR"].values.copy()
        for i in range(len(q) - 2, -1, -1):
            q[i] = min(q[i], q[i + 1])
        df_sorted["q-value"] = q
        df_result = df_sorted.sort_values("orig_idx").reset_index(drop=True)
        return df_result["FDR"], df_result["q-value"]


class IsotonicRegression:
    def __init__(self):
        pass

    def pava_non_decreasing(self, values, counts, min_value=0.0, max_value=1.0):
        if len(values) != len(counts):
            raise ValueError("values and counts must have the same length.")
        n = len(values)
        if n == 0:
            return []
        stack = []
        for i in range(n):
            stack.append({"sum": values[i] * counts[i], "count": counts[i], "avg": values[i]})
            while len(stack) > 1:
                top = stack[-1]
                sec_top = stack[-2]
                if sec_top["avg"] > top["avg"]:
                    merged_sum = sec_top["sum"] + top["sum"]
                    merged_count = sec_top["count"] + top["count"]
                    stack.pop()
                    stack.pop()
                    stack.append({"sum": merged_sum, "count": merged_count, "avg": merged_sum / merged_count})
                else:
                    break
        result = []
        for block in stack:
            clamped_avg = min(max(block["avg"], min_value), max_value)
            result.extend([clamped_avg] * block["count"])
        return result

    def _make_default_knots(self, x_norm, degree=3):
        n = len(x_norm)
        x_sorted = np.sort(x_norm)
        lo = float(x_sorted[0])
        hi = float(x_sorted[-1])
        if hi <= lo:
            hi = lo + 1e-6

        num_internal = min(200, int(np.sqrt(n)))
        internal = []
        prev = lo
        for k in range(1, num_internal + 1):
            q = k / (num_internal + 1)
            idx = min(int(q * (n - 1)), n - 1)
            val = float(x_sorted[idx])
            if val > lo + 1e-12 and val < hi - 1e-12 and val > prev + 1e-12:
                internal.append(val)
                prev = val

        order = degree + 1
        t = [lo] * order + internal + [hi] * order
        return np.array(t, dtype=float)

    def _build_ispline_basis(self, x, t, degree=3, include_intercept=True):
        n_knots = len(t)
        n_bspline = n_knots - degree - 1

        if n_bspline <= 1:
            return np.ones((len(x), 1))

        I = np.zeros((len(x), n_bspline))
        for j in range(n_bspline):
            span = t[j + degree + 1] - t[j]
            if span <= 0:
                continue
            c = np.zeros(n_bspline)
            c[j] = 1.0
            ibspl = BSpline(t, c, degree).antiderivative()
            I[:, j] = np.clip(
                (degree + 1) / span * (ibspl(x) - ibspl(t[j])), 0.0, 1.0
            )

        if include_intercept:
            return np.column_stack([np.ones(len(x)), I])
        return I

    def _build_ispline_basis_deboor(self, x, t, degree=3, include_intercept=True):
        k = degree + 1              # order of M-spline B-splines
        n_knots = len(t)
        n_orig = n_knots - k        # number of order-k B-splines

        # Build augmented knot vector for order k+1: add one extra boundary knot at each end
        t_aug = np.concatenate([[t[0]], t, [t[-1]]])

        k1 = k + 1
        n_high = len(t_aug) - k1    # number of order-k1 B-splines

        if n_orig <= 1 or n_high <= 0:
            return np.ones((len(x), 1))

        # Evaluate order-k1 B-splines on the augmented knot vector
        B1 = np.zeros((len(x), n_high))
        for j in range(n_high):
            c = np.zeros(n_high)
            c[j] = 1.0
            B1[:, j] = BSpline(t_aug, c, degree + 1)(x)

        # Weighted reverse cumulative sum
        I = np.zeros((len(x), n_orig))
        for i in range(n_orig):
            denom = t[i + k] - t[i]
            if denom <= 0:
                continue
            col = np.zeros(len(x))
            for j in range(i, n_high):
                numer = t_aug[j + k] - t_aug[j]
                if numer <= 0:
                    continue
                col += (numer / denom) * B1[:, j]
            I[:, i] = np.clip(col, 0.0, 1.0)

        if include_intercept:
            return np.column_stack([np.ones(len(x)), I])
        return I

    def _fit_ispline(self, x, y, ridge_lambda=DEFAULT_LAMBDA,
                     smooth_lambda=DEFAULT_SMOOTH_LAMBDA,
                     min_value=0.0, max_value=1.0, degree=3):
        n = len(x)
        t = self._make_default_knots(x, degree=degree)
        X = self._build_ispline_basis(x, t, degree=degree, include_intercept=True)
        p = X.shape[1]
        n_ispline = p - 1  # columns 1..p-1 are I-splines; column 0 is intercept

        inv_sqrt_n = 1.0 / np.sqrt(max(1, n))
        Xw = X * inv_sqrt_n
        yw = y * inv_sqrt_n

        X_aug = np.vstack([Xw, np.sqrt(ridge_lambda) * np.eye(p)])
        y_aug = np.concatenate([yw, np.zeros(p)])

        if smooth_lambda > 0.0 and n_ispline > 2:
            nd = n_ispline - 2
            D = np.zeros((nd, p))
            for di in range(nd):
                D[di, 1 + di] = 1.0
                D[di, 1 + di + 1] = -2.0
                D[di, 1 + di + 2] = 1.0
            X_aug = np.vstack([X_aug, np.sqrt(smooth_lambda) * D])
            y_aug = np.concatenate([y_aug, np.zeros(nd)])

        lower_bounds = np.concatenate(([-np.inf], np.zeros(n_ispline)))
        upper_bounds = np.full(p, np.inf)
        res = lsq_linear(X_aug, y_aug, bounds=(lower_bounds, upper_bounds))
        return np.clip(X @ res.x, min_value, max_value)

    def ispline_non_decreasing(self, raw_pep, min_value=0.0, max_value=1.0,
                                ridge_lambda=DEFAULT_LAMBDA,
                                smooth_lambda=DEFAULT_SMOOTH_LAMBDA):
        y = np.asarray(raw_pep, dtype=float)
        N = len(y)
        if N == 0:
            return []
        x = np.linspace(0.0, 1.0, N)
        fitted = self._fit_ispline(x, y, ridge_lambda=ridge_lambda,
                                   smooth_lambda=smooth_lambda,
                                   min_value=min_value, max_value=max_value)
        return fitted.tolist()

    def ispline_non_decreasing_xy(self, x_scores, y, min_value=0.0, max_value=1.0,
                                   ridge_lambda=DEFAULT_LAMBDA,
                                   smooth_lambda=DEFAULT_SMOOTH_LAMBDA):
        x_sc = np.asarray(x_scores, dtype=float)
        y_arr = np.asarray(y, dtype=float)
        N = len(y_arr)
        if N == 0:
            return []
        x_neg = -x_sc
        x_min, x_max = x_neg.min(), x_neg.max()
        if x_max <= x_min:
            return self.ispline_non_decreasing(
                y_arr, min_value=min_value, max_value=max_value,
                ridge_lambda=ridge_lambda, smooth_lambda=smooth_lambda,
            )
        x_norm = (x_neg - x_min) / (x_max - x_min)
        order = np.argsort(x_norm, kind="mergesort")
        x_sorted = x_norm[order]
        y_sorted = y_arr[order]
        fitted_sorted = self._fit_ispline(x_sorted, y_sorted,
                                          ridge_lambda=ridge_lambda,
                                          smooth_lambda=smooth_lambda,
                                          min_value=min_value, max_value=max_value)
        fitted = np.empty(N)
        fitted[order] = fitted_sorted
        return fitted.tolist()


class TDCIsotonicPEP(IsotonicRegression):
    def __init__(self, pava=False, max_iter=DEFAULT_MAX_ITER):
        self.pava = pava
        self.max_iter = max_iter

    def tdc_to_pep(self, df_obs, pava=None, max_iter=None):
        pava = self.pava if pava is None else pava

        df_sorted = df_obs.sort_values(by="score", ascending=False, kind="mergesort").reset_index(drop=True)
        scores = df_sorted["score"].values.astype(float)
        is_decoy = (df_sorted["label"] == 1).astype(float).values

        score_span = scores.max() - scores.min() if len(scores) > 1 else 1.0
        delta = max(score_span * 1e-6, 1e-12)
        scores_aug = np.concatenate([[scores.max() + delta], scores])
        is_decoy_aug = np.concatenate([[0.0], is_decoy])

        if pava:
            fitted = self.pava_non_decreasing(
                list(is_decoy_aug), [1] * len(is_decoy_aug),
                min_value=1e-20, max_value=1.0 - 1e-20,
            )
        else:
            fitted = self.ispline_non_decreasing_xy(
                scores_aug, is_decoy_aug,
                min_value=1e-20, max_value=1.0 - 1e-20,
            )

        fitted = np.array(fitted)
        decoy_rate = fitted[1:]
        with np.errstate(divide="ignore", invalid="ignore"):
            pep = decoy_rate / (1.0 - decoy_rate)
            pep = np.clip(pep, 0.0, 1.0)

        df_sorted["PEP"] = pep
        df_sorted["orig_idx"] = df_sorted["orig_idx"].astype(int)
        df_result = df_sorted.sort_values(by="orig_idx", kind="mergesort").reset_index(drop=True)
        return df_result["PEP"].reset_index(drop=True)


class IsotonicPEP(PreProcessing, TDCIsotonicPEP):
    def __init__(self, pava=False, max_iter=DEFAULT_MAX_ITER):
        PreProcessing.__init__(self)
        TDCIsotonicPEP.__init__(self, pava=pava, max_iter=max_iter)

    def calc_q_from_pep(self, pep_array):
        pep = np.asarray(pep_array, dtype=float)
        q = np.cumsum(pep) / np.arange(1, len(pep) + 1)
        q = np.maximum.accumulate(q)
        return q

    def q_to_pep(self, q_values, scores=None, pava=False, max_iter=DEFAULT_MAX_ITER):
        if not isinstance(q_values, pd.Series):
            q_series = pd.Series(q_values)
        else:
            q_series = q_values.copy()
        q_list = q_series.values.tolist()
        n = len(q_list)

        qn = []
        for i in range(n):
            qn.append(q_list[i] * (i + 1))
            if i < n - 1 and q_list[i] > q_list[i + 1]:
                raise ValueError("q-values must be non-decreasing.")
        raw_pep = [qn[0]] + [qn[i] - qn[i - 1] for i in range(1, n)]

        if n > 0:
            pseudo = 0.5 / n
            raw_pep = [v + pseudo for v in raw_pep]

        if scores is None:
            if pava:
                final_pep = self.pava_non_decreasing(raw_pep, [1] * n, min_value=PEP_CLIP_LO)
            else:
                final_pep = self.ispline_non_decreasing(raw_pep, min_value=PEP_CLIP_LO)
        else:
            sc = np.asarray(scores, dtype=float)
            if len(sc) != n:
                raise ValueError("q_values and scores must have the same length.")
            if pava:
                ord_desc = np.argsort(-sc, kind="mergesort")
                rp_arr = np.array(raw_pep)
                fitted_sorted = self.pava_non_decreasing(rp_arr[ord_desc].tolist(), [1] * n, min_value=PEP_CLIP_LO)
                inv_ord = np.empty(n, dtype=int)
                inv_ord[ord_desc] = np.arange(n)
                final_pep = np.array(fitted_sorted)[inv_ord].tolist()
            else:
                final_pep = self.ispline_non_decreasing_xy(sc, raw_pep, min_value=PEP_CLIP_LO)

        return pd.Series(final_pep, index=q_series.index)

    def pep_regression(self, q_values=None, obs=None, target_scores=None,
                       calc_q_from_fdr=False, calc_q_from_pep=False,
                       method="q2pep", pava=None, max_iter=None):
        pava = self.pava if pava is None else pava

        if method in ("q2pep", "qns2pep"):
            if calc_q_from_fdr:
                if obs is None:
                    raise ValueError("obs must be provided when calc_q_from_fdr is True.")
                fdr_series, q1_series = self.calc_q_from_fdr(obs=obs)
                target_mask = obs[:, 1] == 0
                q_input_raw = q1_series[target_mask].values
                fdr_array = fdr_series[target_mask].values
            else:
                if q_values is None:
                    raise ValueError("Provide q_values or enable calc_q_from_fdr.")
                q_input_raw = np.asarray(q_values, dtype=float)
                fdr_array = None

            scores_for_fit_raw = None
            if method == "qns2pep":
                if target_scores is not None:
                    scores_for_fit_raw = np.asarray(target_scores, dtype=float)
                elif obs is not None:
                    t_mask = obs[:, 1] == 0
                    scores_for_fit_raw = obs[t_mask, 0].astype(float)
                else:
                    raise ValueError("qns2pep requires either target_scores or obs with scores.")
                if len(scores_for_fit_raw) != len(q_input_raw):
                    raise ValueError("Number of target scores does not match q_values length.")

            idx_sorted = np.argsort(q_input_raw, kind="mergesort")
            q1_sorted = q_input_raw[idx_sorted]
            scores_for_fit = scores_for_fit_raw[idx_sorted] if scores_for_fit_raw is not None else None

            pep_sorted = self.q_to_pep(
                q_values=q1_sorted,
                scores=scores_for_fit,
                pava=pava,
            ).values
            q2_sorted = self.calc_q_from_pep(pep_sorted) if calc_q_from_pep else None

            n = len(q_input_raw)
            pep_array = np.empty(n)
            pep_array[idx_sorted] = pep_sorted
            q1_array = np.empty(n)
            q1_array[idx_sorted] = q1_sorted
            q2_array = None
            if q2_sorted is not None:
                q2_array = np.empty(n)
                q2_array[idx_sorted] = q2_sorted
            return fdr_array, q1_array, pep_array, q2_array

        elif method == "tdc2pep":
            if obs is None:
                raise ValueError("tdc2pep requires obs.")
            df_obs = self.process_obs(obs=obs)
            pep_all = self.tdc_to_pep(df_obs, pava=pava)

            target_idx = np.where(obs[:, 1] == 0)[0]
            pep_target_orig = pep_all.values[target_idx]
            scores_target = obs[target_idx, 0].astype(float)
            order_score_desc = np.argsort(-scores_target, kind="mergesort")
            pep_sorted = pep_target_orig[order_score_desc]

            if calc_q_from_fdr:
                fdr_series, q1_series = self.calc_q_from_fdr(obs=obs)
                fdr_sorted = fdr_series.values[target_idx][order_score_desc]
                q1_sorted = q1_series.values[target_idx][order_score_desc]
            else:
                fdr_sorted = None
                q1_sorted = None

            q2_sorted = self.calc_q_from_pep(pep_sorted) if calc_q_from_pep else None

            rev_order = np.empty_like(order_score_desc)
            rev_order[order_score_desc] = np.arange(len(order_score_desc))

            pep_array = pep_sorted[rev_order]
            q2_array = q2_sorted[rev_order] if q2_sorted is not None else None
            fdr_array = fdr_sorted[rev_order] if fdr_sorted is not None else None
            q1_array = q1_sorted[rev_order] if q1_sorted is not None else None
            return fdr_array, q1_array, pep_array, q2_array

        else:
            raise ValueError("Unknown method. Use 'q2pep', 'qns2pep', or 'tdc2pep'.")
