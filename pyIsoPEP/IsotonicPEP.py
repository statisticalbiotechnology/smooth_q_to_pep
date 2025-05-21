import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

# =============================================================================
# Module-Level Default Parameters
# =============================================================================
DEFAULT_MAX_ITER = 5000     # Default maximum iterations for iterative solvers.
DEFAULT_TOL = 1e-6          # Default tolerance for convergence.
DEFAULT_NUM_BINS = 10000    # hard cap after adaptive binning
DEFAULT_MAX_BIN_WIDTH = 20  # hard cap on rows per bin  (early region)
DEFAULT_MIN_DECOYS = 3      # aim for ≥ this many decoys per bin
DEFAULT_SKEW_FACTOR = 0.75  # <1 --> knots left-biased
DEFAULT_LAMBDA = 1e-6       # ridge penalty
EARLY_WEIGHT_FACTOR = 1.0   # 0 = off, 1 = linear boost towards top scores


# =============================================================================
# PreProcessing Class: Handle observation input and optionally compute FDRs and 
# q-values from target and decoy observations.
# =============================================================================
class PreProcessing:
    def process_obs(self, obs):
        """
        Build a DataFrame with columns:
          - score: float
          - label: 0 (target) or 1 (decoy)
          - orig_idx: original row index for order restoration
        """
        if not isinstance(obs, np.ndarray) or obs.ndim != 2 or obs.shape[1] != 2:
            raise ValueError("Input observations obs must be a numpy array with shape (n, 2).")
        scores = obs[:, 0].astype(float)
        labels = obs[:, 1].astype(int)
        if not set(np.unique(labels)).issubset({0, 1}):
            raise ValueError("Labels must be 0 (target) or 1 (decoy)")
        df = pd.DataFrame({"score": scores, "label": labels})
        df["orig_idx"] = np.arange(len(df))
        return df

    def calc_q_from_fdr(self, obs):
        """
        Compute FDRs and q-values from target and decoy observations.
        Returns a tuple of two pandas Series (FDR, q-value) aligned to original order.
        """
        df = self.process_obs(obs)
        df_sorted = df.sort_values(by="score", ascending=False, kind="mergesort").reset_index(drop=True)
        df_sorted["cumulative_target"] = (df_sorted["label"] == 0).cumsum()
        df_sorted["cumulative_decoy"] = (df_sorted["label"] == 1).cumsum()
        df_sorted["FDR"] = (df_sorted["cumulative_decoy"] + 0.5) / df_sorted["cumulative_target"]

        q = df_sorted["FDR"].values.copy()
        for i in range(len(q) - 2, -1, -1):
            q[i] = min(q[i], q[i + 1])
        df_sorted["q-value"] = q
        # align both to original order
        df_result = df_sorted.sort_values("orig_idx").reset_index(drop=True)
        return df_result["FDR"], df_result["q-value"]


# =============================================================================
# Base Class for Isotonic Regression in Real Space
# -----------------------------------------------------------------------------
# This class implements methods for enforcing a non-decreasing constraint 
# on a sequence of values. It provides:
#
# Regression:
#
#   - pava_non_decreasing():
#       Implements the standard Pool-Adjacent-Violators Algorithm (PAVA)
#       for stepwise-constant isotonic regression.
#
#   - ispline_non_decreasing():
#       Performs I-Spline-based monotonic regression on a given sequence.
#
#   - constrained_least_squares():
#       A helper method using FISTA to solve constrained least squares 
#       problems.
#
# Interpolation:
#
#   - pava_non_decreasing_interpolation():
#       A variant that performs I-Spline or Fritsch-Carlson interpolation 
#       between PAVA-derived block centers to produce a smooth, monotonic 
#       cubic spline.
#
#   - ispline_monotonic_interpolate():
#       Performs monotonic interpolation using I-Splines.
#
#   - fritsch_carlson_monotonic_interpolate():
#       Performs monotonic interpolation using Fritsch-Carlson monotonic 
#       cubic interpolation (PCHIP).
#
#   - boundary_derivative_fritsch_carlson():
#       Computes the boundary derivatives used in Fritsch-Carlson 
#       interpolation to maintain monotonicity at endpoints.
# =============================================================================
class IsotonicRegression:
    def __init__(self):
        pass

    def pava_non_decreasing(self, values, counts, min_value=0.0, max_value=1.0):
        """
        Enforces non-decreasing constraint on a sequence using PAVA.

        Each input value is repeated according to its corresponding count.
        After merging adjacent blocks that violate the monotonicity, each block's
        average is clamped to the interval [min_value, max_value].

        Parameters:
            values    : list of floats
                        The data values.
            counts    : list of ints
                        The number of times each value is repeated.
            min_value : float, optional (default=0.0)
                        Lower bound for clamping.
            max_value : float, optional (default=1.0)
                        Upper bound for clamping.

        Returns:
            list of floats:
                The expanded non-decreasing sequence where each block's average 
                (clamped between min_value and max_value) is repeated count times.
        """
        if len(values) != len(counts):
            raise ValueError("values and counts must have the same length.")
        n = len(values)
        if n == 0:
            return []

        # Use a stack to store blocks.
        # Each block is represented as a dict with keys: 'sum', 'count', 'avg'
        stack = []
        for i in range(n):
            sum_i = values[i] * counts[i]
            new_block = {'sum': sum_i, 'count': counts[i], 'avg': values[i]}
            stack.append(new_block)

            # Merge adjacent blocks while the non-decreasing constraint is violated.
            while len(stack) > 1:
                top = stack[-1]
                sec_top = stack[-2]
                if sec_top['avg'] > top['avg']:
                    merged_sum = sec_top['sum'] + top['sum']
                    merged_count = sec_top['count'] + top['count']
                    merged_avg = merged_sum / merged_count
                    stack.pop()
                    stack.pop()
                    stack.append({'sum': merged_sum, 'count': merged_count, 'avg': merged_avg})
                else:
                    break

        # Expand the final solution: repeat each block's clamped average according to its count.
        result = []
        for block in stack:
            # Clamp the block average to [min_value, max_value]
            clamped_avg = min(max(block['avg'], min_value), max_value)
            result.extend([clamped_avg] * block['count'])
        return result
    
    def constrained_least_squares(self, X, y, lower_bounds, upper_bounds, max_iter=DEFAULT_MAX_ITER, tol=DEFAULT_TOL):
        """
        Solves a constrained least squares problem using an accelerated projected
        gradient method (FISTA). The objective is to minimize 0.5 * ||Xc - y||^2 
        subject to lower_bounds <= c <= upper_bounds.

        Parameters:
            X            : numpy.ndarray, shape (n_samples, n_features)
                           The design matrix.
            y            : numpy.ndarray, shape (n_samples,)
                           The response vector.
            lower_bounds : numpy.ndarray, shape (n_features,)
                           The lower bounds for the coefficients.
            upper_bounds : numpy.ndarray, shape (n_features,)
                           The upper bounds for the coefficients.
            max_iter     : int, optional (default=DEFAULT_MAX_ITER)
                           Maximum number of iterations.
            tol          : float, optional (default=DEFAULT_TOL)
                           Tolerance for convergence.

        Returns:
            numpy.ndarray, shape (n_features,):
                The solution vector c that minimizes the constrained least squares.
        """
        # Ensure  inputs are numpy arrays.
        X = np.asarray(X)
        y = np.asarray(y)
        lower_bounds = np.asarray(lower_bounds)
        upper_bounds = np.asarray(upper_bounds)
        
        n_features = X.shape[1]
        
        # Feasible initialization: choose midpoint if both bounds are finite; otherwise, choose a finite bound or zero.
        c = np.empty(n_features)
        for i in range(n_features):
            if np.isfinite(lower_bounds[i]) and np.isfinite(upper_bounds[i]):
                c[i] = (lower_bounds[i] + upper_bounds[i]) / 2.0
            elif np.isfinite(lower_bounds[i]):
                c[i] = lower_bounds[i]
            elif np.isfinite(upper_bounds[i]):
                c[i] = upper_bounds[i]
            else:
                c[i] = 0.0

        # Initialize FISTA parameters.
        z = c.copy()  # The extrapolated point.
        t = 1.0       # Momentum parameter.
        
        # Estimate Lipschitz constant L for the gradient: L = ||X||_2^2.
        L = np.linalg.norm(X, 2)**2
        step_size = 1.0 / L

        for iteration in range(max_iter):
            c_old = c.copy()
            # Compute gradient at the extrapolated point z: grad = X^T(Xz - y)
            grad = X.T @ (X @ z - y)
            # Gradient descent step with projection onto the feasible set.
            c = z - step_size * grad
            c = np.clip(c, lower_bounds, upper_bounds)
            
            # Update momentum term.
            t_new = (1 + np.sqrt(1 + 4 * t * t)) / 2.0
            z = c + ((t - 1) / t_new) * (c - c_old)
            t = t_new
            
            # Check convergence.
            if np.linalg.norm(c - c_old) < tol:
                break
        return c
    
    def bin_data(self, x, y, is_decoy, max_bins=DEFAULT_NUM_BINS, min_decoys=DEFAULT_MIN_DECOYS, max_bin_width=DEFAULT_MAX_BIN_WIDTH):
        """
        Adaptive binning that
          • tries to collect ≥ target_decoys per bin, AND
          • never lets a bin exceed max_bin_width rows, AND
          • never outputs more than max_bins bins.
        
        Parameters:
                x, y          : Original observations (sorted in descending score order)
                is_decoy      : Boolean, 1 indicates the row is a decoy
                min_decoys    : Desired minimum number of decoys per bin
                max_bin_width : Maximum number of observations allowed per bin (to prevent overly large bins in high-score regions)
      
        Returns:
                x_bin, y_bin, w_bin  (all numpy arrays).
        """
        n = len(x)
        if n == 0 or max_bins <= 0:
            return np.array([]), np.array([]), np.array([])

        x_out, y_out, w_out = [], [], []
        start, decoys = 0, 0

        for i in range(n):
            if is_decoy[i]:
                decoys += 1
            
            need_close = (
                decoys >= min_decoys or                         # enough decoys
                (i - start + 1) >= max_bin_width or             # bin too wide
                len(x_out) + 1 >= max_bins or                   # hit bin cap
                i == n - 1                                      # last row
            )

            if need_close:
                end = i + 1  # exclusive
                size = end - start
                x_avg = x[start:end].mean()
                y_avg = y[start:end].mean()

                x_out.append(x_avg)
                y_out.append(y_avg)
                w_out.append(size)

                start, decoys = end, 0

        return (np.asarray(x_out), np.asarray(y_out), np.asarray(w_out, dtype=float))

    def ispline_non_decreasing(self, raw_pep, is_decoy=None, min_value=0.0, max_value=1.0, max_iter=DEFAULT_MAX_ITER, skew_factor=DEFAULT_SKEW_FACTOR, max_bins=DEFAULT_NUM_BINS, ridge_lambda=DEFAULT_LAMBDA, early_weight_factor=EARLY_WEIGHT_FACTOR, min_decoys=DEFAULT_MIN_DECOYS, max_bin_width=DEFAULT_MAX_BIN_WIDTH):
        """
        Performs I-Spline based isotonic regression on raw PEP values.
        
        The model fitted is:
            f(x) = c_0 + c_1 I_1(x) + c_2 I_2(x) + ... + c_m I_m(x)
        where the I-Spline basis functions I_j(x) are constructed on a normalized 
        domain and non-negative constraints on c_1, ..., c_m ensure that f(x) is 
        monotonic.

        Parameters:
            raw_pep   : list or numpy array
                        Raw PEP values.
            min_value : float, optional (default=0.0)
                        Lower bound for clamping the fitted values.
            max_value : float, optional (default=1.0)
                        Upper bound for clamping the fitted values.
            max_iter  : int, optional (default=DEFAULT_MAX_ITER)
                        Maximum iterations for the constrained least squares solver.
            skew_factor : float
                        < 1 pushes knots towards the low-score end; > 1 towards high scores.
            max_bins : int
                        Upper bound on the number of bins produced by the adaptive binning.
            max_knots : int
                        Upper cap on the number of spline intervals (k).
            ridge_lambda : float
                        Small ℓ² regularization term λ.
            early_weight_factor : float in [0, inf)
                        0  : no extra emphasis;  
                        1  : linear down-weighting from first to last bin.

        Returns:
            list of floats:
                The fitted monotone PEP values (clamped between min_value and max_value).
        """
        # Convert raw_pep to numpy array
        y = np.array(raw_pep)
        N = len(y)
        if N == 0:
            return []
        
        # 1. Prepare binned representation (speeds up regression)
        # Normalize x positions to [0, 1] and perform adaptive binning
        x = np.linspace(0, 1, N)
        if is_decoy is not None:
            x_bin, y_bin, w_bin = self.bin_data(
                x, y, is_decoy,
                max_bins=max_bins,
                min_decoys=min_decoys,
                max_bin_width=max_bin_width
            )
        else:   # equal-size bins
            x_bin, y_bin, w_bin = self.bin_data(x, y, max_bins=max_bins)
        n_bin = len(x_bin)

        # Emphasise early bins (optional)
        if early_weight_factor > 0.0 and n_bin > 1:
            scale = 1.0 + early_weight_factor * (n_bin - 1 - np.arange(n_bin)) / (n_bin - 1)
            w_bin *= scale

        # 2. Construct I‑Spline design matrix
        k = int(np.sqrt(n_bin))
        # compute adaptive knots and build knot vector
        # total knots = k + 1  (including the two ends)
        knots = [x_bin[0]]
        for i in range(1, k):
            q = 1.0 - (1.0 - i / k) ** skew_factor      # double q = 1 - pow(...)
            idx = int(round(q * (n_bin - 1)))                 # size_t idx = q*(x.size()-1)
            knots.append(x_bin[idx])
        knots.append(x_bin[-1])
        knots = np.asarray(knots)
        m = len(knots) - 1        # number of basis functions (intervals)

        # Construct the cubic I-Spline basis matrix for binned data.
        # For each interval [knots[j], knots[j+1]], define:
        #    I_j(x) = 0                         if x < knots[j],
        #             3u^2 - 2u^3               if knots[j] <= x < knots[j+1],
        #             1                         if x >= knots[j+1],
        # where u = (x - knots[j]) / (knots[j+1] - knots[j]).
        B = np.zeros((n_bin, m))
        for j in range(m):
            # Compute the normalized coordinate u.
            u = (x_bin - knots[j]) / (knots[j + 1] - knots[j])
            u = np.clip(u, 0, 1)
            B[:, j] = np.where(
                x_bin < knots[j],
                0.0,
                np.where(
                    x_bin >= knots[j+1],
                    1.0,
                    3 * u**2 - 2 * u**3
                )
            )
        # Add an intercept column to the design matrix.
        X = np.column_stack((np.ones(n_bin), B))
        # apply weights  (W½ = sqrt(w))
        w_sqrt = np.sqrt(w_bin)
        Xw = X * w_sqrt[:, None]
        yw = y_bin * w_sqrt

        # Ridge augmentation: add ridge rows  (√λ I, 0)
        p = X.shape[1]            # = m + 1 parameters
        X_aug = np.vstack([Xw, np.sqrt(ridge_lambda) * np.eye(p)])
        y_aug = np.concatenate([yw, np.zeros(p)])
        
        # Solve the constrained least squares problem:
        # minimize ||Xc - y||^2
        # subject to: c[1: ] >= 0 (the intercept c[0] is unconstrained).
        lower_bounds = np.concatenate(([-np.inf], np.zeros(m)))
        # upper_bounds = np.full(m + 1, np.inf)
        upper_bounds = np.full(p, np.inf)
        res = lsq_linear(X_aug, y_aug, bounds=(lower_bounds, upper_bounds))
        c = res.x
        # c = self.constrained_least_squares(X_aug, y_aug, lower_bounds, upper_bounds, max_iter=max_iter)
        
        # 3. Evaluate fitted spline on original x grid
        # build same basis but on x_full
        B_full = np.zeros((N, m))
        for j in range(m):
            u = (x - knots[j]) / (knots[j + 1] - knots[j])
            u = np.clip(u, 0.0, 1.0)
            B_full[:, j] = np.where(
                x < knots[j],
                0.0,
                np.where(x >= knots[j + 1],
                        1.0,
                        3 * u ** 2 - 2 * u ** 3)
            )
        
        X_full = np.column_stack((np.ones(N), B_full))

        # Compute the fitted values and clamp them to [min_value, max_value]
        fitted = X_full @ c
        fitted = np.clip(fitted, min_value, max_value)
        return fitted.tolist()
    
    def ispline_monotonic_interpolate(self, xB, yB, xEval):
        """
        Performs I-Spline monotonic interpolation using a cumulative smooth-step basis.
        
        The interpolant is represented as:
            f(x) = yB[0] + sum_{j=1}^{n-1} (yB[j]-yB[j-1]) * I_j(x),    
        where for each interval j (with 1 <= j < n) the basis function I_j is defined by:
        
            I_j(x) = 0                           if x <= xB[j-1],
                     3u^2 - 2u^3                 if xB[j-1] < x < xB[j],
                     1                           if x >= xB[j],
        
        with u = (x - xB[j-1]) / (xB[j]-xB[j-1]). This construction guarantees f(xB[j]) = yB[j] 
        and yields a smooth (C^1) monotone interpolant.
        
        Parameters:
            xB    : numpy.ndarray
                    Array of block centers (assumed to be strictly increasing).
            yB    : numpy.ndarray
                    Array of block averages (monotonic non-decreasing).
            xEval : numpy.ndarray
                    Array of x-values at which to evaluate the interpolant.

        Returns:
            numpy.ndarray:
                The interpolated values evaluated at xEval.
        """
        n = len(xB)
        if n == 0:
            return np.array([])
        if n == 1:
            return np.full_like(xEval, yB[0])
            
        # Compute the differences (which are nonnegative due to monotonicity)
        d = np.diff(yB)
        
        # Initialize with the base value for interpolation.
        f_interp = np.full_like(xEval, yB[0], dtype=float)
        
        # For each interval between block centers, add the smooth cumulative contribution.
        for j in range(1, n):
            x_left = xB[j-1]
            x_right = xB[j]
            # Compute the normalized variable u; outside the interval u is clipped to 0 or 1.
            u = (xEval - x_left) / (x_right - x_left)
            u = np.clip(u, 0, 1)
            # Smoothstep function: a common cubic that satisfies S(0)=0, S(1)=1, S'(0)=S'(1)=0.
            I_j = 3 * u**2 - 2 * u**3
            f_interp += d[j-1] * I_j
        
        return f_interp

    def boundary_derivative_fritsch_carlson(self, xB, s, idx):
        """
        Compute the boundary derivative d[0] or d[m-1] for Fritsch-Carlson interpolation.
        
        For idx=0, use xB[0..2]; for idx=m-1, use xB[m-3..m-1].
        """
        m = len(xB)
        # If only 2 points in total, slope is constant, but that case is handled outside.
        if m < 3:
            # Fallback if needed, though normally handled above
            return s[0] if idx == 0 else s[-1]

        if idx == 0:
            # Use xB[0], xB[1], xB[2]
            h0 = xB[1] - xB[0]
            h1 = xB[2] - xB[1]
            s0 = s[0]  # slope between xB[0] and xB[1]
            s1 = s[1]  # slope between xB[1] and xB[2]
            d0 = ((2*h0 + h1)*s0 - h0*s1) / (h0 + h1)
            # Monotonicity check
            if d0 * s0 <= 0:
                d0 = 0
            elif abs(d0) > abs(3*s0):
                d0 = 3*s0
            return d0
        else:
            # idx == m-1 => use xB[m-3], xB[m-2], xB[m-1]
            h1 = xB[m-2] - xB[m-3]
            h2 = xB[m-1] - xB[m-2]
            s1 = s[m-2]  # slope between xB[m-2] and xB[m-1]
            s0 = s[m-3]  # slope between xB[m-3] and xB[m-2]
            dN = ((2*h2 + h1)*s1 - h2*s0) / (h1 + h2)
            if dN * s1 <= 0:
                dN = 0
            elif abs(dN) > abs(3*s1):
                dN = 3*s1
            return dN

    def fritsch_carlson_monotonic_interpolate(self, xB, yB, xEval):
        """
        Fritsch-Carlson monotonic cubic interpolation for a strictly increasing
        sequence xB with corresponding yB (non-decreasing).
        
        Returns the interpolated y-values at xEval (array).
        
        Reference:
        F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic Interpolation" 
        SIAM Journal on Numerical Analysis, 1980
        F. N. Fritsch and J. Butland, "A Method for Constructing Local Monotone Piecewise Cubic Interpolants"
        SIAM Journal on Scientific and Statistical Computing, 1984
        C. Moler, "Numerical Computing with Matlab" 2004
        https://epubs.siam.org/doi/epdf/10.1137/1.9780898717952.ch3 Chapter 3.6, Page 14
        
        Parameters:
            xB: numpy array of shape (m,) with strictly increasing x-coordinates (block centers).
            yB: numpy array of shape (m,) with non-decreasing y-coordinates (block averages).
            xEval: numpy array of x-values at which we want to evaluate.

        Returns:
            A numpy array of interpolated values at each xEval.
        """
        m = len(xB)
        if m == 0:
            return np.array([])
        if m == 1:
            # Only one block => constant function
            return np.full_like(xEval, yB[0])

        # 1) Compute slopes between consecutive points
        h = xB[1:] - xB[:-1]          # interval lengths
        s = (yB[1:] - yB[:-1]) / h    # slopes

        # 2) Compute initial derivatives at each xB[i] (d[i])
        d = np.zeros(m, dtype=float)

        # Handle the simple case m=2 directly
        if m == 2:
            d[0] = s[0]
            d[1] = s[0]
        else:
            # More than two points => boundary derivatives
            d[0]   = self.boundary_derivative_fritsch_carlson(xB, s, idx=0)
            d[m-1] = self.boundary_derivative_fritsch_carlson(xB, s, idx=m-1)

            # Interior derivatives
            for i in range(1, m-1):
                if s[i-1] * s[i] <= 0:
                    # If slopes change sign or one is zero, derivative is 0 to preserve monotonicity
                    d[i] = 0.0
                else:
                    # Weighted harmonic mean of s[i-1] and s[i]
                    alpha = 3 * (h[i] + h[i-1])
                    d[i] = alpha / ( ((2*h[i]+h[i-1])/s[i-1]) + ((h[i]+2*h[i-1])/s[i]) )

        # 3) Compute polynomial coefficients for each interval
        # For interval i, define:
        #   c1 = yB[i]
        #   c2 = d[i]
        #   c3 = (3*s[i] - 2*d[i] - d[i+1]) / h[i]
        #   c4 = (d[i] + d[i+1] - 2*s[i]) / (h[i]^2)
        c1 = yB[:-1]
        c2 = d[:-1]
        c3 = np.zeros(m-1, dtype=float)
        c4 = np.zeros(m-1, dtype=float)

        for i in range(m-1):
            d_i   = d[i]
            d_ip1 = d[i+1]
            s_i   = s[i]
            hi    = h[i]
            c3[i] = (3*s_i - 2*d_i - d_ip1) / hi
            c4[i] = (d_i + d_ip1 - 2*s_i) / (hi**2)

        # 4) Evaluate piecewise polynomial at each xEval
        y_out = np.zeros_like(xEval, dtype=float)

        # For each xEval, find interval via binary search
        idxs = np.searchsorted(xB, xEval) - 1
        idxs = np.clip(idxs, 0, m-2)  # valid range for intervals

        for i in range(len(xEval)):
            seg = idxs[i]
            dx = xEval[i] - xB[seg]
            y_out[i] = c1[seg] + c2[seg]*dx + c3[seg]*(dx**2) + c4[seg]*(dx**3)

        return y_out
    
    def pava_non_decreasing_interpolation(self, x, y, ip_algo="ispline", center_method="mean", min_y=0.0, max_y=1.0):
        """
        Computes a smooth, non-decreasing interpolation using a variant of PAVA with
        I-Spline or Fritsch-Carlson monotonic cubic interpolation.
        
        Steps:
        (a) Group the data into extended blocks (each point initially forms a block).
        (b) Merge adjacent blocks via PAVA (each merged block spans indices [startIdx, endIdx]).
        (c) For each final block, compute the x-center based on the specified center_method:
                - "mean": the average of x values in that block.
                - "median": the midpoint between the first and last x in the block, 
                            i.e., (x[startIdx] + x[endIdx]) / 2.
        (d) Collect these (xBlock, yBlock) points and build a monotonic cubic interpolant 
            using I-Spline or Fritsch-Carlson method.
        (e) Evaluate this spline at each original x[i] and clamp the results to [min_y, max_y].
        
        Parameters:
            x             : list of floats
                            Sorted positions.
            y             : list of floats
                            Data values.
            ip_algo       : str, optional (default="ispline")
                            Interpolation algorithm to use on block centers; "ispline" or "pchip".
            center_method : str, optional (default="mean")
                            Method for computing the block center; "mean" or "median".
            min_y         : float, optional (default=0.0)
                            Lower bound for clamping.
            max_y         : float, optional (default=1.0)
                            Upper bound for clamping.

        Returns:
            list of floats:
                The interpolated, non-decreasing values corresponding to the input x.
        """
        if len(x) != len(y):
            raise ValueError("x and y must have the same length.")
        n = len(y)
        if n == 0:
            return []

        # Convert to numpy arrays for convenience
        x_arr = np.array(x, dtype=float)
        y_arr = np.array(y, dtype=float)

        # (a) Build an extended block for each point.
        blocks = []
        for i in range(n):
            block = {'sum': y_arr[i], 'count': 1, 'avg': y_arr[i], 'startIdx': i, 'endIdx': i}
            blocks.append(block)

        # (b) Merge blocks using standard PAVA.
        stack = []
        for block in blocks:
            stack.append(block)
            while len(stack) > 1:
                top = stack[-1]
                sec_top = stack[-2]
                if sec_top['avg'] > top['avg']:
                    new_sum = sec_top['sum'] + top['sum']
                    new_count = sec_top['count'] + top['count']
                    new_avg = new_sum / new_count
                    sIdx = sec_top['startIdx']
                    eIdx = top['endIdx']
                    stack.pop()
                    stack.pop()
                    stack.append({
                        'sum': new_sum,
                        'count': new_count,
                        'avg': new_avg,
                        'startIdx': sIdx,
                        'endIdx': eIdx
                    })
                else:
                    break

        # (c) For each final block, compute the x-center.
        xBlock = []
        yBlock = []
        for blk in stack:
            if center_method == "mean":
                sum_x = np.sum(x_arr[blk['startIdx'] : blk['endIdx'] + 1])
                length = (blk['endIdx'] - blk['startIdx'] + 1)
                center_x = sum_x / length
            elif center_method == "median":
                center_x = 0.5 * (x_arr[blk['startIdx']] + x_arr[blk['endIdx']])
            else:
                raise ValueError("Unknown center_method. Use 'mean' or 'median'.")

            xBlock.append(center_x)
            yBlock.append(blk['avg'])

        xBlock = np.array(xBlock, dtype=float)
        yBlock = np.array(yBlock, dtype=float)

        # (d) Use I-Spline monotonic cubic interpolation on (xBlock, yBlock)
        #     and evaluate at each original x[i].
        if ip_algo == "ispline":
            y_interp = self.ispline_monotonic_interpolate(xBlock, yBlock, x_arr)
        elif ip_algo == "pchip":
            y_interp = self.fritsch_carlson_monotonic_interpolate(xBlock, yBlock, x_arr)
        else:
            raise ValueError("Unknown ip_algo. Use 'ispline' or 'pchip'.")

        # (e) Clamp the result to [min_y, max_y] and return as list.
        y_clamped = np.clip(y_interp, min_y, max_y)
        return y_clamped.tolist()
    

# =============================================================================
# TDCIsotonicPEP Class: Regression for Target-Decoy Competition Data
# =============================================================================
class TDCIsotonicPEP(IsotonicRegression):
    def __init__(self, regression_algo="ispline", max_iter=DEFAULT_MAX_ITER):
        self.regression_algo = regression_algo
        self.max_iter = max_iter

    def tdc_binomial_regression(self, df_obs, regression_algo=None, max_iter=None):
        """
        Compute PEP values from observation data using isotonic regression.
        
        Steps:
          - (a) Sort by score descending.
          - (b) Prepend a pseudo observation with type 0.5.
          - (c) Apply PAVA or I-Spline regression on the binary sequence.
          - (d) Remove the pseudo observation.
          - (e) Compute PEP = decoy_prob/(1 - decoy_prob), clipped to [0,1].
          - (f) Restore original order.
        
        Parameters:
            df_obs         : pandas.DataFrame
                             DataFrame with observation data.
            regression_algo: str, optional (default="ispline")
                             Regression algorithm to use ("PAVA" or "ispline").
            max_iter       : int, optional
                             Maximum iterations for the regression; if None, uses the object's default.

        Returns:
            pandas.Series: A Series of target PEP values aligned with the original order.
        """
        regression_algo = regression_algo or self.regression_algo
        max_iter = max_iter or self.max_iter

        df_sorted = df_obs.sort_values(by="score", ascending=False, kind="mergesort").reset_index(drop=True)
        pseudo = pd.DataFrame({"score": [np.nan], "label": [0.5]})
        df_aug = pd.concat([pseudo, df_sorted], ignore_index=True)
        y_values = df_aug["label"].values
        is_dec = (y_values == 1).astype(bool)             # bool vector, pseudo row = False
        if self.regression_algo == "PAVA":
            fitted = self.pava_non_decreasing(list(y_values), [1] * len(y_values))
        elif self.regression_algo == "ispline":
            fitted = self.ispline_non_decreasing(list(y_values), is_decoy=is_dec, max_iter=max_iter)
        else:
            raise ValueError("Unknown regression_algo. Use 'PAVA' or 'ispline'.")
        
        fitted = np.array(fitted)
        fitted_decoy_prob = fitted[1:]  # remove pseudo observation
        with np.errstate(divide='ignore', invalid='ignore'):
            pep = fitted_decoy_prob / (1 - fitted_decoy_prob)
            pep = np.clip(pep, 0, 1)
        df_sorted["PEP"] = pep
        # Restore original order.
        df_sorted["orig_idx"] = df_sorted["orig_idx"].astype(int)
        df_result = df_sorted.sort_values(by="orig_idx", kind="mergesort").reset_index(drop=True)
        df_target = df_result[df_result["label"] == 0]
        return df_target["PEP"].reset_index(drop=True)
    

# =============================================================================
# IsotonicPEP Class: Unified Interface for PEP Estimation
# =============================================================================
class IsotonicPEP(PreProcessing, TDCIsotonicPEP):
    def __init__(self, regression_algo="ispline", ip_algo="ispline", center_method="mean", max_iter=DEFAULT_MAX_ITER):
        PreProcessing.__init__(self)
        TDCIsotonicPEP.__init__(self, regression_algo=regression_algo, max_iter=max_iter)
        self.ip_algo = ip_algo
        self.center_method = center_method
        
    def calc_q_from_pep(self, pep_array):
        """
        Given a PEP array, compute q-values via:
           q(i) = (1 / i) * cumsum(pep[0]...pep[i]),
        where pep_array is sorted in ascending order.

        Steps:
          - (a) sort pep_array ascending using mergesort,
          - (b) cumsum,
          - (c) divide by rank,
          - (d) restore original order.

        Returns:
            An array of estimated q-values aligned with the input order.
        """
        pep_array = np.array(pep_array, dtype=float)
        idx_sorted = np.argsort(pep_array, kind="mergesort")
        pep_sorted = pep_array[idx_sorted]
        csum = np.cumsum(pep_sorted)
        ranks = np.arange(1, len(pep_sorted) + 1)
        q_sorted = csum / ranks
        # restore
        q_est = np.zeros_like(pep_array)
        for i, idx in enumerate(idx_sorted):
            q_est[idx] = q_sorted[i]
        return q_est
    
    def q_to_pep(self, q_values, regression_algo="ispline", ip=False, ip_algo=None, center_method=None, max_iter=DEFAULT_MAX_ITER):
        """
        Compute smoothed PEP values from q-values (q2pep).

        qn[i] = q_values[i] * (i+1)
        raw_pep[0] = qn[0],
        for i >= 1: raw_pep[i] = qn[i] - qn[i-1].
        
        Parameters:
            q_values      : array-like or pandas.Series
                            Sorted in ascending order.
            regression_algo: str, optional (default="ispline")
                             Regression method ("PAVA" or "ispline").
            ip            : bool, optional (default=False)
                             If True and using PAVA, perform interpolation.
            ip_algo       : str, optional
                             interpolation algorithm to use on block centers;
                             defaults to self.ip_algo.
            center_method : str, optional
                             Method for computing block centers; defaults to self.center_method.
            max_iter      : int, optional (default=DEFAULT_MAX_ITER)
                             Maximum iterations.

        Returns:
            pandas.Series: Smoothed PEP values aligned with the original indices.
        """
        # Convert to Series and record original index.
        if not isinstance(q_values, pd.Series):
            q_series = pd.Series(q_values)
        else:
            q_series = q_values.copy()
        orig_idx = q_series.index.copy()

        # Sort q_values in ascending order.
        q_series_sorted = q_series.sort_values(ascending=True, kind="mergesort")
        q_list = q_series_sorted.values.tolist()
        n = len(q_list)
        qn = []
        for i in range(n):
            qn.append(q_list[i]*(i+1))
            if i < n - 1 and q_list[i] > q_list[i+1]:
                raise AssertionError("q_values must be non-decreasing")
        raw_pep = [qn[0]] + [qn[i] - qn[i-1] for i in range(1, n)]

        if regression_algo == "PAVA":
            if ip:
                x_positions = list(range(len(raw_pep)))
                final_pep = self.pava_non_decreasing_interpolation(x_positions, raw_pep, ip_algo=ip_algo or self.ip_algo, center_method=center_method or self.center_method)
            else:
                final_pep = self.pava_non_decreasing(raw_pep, [1] * len(raw_pep))
        elif regression_algo == "ispline":
            final_pep = self.ispline_non_decreasing(raw_pep, max_iter=max_iter)
        else:
            raise ValueError("Unknown regression_algo. Use 'PAVA' or 'ispline'.")
        
        pep_sorted = pd.Series(final_pep, index=q_series_sorted.index)
        pep_result = pep_sorted.reindex(orig_idx)
        return pep_result
    
    def dprob_to_pep(self, obs, regression_algo=None, max_iter=None):
        """
        Compute PEP values from target-decoy observations (d2pep).
        
        Parameters:
            obs            : Concatenated observation list [score, label_numeric] where
                             label_numeric is 0 for targets and 1 for decoys.
            regression_algo: str, optional (default="ispline")
                             Regression method ("PAVA" or "ispline").
            max_iter       : int, optional (default=DEFAULT_MAX_ITER)
                             Maximum iterations.

        Returns:
            pandas.Series: Target PEP values aligned with the original order.
        """
        if obs is None:
            raise ValueError("For d2pep, provide concatenated observations as obs.")
        df_obs = self.process_obs(obs=obs)
        return self.tdc_binomial_regression(df_obs, regression_algo=regression_algo, max_iter=max_iter)
    
    def pep_regression(self, q_values=None, obs=None, calc_q_from_fdr=False, calc_q_from_pep=False, 
                       method="q2pep", regression_algo=None, max_iter=None, ip=False, ip_algo=None, center_method=None):
        """
        Unified interface for computing PEP values, then optionally computing q-values from PEPs.
        
        Parameters:
            q_values : 1-D array-like, optional
                Target q-values.  Used only when method='q2pep' and
                calc_q_from_fdr is False.
            obs : ndarray (n, 2), optional
                Concatenated observation list [score, label_numeric] where
                label_numeric is 0 for targets and 1 for decoys. Mandatory
                for the d2pep workflow and for q2pep when calc_q_from_fdr
                is True.
            calc_q_from_fdr : bool, default False
                Compute a running FDR curve from obs and turn it into q-values.
                If method='q2pep' and this flag is False, the function expects
                user-supplied q_values.
            calc_q_from_pep : bool, default False
                After estimating the PEPs, derive q-values from them.
            method : {'q2pep', 'd2pep'}, default 'q2pep'
                Workflow selector.
            regression_algo : {'PAVA', 'ispline'}, optional
                Override the algorithm chosen when the object was constructed.
            max_iter : int, optional
                Maximum iterations for the I-Spline solver. Falls back to the
                instance's default when None.
            ip : bool, default False
                When regression_algo=='PAVA' in the q2pep path,
                activate smooth monotone interpolation between PAVA block centres.
            ip_algo : {'ispline', 'pchip'}, optional
                Interpolator backend used when ip is True (defaults to the
                object's self.ip_algo setting).
            center_method : {'mean', 'median'}, optional
                How to place the x-coordinate of each PAVA block when ip is
                True (defaults to self.center_method).

        Returns:
            tuple (fdr_array, q1_array, pep_array, q2_array)
            Arrays that were not requested (or cannot be produced) are returned as None.
            Arrays are restricted to targets and in the original order:
                * fdr_array : running FDR estimates  
                * q1_array  : q-values derived from FDR    
                * pep_array : smoothed PEPs  
                * q2_array  : q-values derived from the PEPs  
        """
        # Set defaults if not provided
        regression_algo = regression_algo or self.regression_algo
        max_iter = max_iter or self.max_iter
        ip_algo = ip_algo or self.ip_algo
        center_method = center_method or self.center_method

        if method == "q2pep":
            if calc_q_from_fdr:
                if obs is None:
                    raise ValueError("For q2pep, obs must be provided when --calc-q-from-fdr is used.")
                fdr_series, q1_series = self.calc_q_from_fdr(obs=obs)
                target_idx = np.where(obs[:, 1] == 0)[0]
                q1_array = q1_series.values[target_idx]
                q_input = q1_array[target_idx]
                fdr_array = fdr_series.values[target_idx]
            else:
                if q_values is None:
                    raise ValueError("Provide q-values or enable --calc-q-from-fdr.")
                q_input = np.array(q_values, dtype=float)
                q1_array = None
                fdr_array = None
            pep_series = self.q_to_pep(q_values=q_input, regression_algo=regression_algo, ip=ip, ip_algo=ip_algo, center_method=center_method, max_iter=max_iter)
            pep_array = pep_series.values
            if calc_q_from_pep:
                q2_array = self.calc_q_from_pep(pep_array)
            else:
                q2_array = None
            
            return fdr_array, q1_array, pep_array, q2_array

        elif method == "d2pep":
            if obs is None:
                raise ValueError("For d2pep, obs must be provided.")
            
            pep_series = self.dprob_to_pep(obs=obs, regression_algo=regression_algo, max_iter=max_iter)
            pep_array = pep_series.values
            if calc_q_from_fdr:
                fdr_series, q1_series = self.calc_q_from_fdr(obs=obs)
                target_idx = np.where(obs[:, 1] == 0)[0]
                q1_array = q1_series.values[target_idx]
                fdr_array = fdr_series.values[target_idx]
            else:
                fdr_array = None
                q1_array = None
            
            if calc_q_from_pep:
                q2_array = self.calc_q_from_pep(pep_array)
            else:
                q2_array = None
            
            return fdr_array, q1_array, pep_array, q2_array

        else:
            raise ValueError("Unknown method. Use 'q2pep' or 'd2pep'.")