import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

# =============================================================================
# Module-Level Default Parameters
# =============================================================================
DEFAULT_MAX_ITER = 5000   # Default maximum iterations for iterative solvers.
DEFAULT_TOL = 1e-6        # Default tolerance for convergence.

# =============================================================================
# Base Class for Isotonic Regression in Real Space
# -----------------------------------------------------------------------------
# This class implements methods for enforcing a non-decreasing constraint 
# on a sequence of values. It provides:
#
#   - pava_non_decreasing(): The standard Pool-Adjacent-Violators Algorithm 
#                            (PAVA) for stepwise constant regression.
#
#   - pava_non_decreasing_interpolation(): A variant that performs I-Spline 
#                            interpolation between PAVA-derived block centers 
#                            to yield a smooth monotonic cubic spline.
#
#   - ispline_non_decreasing(): Performs I-Spline monotonic regression on a 
#                               given sequence.
#
#   - constrained_least_squares(): A helper method using FISTA for solving 
#                                  constrained least squares problems.
#
#   - ispline_monotonic_interpolate(): Implements I-Spline based monotonic 
#                                     interpolation.
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
    
    def ispline_non_decreasing(self, raw_pep, min_value=0.0, max_value=1.0, max_iter=DEFAULT_MAX_ITER):
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

        Returns:
            list of floats:
                The fitted PEP values (clamped between min_value and max_value).
        """
        # Convert raw_pep to numpy array
        y = np.array(raw_pep)
        N = len(y)
        if N == 0:
            return []
        
        # Normalize x positions to [0, 1]
        x = np.linspace(0, 1, N)
        
        # Determine the number of basis functions.
        # Here we choose m as the integer square root of N (with an upper cap if desired).
        m = int(np.sqrt(N))
        # Create equally spaced knots in [0, 1]
        knots = np.linspace(0, 1, m + 1)
        # quantiles = np.linspace(0, 1, m + 1)
        # knots = np.quantile(x, quantiles)

        # Construct the cubic I-Spline basis matrix.
        # For each interval [knots[j], knots[j+1]], define:
        #    I_j(x) = 0                         if x < knots[j],
        #             3u^2 - 2u^3               if knots[j] <= x < knots[j+1],
        #             1                         if x >= knots[j+1],
        # where u = (x - knots[j]) / (knots[j+1] - knots[j]).
        B = np.zeros((N, m))
        for j in range(m):
            # Compute the normalized coordinate u.
            u = (x - knots[j]) / (knots[j+1] - knots[j])
            # Ensure u is in [0, 1].
            u = np.clip(u, 0, 1)
            B[:, j] = np.where(
                x < knots[j],
                0.0,
                np.where(
                    x >= knots[j+1],
                    1.0,
                    3 * u**2 - 2 * u**3
                )
            )
        # Add an intercept column to the design matrix.
        X = np.column_stack((np.ones(N), B))
        
        # Solve the constrained least squares problem:
        #   minimize ||Xc - y||^2
        # subject to: c[1:] >= 0 (the intercept c[0] is unconstrained).
        lower_bounds = np.concatenate(([-np.inf], np.zeros(m)))
        upper_bounds = np.full(m + 1, np.inf)
        res = lsq_linear(X, y, bounds=(lower_bounds, upper_bounds))
        c = res.x
        # c = self.constrained_least_squares(X, y, lower_bounds, upper_bounds, max_iter=max_iter)
        
        # Compute the fitted values and clamp them to [min_value, max_value]
        fitted = X.dot(c)
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
        I-Spline monotonic cubic interpolation.
        
        Steps:
        (a) Group the data into extended blocks (each point initially forms a block).
        (b) Merge adjacent blocks via PAVA (each merged block spans indices [startIdx, endIdx]).
        (c) For each final block, compute the x-center based on the specified center_method:
                - "mean": the average of x values in that block.
                - "median": the midpoint between the first and last x in the block, 
                            i.e., (x[startIdx] + x[endIdx]) / 2.
        (d) Collect these (xBlock, yBlock) points and build a monotonic cubic interpolant 
            using I-Spline method.
        (e) Evaluate this spline at each original x[i] and clamp the results to [min_y, max_y].
        
        Parameters:
            x             : list of floats
                            Sorted positions.
            y             : list of floats
                            Data values.
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
        if ip_algo=="ispline":
            y_interp = self.ispline_monotonic_interpolate(xBlock, yBlock, x_arr)
        elif ip_algo=="pchip":
            y_interp = self.fritsch_carlson_monotonic_interpolate(xBlock, yBlock, x_arr)
        else:
            raise ValueError("Unknown ip_algo. Use 'ispline' or 'pchip'.")

        # (e) Clamp the result to [min_y, max_y] and return as list.
        y_clamped = np.clip(y_interp, min_y, max_y)
        return y_clamped.tolist()

# =============================================================================
# TDCIsotonicPEP Class: Regression for Target-Decoy Competition Data
# -----------------------------------------------------------------------------
# Inherits from IsotonicRegression and implements methods to process binary 
# observations (target vs. decoy) and to compute posterior error probabilities (PEP)
# using isotonic regression.
# =============================================================================
class TDCIsotonicPEP(IsotonicRegression):
    def __init__(self, regression_algo="ispline", max_iter=DEFAULT_MAX_ITER):
        self.regression_algo = regression_algo
        self.max_iter = max_iter
    
    def process_obs(self, target=None, decoy=None, obs=None, target_label="target", decoy_label="decoy"):
        """
        Process observation data for d2pep.

        Acceptable inputs:
          - If obs is provided:
              * a numpy array (2D) with at least 2 columns.
              * a tuple; the two arrays must have the same length.
              * a DataFrame, it must contain the score and type columns.
          - Otherwise, if target and decoy are provided separately, they are used.
        
        Returns:
            A DataFrame with columns 'score' and 'label' (0 for target, 1 for decoy) and an "orig_idx" column preserving the original order.
            In the case of separate inputs, a "group" column is added ("target" or "decoy").
        """
        def convert_label(x, target_label, decoy_label):
            """
            Convert the type values in the observation data to numeric values
            based on the provided target_label and decoy_label.
            """
            try:
                # Convert x to float, then to int, then to string.
                x_str = str(int(float(x)))
            except Exception:
                # Fallback: use the original string representation (trimmed and lowercased)
                x_str = str(x).strip().lower()
            if x_str == str(target_label).strip().lower():
                return 0
            elif x_str == str(decoy_label).strip().lower():
                return 1
            else:
                raise ValueError("Invalid label: " + str(x))
            
        if obs is not None:
            if isinstance(obs, np.ndarray):
                if obs.ndim == 2 and obs.shape[1] >= 2:
                    df = pd.DataFrame(obs, columns=["score", "label"])
                else:
                    raise ValueError("The numpy array for obs must be 2D with at least 2 columns.")
            elif isinstance(obs, (list, tuple)):
                arr1 = np.array(obs[0])
                arr2 = np.array(obs[1])
                if len(arr1) != len(arr2):
                    raise ValueError("For concatenated input, the two arrays must have the same length.")
                df = pd.DataFrame({"score": arr1, "label": arr2})
            elif isinstance(obs, pd.DataFrame):
                df = obs.copy()
                df = df.rename(columns={df.columns[0]: "score", df.columns[1]: "label"})
            else:
                raise ValueError("obs must be a numpy array, tuple, or DataFrame.")
            df["label"] = df["label"].apply(lambda x: convert_label(x, target_label, decoy_label))
        else:
            if target is not None and decoy is not None:
                df_target = pd.DataFrame({"score": np.array(target).astype(float)})
                df_target["label"] = 0
                df_target["group"] = "target"
                df_decoy = pd.DataFrame({"score": np.array(decoy).astype(float)})
                df_decoy["label"] = 1
                df_decoy["group"] = "decoy"
                df = pd.concat([df_target, df_decoy], ignore_index=True)
            else:
                raise ValueError("For obs2pep, provide either obs or both target and decoy.")
            
        # Preserve original order.
        if "group" in df.columns:
            df["orig_idx"] = df.groupby("group").cumcount()
        else:
            df["orig_idx"] = np.arange(len(df))
        return df

    def tdc_binomial_regression(self, df_obs, regression_algo="ispline", max_iter=None):
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
            pandas.Series:
                A Series of target PEP values aligned with the original order.
        """
        if max_iter is None:
            max_iter = self.max_iter
        
        df_sorted = df_obs.sort_values(by="score", ascending=False, kind="mergesort").reset_index(drop=True)
        pseudo = pd.DataFrame({"score": [np.nan], "label": [0.5]})
        df_aug = pd.concat([pseudo, df_sorted], ignore_index=True)
        y_values = df_aug["label"].values
        if regression_algo == "PAVA":
            fitted = self.pava_non_decreasing(list(y_values), [1] * len(y_values))
        elif regression_algo == "ispline":
            fitted = self.ispline_non_decreasing(list(y_values), max_iter=max_iter)
        else:
            raise ValueError("Unknown regression_algo. Use 'PAVA' or 'ispline'.")
        
        fitted = np.array(fitted)
        fitted_decoy_prob = fitted[1:]  # remove pseudo observation
        with np.errstate(divide='ignore', invalid='ignore'):
            pep = fitted_decoy_prob / (1 - fitted_decoy_prob)
            pep = np.clip(pep, 0, 1)
        df_sorted["pep"] = pep
        # Restore original order.
        if "group" in df_obs.columns:
            df_sorted["orig_idx"] = df_sorted["orig_idx"].astype(int)
            df_target = df_sorted[df_sorted["group"]=="target"].sort_values(by="orig_idx", kind="mergesort")
            return df_target["pep"].reset_index(drop=True)
        else:
            df_result = df_sorted.sort_values(by="orig_idx", kind="mergesort").reset_index(drop=True)
            df_target = df_result[df_result["label"]==0]
            return df_target["pep"]

# =============================================================================
# IsotonicPEP Class: Unified Interface for PEP Estimation
# -----------------------------------------------------------------------------
# Provides a unified interface pep_regression() that accepts either:
#   - For q2pep: a Series/array/list of q-values.
#   - For d2pep: either a single DataFrame (or tuple/numpy array) containing score and type, 
#                  or separate target and decoy inputs.
#
# In both cases, the function returns a Series (or a tuple of two Series) of PEP values,
# aligned with the original order so that users can directly assign them as new columns.
# =============================================================================
class IsotonicPEP(TDCIsotonicPEP):
    def __init__(self, regression_algo="ispline", ip_algo="ispline", center_method="mean", max_iter=DEFAULT_MAX_ITER):
        super().__init__(regression_algo=regression_algo, max_iter=max_iter)
        self.center_method = center_method
        self.ip_algo = ip_algo
    
    def q_from_pep(self, pep_array):
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
        ranks = np.arange(1, len(pep_sorted)+1)
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
            center_method : str, optional
                             Method for computing block centers; defaults to self.center_method.
            max_iter      : int, optional (default=DEFAULT_MAX_ITER)
                             Maximum iterations.

        Returns:
            pandas.Series:
                Smoothed PEP values aligned with the original indices.
        """
        if center_method is None:
            center_method = self.center_method
        if ip_algo is None:
            ip_algo = self.ip_algo
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
            if i < n-1 and q_list[i] > q_list[i+1]:
                raise AssertionError("q_values must be non-decreasing")
        raw_pep = [qn[0]] + [qn[i]-qn[i-1] for i in range(1,n)]

        if regression_algo == "PAVA":
            if ip:
                x_positions = list(range(len(raw_pep)))
                final_pep = self.pava_non_decreasing_interpolation(x_positions, raw_pep, ip_algo=ip_algo, center_method=center_method)
            else:
                final_pep = self.pava_non_decreasing(raw_pep, [1] * len(raw_pep))
        elif regression_algo == "ispline":
            final_pep = self.ispline_non_decreasing(raw_pep, max_iter=max_iter)
        else:
            raise ValueError("Unknown regression_algo. Use 'PAVA' or 'ispline'.")
        
        pep_sorted = pd.Series(final_pep, index=q_series_sorted.index)
        pep_result = pep_sorted.reindex(orig_idx)
        return pep_result
    
    def dprob_to_pep(self, obs=None, target=None, decoy=None, target_label="target", decoy_label="decoy", regression_algo="ispline", max_iter=DEFAULT_MAX_ITER):
        """
        Compute PEP values from target-decoy observations (d2pep).
        
        Parameters:
            obs            : DataFrame/tuple/numpy array, optional
                             Contains score and label information.
            target, decoy  : array-like, optional
                             Separate target and decoy scores.
            target_label   : str, optional (default="target")
                             Label for targets.
            decoy_label    : str, optional (default="decoy")
                             Label for decoys.
            regression_algo: str, optional (default="ispline")
                             Regression method ("PAVA" or "ispline").
            max_iter       : int, optional (default=DEFAULT_MAX_ITER)
                             Maximum iterations.

        Returns:
            pandas.Series:
                Target PEP values aligned with the original order.
        """
        if target is not None and decoy is not None:
            df_obs = self.process_obs(target=target, decoy=decoy, obs=None, target_label=target_label, decoy_label=decoy_label)
        else:
            if obs is None:
                raise ValueError("For d2pep, provide either concatenated observations as obs or both target and decoy.")
            df_obs = self.process_obs(obs=obs, target_label=target_label, decoy_label=decoy_label)
        return self.tdc_binomial_regression(df_obs, regression_algo=regression_algo, max_iter=max_iter)

    def pep_regression(self, q_values=None, obs=None, target=None, decoy=None, target_label="target", decoy_label="decoy", 
                       method="q2pep", regression_algo="ispline", max_iter=DEFAULT_MAX_ITER, ip=False, ip_algo=None, center_method=None, calc_q=True):
        """
        Unified interface for computing PEP values,
        then optionally computing q-values from PEPs.
        
        For method "q2pep":
            - q_values: a Series/array/list of q-values,
            - calc_q: Boolean; if True, estimate q-values from calculated PEPs
        
        For method "d2pep":
            - Either provide data as a DataFrame (or tuple) containing score and type,
              or provide target and decoy separately,
                * target_label, decoy_label: Target and decoy labels in column "label" when using concatenated input.
            - calc_q: Boolean; if True, estimate q-values from calculated PEPs.
        
        Parameters:
            q_values       : array-like or pandas.Series, optional
                             Input q-values (for q2pep).
            obs            : DataFrame/tuple/numpy array, optional
                             Observation data for d2pep.
            target, decoy  : array-like, optional
                             Separate target and decoy scores.
            target_label   : str, optional (default="target")
                             Label for targets.
            decoy_label    : str, optional (default="decoy")
                             Label for decoys.
            method         : str, optional (default="q2pep")
                             Either "q2pep" or "d2pep".
            regression_algo: str, optional (default="ispline")
                             Regression method.
            max_iter       : int, optional (default=DEFAULT_MAX_ITER)
                             Maximum iterations.
            ip             : bool, optional (default=False)
                             If True and using PAVA for q2pep, apply interpolation.
            center_method  : str, optional
                             Method for computing block centers; defaults to self.center_method.
            calc_q         : bool, optional (default=True)
                             If True, also compute q-values from the estimated PEPs.

        Returns:
            Depending on the method and calc_q:
              - For q2pep: A tuple (pep_array, q_array) if calc_q is True; otherwise (pep_array, None).
              - For d2pep: Either a Series of PEP values or a tuple (pep_array, q_array).
        """
        if center_method is None:
            center_method = self.center_method
        if ip_algo is None:
            ip_algo = self.ip_algo

        if method == "q2pep":
            if q_values is None:
                raise ValueError("For q2pep, q-values must be provided.")
            pep_series = self.q_to_pep(q_values=q_values, regression_algo=regression_algo, ip=ip, ip_algo=ip_algo, center_method=center_method, max_iter=max_iter)
            pep_array = pep_series.values
            if not calc_q:
                return pep_array, None
            q_array = self.q_from_pep(pep_array)
            return pep_array, q_array
        elif method == "d2pep":
            target_pep = self.dprob_to_pep(obs=obs, target=target, decoy=decoy, target_label=target_label, decoy_label=decoy_label, regression_algo=regression_algo, max_iter=max_iter)
            pep_array = target_pep.values
            if not calc_q:
                return pep_array
            q_array = self.q_from_pep(pep_array)
            return pep_array, q_array
        else:
            raise ValueError("Unknown method. Use 'q2pep' or 'd2pep'.")