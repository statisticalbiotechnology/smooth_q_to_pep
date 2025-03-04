import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# Base class for isotonic regression in real space.
# Provides:
#   - pava_non_decreasing(): Stepwise constant solution.
#   - pava_non_decreasing_interpolation(): Linear interpolation between block centers.
# ----------------------------------------------------------------------
class IsotonicRegression:
    def __init__(self):
        pass

    def pava_non_decreasing(self, values, counts, min_value=0.0, max_value=1.0):
        """
        Perform standard PAVA (Pool-Adjacent-Violators Algorithm) to enforce a 
        non-decreasing sequence.
        
        Each input value is repeated according to its corresponding count.
        After merging, each block's average is clamped to the interval [min_value, max_value].

        
        Parameters:
            values: list of floats, the data values.
            counts: list of ints, the number of times each value is repeated.
            min_value: lower bound for clamping (default 0.0).
            max_value: upper bound for clamping (default 1.0).

        Returns:
            A list (expanded to length = sum(counts)) where each block's average 
            is repeated count times.
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

        # Expand the final solution: repeat each block's average according to its count.
        result = []
        for block in stack:
            # Clamp the block average to [min_value, max_value]
            clamped_avg = min(max(block['avg'], min_value), max_value)
            result.extend([clamped_avg] * block['count'])
        return result

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

    def pava_non_decreasing_interpolation(self, x, y, center_method="mean", min_y=0.0, max_y=1.0):
        """
        'Interpolated' variant of non-decreasing PAVA using Fritsch-Carlson monotonic cubic interpolation.
        
        Instead of returning a stepwise-constant output or simple linear interpolation,
        this version uses Fritsch-Carlson method to achieve a smooth, monotonic cubic spline
        across the block centers.
        
        Each interpolated value is clamped to the interval [min_y, max_y].
        
        Steps:
        (a) Group the data into extended blocks (each point initially forms a block).
        (b) Merge adjacent blocks via PAVA (each merged block spans indices [startIdx, endIdx]).
        (c) For each final block, compute the x-center based on the specified center_method:
                - "mean": the average of x values in that block.
                - "median": the midpoint between the first and last x in the block, 
                            i.e., (x[startIdx] + x[endIdx]) / 2.
        (d) Collect these (xBlock, yBlock) points and build a monotonic cubic spline 
            using Fritsch-Carlson method.
        (e) Evaluate this spline at each original x[i] and clamp the results to [min_y, max_y].
        
        Parameters:
            x: list of floats, sorted positions.
            y: list of floats, the data values.
            center_method: "mean" or "median" to determine block center x.
            min_y: lower bound for clamping (default 0.0).
            max_y: upper bound for clamping (default 1.0).

        Returns:
            A list of floats representing the interpolated, non-decreasing values.
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

        # (d) Use Fritsch-Carlson monotonic cubic interpolation on (xBlock, yBlock)
        #     and evaluate at each original x[i].
        y_interp = self.fritsch_carlson_monotonic_interpolate(xBlock, yBlock, x_arr)

        # (e) Clamp the result to [min_y, max_y] and return as list.
        y_clamped = np.clip(y_interp, min_y, max_y)
        return y_clamped.tolist()

    def pava_non_decreasing_interpolation_linear(self, x, y, center_method="mean", min_y=0.0, max_y=1.0):
        """
        "Interpolated" variant of non-decreasing PAVA.
        
        Instead of returning a stepwise-constant output, this returns a piecewise-linear 
        function across block centers.
        Each interpolated value is clamped to the interval [min_y, max_y].
        
        Steps:
          (a) Group the data into extended blocks (each point initially forms a block).
          (b) Merge adjacent blocks via PAVA (each merged block spans indices [startIdx, endIdx]).
          (c) For each final block, compute the x-center based on the specified center_method:
                - "mean": the average of x values in that block.
                - "median": the midpoint between the first and last x in the block, i.e., (x[startIdx] + x[endIdx]) / 2.
          (d) For each point in the block, linearly interpolate between the block's center and the next block's center.
          (e) Clamp each resulting value to [min_y, max_y].
        
        Parameters:
            x: list of floats, sorted positions.
            y: list of floats, the data values.
            min_y: lower bound for clamping (default 0.0).
            max_y: upper bound for clamping (default 1.0).

        Returns:
            A list of floats representing the interpolated, non-decreasing values.
        """
        if len(x) != len(y):
            raise ValueError("x and y must have the same length.")
        n = len(y)
        if n == 0:
            return []

        # (a) Build an extended block for each point.
        # Each block is a dict with keys: 'sum', 'count', 'avg', 'startIdx', 'endIdx'
        blocks = []
        for i in range(n):
            block = {'sum': y[i], 'count': 1, 'avg': y[i], 'startIdx': i, 'endIdx': i}
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
                    stack.append({'sum': new_sum, 'count': new_count, 'avg': new_avg, 
                                  'startIdx': sIdx, 'endIdx': eIdx})
                else:
                    break

        # (c) For each final block, compute the x-center and store it in the 'sum' field.
        for block in stack:
            if center_method == "mean":
                sum_x = sum(x[i] for i in range(block['startIdx'], block['endIdx'] + 1))
                length = block['endIdx'] - block['startIdx'] + 1
                center = sum_x / length
            elif center_method == "median":
                # Compute the median: (first x + last x) / 2.
                center = (x[block['startIdx']] + x[block['endIdx']]) / 2
            else:
                raise ValueError("Unknown center_method. Use 'mean' or 'median'.")
            # Store the computed center in the 'sum' field.
            block['sum'] = center

        # (d) Interpolate each final block's points.
        result = [0.0] * n
        for b in range(len(stack)):
            curBlk = stack[b]
            curAvg = curBlk['avg']
            curXc = curBlk['sum']  # current block x-center

            # Determine next block's average and x-center if it exists.
            if b < len(stack) - 1:
                nextBlk = stack[b+1]
                nextAvg = nextBlk['avg']
                nextXc = nextBlk['sum']
            else:
                nextAvg = curAvg
                nextXc = curXc

            for i in range(curBlk['startIdx'], curBlk['endIdx'] + 1):
                # For the last block or if the x-center difference is almost zero, use constant value.
                if b == len(stack)-1 or abs(nextXc - curXc) < 1e-15:
                    result[i] = curAvg
                else:
                    t = (x[i] - curXc) / (nextXc - curXc)
                    t = max(0.0, min(1.0, t))
                    result[i] = curAvg + t * (nextAvg - curAvg)
                # (e) Clamp the result to [min_y, max_y]
                result[i] = min(max(result[i], min_y), max_y)
        return result

# ----------------------------------------------------------------------
# Inherits from the IsotonicRegression base class to perform binomial regression on
# a stream of binary observations yielded by target-decoy competition (TDC) method.
# ----------------------------------------------------------------------
class TDCIsotonicPEP(IsotonicRegression):
    def __init__(self):
        pass
    
    def process_obs(self, target=None, decoy=None, obs=None, target_label="target", decoy_label="decoy"):
        """
        Process observation data for obs2pep.

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

    def tdc_binomial_regression(self, df_obs):
        """
        Compute PEP values from observation data using isotonic regression.
        
        Steps:
          - (a) Sort by score descending.
          - (b) Prepend a pseudo observation with type 0.5.
          - (c) Apply PAVA on the binary sequence.
          - (d) Remove the pseudo observation.
          - (e) Compute PEP = decoy_prob/(1 - decoy_prob), clipped to [0,1].
          - (f) Restore original order.
        
        Returns:
            A Series of target PEP values aligned with the original order.
        """
        df_sorted = df_obs.sort_values(by="score", ascending=False, kind="mergesort").reset_index(drop=True)
        pseudo = pd.DataFrame({"score": [np.nan], "label": [0.5]})
        df_aug = pd.concat([pseudo, df_sorted], ignore_index=True)
        y_values = df_aug["label"].values
        fitted = self.pava_non_decreasing(list(y_values), [1] * len(y_values))
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

# ----------------------------------------------------------------------
# Unified class for isotonic PEP estimation.
# Provides a unified interface pep_regression() that accepts either:
#   - For q2pep: a Series/array/list of q-values.
#   - For obs2pep: either a single DataFrame (or tuple/numpy array) containing score and type, 
#                  or separate target and decoy inputs.
#
# In both cases, the function returns a Series (or a tuple of two Series) of PEP values,
# aligned with the original order so that users can directly assign them as new columns.
# ----------------------------------------------------------------------
class IsotonicPEP(TDCIsotonicPEP):
    def __init__(self, pava_method="basic", center_method="mean"):
        self.pava_method = pava_method
        self.center_method = center_method

    def create_blocks_in_unit_interval_and_unfold(self, y):
        """
        Merge consecutive y[i] so that each block's average is in [0,1] and then unfold the blocks.
        
        For each value in y, accumulate the sum and count. Compute the current average.
        If the average is between 0 and 1 (inclusive), finalize the current block by creating
        a list of length equal to the current count with every element equal to the average.
        Then, reset the accumulator and continue. After processing all values, if there is a leftover,
        clip its average to be within [0,1] and create a block.
        Finally, unfold (concatenate) all blocks into one list.
        
        Parameters:
            y: list of floats
            
        Returns:
            result: a flattened list after merging (unfolded blocks)
        """
        blocks = [] # This will store the finalized blocks (each block is a list)
        current_sum = 0.0
        current_count = 0

        for val in y:
            current_sum += val
            current_count += 1
            avg = current_sum / current_count
            if avg >= 0.0 and avg <= 1.0:
                # Finalize this block: create a block of length current_count filled with avg.
                block = [avg] * current_count
                blocks.append(block)
                # Reset the accumulator for the next block.
                current_sum = 0.0
                current_count = 0
        # Process any leftover values.
        if current_count > 0:
            avg = current_sum / current_count
            # If the average is out of [0,1], clip it to a near-boundary value.
            if avg < 0.0:
                avg = 0
            if avg > 1.0:
                avg = 1.0
            block = [avg] * current_count
            blocks.append(block)
        # "Unfold" all blocks into a single list by concatenating them.
        result = []
        for block in blocks:
            result.extend(block)
        return result
    
    def q_from_pep(self, pep_array):
        """
        Given a PEP array, compute q-values via:
           q(i) = (1 / i) * cumsum(pep up to i),
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

    def q_to_pep(self, q_values, smooth=False, pava_method=None, center_method=None):
        """
        Compute smoothed PEP values from q-values (q2pep).

        qn[i] = q_values[i] * (i+1)
        raw_pep[0] = qn[0],
        for i >= 1: raw_pep[i] = qn[i] - qn[i-1].
        
        Parameters:
            q_values: a Series/array/list of q-values.
            smooth: Boolean; if True, apply block-merge pre-processing.
            pava_method: "basic" or "ip". Defaults to self.pava_method.
            center_method: "mean" or "median". Defaults to self.center_method.
        
        Returns:
            A Series of PEP values aligned with the original index.
        """
        if pava_method is None:
            pava_method = self.pava_method
        if center_method is None:
            center_method = self.center_method
        
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

        if smooth:
            processed = self.create_blocks_in_unit_interval_and_unfold(raw_pep)
        else:
            processed = raw_pep
        
        if pava_method == "basic":
            final_pep = self.pava_non_decreasing(processed, [1]*len(processed))
        elif pava_method == "ip":
            x_positions = list(range(len(processed)))
            final_pep = self.pava_non_decreasing_interpolation(x_positions, processed, center_method=center_method)
        else:
            raise ValueError("Unknown pava_method. Use 'basic' or 'ip'.")
        
        pep_sorted = pd.Series(final_pep, index=q_series_sorted.index)
        pep_result = pep_sorted.reindex(orig_idx)
        return pep_result
    
    def obs_to_pep(self, obs=None, target=None, decoy=None, target_label="target", decoy_label="decoy"):
        """
        Compute PEP values from target-decoy observations (obs2pep).
        
        Parameters:
            Either:
              - obs: a DataFrame (or tuple or numpy array) containing score and type.
              - or target and decoy: separate inputs (e.g., arrays or Series) for target and decoy scores.
            target_label, decoy_label: Target and decoy labels in column "label" when using concatenated input.
        
        Returns:
            a Series of target PEP values aligned with the original order.
        """
        if target is not None and decoy is not None:
            df_obs = self.process_obs(target=target, decoy=decoy, obs=None, target_label=target_label, decoy_label=decoy_label)
        else:
            if obs is None:
                raise ValueError("For obs2pep, provide either concatenated observations as obs or both target and decoy.")
            df_obs = self.process_obs(obs=obs, target_label=target_label, decoy_label=decoy_label)
        return self.tdc_binomial_regression(df_obs)

    def pep_regression(self, q_values=None, obs=None, target=None, decoy=None, method="q2pep", target_label="target",
                       decoy_label="decoy", smooth=False, pava_method=None, center_method=None, calc_q=True):
        """
        Unified interface for computing PEP values,
        then optionally computing q-values from PEPs.
        
        For method "q2pep":
            - q_values: a Series/array/list of q-values,
            - calc_q: Boolean; if True, estimate q-values from calculated PEPs,
            - smooth: Boolean; if True, apply block-merge pre-processing.
        
        For method "obs2pep":
            - Either provide data as a DataFrame (or tuple) containing score and type,
              or provide target and decoy separately,
                * target_label, decoy_label: Target and decoy labels in column "label" when using concatenated input.
            - calc_q: Boolean; if True, estimate q-values from calculated PEPs.
        
        Returns:
            an array of PEP values aligned with the input index,
            if calc_q is True, a tuple of two arrays (pep_array, q_array).

        By default, calc_q=True.
        """
        if method == "q2pep":
            if q_values is None:
                raise ValueError("For q2pep, q-values must be provided.")
            pep_series = self.q_to_pep(q_values=q_values, smooth=smooth, pava_method=pava_method, center_method=center_method)
            pep_array = pep_series.values
            if not calc_q:
                return pep_array, None
            q_array = self.q_from_pep(pep_array)
            return pep_array, q_array
        elif method == "obs2pep":
            target_pep = self.obs_to_pep(obs=obs, target=target, decoy=decoy, target_label=target_label, decoy_label=decoy_label)
            pep_array = target_pep.values
            if not calc_q:
                return pep_array
            q_array = self.q_from_pep(pep_array)
            return pep_array, q_array
        else:
            raise ValueError("Unknown method. Use 'q2pep' or 'obs2pep'.")