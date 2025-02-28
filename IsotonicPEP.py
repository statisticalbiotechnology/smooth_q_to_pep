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

    def pava_non_decreasing_interpolation(self, x, y, center_method="mean", min_y=0.0, max_y=1.0):
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
    
    def process_obs(self, target=None, decoy=None, obs=None, score_col="score", type_col="type"):

        """
        Process observation data for obs2pep.

        Acceptable inputs:
          - If obs is provided:
              * a numpy array (2D) with at least 2 columns.
              * a tuple; the two arrays must have the same length.
              * a DataFrame, it must contain the score and type columns.
          - Otherwise, if target and decoy are provided separately, they are used.
        
        Returns:
            A DataFrame with columns 'score' and 'type' (0 for target, 1 for decoy) and an "orig_idx" column preserving the original order.
            In the case of separate inputs, a "group" column is added ("target" or "decoy").
        """
        if obs is not None:
            if isinstance(obs, np.ndarray):
                if obs.ndim == 2 and obs.shape[1] >= 2:
                    df = pd.DataFrame(obs, columns=[score_col, type_col])
                else:
                    raise ValueError("The numpy array for obs must be 2D with at least 2 columns.")
            elif isinstance(obs, (list, tuple)):
                arr1 = np.array(obs[0])
                arr2 = np.array(obs[1])
                if len(arr1) != len(arr2):
                    raise ValueError("For concatenated input, the two arrays must have the same length.")
                df = pd.DataFrame({score_col: arr1, type_col: arr2})
            elif isinstance(obs, pd.DataFrame):
                df = obs.copy()
                df = df.rename(columns={df.columns[0]: score_col, df.columns[1]: type_col})
            else:
                raise ValueError("obs must be a numpy array, tuple, or DataFrame.")
        else:
            if target is not None and decoy is not None:
                df_target = pd.DataFrame({score_col: np.array(target).astype(float)})
                df_target[type_col] = 0
                df_target["group"] = "target"
                df_decoy = pd.DataFrame({score_col: np.array(decoy).astype(float)})
                df_decoy[type_col] = 1
                df_decoy["group"] = "decoy"
                df = pd.concat([df_target, df_decoy], ignore_index=True)
            else:
                raise ValueError("For obs2pep, provide either obs or both target and decoy.")
        # Convert type values to numeric.
        def convert_label(x):
            try:
                num = float(x)
                if num in [0, 1]:
                    return int(num)
            except:
                pass
            x_str = str(x).strip().lower()
            if x_str == "target":
                return 0
            elif x_str == "decoy":
                return 1
            else:
                raise ValueError("Invalid label: " + str(x))
        df[type_col] = df[type_col].apply(convert_label)
        # Preserve original order.
        if "group" in df.columns:
            df["orig_idx"] = df.groupby("group").cumcount()
        else:
            df["orig_idx"] = np.arange(len(df))
        return df

    def tdc_binomial_regression(self, df_obs, score_col="score", type_col="type"):
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
            If a "group" column exists (i.e. separate target and decoy inputs), returns a tuple of two Series:
                (target_pep, decoy_pep);
            otherwise, returns a Series of PEP values aligned with the original order.
        """
        df_sorted = df_obs.sort_values(by=score_col, ascending=False).reset_index(drop=True)
        pseudo = pd.DataFrame({score_col: [np.nan], type_col: [0.5]})
        df_aug = pd.concat([pseudo, df_sorted], ignore_index=True)
        y_values = df_aug[type_col].values
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
            df_target = df_sorted[df_sorted["group"]=="target"].sort_values(by="orig_idx")
            df_decoy = df_sorted[df_sorted["group"]=="decoy"].sort_values(by="orig_idx")
            return df_target["pep"].reset_index(drop=True), df_decoy["pep"].reset_index(drop=True)
        else:
            df_result = df_sorted.sort_values(by="orig_idx").reset_index(drop=True)
            return df_result["pep"]

# ----------------------------------------------------------------------
# Unified class for isotonic PEP estimation.
# Provides a unified interface pep_regression() that accepts either:
#   - For q2pep: a Series/array/list of q-values.
#   - For obs2pep: either a single DataFrame (or tuple/numpy array) containing score and type,
#                or separate target and decoy inputs.
#
# In both cases, the function returns a Series (or tuple of two Series) of PEP values,
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
        q_series_sorted = q_series.sort_values(ascending=True)
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
    
    def obs_to_pep(self, obs=None, target=None, decoy=None, score_col="score", type_col="type"):
        """
        Compute PEP values from target-decoy observations (obs2pep).
        
        Parameters:
            Either:
              - data: a DataFrame (or tuple or numpy array) containing score and type.
              - or target and decoy: separate inputs (e.g., arrays or Series) for target and decoy scores.
            score_col, type_col: Column names when data is a DataFrame.
        
        Returns:
            If separate inputs are provided, returns a tuple of two Series (target_pep, decoy_pep);
            otherwise, returns a Series of PEP values aligned with the original order.
        """
        if target is not None and decoy is not None:
            df_obs = self.process_obs(target=target, decoy=decoy, obs=None, score_col=score_col, type_col=type_col)
        else:
            if obs is None:
                raise ValueError("For obs2pep, provide either concatenated observations as obs or both target and decoy.")
            df_obs = self.process_obs(obs=obs, score_col=score_col, type_col=type_col)
        return self.tdc_binomial_regression(df_obs, score_col, type_col)

    def pep_regression(self, q_values=None, obs=None, target=None, decoy=None, method="q2pep",
                    score_col="score", type_col="type", smooth=False, pava_method=None, center_method=None):
        """
        Unified interface for computing PEP values.
        
        For method "q2pep":
            - q_values: a Series/array/list of q-values.
            - smooth: Boolean; if True, apply block-merge pre-processing.
        
        For method "obs2pep":
            - Either provide data as a DataFrame (or tuple) containing score and type,
              or provide target and decoy separately.
        
        Returns:
            For q2pep: a Series of PEP values aligned with the input index.
            For obs2pep: if separate inputs are provided, a tuple (target_pep, decoy_pep);
                       otherwise, a Series of PEP values aligned with the input index.
        """
        if method == "q2pep":
            if q_values is None:
                raise ValueError("For q2pep, q-values must be provided.")
            return self.q_to_pep(q_values=q_values, smooth=smooth, pava_method=pava_method, center_method=center_method)
        elif method == "obs2pep":
            return self.obs_to_pep(obs=obs, target=target, decoy=decoy, score_col=score_col, type_col=type_col)
        else:
            raise ValueError("Unknown method. Use 'q2pep' or 'obs2pep'.")