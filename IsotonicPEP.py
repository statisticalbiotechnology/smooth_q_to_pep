import math

# ----------------------------------------------------------------------
# Base class for isotonic regression in real space (no logistic transformations).
# Provides:
#   - pava_non_decreasing(): Stepwise constant solution.
#   - pava_non_decreasing_interpolation(): Linear interpolation between block centers.
# ----------------------------------------------------------------------
class IsotonicRegression:
    def __init__(self):
        pass

    def pava_non_decreasing(self, values, counts):
        """
        Perform standard PAVA (Pool-Adjacent-Violators Algorithm) to enforce a 
        non-decreasing sequence.
        
        Each input value is repeated according to its corresponding count.
        
        Parameters:
            values: list of floats, the data values.
            counts: list of ints, the number of times each value is repeated.
            
        Returns:
            A list (expanded to length = sum(counts)) where each block's average 
            is repeated count times.
        """
        if len(values) != len(counts):
            raise ValueError("values and counts must have the same length")
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
            result.extend([block['avg']] * block['count'])
        return result

    def pava_non_decreasing_interpolation(self, x, y):
        """
        "Interpolated" variant of non-decreasing PAVA.
        
        Instead of returning a stepwise-constant output, this returns a piecewise-linear 
        function across block centers.
        
        Steps:
          (a) Group the data into extended blocks (each point initially forms a block).
          (b) Merge adjacent blocks via PAVA (each merged block spans indices [startIdx, endIdx]).
          (c) For each final block, compute the x-center (e.g., the average of x values in that block).
          (d) For each point in the block, linearly interpolate between the block’s center and the next block’s center.
        
        Parameters:
            x: list of floats, sorted positions.
            y: list of floats, the data values.
            
        Returns:
            A list of floats representing the interpolated, non-decreasing values.
        """
        if len(x) != len(y):
            raise ValueError("x and y must have the same length")
        n = len(y)
        if n == 0:
            return []

        # (1) Build an extended block for each point.
        # Each block is a dict with keys: 'sum', 'count', 'avg', 'startIdx', 'endIdx'
        blocks = []
        for i in range(n):
            block = {'sum': y[i], 'count': 1, 'avg': y[i], 'startIdx': i, 'endIdx': i}
            blocks.append(block)

        # (2) Merge blocks using standard PAVA.
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

        # (3) For each final block, compute the x-center and store it in the 'sum' field.
        for block in stack:
            sum_x = sum(x[i] for i in range(block['startIdx'], block['endIdx'] + 1))
            length = block['endIdx'] - block['startIdx'] + 1
            center = sum_x / length
            block['sum'] = center

        # (4) Interpolate each final block's points.
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
        return result

# ----------------------------------------------------------------------
# Derived class that applies logistic transformations around the base PAVA routines.
# Provides:
#   - clamp_probability(), logistic(), logit()
#   - logistic_isotonic_regression(): stepwise-constant in logit space.
#   - logistic_isotonic_interpolation(): piecewise-linear in logit space.
# ----------------------------------------------------------------------
class LogisticIsotonicRegression(IsotonicRegression):
    @staticmethod
    def clamp_probability(p, eps=1e-12):
        """
        Clamp probability p to the interval [eps, 1-eps] to avoid numerical issues.
        """
        if p < eps:
            return eps
        if p > 1.0 - eps:
            return 1.0 - eps
        return p

    @staticmethod
    def logistic(x):
        """
        Logistic (sigmoid) function: 1/(1+exp(-x)).
        """
        return 1.0 / (1.0 + math.exp(-x))

    @staticmethod
    def logit(p):
        """
        Logit function: log(p/(1-p)). Uses clamping to ensure p is within (0,1).
        """
        p = LogisticIsotonicRegression.clamp_probability(p)
        return math.log(p/(1.0-p))

    def logistic_isotonic_regression(self, y, counts=None):
        """
        Perform logistic isotonic regression:
         1. Convert each y value to logit space.
         2. Run PAVA (non-decreasing) in logit space.
         3. Convert the result back to probability space using logistic().
         
        Parameters:
            y: list of floats, input values (which might be outside [0,1]).
            counts: optional list of ints; if omitted, each point is assumed to have count=1.
            
        Returns:
            A list of calibrated probability values.
        """
        n = len(y)
        if n == 0:
            return []
        if counts is None or len(counts) == 0:
            counts = [1] * n

        # Convert y to logit space.
        logit_vals = [self.logit(self.clamp_probability(val)) for val in y]

        # Perform PAVA in logit space.
        merged_logits = self.pava_non_decreasing(logit_vals, counts)

        # Convert back to probability space.
        result = [self.logistic(val) for val in merged_logits]
        return result

    def logistic_isotonic_interpolation(self, x, y):
        """
        Perform logistic isotonic interpolation:
         1. Convert y values to logit space.
         2. Apply the interpolated PAVA in logit space.
         3. Convert the result back to probability space.
         
        Parameters:
            x: list of floats, positions.
            y: list of floats, input values.
            
        Returns:
            A list of calibrated probability values.
        """
        n = len(y)
        if n == 0:
            return []
        if len(x) != n:
            raise ValueError("x and y must have the same length")

        logit_vals = [self.logit(self.clamp_probability(val)) for val in y]
        logit_interp = self.pava_non_decreasing_interpolation(x, logit_vals)
        result = [self.logistic(val) for val in logit_interp]
        return result

# ----------------------------------------------------------------------
# Derived class that includes PEP pre-processing.
# The q_to_pep function computes raw PEPs from q_values as follows:
#   - Compute qn values: qn[i] = q_values[i] * (i+1)
#   - Compute raw_pep: raw_pep[0] = qn[0], raw_pep[i] = qn[i] - qn[i-1] for i>=1
#
# Then, raw_pep is processed into [0,1] using one of two methods:
#   1. "clip": Directly clip values to [0,1].
#   2. "block_merge": Merge consecutive y[i] so that each partial block average is in [0,1]. 
#       Then "unfold" back into a single vector.
#
# Next, we choose the processing space:
#   - "real": work directly in real space.
#   - "logit": convert data into logit space, process, then convert back via logistic().
#
# Finally, we choose which PAVA version to use:
#   - "basic": basic PAVA (stepwise constant).
#   - "interp": interpolated PAVA (piecewise-linear).
#
# The final output is a set of smoothed PEP values.
# ----------------------------------------------------------------------
class IsotonicPEP(LogisticIsotonicRegression):
    def __init__(self):
        self.qs = []
        self.pep_iso = []

    def create_blocks_in_unit_interval_and_unfold(self, y):
        """
        Merge consecutive y[i] values so that each block's average is in [0,1],
        then flatten the blocks into a single list.
        
        The approach:
          - Accumulate values until the average falls in (0,1), then finalize the block.
          - For any leftover block, clip its average to [0,1] and fill.
        
        Parameters:
            y: list of floats.
            
        Returns:
            A flattened list after merging.
        """
        blocks = []
        current_sum = 0.0
        current_count = 0
        current_block = []
        for val in y:
            current_sum += val
            current_count += 1
            current_block.append(val)
            avg = current_sum / current_count
            if avg > 0.0 and avg < 1.0:
                blocks.append(current_block.copy())
                current_sum = 0.0
                current_count = 0
                current_block = []
        if current_block:
            avg = current_sum / current_count
            avg = max(0.0, min(1.0, avg))
            clipped_block = [avg] * current_count
            blocks.append(clipped_block)
        result = []
        for block in blocks:
            result.extend(block)
        return result

    def q_to_pep(self, q_values):
        """
        Convert q_values to raw PEP values and then process them.
        
        This function:
          1. Computes qn values: qn[i] = q_values[i]*(i+1)
          2. Computes raw_pep: raw_pep[0] = qn[0], for i>=1: raw_pep[i] = qn[i] - qn[i-1]
          3. Returns raw_pep; further processing is done in pep_regression.
        
        Parameters:
            q_values: list of floats (must be non-decreasing).
            
        Returns:
            A list of raw PEP values.
        """
        self.qs = q_values[:]  # store a copy of q_values
        n = len(q_values)
        qn = []
        for i in range(n):
            qn.append(q_values[i] * (i + 1))
            if i < n-1 and q_values[i] > q_values[i+1]:
                raise AssertionError("q_values must be non-decreasing")
        raw_pep = [qn[0]] + [qn[i] - qn[i-1] for i in range(1, n)]
        return raw_pep

    def pep_regression(self, q_values, raw_process_method="block_merge", space="real", pava_method="interp"):
        """
        Compute smoothed PEP values from q_values with different options.
        
        Steps:
          1. Calculate qn values and compute raw_pep from q_values.
          2. Process raw_pep into [0,1]:
                 - "clip": Directly clip values to [0,1].
                 - "block_merge": Merge consecutive points so that each block's average is in [0,1].
          3. Choose processing space and apply PAVA:
                 - "real": Work in real space using IsotonicRegression functions.
                    * "basic": Use basic PAVA (pava_non_decreasing).
                    * "interp": Use interpolated PAVA (pava_non_decreasing_interpolation).
                 - "logit": Work in logit space using LogisticIsotonicRegression functions.
                    * "basic": Use logistic basic PAVA (logistic_isotonic_regression).
                    * "interp": Use logistic interpolated PAVA (logistic_isotonic_interpolation).
        
        Parameters:
            q_values: list of floats, the q_values (assumed non-decreasing).
            raw_process_method: "clip" or "block_merge" (default "block_merge").
            space: "real" or "logit" (default "real").
            pava_method: "basic" or "interp" (default "interp").
            
        Returns:
            A list of final smoothed PEP values.
        """
        # Step 1: Compute raw_pep from q_values.
        raw_pep = self.q_to_pep(q_values)
        
        # Step 2: Process raw_pep into [0,1].
        if raw_process_method == "clip":
            processed = [max(0.0, min(1.0, val)) for val in raw_pep]
        elif raw_process_method == "block_merge":
            processed = self.create_blocks_in_unit_interval_and_unfold(raw_pep)
        else:
            raise ValueError("Unknown raw_process_method. Use 'clip' or 'block_merge'.")

        # Step 3: Choose processing space and apply PAVA.
        if space == "real":
            # Use IsotonicRegression functions in real space.
            if pava_method == "basic":
                pep_after_pava = self.pava_non_decreasing(processed, [1] * len(processed))
            elif pava_method == "interp":
                x_positions = list(range(len(processed)))
                pep_after_pava = self.pava_non_decreasing_interpolation(x_positions, processed)
            else:
                raise ValueError("Unknown pava_method. Use 'basic' or 'interp'.")
        elif space == "logit":
            # Use LogisticIsotonicRegression functions in logit space.
            if pava_method == "basic":
                pep_after_pava = self.logistic_isotonic_regression(processed)
            elif pava_method == "interp":
                x_positions = list(range(len(processed)))
                pep_after_pava = self.logistic_isotonic_interpolation(x_positions, processed)
            else:
                raise ValueError("Unknown pava_method. Use 'basic' or 'interp'.")
        else:
            raise ValueError("Unknown space. Use 'real' or 'logit'.")
        
        final_pep = pep_after_pava
        return final_pep
