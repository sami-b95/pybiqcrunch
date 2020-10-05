import os
import subprocess
import tempfile
from .binary_encoder import BinaryEncoder


class BiqCrunchOptimizer:
    """
    A wrapper to solve quadratic constrained binary optimization problems using BiqCrunch
    (based on a branch-and-bound algorithm).
    """
    def __init__(self, variables, objective, equality_constraints, inequality_constraints, python_bin="python3.7"):
        """
        Create a BiqCrunch optimizer instance.

        Parameters
        ----------
        variables: dict
            A dictionary describing the variables of the cost function. Each
            (key, value) pair corresponds to a variable; the key gives the name
            of the variable and the value its description. This description
            is a dictionary with 3 keys: "min", "max", "bits". "min" and "max"
            stand for the minimum and maximum values the discretized version
            of the variable can take; "bits" is a list containing the indices
            of the bits in which the variable is to be encoded.
        objective: dict
            A dictionary describing the polynomial objective function to minimize. A
            monomial $c\cdot x_1\dots x_k$ (with $c$ a scalar) is represented by the
            (key, value) pair (("x_1", ..., "x_k"), c).
        equality_constraints: list
            A list of dictionaries describing the equality constraints. Each equality
            constraint is expressed as a polynomial $P$ (represented in the same format
            as the 'objective' argument) such that the constraint reads $P = 0$.
        inequality_constraints: list
            A list of dictionaries describing the inequality constraints. Each inequality
            constraint is expressed as a polynomial $P$ (represented in the same format
            as the 'objective' argument) such that the constraint reads $P \leq 0$.
        python_bin: str
            An optional argument to specify the python binary with which to execute the
            BiqCrunch python scripts.
        """
        self.variables = variables
        self.check_quadratic(objective)
        self.objective = objective
        for inequality_constraint in inequality_constraints:
            self.check_quadratic(inequality_constraint)
        self.inequality_constraints = inequality_constraints
        for equality_constraint in equality_constraints:
            self.check_quadratic(equality_constraint)
        self.equality_constraints = equality_constraints
        self.biqcrunch_dir = os.path.join(os.path.dirname(__file__), "biqcrunch/BiqCrunch_second_release/BiqCrunch")
        self.python_bin = python_bin

    @classmethod
    def check_quadratic(cls, couplings):
        for indices, coupling in couplings.items():
            if len(indices) > 2:
                raise ValueError("BiqCrunch solver can only handle objectives/constraints of degree at most 2")

    def couplings_to_string(self, couplings, coupling_rescaling=None):
        binary_encoder = BinaryEncoder(self.variables, couplings)
        bits_couplings = binary_encoder.to_binary()
        string = " "
        constant = 0
        for indices, coupling in bits_couplings.items():
            if coupling_rescaling:
                coupling = int(coupling * coupling_rescaling)
            if len(indices) == 0:
                constant = coupling
                continue
            if coupling == 0:
                continue
            elif coupling > 0:
                string += f" - {coupling} "
            else:
                string += f" + {-coupling} "
            string += "*".join([f"x{index}" for index in indices])
        return string, constant

    def optimize(self, objective_multiplier):
        """
        Run BiqCrunch on the optimization problem.

        Parameters
        ----------
        objective_multiplier: float
            A number to multiply the coefficients of the objective function with before
            rounding them (BiqCrunch only handles integer coefficients indeed).

        Returns
        -------
        tuple
            A pair whose first element is the achieved objective (on the objective with
            rounded coefficients, see above) and second element is a dictionary mapping
            each variable to its optimal value.
        """
        binary_encoder = BinaryEncoder(self.variables, self.objective)
        num_bits = binary_encoder.num_spins()
        with tempfile.TemporaryDirectory() as tmp_dir_path:
            # Generate LP file
            lp_filename = os.path.join(tmp_dir_path, "example.lp")
            with open(lp_filename, "w") as lp_file:
                # Objective
                lp_file.write("Maximize\n")
                objective_string, objective_constant = self.couplings_to_string(self.objective, coupling_rescaling=objective_multiplier)
                lp_file.write(objective_string)
                lp_file.write("\n\nSubject to\n")
                # Equality constraints
                for equality_constraint in self.equality_constraints:
                    equality_constraint_string, equality_constraint_constant = self.couplings_to_string(equality_constraint)
                    lp_file.write(equality_constraint_string + f" = {-equality_constraint_constant}\n")
                lp_file.write("\n")
                # Inequality constraints
                for inequality_constraint in self.inequality_constraints:
                    inequality_constraint_string, inequality_constraint_constant = self.couplings_to_string(inequality_constraint)
                    lp_file.write(inequality_constraint_string + f" <= {-inequality_constraint_constant}\n")
                lp_file.write(f"\nBinary\n  {' '.join([f'x{index}' for index in range(num_bits)])}\n\nEnd")
            # Convert LP to BC and run BiqCrunch
            bc_filename = os.path.join(tmp_dir_path, "example.bc")
            with open(bc_filename, "wb") as bc_file:
                subprocess.call([self.python_bin, os.path.join(self.biqcrunch_dir, "tools", "lp2bc.py"), lp_filename], stdout=bc_file)
            # Run BiqCrunch
            output = subprocess.check_output([os.path.join(self.biqcrunch_dir, "problems", "generic", "biqcrunch"), bc_filename, os.path.join(self.biqcrunch_dir, "biq_crunch.param")])
            # Process result
            if "Problem is infeasible" in str(output):
                raise Exception("Problem is infeasible.")
            bits = {bit: 0 for bit in range(num_bits)}
            ###print(output)
            for bit_idx in [int(bit_idx_str) for bit_idx_str in str(output).split("Solution = ")[1].split("\\n")[0][2:-2].split(" ") if len(bit_idx_str)]:
                bits[bit_idx - 1] = 1
            ###print("Maximum value: {}".format(int(str(output).split("Maximum value = ")[1].split("\\n")[0])))
            return (-int(str(output).split("Maximum value = ")[1].split("\\n")[0]) + objective_constant) / objective_multiplier, binary_encoder.bits_to_variables(bits)
