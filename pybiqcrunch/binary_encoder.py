from copy import deepcopy
from functools import reduce
from itertools import product, combinations
from numpy import round


class BinaryEncoder:
    """
    A binary encoder allows to encode a cost function as
    a polynomial in binary variables, which may either be
    bits (0/1-valued) or spins (-1/+1v+___.
    """

    def __init__(self, variables, couplings):
        r"""
        Create a binary encoder. The variables of the function being
        discretized are encoded as bitstrings. The cost function has to be
        polynomial. Besides, it allows to decode the latter bitstrings to their
        corresponding variables.

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
        couplings: dict
            A dictionary describing the monomial terms in the cost function.
            A term $c\cdot x_1\dots x_k$ (with $c$ a scalar) is represented by the
            (key, value) pair (("x_1", ..., "x_k"), c).
        """

        self.variables = variables
        self.couplings = deepcopy(couplings)

    def update_couplings(self, new_couplings):
        """Add monomial terms to the cost function.

        Parameters
        ----------
        new_couplings: dict
            Similar to the class constructor, a dictionary describing the
            monomial terms to add to the cost function.
        """

        for index, value in new_couplings.items():
            if index not in self.couplings:
                self.couplings[index] = 0

            self.couplings[index] += value

    def num_spins(self):
        return sum(len(variable["bits"]) for variable in self.variables.values())

    def to_binary(self):
        r"""Generate the polynomial in binary operators encoding the cost function.

        Returns
        -------
        dict
            A dictionary describing the monomials of the polynomial in bits
            that encodes the cost function. Precisely, a monomial
            $c\cdot b_{i_1}\cdots b_{i_k}$ is represented by the (key, value) pair
            ((i_1, ..., i_k), c).
        """
        bit_couplings = {}

        # Iterate over monomial terms.
        for index, value in self.couplings.items():

            # The monomial term will be represented by a linear combination of monomials of bits.
            # Each monomial is (by definition) the product of a certain number of bits; iterate
            # over this number.
            for num_bit_factors in range(len(index) + 1):

                # Iterate over the combinations of num_bit_factors variables to which the
                # num_bit_factors bits relate.
                for bit_factors_positions in combinations(range(len(index)), num_bit_factors):
                    variables_at_positions = [index[position] for position in bit_factors_positions]

                    # For each variable related to a bit, iterate over the bits representing
                    # the variable.
                    for bits_indices in product(*[
                        range(len(self.variables[variable_index]["bits"]))
                        for variable_index in variables_at_positions
                    ]):

                        bit_couplings_index = tuple(sorted(self.variables[variable_index]["bits"][bits_indices[i]] for i, variable_index in enumerate(variables_at_positions)))

                        # If a bit is repeated in the bits monomial, the latter can be simplified (b^2 = b). This
                        # simplification is performed by a conversion to a set.
                        bit_couplings_index = tuple(sorted(list(set(bit_couplings_index))))

                        a = [(self.variables[variable_index]["max"] - self.variables[variable_index]["min"]) / (2 ** len(self.variables[variable_index]["bits"]) - 1) * 2 ** bits_indices[i] for i, variable_index in enumerate(variables_at_positions)]

                        b = [self.variables[variable_index]["min"] for i, variable_index in enumerate(index) if i not in bit_factors_positions]

                        bit_couplings_value = reduce(lambda x, y: x * y, a + b, 1) * value

                        if bit_couplings_index not in bit_couplings:
                            bit_couplings[bit_couplings_index] = 0

                        bit_couplings[bit_couplings_index] += bit_couplings_value

        return bit_couplings

    def to_spins(self):
        r"""Output the polynomial in spin variables encoding the cost function.

        Returns
        -------
        dict
            A dictionary describing the monomials of the polynomial in spin
            variables that encodes the cost function. Precisely, a monomial
            $c\cdot Z_{i_1}\cdots Z_{i_k}$ is represented by the (key, value) pair
            ((i_1, ..., i_k), c).
        """

        spin_couplings = {}

        # Iterate over monomial terms.
        for index, value in self.couplings.items():

            # The monomial term will be represented by a linear combination of monomials in
            # spin variables. Each monomial is the product of a certain number of spin
            # variables; iterate over this number.
            for num_sigma_factor in range(len(index) + 1):

                # Iterate over the combinations of num_sigma_factor variables to which the
                # num_sigma_factor spin variables relate.
                for sigma_factors_positions in combinations(range(len(index)), num_sigma_factor):
                    variables_at_positions = [index[position] for position in sigma_factors_positions]

                    # For each variable related to a spin variable, iterate over the bits representing
                    # the variable.
                    for bits_indices in product(*[
                        range(len(self.variables[variable_index]["bits"]))
                        for variable_index in variables_at_positions
                    ]):

                        spin_couplings_index = tuple(sorted(self.variables[variable_index]["bits"][bits_indices[i]] for i, variable_index in enumerate(variables_at_positions)))

                        # If a spin variable is repeated in the monomial, the latter can be simplified (as sigma^2 = 1).
                        # Perform this simplification, storing the result to reduced_spin_couplings_index, then to
                        # spin_couplings_index.
                        reduced_spin_couplings_index = []
                        identical_indices = 1

                        if len(spin_couplings_index):
                            for i in range(1, len(spin_couplings_index)):
                                if spin_couplings_index[i] == spin_couplings_index[i - 1]:
                                    identical_indices += 1
                                else:
                                    if identical_indices % 2:
                                        reduced_spin_couplings_index.append(spin_couplings_index[i - 1])
                                    identical_indices = 1

                            if identical_indices % 2:
                                reduced_spin_couplings_index.append(spin_couplings_index[len(spin_couplings_index) - 1])

                        spin_couplings_index = tuple(reduced_spin_couplings_index)

                        a = [(self.variables[variable_index]["max"] - self.variables[variable_index]["min"]) / (2 * (2 ** len(self.variables[variable_index]["bits"]) - 1)) * 2 ** bits_indices[i] for i, variable_index in enumerate(variables_at_positions)]

                        b = [(self.variables[variable_index]["max"] + self.variables[variable_index]["min"]) / 2 for i, variable_index in enumerate(index) if i not in sigma_factors_positions]

                        spin_couplings_value = reduce(lambda x, y: x * y, a + b, 1) * value

                        if spin_couplings_index not in spin_couplings:
                            spin_couplings[spin_couplings_index] = 0

                        spin_couplings[spin_couplings_index] += spin_couplings_value

        return spin_couplings

    def bits_to_variables(self, bits):
        """
        Convert an iterable sequence of bits (0/1) to the values of the
        variables they encode.

        Parameters
        ----------
        bits: iterable
            The bits encoding the values of the variables.


        Returns
        -------
        dict
            A dictionary mapping each variable name to its value encoded in the
            bits.
        """

        d = {}

        for idx, var in self.variables.items():
            d[idx] = var["min"] + (var["max"] - var["min"]) / (2 ** len(var["bits"]) - 1) * sum([bits[bit] * 2 ** bit_pos for bit_pos, bit in enumerate(var["bits"])])

        return d

    def spins_to_variables(self, spins):
        """
        Convert an iterable sequence of spins(-1/1) to the values of the
        variables they encode.

        Parameters
        ----------
        spins: iterable
            The spins encoding the values of the variables.

        Returns
        -------
        dict
            A dictionary mapping each variable name to its value encoded in the spins.
        """

        return self.bits_to_variables(tuple((1 + spin) / 2 for spin in spins))

    def variables_to_bits(self, variables):
        """
        Convert a dictionary variables to the corresponding sequence of
        encoding bits.

        Parameters
        ----------
        variables: dict
            A dictionary mapping each variable name to its value.

        Returns
        -------
        list
            A sequence of bits encoding the variables passed in argument.
        """
        bits = [0] * self.num_spins()
        for variable_name, variable_encoding in self.variables.items():
            num_variable_bits = len(variable_encoding["bits"])
            variable_integer_value = int(
                round(
                    (2 ** num_variable_bits - 1)
                    * (variables[variable_name] - variable_encoding["min"]) / (variable_encoding["max"] - variable_encoding["min"])
                )
            )
            variable_bits = [int(bit) for bit in f"{{:0{num_variable_bits}b}}".format(variable_integer_value)[::-1]]
            for idx, variable_bit in enumerate(variable_bits):
                bits[variable_encoding["bits"][idx]] = variable_bit
        return bits

    def variables_to_spins(self, variables):
        """
        Convert a dictionary variables to the corresponding sequence of
        encoding spins.

        Parameters
        ----------
        variables: dict
            A dictionary mapping each variable name to its value.

        Returns
        -------
        list
            A sequence of spins encoding the variables passed in argument.
        """
        return [2 * bit - 1 for bit in self.variables_to_bits(variables)]

    def bits_cost(self, bits):
        """Compute the cost function for the values of the variables specified by bits.

        Parameters
        ----------
        bits: iterable
            The bits encoding the values of the variables.

        Returns
        -------
        float
            The value of the cost function.
        """
        spin_couplings = self.to_spins()
        cost = 0

        for index, value in spin_couplings.items():
            cost += value * reduce(lambda x, y: x * y, [2 * bits[spin_index] - 1 for spin_index in index], 1)

        return cost

    def spins_cost(self, spins):
        """Compute the cost function for the values of the variables specified by spins.

        Parameters
        ----------
        spins: iterable
            The spins (+/-1 integers) encoding the values of the variables.

        Returns
        -------
        float
            The value of the cost function.
        """

        return self.bits_cost(tuple((1 + spin) / 2 for spin in spins))

    def variables_cost(self, variables):
        """Compute the cost function for some assignment of the variables.

        Parameters
        ----------
        variables: dict
            A dictionary mapping each variable name to its value.

        Returns
        -------
            The value of the cost function.
        """

        cost = 0

        for index, value in self.couplings.items():
            cost += value * reduce(lambda x, y: x * y, [variables[variable_index] for variable_index in index], 1)

        return cost

    def objectives_enumerate(self):
        """
        Enumerate the values taken by the cost function together with the
        corresponding variable assignments.

        Returns
        -------
        list
            A list of (cost, assignment) pairs where 'assignment' is a
            dictionary mapping each variable name to the value assigned to it
            and 'cost' is the value of the cost function for this assignment.
        """

        objectives = []

        a = product(*[[variable["min"] + (variable["max"] - variable["min"]) * discrete_value / (2 ** len(variable["bits"]) - 1) for discrete_value in range(2 ** len(variable["bits"]))] for variable in self.variables.values()])

        for values_flat in a:

            values = dict(zip(self.variables.keys(), values_flat))
            objective = 0

            for index, coupling in self.couplings.items():
                objective += coupling * reduce(lambda x, y: x * y, [values[variable_index] for variable_index in index], 1)

            objectives.append((objective, values))

        return objectives

    def objectives_encoded_enumerate(self):
        """
        Enumerate the values taken by the cost function together with the
        corresponding variable assignment, the values of the variables being
        encoded in bits.

        Returns
        -------
        list
            A list of (cost, spins) pairs where 'spins' is a tuple of spins
            (-1/+1 integers) describing the assignment of the variables as
            encoded in binary and 'cost' is the value of the cost function for
            this assignment.
        """

        encoded_couplings = self.to_spins()

        max_spin = max(max(variable["bits"]) for variable in self.variables.values())
        objectives = []

        for values in product([-1, 1], repeat=(max_spin + 1)):
            objective = 0

            for index, coupling in encoded_couplings.items():
                objective += coupling * reduce(lambda x, y: x * y, [values[spin_index] for spin_index in index], 1)

            objectives.append((objective, values))

        return objectives
