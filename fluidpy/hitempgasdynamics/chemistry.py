# The below is to ignore numpy errors from pylint:
#     pylint: disable=E1101

"""
This module has some classes used to represent equilibrium air.
It only works for thermally perfect gases, i.e., internal energy,
enthalpy, and specific heats are functions of temperature alone.
"""

import sys
import numpy as np

class Species:
    """Represents a chemical species.

    This class represents a chemical species.  It can be thought of
    as a C-style structure which contains data for the species and
    some methods.

    :param molecular_weight: molecular weight of the given species in kg/mol.
    :param enthalpy_of_formation: enthalpy of formation of the species in J/kg.
    :param name: name of species, e.g., N, O, N2, O2, NO.
    :param stoichiometric_coeff: stoichiometric coeffient for species in reaction.
    :param molecule_type: monatomic or diatomic.
    :param theta_v: constant used to compute internal energy and enthalpy for diatomic molecules.
    :param temperature_ref: reference temperature used for enthalpy of formation.
    """
    def __init__(self, molecular_weight, enthalpy_of_formation, name, stoichiometric_coeff,
                 molecule_type, theta_v=None, temperature_ref=5.0):
        self.r_u = 8.314                                   # J/(mol K)
        self.molecular_weight = molecular_weight           # kg/mol
        self.gas_constant = self.r_u/self.molecular_weight # J/(kg K)
        self.h_f = enthalpy_of_formation                   # J/kg
        self.molecule_type = molecule_type                 # e.g., diatomic
        self.theta_v = theta_v                             # constant used to compute e,h
        self.name = name                                   # Name of the species
        self.pressure = None                               # Partial pressure of the mixture
        self.stoichiometric_coeff = stoichiometric_coeff   # stoichiometric coefficient
        self.temperature_ref = temperature_ref             # Reference temp for h_f
        self.k = 10                                        # number used to create temperature_vec
        self.mole_fraction = None                          # mole fraction
        self.mass_fraction = None                          # mass fraction

    def c_v(self, temperature):
        """Specific heat at constant volume.

        This method computes the specific heat at constant volume when provided
        with the temperature.  It only works for monatomic and diatomic molecules.
        This is done using kinetic theory and the equipartition theorem.

        :param temperature: temperature in Kelvin.
        """
        if self.molecule_type == "monatomic":
            return 3.0/2.0*self.gas_constant + 0.0*temperature
        elif self.molecule_type == "diatomic":
            return 3.0/2.0*self.gas_constant + self.gas_constant + (self.theta_v/temperature)**2\
                   *np.exp(self.theta_v/temperature)/(np.exp(self.theta_v/temperature)-1.0)**2\
                   *self.gas_constant

    def c_p(self, temperature):
        """Specific heat at constant pressure.

        This method computes the specific heat at constant pressure when provided
        with the temperature.  It only works for monatomic and diatomic molecules.
        This is accomplished using kinetic theory and the equipartition theorem.

        :param temperature: temperature in Kelvin.
        """
        return self.c_v(temperature) + self.gas_constant

    def internal_energy(self, temperature):
        """Computes internal energy per unit mass.

        This method computes the internal energy of the species when provided
        with the temperature.

        :param temperature: temperature in Kelvin.
        """
        if self.molecule_type == "monatomic":
            return 3.0/2.0*self.gas_constant*temperature + self.h_f
        elif self.molecule_type == "diatomic":
            return 3.0/2.0*self.gas_constant*temperature + self.gas_constant*temperature \
              + self.theta_v*self.gas_constant/(np.exp(self.theta_v/temperature)-1.0) + self.h_f

    def enthalpy(self, temperature):
        """Computes enthalpy per unit mass.

        This method computes the enthalpy of the species when provided with
        the temperature.

        :param temperature: temperature in Kelvin.
        """
        return self.internal_energy(temperature) + self.gas_constant*temperature

class Reaction:
    """Represents a chemical reaction.

    This class represents a chemical reaction.  It contains data for
    the reaction and lists of species objects for reactants and products.
    """

    def __init__(self, reactants, products, eqn, c_c=None, eta_c=None, theta=None):
        """Constructor for `Reaction` class.

        :param reactants: reactants in chemical reaction.
        :param products: products in the chemical reaction.
        :param eqn: equation for chemical reaction (e.g., 'N2 -> 2N').
        :param c_c: constant used to compute equilibrium constant.
        :param eta_c: constant used to compute equilibrium constant.
        :param theta: constant used to compute equilibrium constant.
        """
        self.reactants = reactants
        self.products = products
        self.c_c = c_c
        self.eta_c = eta_c
        self.theta = theta
        self.r_u = 8.314      # J/(mol K)
        self.eqn = eqn

    def get_equilibrium_constant(self, temperature):
        """This method returns the equilibrium constant in Pascals.

        :param temperature: temperature in Kelvin.
        """
        k_c = self.c_c*temperature**self.eta_c*np.exp(-self.theta/temperature)
        nu_sum = 0.0
        for reactant in self.reactants:
            nu_sum -= reactant.stoichiometric_coeff
        for product in self.products:
            nu_sum += product.stoichiometric_coeff

        return (self.r_u*temperature)**nu_sum*k_c*1.0e3

class ReactingMixture:
    """Represents a reacting mixture of gases.

    This class represents a reacting mixture.  It holds a list of
    reactions.  It's main purpose is to find equilibrium.

    Methods:

      __init__

        Arguments:
          - reactions: a list of reaction objects.
          - pressure : the pressure of the mixture.
          - eqn      : formula for reaction (e.g., O2 -> 2O)

    """

    def __init__(self, reactions, pressure=None, tol=1.0e-3, max_iter=5000):
        """Constructor for `ReactingMixture` class.

        :param reactions: list of reaction objects.
        :param pressure: the pressure of the mixture.
        :param tol: tolerance used in iterative procedure.
        :param max_iter: max number of iterations.
        """
        self.reactions = reactions
        self.p_m = pressure
        self.tol = tol
        self.hist = []
        self.max_iter = max_iter

    def find_equilibrium(self):
        """Method for finding chemical equilibrium.

        This is not done for the general case yet.
        """
        print("Not done yet.")
        sys.exit()

#####################################################################

class Air(ReactingMixture):
    """
    This class is derived from the reacting mixture.  It is only for air.
    The difference is that it uses a special iterative method to solve
    the nonlinear system of equations.  The iterative method is different
    because specific assumptions about dissociation, etc, come into play.

    .. note:: Air only works for pressure in atm.
    """

    def __init__(self, reactions, tol=1.0e-6, molar_ratio=3.76, max_iter=5000):
        self.molar_ratio = molar_ratio   # Ratio of N2 to O2
        super().__init__(reactions, tol=tol, max_iter=max_iter)

    def find_equilibrium(self, p_m, T):
        """Method for finding the equilibrium composition of air given the state.

        This method uses a few assumptions along with what amounts to a nonlinear
        Gauss-Seidel procedure to determine the chemical composition for air in
        equilibrium.

        :param p_m: pressure of the mixture in atmospheres (i.e., relative to 101.325e3 Pa).
        :param T: temperature in Kelvin.
        """

        # Sorting reactions
        use_no = False
        for reaction in self.reactions:
            if reaction.eqn == "O2 -> 2O":
                o2_dissociation = reaction
            elif reaction.eqn == "N2 -> 2N":
                n2_dissociation = reaction
            elif reaction.eqn == "NO -> N + O":
                print("Using NO synthesis.")
                no_synthesis = reaction
                use_no = True

        # Sorting species
        for reaction in self.reactions:
            for species in reaction.reactants:
                if species.name == "N2":
                    n2 = species
                elif species.name == "O2":
                    o2 = species
                elif species.name == "NO":
                    no = species
            for species in reaction.products:
                if species.name == "N":
                    n = species
                elif species.name == "O":
                    o = species

        # Checking to see if NO is being produced
        if use_no:

            # Initializing variables
            resid = 1.0
            nhist = []
            ohist = []
            n2hist = []
            o2hist = []
            nohist = []

            # Setting up initial guess
            equilibrium_constant1 = o2_dissociation.get_equilibrium_constant(T)/101.325e3
            equilibrium_constant2 = n2_dissociation.get_equilibrium_constant(T)/101.325e3
            equilibrium_constant3 = no_synthesis.get_equilibrium_constant(T)/101.325e3
            print("\nKp1 = {:1.4e} atm Kp2 = {:1.4e} atm Kp3 = {:1.4e} atm\n"\
                  .format(equilibrium_constant1, equilibrium_constant2, equilibrium_constant3))
            o.pressure = 2.0*p_m/(2.0 + self.molar_ratio)
            n2.pressure = o.pressure/2.0*self.molar_ratio
            n.pressure = np.sqrt(n2.pressure*equilibrium_constant2)
            o2.pressure = o.pressure**2/equilibrium_constant1
            no.pressure = o.pressure*n.pressure/equilibrium_constant3
            nhist.append(n.pressure)
            ohist.append(o.pressure)
            n2hist.append(n2.pressure)
            o2hist.append(o2.pressure)
            nohist.append(no.pressure)

            # Iterating (after this loop is done, each species object has
            # the correct partial pressure.)
            ctr = 1
            while resid > self.tol and ctr < self.max_iter:

                o2.pressure = o.pressure**2/equilibrium_constant1
                n.pressure = np.sqrt(equilibrium_constant2*n2.pressure)
                no.pressure = o.pressure*n.pressure/equilibrium_constant3
                n2.pressure = 0.5*(self.molar_ratio*(2.0*o2.pressure + o.pressure + no.pressure) \
                              - n.pressure - no.pressure)
                o.pressure = p_m - o2.pressure - n.pressure - no.pressure - n2.pressure

                nhist.append(n.pressure)
                ohist.append(o.pressure)
                n2hist.append(n2.pressure)
                o2hist.append(o2.pressure)
                nohist.append(no.pressure)
                resid = np.fabs(n2hist[n] - n2hist[n-1])

                print("Iteration: {0:>5d} Residual: {1:>1.4e}".format(ctr, resid))

                ctr += 1
            print("\nDone.\n")
            self.hist = [ohist, nhist, o2hist, n2hist, nohist]

        else:

            # Initializing variables
            resid = 1.0
            nhist = []
            ohist = []
            n2hist = []
            o2hist = []

            # Setting up initial guess
            equilibrium_constant1 = o2_dissociation.get_equilibrium_constant(T)/101.325e3
            equilibrium_constant2 = n2_dissociation.get_equilibrium_constant(T)/101.325e3
            print("\nKp1 = {:1.4e} atm Kp2 = {:1.4e} atm\n"\
                  .format(equilibrium_constant1, equilibrium_constant2))
            o.pressure = 2.0*p_m/(2.0 + self.molar_ratio)
            n2.pressure = o.pressure/2.0*self.molar_ratio
            n.pressure = np.sqrt(n2.pressure*equilibrium_constant2)
            o2.pressure = o.pressure**2/equilibrium_constant1
            nhist.append(n.pressure)
            ohist.append(o.pressure)
            n2hist.append(n2.pressure)
            o2hist.append(o2.pressure)

            # Iterating (after this loop is done, each species object has
            # the correct partial pressure.)
            ctr = 1
            while resid > self.tol and ctr < self.max_iter:

                o2.pressure = o.pressure**2/equilibrium_constant1
                n.pressure = np.sqrt(equilibrium_constant2*n2.pressure)
                n2.pressure = 0.5*(self.molar_ratio*(2.0*o2.pressure + o.pressure) - n.pressure)
                o.pressure = (2.0*p_m - 2.0*o2.pressure*(1.0 + self.molar_ratio) - n.pressure)/\
                             (2.0 + self.molar_ratio)

                nhist.append(n.pressure)
                ohist.append(o.pressure)
                n2hist.append(n2.pressure)
                o2hist.append(o2.pressure)
                resid = np.fabs(n2hist[ctr] - n2hist[ctr-1])

                print("Iteration: {0:>5d} Residual: {1:>1.4e}".format(n, resid))

                ctr += 1
            print("\nDone.\n")
            self.hist = [ohist, nhist, o2hist, n2hist]
