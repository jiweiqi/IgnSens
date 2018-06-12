# IgnitionDelaySensitivity
Efficient approach for evaluating the sensitivity of ignition delay time. Implement with Cantera/Python

# Instruction

# API
+ idtSens(gas, T, P, X, Type = 'UV', factor = 0.05 )

    + gas: The object of IdealGasMixture
    + T: [K]
    + P: [atm]
    + X: Mixture of Fuel and Air, e.g., 'H2:1.0,O2:1.0,N2:3.76'
    + Type: 'UV': Constant Vloume Reactor, 'HP': Constant Pressure Reactor
    + factor: the pertubation factor for the finite difference (brute force) approach, default value is 0.05

Alternatively, the reaction mechanism .xml file can be passed. [Proposal]

# Demo
GRImech 3.0 which is included in the Cantera package