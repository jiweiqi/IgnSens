# IgnitionDelaySensitivity

Efficient approach for evaluating the sensitivity of ignition delay time. Implement with Cantera/Python

# Instruction

+ If you need the relative sensitivities, you can simply use the local temeprature sensitvity. But take care of the positive/negative sign. You can figure out the sign based on intution, or compute one of the sensitvity with brute force approach.

+ If you need the absolute value of the sensitvities. You can follow the paper and compute two top sensitive reactions.

# API 【Initially, I planed to write an api, later on , I realized that it is unnessary ... But if you find that an api might be useful, you can start from the proposal bellow. 

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

# Cite

Weiqi Ji, Zhuyin Ren, Chung K. Law, Evolution of Sensitivity Directions during Autoignition, Proceedings of the Combustion Instituts, Accepted Manuscript. 2018