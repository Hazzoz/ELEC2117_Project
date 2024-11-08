module ELEC2117_Project

    # Necessary Packages
    using DifferentialEquations
    using Plots

    # SIR Functions
    include("mycode.jl")
    export solve_SIR, plot_SIR, plot_infected, error_beta, plot_intervention_with_error, compare_intervention, plot_infected_intervention, error_coverage, plot_coverage
    export SIRFoI, SIRFoIIntervention

end
