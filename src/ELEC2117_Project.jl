module ELEC2117_Project

    # Necessary Packages
    using DifferentialEquations
    using Plots

    # SIR Functions
    include("mycode.jl")
    export solve_SIR, plot_SIR, plot_infected, error_beta, plot_intervention_with_error, compare_intervention
    export plot_infected_intervention, error_coverage, plot_coverage, compare_coverage, plot_coverage_with_error
    export best_beta_day_plot
    export SIRFoI, SIRFoIIntervention

end
