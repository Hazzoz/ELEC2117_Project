module ELEC2117_Project

    # Necessary Packages
    using DifferentialEquations
    using Plots

    # SIR Functions
    include("mycode.jl")
    export SIRFoI, SIRFoIIntervention
    export solve_SIR, plot_SIR
    export Level2_plot_infected
    export Level3_plot_infected, Level3_error_beta
    export Level4_1_plot_intervention_with_error, Level4_1_compare_intervention
    export Level4_2_plot_coverage, Level4_2_error_coverage, Level4_2_plot_coverage_with_error, Level4_2_compare_intervention_coverage
    export Level5_best_beta_day_plot, Level5_plot_second_town

end
