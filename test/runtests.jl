using ELEC2117_Project
using Test

@testset "Level 2" begin
    # Dummy values for parameters
    S0 = 5999        # Initial susceptible population
    I0 = 1           # Initial infected population
    SI0 = 0          # Initial seriously infected population
    R0 = 0           # Initial recovered population
    days = 25        # Number of days to simulate

    # Example parameters structure, replace with actual SIR parameters as needed
    params = SIRFoI(0.036, 1/7, 8, 0.15, 1/14, 1/28)

    # Define the ratio range for testing
    ratio_range = (0.15, 0.25)

    # Test with infect = 1 for infected population and level = 2
    level = 2
    infect = 1

    # Run the function
    R0 = plot_infected(S0, I0, SI0, R0, days, params, ratio_range, infect, level)

    # Check R0 is computed and output
    @test R0 == (8*0.036/(1/7))
end

@testset "Level 3" begin
    # Dummy values for parameters
    S0 = 5999        # Initial susceptible population
    I0 = 1           # Initial infected population
    SI0 = 0          # Initial seriously infected population
    R0 = 0           # Initial recovered population
    days = 30        # Number of days to simulate

    # Example parameters structure, replace with actual SIR parameters as needed
    params = SIRFoI(0.036, 1/7, 8, 0.15, 1/14, 1/28)

    # Define the ratio range for testing
    ratio_range = (0.15, 0.25)

    # Test with infect = 1 for infected population and level = 2
    level = 2
    infect = 1

    # Run the function
    R0 = plot_infected(S0, I0, SI0, R0, days, params, ratio_range, infect, level)

    # Check R0 is computed and output
    @test R0 == (8*0.036/(1/7))

    # Define the ratio and beta ranges for testing
    beta_range = (0.025, 0.042)       # Example range for beta

    # Run the Level3_error_beta function
    beta_min = Level3_error_beta(S0, I0, SI0, R0, days, params, beta_range, ratio_range)

    @test beta_min ≈ 0.0347 atol=0.01
end

@testset "Level 4" begin
    # Dummy values for parameters
    S0 = 5999        # Initial susceptible population
    I0 = 1           # Initial infected population
    SI0 = 0          # Initial seriously infected population
    R0 = 0           # Initial recovered population
    days = (30,120)  # Days for each phase: [pre-intervention, post-intervention]

    # Example parameters structure, replace with actual SIR parameters as needed
    params = SIRFoI(0.036, 1/7, 8, 0.15, 1/14, 1/28)  # Adjust based on your model
    params2 = SIRFoIIntervention(0.036, 1/7, 8, 0.15, 1/14, 1/28, 0.3,0.8) # Adjust post-intervention parameters

    # Define the ratio and beta ranges for testing
    ratio_range = (0.15, 0.25)      # Example range for SIratio
    beta_range = (0.0335, 0.0355)       # Example range for beta

    # Run the Level4_1_plot_intervention_with_error function
    max_infected_peak, min_infected_peak, max_seriously_infected_peak, min_seriously_infected_peak = Level4_1_plot_intervention_with_error(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range)

    @test min_infected_peak > 0
    @test max_infected_peak > 0
    @test min_seriously_infected_peak > 0
    @test max_seriously_infected_peak > 0
    @test max_infected_peak > min_infected_peak  # Max peak should be higher than min peak
    @test max_seriously_infected_peak > min_seriously_infected_peak  # Same for seriously infected

    @test (max_infected_peak+min_infected_peak)/2 > (max_seriously_infected_peak+min_seriously_infected_peak)/2

    max_infected_peak1, min_infected_peak1, max_seriously_infected_peak1, min_seriously_infected_peak1, max_infected_peak_intervention1, min_infected_peak_intervention1, max_seriously_infected_peak_intervention1, min_seriously_infected_peak_intervention1 = Level4_1_compare_intervention(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, 1, 0)
    max_infected_peak2, min_infected_peak2, max_seriously_infected_peak2, min_seriously_infected_peak2, max_infected_peak_intervention2, min_infected_peak_intervention2, max_seriously_infected_peak_intervention2, min_seriously_infected_peak_intervention2 = Level4_1_compare_intervention(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, 1, 1)

    @test (max_infected_peak1+min_infected_peak1)/2 > (max_infected_peak_intervention2+min_infected_peak_intervention2)/2

    coverage_range = (0,1)

    coverage_min = Level4_2_error_coverage(S0, I0, SI0, R0, (30,25), params, beta_range, params2, ratio_range, coverage_range, 0)

    @test coverage_min ≈ 0.592 atol=0.01
end

@testset "Level 5" begin
    # Dummy values for parameters
    S0 = 10000        # Initial susceptible population
    I0 = 1           # Initial infected population
    SI0 = 0          # Initial seriously infected population
    R0 = 0           # Initial recovered population
    days = (35,45)  # Days for each phase: [pre-intervention, post-intervention]

    # Example parameters structure, replace with actual SIR parameters as needed
    params = SIRFoI(0.036, 1/7, 8, 0.15, 1/14, 1/28)  # Adjust based on your model
    params2 = SIRFoIIntervention(0.036, 1/7, 8, 0.15, 1/14, 1/28, 0.3,0.8) # Adjust post-intervention parameters

    # Define the ratio and beta ranges for testing
    ratio_range = (0.15, 0.25)      # Example range for SIratio
    beta_range = (0.03, 0.05)       # Example range for beta
    coverage_range = (0.5,0.7)

    best_beta, best_day = Level5_best_beta_day_plot(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, 2, (1,15))

    @test best_day == 13
    @test best_beta ≈ 0.0403 atol=0.01
end