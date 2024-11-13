# Structs for different versions of the SIR models
mutable struct SIRFoI
    beta::Float64  # Transmission chance of interaction
    gamma::Float64  # Recovery rate
    contacts::Float64  # Number of daily contacts
    SIratio::Float64  # Ratio of people becoming seriously infected
    delta::Float64  # Seriously infected recovery rate
    alpha::Float64  # Re-susceptible rate
end

mutable struct SIRFoIIntervention
    beta::Float64  # Transmission chance of interaction
    gamma::Float64  # Recovery rate
    contacts::Float64  # Number of daily contacts
    SIratio::Float64  # Ratio of people becoming seriously infected
    delta::Float64  # Seriously infected recovery rate
    alpha::Float64  # Re-susceptible rate
    epsilon::Float64  # Intervention Efficacy Rate
    p::Float64 # Coverage rate
end

"""
This SIR! is a function that represents the differential
equations that defines a more detailed force of infection
version of the SIR model.

Inputs:
- dP = array of gradients for each population
- P = array of current values for each population
- params = array of other necessary parameters, defined by the structs
- t = timespan
"""
function SIR!(dP, P, params::SIRFoI, t)
    N = P[1] + P[2] + P[3] + P[4] # Total Population
    lambda = P[2]/N*params.beta*params.contacts # Force of Infection
    dP[1] = params.alpha*P[4] - lambda*P[1] # Change in Susceptible Population
    dP[2] = lambda*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*params.SIratio*P[2] - params.delta*P[3] # Change in Seriously Infected Population
    dP[4] = params.gamma*(1-params.SIratio)*P[2] + params.delta*P[3] - params.alpha*P[4] # Change in Recovered Population
end

"""
This SIR! is a function that represents the differential
equations that defines a more detailed force of infection
version of the SIR model.

Inputs:
- dP = array of gradients for each population
- P = array of current values for each population
- params = array of other necessary parameters
- t = timespan
"""
function SIR!(dP, P, params::SIRFoIIntervention, t)
    N = P[1] + P[2] + P[3] + P[4] # Total Population
    lambda = P[2]/N*params.beta*params.contacts*(1-params.epsilon*params.p) # Force of Infection
    dP[1] = params.alpha*P[4] - lambda*P[1] # Change in Susceptible Population
    dP[2] = lambda*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*params.SIratio*P[2] - params.delta*P[3] # Change in Seriously Infected Population
    dP[4] = params.gamma*(1-params.SIratio)*P[2] + params.delta*P[3] - params.alpha*P[4] # Change in Recovered Population
end

"""
solve_SIR is a driver function that chooses the required SIR model.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
"""
function solve_SIR(S0, I0, SI0, R0, days, params)
    P0 = [S0, I0, SI0, R0] # Initial populations vector

    tspan = (0, days) # Time span tuple

    solution = solve(ODEProblem(SIR!, P0, tspan, params)) # Solve the ODE with given parameters and timespan

    return solution # Return values
end

"""
plot_SIR is a driver function that runs and then plots an SIR model.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
"""
function plot_SIR(S0, I0, SI0, R0, days, params)
    solution = solve_SIR(S0, I0, SI0, R0, days, params) # Solve the SIR model

    plot(solution, xlabel="Time", ylabel="Population", title="Solution", labels=["Susceptible" "Infected" "Seriously Infected" "Recovered"]) # Plot the model
end

"""
Level2_plot_infected is a driver function that plots either the infected or
seriously infected data against the true data to allow for visual
inspection of the model's fit. This function is for plotting the level 2 data.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- infect = determines whether to plot infected or seriously infected data (1 = infected, else = seriously infected)
"""
function Level2_plot_infected(S0, I0, SI0, R0, days, params, ratio_range, infect)
    plot_infected(S0, I0, SI0, R0, days, params, infect, ratio_range, 2)  # Call the plotting function for level 2
end

"""
Level3_plot_infected is a driver function that plots either the infected or
seriously infected data against the true data to allow for visual
inspection of the model's fit. This function is for plotting the level 3 data.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- infect = determines whether to plot infected or seriously infected data (1 = infected, else = seriously infected)
"""
function Level3_plot_infected(S0, I0, SI0, R0, days, params, ratio_range, infect)
    plot_infected(S0, I0, SI0, R0, days, params, infect, ratio_range, 3)  # Call the plotting function for level 3
end

"""
plot_infected is a function that plots either the infected or
seriously infected data against the true data to allow for visual
inspection of the model's fit.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- infect = determines whether to plot infected or seriously infected data (1 = infected, else = seriously infected)
- level = which TLT level is being addressed
"""
function plot_infected(S0, I0, SI0, R0, days, params, ratio_range, infect, level)
    # Set actual data and corresponding time points based on the level
    actual_infected, ti, actual_seriously_infected, tsi = if level == 3
        ([11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55],
         [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30],
         [0,0,1,2,5,5,5,2,9,4],
         [21,22,23,24,25,26,27,28,29,30])
    else
        ([11,7,20,3,29,14,11,12,16,10,58],
         [15,16,17,18,19,20,21,22,23,24,25],
         [0,0,1,2,5],
         [21,22,23,24,25])
    end

    # Create a base plot
    plot([0],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing)

    # Generate ratios to test
    ratios = range(ratio_range[1], ratio_range[2], 5)

    # Initialize arrays
    infected, seriously_infected, time = Float64[], Float64[], Float64[]
    previous_infected, previous_seriously_infected = Float64[], Float64[]

    # Iterate through each ratio
    for r in ratios
        params.SIratio = r  # Update the seriously infected ratio
        solution = solve_SIR(S0, I0, SI0, R0, days, params)  # Solve the SIR model

        # Reset arrays for current solution
        infected, seriously_infected, time = Float64[], Float64[], Float64[]

        # Extract infected, seriously infected, and time data
        for i = 1:length(solution.t)
            push!(infected, solution.u[i][2])
            push!(seriously_infected, solution.u[i][3])
            push!(time, solution.t[i])
        end

        # Ensure previous arrays are initialized and match lengths
        if isempty(previous_infected)
            previous_infected = deepcopy(infected)
        end
        if isempty(previous_seriously_infected)
            previous_seriously_infected = deepcopy(seriously_infected)
        end

        while length(previous_infected) < length(infected)
            push!(previous_infected, previous_infected[end])
        end
        while length(previous_seriously_infected) < length(seriously_infected)
            push!(previous_seriously_infected, previous_seriously_infected[end])
        end

        # Plot based on the type of infection data
        if infect == 1
            plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels=nothing, fillrange=previous_infected, colour=:blue)
        else
            plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels=nothing, fillrange=previous_seriously_infected, colour=:blue)
        end

        # Update previous arrays for the next iteration
        previous_infected, previous_seriously_infected = deepcopy(infected), deepcopy(seriously_infected)
    end

    # Compute and display R0 if level is 2
    if level == 2
        R0 = params.contacts * params.beta / params.gamma
        println("R0: ", round(R0, digits=3))
    end

    # Overlay actual data on the plot
    if infect == 1
        plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels="Infected", colour=:blue)
        plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected", colour=:red)
    else
        plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected", colour=:blue)
        plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected", colour=:red)
    end
end

"""
Level3_error_beta is a function that plots the RMSE of a range of beta values
using the predicted and actual data provided.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
- beta_range = range of betas that are to be tested
- ratio_range = range of proportions of people that become seriously infected
"""
function Level3_error_beta(S0, I0, SI0, R0, days, params, beta_range, ratio_range)
    # Actual data for infected population and corresponding time points
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]

    # Actual data for seriously infected population and corresponding time points
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4]
    tsi = [21,22,23,24,25,26,27,28,29,30]

    # Setup empty plot for beta vs error
    plot([beta_range[1]],[0], xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing)

    # Initialize metrics to track minimum error and associated beta value
    beta_min = 1
    min = S0

    # Create arrays for storing error values and range of betas and ratios to test
    previous_error_vals = error_vals = Float64[]
    betas, ratios= range(beta_range[1], beta_range[2], length=50), range(ratio_range[1], ratio_range[2], 5)

    # Loop through ratios
    for r in ratios
        params.SIratio = r
        error_vals = Float64[]
        # Loop through range of betas
        for b in betas
            params.beta = b  # Update beta in parameters
            solution = solve_SIR(S0, I0, SI0, R0, days, params)  # Solve the model

            # Initialize error metrics for current beta
            current_infected_error = 0
            current_seriously_infected_error = 0

            # Calculate errors for infected population
            for j = 1:length(solution.t)
                for k = 1:length(ti)
                    if ti[k] ≈ solution.t[j] atol=0.01  # Check if time values are approximately equal
                        current_infected_error += (solution[j][2] - actual_infected[k])^2  # Add squared error
                    end
                end

                # Calculate errors for seriously infected population
                for k = 1:length(tsi)
                    if tsi[k] ≈ solution.t[j] atol=0.01  # Check if time values are approximately equal
                        current_seriously_infected_error += (solution[j][3] - actual_seriously_infected[k])^2
                    end
                end
            end

            # Compute RMSE and save it to error array
            push!(error_vals, sqrt(current_infected_error + current_seriously_infected_error))
        end

        # Check if the current ratio's error curve has a lower minimum error
        if minimum(error_vals) < min
            min = minimum(error_vals)
            beta_min = betas[argmin(error_vals)]
        end

        # Ensure previous error values are not empty for fill range
        if isempty(previous_error_vals)
            previous_error_vals = deepcopy(error_vals)
        end

        # Match the length of previous error values for filling
        while length(previous_error_vals) < length(error_vals)
            push!(previous_error_vals, previous_error_vals[end])
        end

        # Plot current error curve with shaded region
        plot!(betas, error_vals, xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing, fillrange=previous_error_vals, colour=:blue)

        # Save current error curve for the next iteration
        previous_error_vals = deepcopy(error_vals)
    end

    # Output the beta value associated with the minimum error
    println("Beta min: ", beta_min)

    # Final plot of the error curve
    plot!(betas, error_vals, xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing, colour=:blue)
end


"""
Level4_1_plot_intervention_with_error is a function that plots the behaviour of the infected and 
seriously infected curves from the model with the effect of the intervention, with 
associated uncertainties.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- beta_range = range of betas that are to be tested
- params2 = array of other necessary parameters for the model after intervention, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
"""
function Level4_1_plot_intervention_with_error(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range)
    # Initialize variables to track maximum and minimum peaks
    max_infected_peak = min_infected_peak = 0
    max_seriously_infected_peak = min_seriously_infected_peak = 0

    # Plot initial vertical line indicating intervention time
    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, colour=:black)

    # Define ranges for ratio and beta
    ratios, betas = range(ratio_range[1], ratio_range[2], 20), range(beta_range[1], beta_range[2], length=20)

    # Initialize arrays for storing results
    infected, seriously_infected, time = Float64[], Float64[], Float64[]

    # Loop through the range of ratios
    for r in ratios
        params.SIratio = r
        params2.SIratio = r
        previous_infected, previous_seriously_infected = Float64[], Float64[]

        # Loop through range of betas
        for b in betas
            params.beta = b  # Update beta in parameters
            solution = solve_SIR(S0, I0, SI0, R0, days[1], params)

            # Reset arrays for infected and seriously infected
            infected, seriously_infected, time = Float64[], Float64[], Float64[]

            # Populate arrays with data from the first simulation phase
            for i = 1:length(solution.t)
                push!(infected, solution.u[i][2])
                push!(seriously_infected, solution.u[i][3])
                push!(time, solution.t[i])
            end

            # Simulate second phase with post-intervention parameters
            params2.beta = b
            solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

            # Append post-intervention data
            for i = 1:length(solution.t)
                push!(infected, solution.u[i][2])
                push!(seriously_infected, solution.u[i][3])
                push!(time, solution.t[i] + days[1])
            end

            # Handle previous data to plot shaded error regions
            if isempty(previous_infected)
                previous_infected = deepcopy(infected)
            end
            if isempty(previous_seriously_infected)
                previous_seriously_infected = deepcopy(seriously_infected)
            end

            # Ensure lengths of previous and current arrays match
            while length(previous_infected) < length(infected)
                push!(previous_infected, previous_infected[end])
            end
            while length(previous_seriously_infected) < length(seriously_infected)
                push!(previous_seriously_infected, previous_seriously_infected[end])
            end

            # Plot current results with shaded regions
            plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, color=:blue, fillrange=previous_infected)
            plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, color=:red, fillrange=previous_seriously_infected)

            # Update peak values for infected and seriously infected
            if max_infected_peak < maximum(infected)
                max_infected_peak = maximum(infected)
            end
            if max_seriously_infected_peak < maximum(seriously_infected)
                max_seriously_infected_peak = maximum(seriously_infected)
            end
            if min_infected_peak > maximum(infected) || min_infected_peak == 0
                min_infected_peak = maximum(infected)
            end
            if min_seriously_infected_peak > maximum(seriously_infected) || min_seriously_infected_peak == 0
                min_seriously_infected_peak = maximum(seriously_infected)
            end

            # Update previous arrays for the next iteration
            previous_infected, previous_seriously_infected = deepcopy(infected), deepcopy(seriously_infected)
        end
    end

    # Display peak values
    println("Peak Infected: ", round(min_infected_peak, digits=0), " - ", round(max_infected_peak, digits=0))
    println("Peak Seriously Infected: ", round(min_seriously_infected_peak, digits=0), " - ", round(max_seriously_infected_peak, digits=0))

    # Final plot with labels
    plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels="Infected", color=:blue)
    plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels="Seriously Infected", color=:red)
end


"""
Level4_1_compare_intervention is a function that plots the behaviour of infected or 
seriously infected curves from the model with or without the effect of intervention, 
along with associated uncertainties.

Inputs:
- S0: Initial Susceptible Population
- I0: Initial Infected Population
- SI0: Initial Seriously Infected Population
- R0: Initial Recovered Population
- days: Number of days to model 
- params: Parameters for the model before intervention, defined by the structs
- beta_range: Range of beta values to test
- params2: Parameters for the model after intervention, defined by the structs
- ratio_range: Range of proportions of people that become seriously infected
- infect: Determines whether to plot infected or seriously infected data (1 = infected, otherwise seriously infected)
- intervene: Determines whether to plot the intervention or no intervention curve (1 = intervention, otherwise no intervention)
"""
function Level4_1_compare_intervention(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, infect, intervene)
    # Initialize variables to track peaks
    max_infected_peak = min_infected_peak = max_seriously_infected_peak = min_seriously_infected_peak = 0
    max_infected_peak_intervention = min_infected_peak_intervention = max_seriously_infected_peak_intervention = min_seriously_infected_peak_intervention = 0

    # Generate ranges for ratios and betas
    ratios, betas = range(ratio_range[1], ratio_range[2], length=30), range(beta_range[1], beta_range[2], length=20)
    
    # Initialize plot
    plot([31, 31], [0, 250], xlabel="Time (Days)", ylabel="Population", title=" ", labels=nothing, colour=:black)

    if intervene != 1
        for r in ratios
            params.SIratio = r
            previous_infected = previous_seriously_infected = Float64[]

            # Loop through range of betas
            for b in betas
                params.beta = b # Update beta value
                solution = solve_SIR(S0, I0, SI0, R0, (days[1] + days[2]), params)

                # Extract data for infected and seriously infected
                infected, seriously_infected = Float64[], Float64[]
                for i = 1:length(solution.t)
                    push!(infected, solution.u[i][2])
                    push!(seriously_infected, solution.u[i][3])
                end

                # Ensure previous data arrays are initialized and match length
                if isempty(previous_infected)
                    previous_infected = deepcopy(infected)
                end
                if isempty(previous_seriously_infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end

                # Ensure lengths of previous and current arrays match
                while length(previous_infected) < length(infected)
                    push!(previous_infected, previous_infected[end])
                end
                while length(previous_seriously_infected) < length(seriously_infected)
                    push!(previous_seriously_infected, previous_seriously_infected[end])
                end

                # Plot curves
                if infect == 1
                    plot!(solution.t, infected, xlabel="Time (Days)", ylabel="Population", title="Infected No Intervention", labels=nothing, color=:red, fillrange=previous_infected)
                else
                    plot!(solution.t, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected No Intervention", labels=nothing, color=:red, fillrange=previous_seriously_infected)
                end

                # Update peak values
                if max_infected_peak < maximum(infected)
                    max_infected_peak = maximum(infected)
                end
                if max_seriously_infected_peak < maximum(seriously_infected)
                    max_seriously_infected_peak = maximum(seriously_infected)
                end
                if min_infected_peak > maximum(infected) || min_infected_peak == 0
                    min_infected_peak = maximum(infected)
                end
                if min_seriously_infected_peak > maximum(seriously_infected) || min_seriously_infected_peak == 0
                    min_seriously_infected_peak = maximum(seriously_infected)
                end

                # Update previous arrays
                previous_infected, previous_seriously_infected = deepcopy(infected), deepcopy(seriously_infected)
            end
        end
    else
        for r in ratios
            params.SIratio, params2.SIratio = r, r
            previous_infected = previous_seriously_infected = Float64[]

            # Loop through range of betas
            for b in betas
                params.beta = b
                solution = solve_SIR(S0, I0, SI0, R0, days[1], params)

                # Extract data
                infected, seriously_infected, time = Float64[], Float64[], Float64[]
                for i = 1:length(solution.t)
                    push!(infected, solution.u[i][2])
                    push!(seriously_infected, solution.u[i][3])
                    push!(time, solution.t[i])
                end

                params2.beta = b
                solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

                for i = 1:length(solution.t)
                    push!(infected, solution.u[i][2])
                    push!(seriously_infected, solution.u[i][3])
                    push!(time, solution.t[i] + days[1])
                end

                # Ensure previous data arrays are initialized and match length
                if isempty(previous_infected)
                    previous_infected = deepcopy(infected)
                end
                if isempty(previous_seriously_infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end

                # Ensure lengths of previous and current arrays match
                while length(previous_infected) < length(infected)
                    push!(previous_infected, previous_infected[end])
                end
                while length(previous_seriously_infected) < length(seriously_infected)
                    push!(previous_seriously_infected, previous_seriously_infected[end])
                end

                # Plot curves
                if infect == 1
                    plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected Intervention", labels=nothing, color=:blue, fillrange=previous_infected)
                else
                    plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected Intervention", labels=nothing, color=:blue, fillrange=previous_seriously_infected)
                end

                # Update peak values
                if max_infected_peak_intervention < maximum(infected)
                    max_infected_peak_intervention = maximum(infected)
                end
                if max_seriously_infected_peak_intervention < maximum(seriously_infected)
                    max_seriously_infected_peak_intervention = maximum(seriously_infected)
                end
                if min_infected_peak_intervention > maximum(infected) || min_infected_peak_intervention == 0
                    min_infected_peak_intervention = maximum(infected)
                end
                if min_seriously_infected_peak_intervention > maximum(seriously_infected) || min_seriously_infected_peak_intervention == 0
                    min_seriously_infected_peak_intervention = maximum(seriously_infected)
                end

                # Update previous arrays
                previous_infected, previous_seriously_infected = deepcopy(infected), deepcopy(seriously_infected)
            end
        end
    end

    # Print peak statistics
    if infect == 1
        if intervene != 1
            println("Peak Infected No Intervention: ", round(min_infected_peak, digits=0), " - ", round(max_infected_peak, digits=0))
        else
            println("Peak Infected Intervention: ", round(min_infected_peak_intervention, digits=0), " - ", round(max_infected_peak_intervention, digits=0))
        end
    else
        if intervene != 1
            println("Peak Seriously Infected No Intervention: ", round(min_seriously_infected_peak, digits=0), " - ", round(max_seriously_infected_peak, digits=0))
        else
            println("Peak Seriously Infected Intervention: ", round(min_seriously_infected_peak_intervention, digits=0), " - ", round(max_seriously_infected_peak_intervention, digits=0))
        end
    end

    plot!([0], [0], labels=nothing)
end

"""
Level4_2_plot_coverage is a function that plots the infected or seriously infected model curves 
with new coverage against the actual data. It incorporates the uncertainty of beta 
and ratio.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modeled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- beta_range = range of transmission rates of disease
- params2 = array of other necessary parameters for the model after intervention, defined by the structs, with new coverage
- ratio_range = range of proportions of people that become seriously infected
- infect = determines whether to plot infected or seriously infected data (1 = infected, else = seriously infected)
"""
function Level4_2_plot_coverage(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, infect)
    # Actual data for comparison (infected and seriously infected populations)
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119]
    tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]

    # Initialize the plot
    plot([0],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing)

    # Define ranges for ratios and beta values
    ratios = range(ratio_range[1], ratio_range[2], 20)
    betas = range(beta_range[1], beta_range[2], 50)

    # Initialize variables for storing simulation data
    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[]
    previous_infected = Float64[]
    previous_seriously_infected = Float64[]

    # Loop through the range of SI ratios
    for r in ratios
        params.SIratio = r  # Update SI ratio in pre-intervention parameters
        params2.SIratio = r # Update SI ratio in post-intervention parameters

        # Loop through the range of beta values
        for b in betas
            params.beta = b # Update beta value in pre-intervention parameters
            solution = solve_SIR(S0, I0, SI0, R0, days[1], params) # Solve the SIR model

            # Clear data arrays for current beta
            infected = Float64[]
            seriously_infected = Float64[]
            time = Float64[]

            # Extract data from the pre-intervention solution
            for i = 1:length(solution.t)
                push!(infected, solution.u[i][2])
                push!(seriously_infected, solution.u[i][3])
                push!(time, solution.t[i])
            end

            # Update post-intervention parameters and solve again
            params2.beta = b
            solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

            # Extract data from the post-intervention solution
            for i = 1:length(solution.t)
                push!(infected, solution.u[i][2])
                push!(seriously_infected, solution.u[i][3])
                push!(time, solution.t[i] + days[1])
            end

            # Initialize previous data arrays if empty
            if isempty(previous_infected)
                previous_infected = deepcopy(infected)
            end
            if isempty(previous_seriously_infected)
                previous_seriously_infected = deepcopy(seriously_infected)
            end

            # Match lengths of previous and current data arrays
            while length(previous_infected) < length(infected)
                push!(previous_infected, previous_infected[end])
            end
            while length(previous_seriously_infected) < length(seriously_infected)
                push!(previous_seriously_infected, previous_seriously_infected[end])
            end

            # Plot the curves
            if infect == 1
                plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels=nothing, fillrange=previous_infected, colour=:blue)
            else
                plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels=nothing, fillrange=previous_seriously_infected, colour=:blue)
            end

            # Update previous data arrays
            previous_seriously_infected = deepcopy(seriously_infected)
            previous_infected = deepcopy(infected)
        end
    end

    # Plot the actual data for comparison
    if infect == 1
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Infected", labels="Infected", colour=:blue)
        plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected", colour=:red)
    else
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected", colour=:blue)
        plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected", colour=:red)
    end
end


"""
Level4_2_error_coverage is a function that plots the RMSE of a range of coverage values
using the predicted and actual data provided. It incorporates the uncertainty
of beta and ratio, and can be plotted logarithmically to better highlight errors.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modeled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- coverage_range = range of coverage values to be tested
- params2 = array of other necessary parameters for the model after intervention, defined by the structs
- beta_range = range of transmission rates of the disease
- ratio_range = range of proportions of people that become seriously infected
- logarithm = plots the error logarithmically (1 = log, else = not log)
"""
function Level4_2_error_coverage(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, logarithm)
    # Actual infected and seriously infected data for comparison
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119]
    tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]

    # Define ranges for coverage, beta, and SI ratios
    coverages = range(coverage_range[1], coverage_range[2], length=50)
    betas = range(beta_range[1], beta_range[2], length=40)
    ratios = range(ratio_range[1], ratio_range[2], length=10)

    # Initialize plot for error visualization
    plot([coverage_range[1]],[0], xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing)

    # Variables to track the minimum error and corresponding coverage
    coverage_min = 1
    min = S0

    # Loop through SI ratios
    for r in ratios
        params.SIratio = r        # Update SI ratio for pre-intervention parameters
        params2.SIratio = r       # Update SI ratio for post-intervention parameters
        previous_error_vals = Float64[]  # Array for previous error values
        error_vals = Float64[]           # Array for current error values

        # Loop through beta values
        for b in betas
            params.beta = b       # Update beta for pre-intervention parameters
            params2.beta = b      # Update beta for post-intervention parameters
            error_vals = Float64[]  # Reset error values for current beta

            # Loop through coverage values
            for j in coverages
                # Solve the pre-intervention SIR model
                solution = solve_SIR(S0, I0, SI0, R0, days[1], params)
                infected = Float64[]                # Store infected values
                seriously_infected = Float64[]      # Store seriously infected values
                time = Float64[]                    # Store time values

                # Extract results from pre-intervention solution
                for i = 1:length(solution.t)
                    push!(infected, solution.u[i][2])
                    push!(seriously_infected, solution.u[i][3])
                    push!(time, solution.t[i])
                end

                # Solve the post-intervention SIR model
                params2.p = j  # Update coverage value
                solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

                # Extract results from post-intervention solution
                for i = 1:length(solution.t)
                    push!(infected, solution.u[i][2])
                    push!(seriously_infected, solution.u[i][3])
                    push!(time, solution.t[i] + days[1])
                end

                # Initialize error calculations
                current_infected_error = 0
                current_seriously_infected_error = 0

                # Calculate errors by comparing model predictions to actual data
                for j = 1:length(time)
                    for k = 1:length(ti)
                        if ti[k] ≈ time[j] atol=0.01  # Approximate time match for infected
                            current_infected_error += (infected[j] - actual_infected[k])^2
                        end
                    end
                    for k = 1:length(tsi)
                        if tsi[k] ≈ time[j] atol=0.01  # Approximate time match for seriously infected
                            current_seriously_infected_error += (seriously_infected[j] - actual_seriously_infected[k])^2
                        end
                    end
                end

                # Calculate RMSE for current coverage
                if logarithm == 1
                    push!(error_vals, log(sqrt(current_infected_error + current_seriously_infected_error)))
                else
                    push!(error_vals, sqrt(current_infected_error + current_seriously_infected_error))
                end
            end

            # Update minimum error and corresponding coverage
            if minimum(error_vals) < min
                min = minimum(error_vals)
                coverage_min = coverages[argmin(error_vals)]
            end

            # Handle initialization of previous_error_vals
            if isempty(previous_error_vals)
                previous_error_vals = deepcopy(error_vals)
            end

            # Ensure previous_error_vals matches current error length
            while length(previous_error_vals) < length(error_vals)
                push!(previous_error_vals, previous_error_vals[end])
            end

            # Plot error values for current beta
            plot!(coverages, error_vals, xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing, fillrange=previous_error_vals, colour=:blue)

            # Update previous_error_vals
            previous_error_vals = deepcopy(error_vals)
        end
    end

    # Output minimum coverage value and plot final error curve
    println("Coverage min: ", coverage_min)
    plot!([0], [0], xlabel="Coverage", ylabel="Error", title="Coverage vs Error", labels=nothing)
end

"""
Level4_2_plot_coverage_with_error is a function to plot the effect of 
different coverage values on the infected and seriously infected populations,
including error margins based on varying beta and ratio values.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modeled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- beta_range = range of transmission rates of the disease
- params2 = array of other necessary parameters for the model after intervention, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- coverage_range = range of coverage values to be tested
"""
function Level4_2_plot_coverage_with_error(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range)

    # Initialize variables to track the peak values for infected and seriously infected populations
    max_infected_peak = 0
    min_infected_peak = 0
    max_seriously_infected_peak = 0
    min_seriously_infected_peak = 0

    # Initial plot setup with black color and axis labels
    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, colour=:black)

    # Define ranges for ratio, beta, and coverage values to iterate over
    ratios = range(ratio_range[1], ratio_range[2], 10)
    betas = range(beta_range[1], beta_range[2], length=30)
    coverages = range(coverage_range[1], coverage_range[2], length=30)

    # Initialize empty arrays to store infected, seriously infected populations, and time
    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[] 

    # Loop over different coverage values
    for c in coverages
        params2.p = c  # Update coverage parameter in params2
        # Loop over different ratios
        for r in ratios
            params.SIratio = r  # Set the SI ratio for params
            params2.SIratio = r  # Set the SI ratio for params2
            previous_infected = Float64[]  # Array to store the previous infected population
            previous_seriously_infected = Float64[]  # Array to store the previous seriously infected population

            # Loop over different beta values
            for b in betas
                params.beta = b  # Update the current beta value

                # Solve SIR model for the first phase (before intervention)
                solution = solve_SIR(S0, I0, SI0, R0, days[1], params)

                # Initialize empty arrays for infected, seriously infected, and time data
                infected = Float64[]
                seriously_infected = Float64[]
                time = Float64[]

                # Store the infected and seriously infected populations over time for the first phase
                for i = 1:length(solution.t)
                    push!(infected, solution.u[i][2])  # Store infected population
                    push!(seriously_infected, solution.u[i][3])  # Store seriously infected population
                    push!(time, solution.t[i])  # Store time
                end

                # Solve SIR model for the second phase (after intervention)
                params2.beta = b  # Update the beta value in params2
                solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

                # Store the infected and seriously infected populations over time for the second phase
                for i = 1:length(solution.t)
                    push!(infected, solution.u[i][2])  # Store infected population
                    push!(seriously_infected, solution.u[i][3])  # Store seriously infected population
                    push!(time, solution.t[i] + days[1])  # Adjust time for second phase
                end

                # If previous data arrays are empty, initialize them with current data
                if isempty(previous_infected)
                    previous_infected = deepcopy(infected)
                end
                if isempty(previous_seriously_infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end

                # Extend previous data arrays to match the current data length for smooth plotting
                while length(previous_infected) < length(infected)
                    push!(previous_infected, previous_infected[end])
                end
                while length(previous_seriously_infected) < length(seriously_infected)
                    push!(previous_seriously_infected, previous_seriously_infected[end])
                end

                # Plot the infected and seriously infected populations with error range (filled area)
                plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, color=:blue, fillrange=previous_infected)
                plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, color=:red, fillrange=previous_seriously_infected)

                # Track peak values for infected and seriously infected populations
                if max_infected_peak < maximum(infected)
                    max_infected_peak = maximum(infected)
                end
                if max_seriously_infected_peak < maximum(seriously_infected)
                    max_seriously_infected_peak = maximum(seriously_infected)
                end
                if min_infected_peak > maximum(infected) || min_infected_peak == 0
                    min_infected_peak = maximum(infected)
                end
                if min_seriously_infected_peak > maximum(seriously_infected) || min_seriously_infected_peak == 0
                    min_seriously_infected_peak = maximum(seriously_infected)
                end

                # Update previous data for the next iteration
                previous_infected = deepcopy(infected)
                previous_seriously_infected = deepcopy(seriously_infected)
            end
        end
    end

    # Print the range of peak values for infected and seriously infected populations
    println("Peak Infected: ", round(min_infected_peak, digits=0), " - ", round(max_infected_peak, digits=0))
    println("Peak Seriously Infected: ", round(min_seriously_infected_peak, digits=0), " - ", round(max_seriously_infected_peak, digits=0))

    # Final plot with labels for infected and seriously infected populations
    plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels="Infected", color=:blue)
    plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels="Seriously Infected", color=:red)
end


"""
Level4_2_compare_intervention_coverage is a function to compare the effect of intervention 
and coverage on the infected and seriously infected populations. It simulates and plots 
infection dynamics with and without intervention, across different parameter ranges.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modeled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- beta_range = range of transmission rates of the disease
- params2 = array of other necessary parameters for the model after intervention, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- coverage_range = range of coverage values to be tested
- infect: Determines whether to plot infected or seriously infected data (1 = infected, otherwise seriously infected)
- intervene: Determines whether to plot the intervention or no intervention curve (1 = intervention, otherwise no intervention)
"""
function Level4_2_compare_intervention_coverage(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, infect, intervene)
    
    # Initialize variables to track peak values for infected and seriously infected populations
    max_infected_peak = 0
    min_infected_peak = 0
    max_seriously_infected_peak = 0
    min_seriously_infected_peak = 0
    max_infected_peak_intervention = 0
    min_infected_peak_intervention = 0
    max_seriously_infected_peak_intervention = 0
    min_seriously_infected_peak_intervention = 0

    # Define ranges for ratio, beta, and coverage values to iterate over
    ratios = range(ratio_range[1], ratio_range[2], length=20)
    betas = range(beta_range[1], beta_range[2], length=20)
    coverages = range(coverage_range[1], coverage_range[2], length=15)
    
    # Initialize an empty plot to hold the results
    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title=" ", labels=nothing, colour=:black)

    # Case when intervention is NOT applied
    if intervene != 1
        # Loop through different coverage values
        for c in coverages
            params2.p = c  # Update coverage value in params
            # Loop through different ratios
            for r in ratios
                params.SIratio = r  # Update SI ratio in params
                # Initialize empty arrays for storing previous infected and seriously infected data for smoothing
                previous_infected = Float64[]
                previous_seriously_infected = Float64[]
                
                # Loop through different beta values
                for b in betas
                    params.beta = b  # Set the current beta value
                    solution = solve_SIR(S0, I0, SI0, R0, days[1]+days[2], params)  # Solve SIR model with no intervention

                    # Arrays to store infected and seriously infected populations over time
                    infected = Float64[]
                    seriously_infected = Float64[]

                    # Store the infected and seriously infected data over time
                    for i = 1:length(solution.t)
                        push!(infected, solution.u[i][2])  # Store infected population
                        push!(seriously_infected, solution.u[i][3])  # Store seriously infected population
                    end

                    # Initialize previous data if it's empty
                    if isempty(previous_infected)
                        previous_infected = deepcopy(infected)
                    end
                    if isempty(previous_seriously_infected)
                        previous_seriously_infected = deepcopy(seriously_infected)
                    end

                    # Extend previous data arrays to match current data length for smooth plotting
                    while length(previous_infected) < length(infected)
                        push!(previous_infected, previous_infected[end])
                    end
                    while length(previous_seriously_infected) < length(seriously_infected)
                        push!(previous_seriously_infected, previous_seriously_infected[end])
                    end

                    # Plot the infected or seriously infected data
                    if infect == 1
                        plot!(solution.t, infected, xlabel="Time (Days)", ylabel="Population", title="Infected No Intervention", labels=nothing, color=:red, fillrange=previous_infected)
                    else
                        plot!(solution.t, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected No Intervention", labels=nothing, color=:red, fillrange=previous_seriously_infected)
                    end

                    # Track peak values for infected and seriously infected populations
                    if max_infected_peak < maximum(infected)
                        max_infected_peak = maximum(infected)
                    end
                    if max_seriously_infected_peak < maximum(seriously_infected)
                        max_seriously_infected_peak = maximum(seriously_infected)
                    end
                    if min_infected_peak > maximum(infected) || min_infected_peak == 0
                        min_infected_peak = maximum(infected)
                    end
                    if min_seriously_infected_peak > maximum(seriously_infected) || min_seriously_infected_peak == 0
                        min_seriously_infected_peak = maximum(seriously_infected)
                    end

                    # Store the current infected and seriously infected data for the next iteration
                    previous_infected = deepcopy(infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end
            end
        end

    # Case when intervention IS applied
    else
        # Loop through different coverage values
        for c in coverages
            params2.p = c  # Update coverage value in params
            # Loop through different ratios
            for r in ratios
                params.SIratio = r  # Update SI ratio in params
                params2.SIratio = r  # Update SI ratio in params2
                previous_infected = Float64[]  # To store previous infected data for smoothing
                previous_seriously_infected = Float64[]  # To store previous seriously infected data for smoothing

                # Loop through different beta values
                for b in betas
                    params.beta = b  # Set the current beta value
                    solution = solve_SIR(S0, I0, SI0, R0, days[1], params)  # Solve SIR model for the first phase of intervention

                    # Arrays to store infected, seriously infected, and time data
                    infected = Float64[]
                    seriously_infected = Float64[]
                    time = Float64[]

                    # Store the infected and seriously infected data over time for the first phase
                    for i = 1:length(solution.t)
                        push!(infected, solution.u[i][2])  # Store infected population
                        push!(seriously_infected, solution.u[i][3])  # Store seriously infected population
                        push!(time, solution.t[i])  # Store time
                    end

                    # Solve for the second phase of intervention and update params2
                    params2.beta = b
                    solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

                    # Store the data for the second phase of the model
                    for i = 1:length(solution.t)
                        push!(infected, solution.u[i][2])  # Store infected population
                        push!(seriously_infected, solution.u[i][3])  # Store seriously infected population
                        push!(time, solution.t[i] + days[1])  # Adjust time for the second phase
                    end

                    # Initialize previous data if it's empty
                    if isempty(previous_infected)
                        previous_infected = deepcopy(infected)
                    end
                    if isempty(previous_seriously_infected)
                        previous_seriously_infected = deepcopy(seriously_infected)
                    end

                    # Extend previous data arrays to match current data length for smooth plotting
                    while length(previous_infected) < length(infected)
                        push!(previous_infected, previous_infected[end])
                    end
                    while length(previous_seriously_infected) < length(seriously_infected)
                        push!(previous_seriously_infected, previous_seriously_infected[end])
                    end

                    # Plot the infected or seriously infected data
                    if infect == 1
                        plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected Intervention", labels=nothing, color=:blue, fillrange=previous_infected)
                    else
                        plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected Intervention", labels=nothing, color=:blue, fillrange=previous_seriously_infected)
                    end

                    # Track peak values for infected and seriously infected populations with intervention
                    if max_infected_peak_intervention < maximum(infected)
                        max_infected_peak_intervention = maximum(infected)
                    end
                    if max_seriously_infected_peak_intervention < maximum(seriously_infected)
                        max_seriously_infected_peak_intervention = maximum(seriously_infected)
                    end
                    if min_infected_peak_intervention > maximum(infected) || min_infected_peak_intervention == 0
                        min_infected_peak_intervention = maximum(infected)
                    end
                    if min_seriously_infected_peak_intervention > maximum(seriously_infected) || min_seriously_infected_peak_intervention == 0
                        min_seriously_infected_peak_intervention = maximum(seriously_infected)
                    end

                    # Store the current infected and seriously infected data for the next iteration
                    previous_infected = deepcopy(infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end
            end
        end
    end

    # Print the peak values for infected and seriously infected populations based on the intervention status
    if infect == 1
        if intervene != 1
            println("Peak Infected No Intervention: ", round(min_infected_peak, digits=0), " - ", round(max_infected_peak, digits=0))
        else
            println("Peak Infected Intervention: ", round(min_infected_peak_intervention, digits=0), " - ", round(max_infected_peak_intervention, digits=0))
        end
    else
        if intervene != 1
            println("Peak Seriously Infected No Intervention: ", round(min_seriously_infected_peak, digits=0), " - ", round(max_seriously_infected_peak, digits=0))
        else
            println("Peak Seriously Infected Intervention: ", round(min_seriously_infected_peak_intervention, digits=0), " - ", round(max_seriously_infected_peak_intervention, digits=0))
        end
    end

    # Finalize the plot with an empty series to complete rendering
    plot!([0], [0], labels=nothing)
end

"""
Level5_best_beta_day_plot is a find the beta and start day by minimising the erro of test values.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modeled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- beta_range = range of transmission rates of the disease
- params2 = array of other necessary parameters for the model after intervention, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- coverage_range = range of coverage values to be tested
- town = which town is being worked on
- day_range = the range of start days to test
"""
function Level5_best_beta_day_plot(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, town, day_range)
    # Actual data up to day 80
    actual_infected = [21, 29, 25, 30, 28, 34, 28, 54, 57,92,73,80,109,102,128,135,163,150,211,196,233,247,283,286,332,371,390,404,467,529,598,
    641,704,702,788,856,854,955,995,1065,1106,1159,1217,1269,1298,1328,1339,1383,1431,1422,1414,1485,1464,1480]
    ti = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    actual_seriously_infected = [3, 3, 4, 7, 3, 8, 7, 5, 9,13,15,3,20,13,11,20,16,11,15,18,27,24,28,36,41,35,41,55,63,66,72,80,90,104,109,
    115,127,135,147,162,163,186,194,200,216,223,241,249,258,275,277,299,302,300]
    tsi = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    
    if town == 1
        # Actual data up to day 55 for town 1
        actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352,385,221,420,544,329,440,427,369,606,416,546,475,617,593,352,337,473,673,653,523,602,551,686,556,600]
        ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
        actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119,76,95,28,109,136,119,104,121,93,147,129,130,161,133,136,138,139,181,181,218,183,167,164,219,220]
        tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    end

    # Initialize ranges for beta, days, coverage, and ratio
    betas = range(beta_range[1], beta_range[2], length=40)   # Beta values to test
    start_days = range(day_range[1], day_range[2], length=day_range[2]-day_range[1]+1)  # Days to start intervention
    coverages = range(coverage_range[1], coverage_range[2], length=20)   # Coverage levels
    ratios = range(ratio_range[1], ratio_range[2], length=20)   # Ratio range for SI (Seriously Infected)

    plot([beta_range[1]],[0], labels=nothing)   # Initial empty plot

    # Track minimum error and corresponding beta and day
    min_error = Inf   # Start with infinite error
    best_beta = 0
    best_day = 0

    # Initialize array to store error for plotting
    beta_errors = Float64[Inf]  # Track the errors for different beta values

    # Loop through each coverage and ratio to test
    for c in coverages
        params2.p = c  # Set coverage value for the second set of parameters
        for r in ratios
            params.SIratio = r  # Set SI ratio for the first set of parameters
            params2.SIratio = r  # Set SI ratio for the second set of parameters
            for start_day in start_days
                error_vals = Float64[]  # Store the error values for each beta
                for b in betas
                    params.beta = b  # Set current beta value
                    params2.beta = b  # Set beta value for the second set of parameters

                    # Solve SIR model for current beta, start day, coverage, and ratio
                    solution = solve_SIR(S0, I0, SI0, R0, days[1] - (start_day), params)

                    infected = Float64[]  # Track infected population
                    seriously_infected = Float64[]  # Track seriously infected population
                    time = Float64[]  # Track the time steps

                    # Collect results from the first phase
                    for i = 1:length(solution.t)
                        push!(infected, solution.u[i][2])
                        push!(seriously_infected, solution.u[i][3])
                        push!(time, solution.t[i]+start_day)
                    end

                    # Re-solve for second phase (post-intervention)
                    solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

                    for i = 1:length(solution.t)
                        push!(infected, solution.u[i][2])
                        push!(seriously_infected, solution.u[i][3])
                        push!(time, solution.t[i] + days[1])
                    end

                    # Calculate error (RMSE) for the model
                    current_infected_error = 0
                    current_seriously_infected_error = 0

                    for j = 1:length(time)
                        for k = 1:length(ti)  # Compare infected with actual data
                            if abs(ti[k] - time[j]) < 0.01
                                current_infected_error += (infected[j] - actual_infected[k])^2
                            end
                        end
                        for k = 1:length(tsi)  # Compare seriously infected with actual data
                            if abs(tsi[k] - time[j]) < 0.01
                                current_seriously_infected_error += (seriously_infected[j] - actual_seriously_infected[k])^2
                            end
                        end
                    end

                    push!(error_vals, sqrt(current_infected_error + current_seriously_infected_error))  # Calculate RMSE at each coverage value

                    # Track minimum error across beta and days
                    if sqrt(current_infected_error + current_seriously_infected_error) < min_error
                        min_error = sqrt(current_infected_error + current_seriously_infected_error)
                        best_beta = b
                        best_day = start_day
                    end
                end
                plot!([betas], [error_vals], labels=nothing, colour=:blue)  # Plot error for each beta value
                if minimum(error_vals) < minimum(beta_errors)
                    beta_errors = error_vals
                end
            end
        end
    end

    # Plot error vs beta
    println("Best Beta: ", best_beta)
    println("Start Day: ", round(best_day,digits=0))  # Output best beta and day
    plot!([betas], [beta_errors], xlabel="Beta", ylabel="Error", title="Error vs Beta across Days", labels=nothing, colour=:red)  # Final plot showing error vs beta
end

"""
Level5_plot_second_town is a function to plot the behaviourof the system compared to actual data
and plot hypotheticals about the town.

Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modeled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- beta_range = range of transmission rates of the disease
- params2 = array of other necessary parameters for the model after intervention, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- coverage_range = range of coverage values to be tested
- start_day = day to start plotting from
- infect: Determines whether to plot infected or seriously infected data (1 = infected, otherwise seriously infected)
- plot_actual: plots the actual data (1 = plot, else = don't plot)
"""
function Level5_plot_second_town(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, start_day, town, infect, plot_actual)
    # Define actual data for infected population and seriously infected population up to day 80
    actual_infected = [21, 29, 25, 30, 28, 34, 28, 54, 57,92,73,80,109,102,128,135,163,150,211,196,233,247,283,286,332,371,390,404,467,529,598,
    641,704,702,788,856,854,955,995,1065,1106,1159,1217,1269,1298,1328,1339,1383,1431,1422,1414,1485,1464,1480]
    ti = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    actual_seriously_infected = [3, 3, 4, 7, 3, 8, 7, 5, 9,13,15,3,20,13,11,20,16,11,15,18,27,24,28,36,41,35,41,55,63,66,72,80,90,104,109,
    115,127,135,147,162,163,186,194,200,216,223,241,249,258,275,277,299,302,300]
    tsi = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]

    # Adjust the data depending on the town number
    if town == 1
        # Actual data for the first town up to day 55
        actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352,385,221,420,544,329,440,427,369,606,416,546,475,617,593,352,337,473,673,653,523,602,551,686,556,600]
        ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
        actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119,76,95,28,109,136,119,104,121,93,147,129,130,161,133,136,138,139,181,181,218,183,167,164,219,220]
        tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    end

    # Initialize the plot with an empty starting point
    plot([start_day],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing)

    # Define parameter ranges for ratios, beta values, and coverage
    ratios = range(ratio_range[1], ratio_range[2], 15)
    betas = range(beta_range[1], beta_range[2], 40)
    coverages = range(coverage_range[1], coverage_range[2],15)

    # Initialize arrays for the model output
    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[]

    # Loop through different coverage values
    for c in coverages
        params2.p = c
        # Loop through different ratio values
        for r in ratios
            params.SIratio = r
            params2.SIratio = r
            previous_infected = Float64[]
            previous_seriously_infected = Float64[]

            # Loop through different beta values
            for b in betas
                params.beta = b
                # Solve the SIR model for the first period
                solution = solve_SIR(S0, I0, SI0, R0, days[1]-(start_day), params)

                # Reset arrays for storing results
                infected = Float64[]
                seriously_infected = Float64[]
                time = Float64[]

                # Store the results from the first period of simulation
                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i]+(start_day))
                end

                # Solve the SIR model for the second period
                params2.beta = b
                solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

                # Store the results from the second period of simulation
                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i]+days[1])
                end

                # Ensure previous data arrays are aligned with current results
                if isempty(previous_infected)
                    previous_infected = deepcopy(infected)
                end
                if isempty(previous_seriously_infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end

                while length(previous_infected) < length(infected)
                    push!(previous_infected,previous_infected[end])
                end
                while length(previous_seriously_infected) < length(seriously_infected)
                    push!(previous_seriously_infected,previous_seriously_infected[end])
                end

                # Plot the results for infected or seriously infected population depending on the flag
                if infect == 1
                    plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels=nothing, fillrange=previous_infected, colour=:blue)
                else
                    plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels=nothing, fillrange=previous_seriously_infected, colour=:blue)
                end

                # Update the previous values for the next iteration
                previous_seriously_infected = deepcopy(seriously_infected)
                previous_infected = deepcopy(infected)
            end
        end
    end

    # Plot the actual data if requested
    if infect == 1
        if plot_actual == 1
            plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected", colour=:red)
        end
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Infected", labels="Infected", colour=:blue)
    else
        if plot_actual == 1
            plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected", colour=:red)
        end
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected", colour=:blue)
    end
end


"""
plot_SIR is a driver function that runs and then plots an SIR model with uncertainties on variables.
Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
- beta_range = range of transmission rates of disease
- ratio_range = range of proportions of people that become seriously infected
"""
function plot_SIR_with_uncertainties(S0, I0, SI0, R0, days, params, beta_range, ratio_range)
    plot([0],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing) # Open a base plot

    ratios = range(ratio_range[1], ratio_range[2], 10) # Create a range of values defined by the ratio range
    betas = range(beta_range[1], beta_range[2], 20) # Create a range of values defined by the beta range

    # Loop through ratios
    for r in ratios
        params.SIratio = r # Set new SIratio value
        previous_infected = Float64[]
        previous_seriously_infected = Float64[]
        previous_susceptible = Float64[]
        previous_recovered = Float64[]
        for b in betas
            params.beta = b
            solution = solve_SIR(S0, I0, SI0, R0, days, params) # Solve the SIR model

            # Ensure arrays are empty
            infected = Float64[]
            seriously_infected = Float64[]
            time = Float64[]
            susceptible = Float64[]
            recovered = Float64[]

            # Extract infected, seriously infected, and time graphs
            for i = 1:length(solution.t)
                push!(infected,solution.u[i][2])
                push!(seriously_infected,solution.u[i][3])
                push!(time, solution.t[i])
                push!(susceptible,solution.u[i][1])
                push!(recovered,solution.u[i][4])
            end

            # Ensure if fill range arrays aren't empty to avoid erros
            if isempty(previous_infected)
                previous_infected = deepcopy(infected)
            end
            if isempty(previous_seriously_infected)
                previous_seriously_infected = deepcopy(seriously_infected)
            end
            if isempty(susceptible)
                previous_susceptible = deepcopy(susceptible)
            end
            if isempty(recovered)
                previous_recovered = deepcopy(recovered)
            end

            # Ensure fill range arrays are long enough to fill total
            while length(previous_infected) < length(infected)
                push!(previous_infected,previous_infected[end])
            end
            while length(previous_seriously_infected) < length(seriously_infected)
                push!(previous_seriously_infected,previous_seriously_infected[end])
            end
            while length(previous_susceptible) < length(susceptible)
                push!(previous_susceptible,previous_susceptible[end])
            end
            while length(previous_seriously_infected) < length(recovered)
                push!(previous_recovered,previous_recovered[end])
            end

            # Determines whether to plot infected or seriously infected graph
            plot!(time, infected, xlabel="Time", ylabel="Population", labels=nothing, fillrange=previous_infected, colour=:yellow) # Added curve to plot 
            plot!(time, seriously_infected, xlabel="Time", ylabel="Population", labels=nothing, fillrange=previous_seriously_infected, colour=:red) # Added curve to plot
            plot!(time, susceptible, xlabel="Time", ylabel="Population", labels=nothing, fillrange=previous_seriously_infected, colour=:blue) # Added curve to plot
            plot!(time, recovered, xlabel="Time", ylabel="Population", labels=nothing, fillrange=previous_seriously_infected, colour=:green) # Added curve to plot 

            # Keep previous curve
            previous_seriously_infected = deepcopy(seriously_infected)
            previous_infected = deepcopy(infected)
            previous_susceptible = deepcopy(susceptible)
            previous_recovered = deepcopy(recovered)
        end
    end

    # Plot graph
    plot!([0], [0], xlabel="Time", ylabel="Population", title="Infection Behaviour", labels="Infected", colour=:yellow) # Plot the model
    plot!([0], [0], xlabel="Time", ylabel="Population", labels="Seriously Infected", colour=:red) # Plot the model
    plot!([0], [0], xlabel="Time", ylabel="Population", labels="Susceptible", colour=:blue) # Plot the model
    plot!([0], [0], xlabel="Time", ylabel="Population", labels="Recovered", colour=:green) # Plot the model
end

"""
plot_SIR_intervention_with_uncertainties is a function that runs and then plots an 
SIR model with an intervention and uncertainties.
Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
- params2 = array of other necessary parameters, defined by the structs
- beta_range = range of transmission rates of disease
- ratio_range = range of proportions of people that become seriously infected
- coverage_range = range of proportions of people that followed the intervention procedure
"""
function plot_SIR_intervention_with_uncertainties(S0, I0, SI0, R0, days, params, params2, beta_range, ratio_range, coverage_range)
    plot([0],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing) # Open a base plot

    ratios = range(ratio_range[1], ratio_range[2], 10) # Create a range of values defined by the ratio range
    betas = range(beta_range[1], beta_range[2], 20) # Create a range of values defined by the beta range
    coverages = range(coverage_range[1], coverage_range[2], 20) # Create a range of values defined by the beta range

    # Loop through ratios
    for c in coverages
        params2.p = c
        for r in ratios
            params.SIratio = r # Set new SIratio value
            params2.SIratio = r
            previous_infected = Float64[]
            previous_seriously_infected = Float64[]
            previous_susceptible = Float64[]
            previous_recovered = Float64[]
            for b in betas
                params.beta = b
                solution = solve_SIR(S0, I0, SI0, R0, days[1], params) # Solve the SIR model

                # Ensure arrays are empty
                infected = Float64[]
                seriously_infected = Float64[]
                time = Float64[]
                susceptible = Float64[]
                recovered = Float64[]

                # Extract infected, seriously infected, and time graphs
                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i])
                    push!(susceptible,solution.u[i][1])
                    push!(recovered,solution.u[i][4])
                end

                params2.beta = b
                solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2) # Solve the SIR model

                # Extract infected, seriously infected, and time graphs
                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i]+days[1])
                    push!(susceptible,solution.u[i][1])
                    push!(recovered,solution.u[i][4])
                end

                # Ensure if fill range arrays aren't empty to avoid erros
                if isempty(previous_infected)
                    previous_infected = deepcopy(infected)
                end
                if isempty(previous_seriously_infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end
                if isempty(susceptible)
                    previous_susceptible = deepcopy(susceptible)
                end
                if isempty(recovered)
                    previous_recovered = deepcopy(recovered)
                end

                # Ensure fill range arrays are long enough to fill total
                while length(previous_infected) < length(infected)
                    push!(previous_infected,previous_infected[end])
                end
                while length(previous_seriously_infected) < length(seriously_infected)
                    push!(previous_seriously_infected,previous_seriously_infected[end])
                end
                while length(previous_susceptible) < length(susceptible)
                    push!(previous_susceptible,previous_susceptible[end])
                end
                while length(previous_seriously_infected) < length(recovered)
                    push!(previous_recovered,previous_recovered[end])
                end

                # Determines whether to plot infected or seriously infected graph
                plot!(time, infected, xlabel="Time", ylabel="Population", labels=nothing, fillrange=previous_infected, colour=:yellow) # Added curve to plot 
                plot!(time, seriously_infected, xlabel="Time", ylabel="Population", labels=nothing, fillrange=previous_seriously_infected, colour=:red) # Added curve to plot
                plot!(time, susceptible, xlabel="Time", ylabel="Population", labels=nothing, fillrange=previous_seriously_infected, colour=:blue) # Added curve to plot
                plot!(time, recovered, xlabel="Time", ylabel="Population", labels=nothing, fillrange=previous_seriously_infected, colour=:green) # Added curve to plot 

                # Keep previous curve
                previous_seriously_infected = deepcopy(seriously_infected)
                previous_infected = deepcopy(infected)
                previous_susceptible = deepcopy(susceptible)
                previous_recovered = deepcopy(recovered)
            end
        end
    end

    # Plot graph
    plot!([0], [0], xlabel="Time", ylabel="Population", title="Infection Behaviour", labels="Infected", colour=:yellow) # Plot the model
    plot!([0], [0], xlabel="Time", ylabel="Population", labels="Seriously Infected", colour=:red) # Plot the model
    plot!([0], [0], xlabel="Time", ylabel="Population", labels="Susceptible", colour=:blue) # Plot the model
    plot!([0], [0], xlabel="Time", ylabel="Population", labels="Recovered", colour=:green) # Plot the model
end