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
version of the SIR model
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
version of the SIR model
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
inspection of the models fit. This function is for plotting the level 2 data.
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
    plot_infected(S0, I0, SI0, R0, days, params, infect, ratio_range, 2) # Use plotting function
end

"""
Level3_plot_infected is a driver function that plots either the infected or
seriously infected data against the true data to allow for visual
inspection of the models fit. This function is for plotting the level 3 data.
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
    plot_infected(S0, I0, SI0, R0, days, params, infect, ratio_range, 3) # Use Plotting Function
end

"""
plot_infected is a function that plots either the infected or
seriously infected data against the true data to allow for visual
inspection of the models fit.
Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- infect = determines whether to plot infected or seriously infected data (0 = infected, else = seriously infected)
- level = which TLT level is being addressed
"""
function plot_infected(S0, I0, SI0, R0, days, params, ratio_range, infect, level)
    # Actual data up to day 25
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58]
    ti = [15,16,17,18,19,20,21,22,23,24,25]
    actual_seriously_infected = [0,0,1,2,5]
    tsi = [21,22,23,24,25]
    if level == 3 # Select which TLT Level
        # Actual data up to day 30
        actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55]
        ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
        actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4]
        tsi = [21,22,23,24,25,26,27,28,29,30]
    end

    plot([0],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing) # Open a base plot

    ratios = range(ratio_range[1], ratio_range[2], 5) # Create a range of values defined by the ratio range

    # Setup arrays
    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[]
    previous_infected = Float64[]
    previous_seriously_infected = Float64[]

    # Loop through ratios
    for j in ratios
        params.SIratio = j # Set new SIratio value
        solution = solve_SIR(S0, I0, SI0, R0, days, params) # Solve the SIR model

        # Ensure arrays are empty
        infected = Float64[]
        seriously_infected = Float64[]
        time = Float64[]

        # Extract infected, seriously infected, and time graphs
        for i = 1:length(solution.t)
            push!(infected,solution.u[i][2])
            push!(seriously_infected,solution.u[i][3])
            push!(time, solution.t[i])
        end

        # Ensure if fill range arrays aren't empty to avoid erros
        if isempty(previous_infected)
            previous_infected = deepcopy(infected)
        end
        if isempty(previous_seriously_infected)
            previous_seriously_infected = deepcopy(seriously_infected)
        end

        # Ensure fill range array is long enough to fill total
        while length(previous_infected) < length(infected)
            push!(previous_infected,previous_infected[end])
        end
        while length(previous_seriously_infected) < length(seriously_infected)
            push!(previous_seriously_infected,previous_seriously_infected[end])
        end

        # Determines whether to plot infected or seriously infected graph
        if infect == 1
            plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels=nothing, fillrange=previous_infected, colour=:blue) # Added curve to plot 
        else
            plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels=nothing, fillrange=previous_seriously_infected, colour=:blue) # Added curve to plot 
        end

        # Keep previous curve
        previous_seriously_infected = deepcopy(seriously_infected)
        previous_infected = deepcopy(infected)
    end

    # Get and print reproduction number if level 2
    if level == 2
        # Print associated R0 value
        R0 = params.contacts*params.beta/params.gamma
        println("R0: ", round(R0, digits=3))
    end

    # Plot either infected or seriously infected graphs
    if infect == 1
        plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels="Infected", colour=:blue) # Plot the model
        plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected", colour=:red) # Plot the model
    else
        plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected", colour=:blue) # Plot the model
        plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected", colour=:red) # Plot the model
    end
end

"""
error_beta is a function that plots the RSME of a range of beta values
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
    # Actual data up to day 30
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4]
    tsi = [21,22,23,24,25,26,27,28,29,30]

    plot([beta_range[1]],[0], xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing)
    beta_min = 1
    min = S0

    # Create error array and range of betas to test
    previous_error_vals = Float64[]
    error_vals = Float64[]
    betas = range(beta_range[1], beta_range[2], length=50)
    ratios = range(ratio_range[1], ratio_range[2], 5)

    for l in ratios
        params.SIratio = l
        error_vals = Float64[]
        # Loops through range of betas
        for i in betas
            params.beta = i # Changes value to new beta
            solution = solve_SIR(S0, I0, SI0, R0, days, params)

            current_infected_error = 0
            current_seriously_infected_error = 0

            # Loop through each solution value
            for j = 1:length(solution.t)
                for k = 1:length(ti)
                    if ti[k] ≈ solution.t[j] atol=0.01 # Check whether time value as approx equal
                        current_infected_error += (solution[j][2] - actual_infected[k])^2 # Calculate error
                    end

                end

                for k = 1:length(tsi)
                    if tsi[k] ≈ solution.t[j] atol=0.01
                        current_seriously_infected_error += (solution[j][3] - actual_seriously_infected[k])^2
                    end
                end
            end
            push!(error_vals, sqrt(current_infected_error + current_seriously_infected_error)) # Calculate RSME at beta value
        end

        if minimum(error_vals) < min
            min = minimum(error_vals)
            beta_min = betas[argmin(error_vals)]
        end
        if isempty(previous_error_vals)
            previous_error_vals = deepcopy(error_vals)
        end

        while length(previous_error_vals) < length(error_vals)
            push!(previous_error_vals,previous_error_vals[end])
        end

        plot!(betas, error_vals, xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing, fillrange=previous_error_vals, colour=:blue)

        previous_error_vals = deepcopy(error_vals)
    end

    println("Beta min: ", beta_min)
    plot!(betas, error_vals, xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing, colour=:blue)
end

"""
plot_intervention_with_error is a function that plots the behaviour of the infected and 
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
- params2 = array of other necessary parametersfor the model after intervention, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
"""
function Level4_1_plot_intervention_with_error(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range)
    max_infected_peak = 0
    min_infected_peak = 0
    max_seriously_infected_peak = 0
    min_seriously_infected_peak = 0

    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, colour=:black)

    ratios = range(ratio_range[1], ratio_range[2], 20)
    betas = range(beta_range[1], beta_range[2], length=20)

    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[] 

    for r in ratios
        params.SIratio = r
        params2.SIratio = r
        previous_infected = Float64[]
        previous_seriously_infected = Float64[]

        # Loops through range of betas
        for j in betas
            params.beta = j # Changes value to new beta
            solution = solve_SIR(S0, I0, SI0, R0, days[1], params)

            infected = Float64[]
            seriously_infected = Float64[]
            time = Float64[]

            for i = 1:length(solution.t)
                push!(infected,solution.u[i][2])
                push!(seriously_infected,solution.u[i][3])
                push!(time, solution.t[i])
            end

            params2.beta = j
            solution = solve_SIR(solution.u[end][1],solution.u[end][2],solution.u[end][3],solution.u[end][4],days[2], params2)

            for i = 1:length(solution.t)
                push!(infected,solution.u[i][2])
                push!(seriously_infected,solution.u[i][3])
                push!(time, solution.t[i]+days[1])
            end

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

            plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, color=:blue, fillrange=previous_infected)
            plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, color=:red, fillrange=previous_seriously_infected)

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

            previous_infected = deepcopy(infected)
            previous_seriously_infected = deepcopy(seriously_infected)
        end
    end

    println("Peak Infected: ", round(min_infected_peak,digits=0), " - ", round(max_infected_peak,digits=0))
    println("Peak Seriously Infected: ", round(min_seriously_infected_peak,digits=0), " - ", round(max_seriously_infected_peak,digits=0))

    plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels="Infected", color=:blue)
    plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels="Seriously Infected", color=:red)
end

"""
compare_intervention is a function that plots the behaviour of the infected or 
seriously infected curves from the model with or without the effect of the intervention,
with associated uncertainties.
Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- beta_range = range of betas that are to be tested
- params2 = array of other necessary parametersfor the model after intervention, defined by the structs
- ratio_range = range of proportions of people that become seriously infected
- infect = determines whether to plot infected or seriously infected data (1 = infected, else = seriously infected)
- intervent = determines whether to plot the intervention or no intervention curve (1 = intervention, else = no intervention)
"""
function Level4_1_compare_intervention(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, infect, intervene)
    max_infected_peak = 0
    min_infected_peak = 0
    max_seriously_infected_peak = 0
    min_seriously_infected_peak = 0
    max_infected_peak_intervention = 0
    min_infected_peak_intervention = 0
    max_seriously_infected_peak_intervention = 0
    min_seriously_infected_peak_intervention = 0

    ratios = range(ratio_range[1], ratio_range[2], length=30)
    betas = range(beta_range[1], beta_range[2], length=20)
    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title=" ", labels=nothing, colour=:black)

    if intervene != 1
        for r in ratios
            params.SIratio = r
            previous_infected = Float64[]
            previous_seriously_infected = Float64[]

            # Loops through range of betas
            for j in betas
                params.beta = j # Changes value to new beta
                solution = solve_SIR(S0, I0, SI0, R0, (days[1]+days[2]), params)

                infected = Float64[]
                seriously_infected = Float64[]

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                end

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

                if infect == 1
                    plot!(solution.t, infected, xlabel="Time (Days)", ylabel="Population", title="Infected No Intervention", labels=nothing, color=:red, fillrange=previous_infected)
                else
                    plot!(solution.t, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected No Intervention", labels=nothing, color=:red, fillrange=previous_seriously_infected)
                end

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

                previous_infected = deepcopy(infected)
                previous_seriously_infected = deepcopy(seriously_infected)
            end
        end
    else
        for r in ratios
            params.SIratio = r
            params2.SIratio = r
            previous_infected = Float64[]
            previous_seriously_infected = Float64[]

            # Loops through range of betas
            for j in betas
                params.beta = j # Changes value to new beta
                solution = solve_SIR(S0, I0, SI0, R0, days[1], params)

                infected = Float64[]
                seriously_infected = Float64[]
                time = Float64[]

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i])
                end

                params2.beta = j
                solution = solve_SIR(solution.u[end][1],solution.u[end][2],solution.u[end][3],solution.u[end][4],days[2], params2)

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i]+days[1])
                end

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

                if infect == 1
                    plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected Intervention", labels=nothing, color=:blue, fillrange=previous_infected)
                else
                    plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected Intervention", labels=nothing, color=:blue, fillrange=previous_seriously_infected)
                end

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

                previous_infected = deepcopy(infected)
                previous_seriously_infected = deepcopy(seriously_infected)
            end
        end
    end

    if infect == 1
        if intervene != 1
            println("Peak Infected No Intervention: ", round(min_infected_peak,digits=0), " - ", round(max_infected_peak,digits=0))
        else
            println("Peak Infected Intervention: ", round(min_infected_peak_intervention,digits=0), " - ", round(max_infected_peak_intervention,digits=0))
        end
    else
        if intervene != 1
            println("Peak Seriously Infected No Intervention: ", round(min_seriously_infected_peak,digits=0), " - ", round(max_seriously_infected_peak,digits=0))
        else
            println("Peak Seriously Infected Intervention: ", round(min_seriously_infected_peak_intervention,digits=0), " - ", round(max_seriously_infected_peak_intervention,digits=0))
        end
    end
    plot!([0],[0],labels=nothing)
end

"""
plot_coverage is a function that plots the infected or seriously infected model curves 
with new coverage against the actual data. It is incorporates the uncertainty of beta 
and ratio.
Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- infect = determines whether to plot infected or seriously infected data (1 = infected, else = seriously infected)
- params2 = array of other necessary parametersfor the model after intervention, defined by the structs, with new coverage
- ratio_range = range of proportions of people that become seriously infected
- beta_range = range of transmission rates of disease
"""
function Level4_2_plot_coverage(S0, I0, SI0, R0, days, params,  beta_range, params2, ratio_range, infect)
    # Actual data up to day 55
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119]
    tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]

    plot([0],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing)

    ratios = range(ratio_range[1], ratio_range[2], 20)
    betas = range(beta_range[1], beta_range[2], 50)

    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[]
    previous_infected = Float64[]
    previous_seriously_infected = Float64[]

    for r in ratios
        params.SIratio = r
        params2.SIratio = r

        for b in betas
            params.beta = b
            solution = solve_SIR(S0, I0, SI0, R0, days[1], params) # Solve the SIR model

            infected = Float64[]
            seriously_infected = Float64[]
            time = Float64[]

            for i = 1:length(solution.t)
                push!(infected,solution.u[i][2])
                push!(seriously_infected,solution.u[i][3])
                push!(time, solution.t[i])
            end

            params2.beta = b
            solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2) # Solve the SIR model

            for i = 1:length(solution.t)
                push!(infected,solution.u[i][2])
                push!(seriously_infected,solution.u[i][3])
                push!(time, solution.t[i]+days[1])
            end

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

            # Determines whether to plot infected or seriously infected graph
            if infect == 1
                plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels=nothing, fillrange=previous_infected, colour=:blue) # Plot the model
            else
                plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels=nothing, fillrange=previous_seriously_infected, colour=:blue) # Plot the model
            end

            previous_seriously_infected = deepcopy(seriously_infected)
            previous_infected = deepcopy(infected)
        end
    end

    if infect == 1
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Infected", labels="Infected", colour=:blue) # Plot the model
        plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected", colour=:red) # Plot the model
    else
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected", colour=:blue) # Plot the model
        plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected", colour=:red) # Plot the model
    end
end

"""
error_coverage is a function that plots the RSME of a range of coverage values
using the predicted and actual data provided. It is incorporates the uncertainty
of beta and ratio, and can be plotted logarithmically to better highlight errors.
Inputs:
- S0 = Initial Susceptible Population
- I0 = Initial Infected Population
- SI0 = Initial Seriously Infected Population
- R0 = Initial Recovered Population
- days = No. of days modelled 
- params = array of other necessary parameters for the model before intervention, defined by the structs
- coverage_range = range of coverages that are to be tested
- params2 = array of other necessary parametersfor the model after intervention, defined by the structs
- beta_range = range of transmission rates of disease
- ratio_range = range of proportions of people that become seriously infected
- logarithm = plots the error logarithmically (1 = log, else = not log)
"""
function Level4_2_error_coverage(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, logarithm)
    # Actual data up to day 55
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119]
    tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]

    coverages = range(coverage_range[1], coverage_range[2], length=50)
    betas = range(beta_range[1], beta_range[2], length=40)
    ratios = range(ratio_range[1], ratio_range[2], length=10)

    plot([coverage_range[1]],[0], xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing)
    coverage_min = 1
    min = S0

    for r in ratios
        params.SIratio = r
        params2.SIratio = r
        previous_error_vals = Float64[]
        error_vals = Float64[]

        for b in betas
            params.beta = b
            params2.beta = b
            error_vals = Float64[]
            # Loops through range of coverages
            for j in coverages
                solution = solve_SIR(S0, I0, SI0, R0, days[1], params)

                infected = Float64[]
                seriously_infected = Float64[]
                time = Float64[]

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i])
                end

                params2.p = j
                solution = solve_SIR(solution.u[end][1],solution.u[end][2],solution.u[end][3],solution.u[end][4],days[2], params2)

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i]+days[1])
                end

                current_infected_error = 0
                current_seriously_infected_error = 0

                # Loop through each solution value
                for j = 1:length(time)
                    for k = 1:length(ti)
                        if ti[k] ≈ time[j] atol=0.01 # Check whether time value as approx equal
                            current_infected_error += (infected[j] - actual_infected[k])^2 # Calculate error
                        end
                    end

                    for k = 1:length(tsi)
                        if tsi[k] ≈ time[j] atol=0.01
                            current_seriously_infected_error += (seriously_infected[j] - actual_seriously_infected[k])^2
                        end
                    end
                end
                if logarithm == 1
                    push!(error_vals, log(sqrt(current_infected_error + current_seriously_infected_error))) # Calculate RSME at coverage value
                else
                    push!(error_vals, sqrt(current_infected_error + current_seriously_infected_error)) # Calculate RSME at coverage value
                end
            end

            if minimum(error_vals) < min
                min = minimum(error_vals)
                coverage_min = coverages[argmin(error_vals)]
            end
            if isempty(previous_error_vals)
                previous_error_vals = deepcopy(error_vals)
            end
    
            while length(previous_error_vals) < length(error_vals)
                push!(previous_error_vals,previous_error_vals[end])
            end
    
            plot!(coverages, error_vals, xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing, fillrange=previous_error_vals, colour=:blue)
    
            previous_error_vals = deepcopy(error_vals)
        end

    end

    println("Coverage min: ", coverage_min)
    plot!([0], [0], xlabel="Coverage", ylabel="Error", title="Coverage vs Error", labels=nothing)
end

function Level4_2_plot_coverage_with_error(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range)
    max_infected_peak = 0
    min_infected_peak = 0
    max_seriously_infected_peak = 0
    min_seriously_infected_peak = 0

    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, colour=:black)

    ratios = range(ratio_range[1], ratio_range[2], 10)
    betas = range(beta_range[1], beta_range[2], length=30)
    coverages = range(coverage_range[1], coverage_range[2], length=30)

    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[] 

    for c in coverages
        params2.p = c
        for r in ratios
            params.SIratio = r
            params2.SIratio = r
            previous_infected = Float64[]
            previous_seriously_infected = Float64[]

            # Loops through range of betas
            for j in betas
                params.beta = j # Changes value to new beta
                solution = solve_SIR(S0, I0, SI0, R0, days[1], params)

                infected = Float64[]
                seriously_infected = Float64[]
                time = Float64[]

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i])
                end

                params2.beta = j
                solution = solve_SIR(solution.u[end][1],solution.u[end][2],solution.u[end][3],solution.u[end][4],days[2], params2)

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i]+days[1])
                end

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

                plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, color=:blue, fillrange=previous_infected)
                plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, color=:red, fillrange=previous_seriously_infected)

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

                previous_infected = deepcopy(infected)
                previous_seriously_infected = deepcopy(seriously_infected)
            end
        end
    end

    println("Peak Infected: ", round(min_infected_peak,digits=0), " - ", round(max_infected_peak,digits=0))
    println("Peak Seriously Infected: ", round(min_seriously_infected_peak,digits=0), " - ", round(max_seriously_infected_peak,digits=0))

    plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels="Infected", color=:blue)
    plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels="Seriously Infected", color=:red)
end


function Level4_2_compare_intervention_coverage(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, infect, intervene)
    max_infected_peak = 0
    min_infected_peak = 0
    max_seriously_infected_peak = 0
    min_seriously_infected_peak = 0
    max_infected_peak_intervention = 0
    min_infected_peak_intervention = 0
    max_seriously_infected_peak_intervention = 0
    min_seriously_infected_peak_intervention = 0

    ratios = range(ratio_range[1], ratio_range[2], length=20)
    betas = range(beta_range[1], beta_range[2], length=20)
    coverages = range(coverage_range[1], coverage_range[2], length=15)
    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title=" ", labels=nothing, colour=:black)

    if intervene != 1
        for c in coverages
            params2.p = c
            for r in ratios
                params.SIratio = r
                # Loops through range of betas
                previous_infected = Float64[]
                previous_seriously_infected = Float64[]
                for j in betas
                    params.beta = j # Changes value to new beta
                    solution = solve_SIR(S0, I0, SI0, R0, days[1]+days[2], params)

                    infected = Float64[]
                    seriously_infected = Float64[]

                    for i = 1:length(solution.t)
                        push!(infected,solution.u[i][2])
                        push!(seriously_infected,solution.u[i][3])
                    end

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

                    if infect == 1
                        plot!(solution.t, infected, xlabel="Time (Days)", ylabel="Population", title="Infected No Intervention", labels=nothing, color=:red, fillrange=previous_infected)
                    else
                        plot!(solution.t, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected No Intervention", labels=nothing, color=:red, fillrange=previous_seriously_infected)
                    end

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

                    previous_infected = deepcopy(infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end
            end
        end

    else 
        for c in coverages
            params2.p = c
            for r in ratios
                params.SIratio = r
                params2.SIratio = r
                previous_infected = Float64[]
                previous_seriously_infected = Float64[]
                # Loops through range of betas
                for j in betas
                    params.beta = j # Changes value to new beta
                    solution = solve_SIR(S0, I0, SI0, R0, days[1], params)

                    infected = Float64[]
                    seriously_infected = Float64[]
                    time = Float64[]

                    for i = 1:length(solution.t)
                        push!(infected,solution.u[i][2])
                        push!(seriously_infected,solution.u[i][3])
                        push!(time, solution.t[i])
                    end

                    params2.beta = j
                    solution = solve_SIR(solution.u[end][1],solution.u[end][2],solution.u[end][3],solution.u[end][4],days[2], params2)

                    for i = 1:length(solution.t)
                        push!(infected,solution.u[i][2])
                        push!(seriously_infected,solution.u[i][3])
                        push!(time, solution.t[i]+days[1])
                    end

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

                    if infect == 1
                        plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected Intervention", labels=nothing, color=:blue, fillrange=previous_infected)
                    else
                        plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected No Intervention", labels=nothing, color=:blue, fillrange=previous_seriously_infected)
                    end

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

                    previous_infected = deepcopy(infected)
                    previous_seriously_infected = deepcopy(seriously_infected)
                end
            end
        end
    end

    if infect == 1
        if intervene != 1
            println("Peak Infected No Intervention: ", round(min_infected_peak,digits=0), " - ", round(max_infected_peak,digits=0))
        else
            println("Peak Infected Intervention: ", round(min_infected_peak_intervention,digits=0), " - ", round(max_infected_peak_intervention,digits=0))
        end
    else
        if intervene != 1
            println("Peak Seriously Infected No Intervention: ", round(min_seriously_infected_peak,digits=0), " - ", round(max_seriously_infected_peak,digits=0))
        else
            println("Peak Seriously Infected Intervention: ", round(min_seriously_infected_peak_intervention,digits=0), " - ", round(max_seriously_infected_peak_intervention,digits=0))
        end
    end
    plot!([0],[0], labels=nothing)
end

function Level5_best_beta_day_plot(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, town, day_range)
    # Actual data up to day 80
    actual_infected = [21, 29, 25, 30, 28, 34, 28, 54, 57,92,73,80,109,102,128,135,163,150,211,196,233,247,283,286,332,371,390,404,467,529,598,
    641,704,702,788,856,854,955,995,1065,1106,1159,1217,1269,1298,1328,1339,1383,1431,1422,1414,1485,1464,1480]
    ti = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    actual_seriously_infected = [3, 3, 4, 7, 3, 8, 7, 5, 9,13,15,3,20,13,11,20,16,11,15,18,27,24,28,36,41,35,41,55,63,66,72,80,90,104,109,
    115,127,135,147,162,163,186,194,200,216,223,241,249,258,275,277,299,302,300]
    tsi = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    if town == 1
        # Actual data up to day 55
        actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352,385,221,420,544,329,440,427,369,606,416,546,475,617,593,352,337,473,673,653,523,602,551,686,556,600]
        ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
        actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119,76,95,28,109,136,119,104,121,93,147,129,130,161,133,136,138,139,181,181,218,183,167,164,219,220]
        tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    end

    # Initialize ranges for beta, days, coverage, and ratio
    betas = range(beta_range[1], beta_range[2], length=40)
    start_days = range(day_range[1], day_range[2], length=day_range[2]-day_range[1]+1)
    coverages = range(coverage_range[1], coverage_range[2], length=20)
    ratios = range(ratio_range[1], ratio_range[2], length=20)

    plot([beta_range[1]],[0], labels=nothing)

    # Track minimum error and corresponding beta and day
    min_error = Inf
    best_beta = 0
    best_day = 0

    # Initialize array to store error for plotting
    beta_errors = Float64[Inf]

    # Loop through each beta
    for c in coverages
        params2.p = c
        for r in ratios
            params.SIratio = r
            params2.SIratio = r
            for start_day in start_days
                error_vals = Float64[]
                for b in betas
                    params.beta = b
                    params2.beta = b

                    # Solve SIR model for current beta, start day, coverage, and ratio
                    solution = solve_SIR(S0, I0, SI0, R0, days[1] - (start_day), params)

                    infected = Float64[]
                    seriously_infected = Float64[]
                    time = Float64[]

                    # Collect results from the first phase
                    for i = 1:length(solution.t)
                        push!(infected, solution.u[i][2])
                        push!(seriously_infected, solution.u[i][3])
                        push!(time, solution.t[i]+start_day)
                    end

                    # Re-solve for second phase
                    solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2)

                    for i = 1:length(solution.t)
                        push!(infected, solution.u[i][2])
                        push!(seriously_infected, solution.u[i][3])
                        push!(time, solution.t[i] + days[1])
                    end

                    # Calculate error (RMSE)
                    current_infected_error = 0
                    current_seriously_infected_error = 0

                    for j = 1:length(time)
                        for k = 1:length(ti)
                            if abs(ti[k] - time[j]) < 0.01
                                current_infected_error += (infected[j] - actual_infected[k])^2
                            end
                        end
                        for k = 1:length(tsi)
                            if abs(tsi[k] - time[j]) < 0.01
                                current_seriously_infected_error += (seriously_infected[j] - actual_seriously_infected[k])^2
                            end
                        end
                    end

                    push!(error_vals, sqrt(current_infected_error + current_seriously_infected_error)) # Calculate RSME at coverage value
                    # Track minimum error across beta and days
                    if sqrt(current_infected_error + current_seriously_infected_error) < min_error
                        min_error = sqrt(current_infected_error + current_seriously_infected_error)
                        best_beta = b
                        best_day = start_day
                    end
                end
                plot!([betas], [error_vals], labels=nothing, colour=:blue)
                if minimum(error_vals) < minimum(beta_errors)
                    beta_errors = error_vals
                end

            end
        end
    end

    # Plot error vs beta
    println("Best Beta: ", best_beta)
    println("Start Day: ", round(best_day,digits=0))
    plot!([betas], [beta_errors], xlabel="Beta", ylabel="Error", title="Error vs Beta across Days", labels=nothing, colour=:red)
end

function Level5_plot_second_town(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range, coverage_range, start_day, town, infect, plot_actual)
    # Actual data up to day 80
    actual_infected = [21, 29, 25, 30, 28, 34, 28, 54, 57,92,73,80,109,102,128,135,163,150,211,196,233,247,283,286,332,371,390,404,467,529,598,
    641,704,702,788,856,854,955,995,1065,1106,1159,1217,1269,1298,1328,1339,1383,1431,1422,1414,1485,1464,1480]
    ti = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    actual_seriously_infected = [3, 3, 4, 7, 3, 8, 7, 5, 9,13,15,3,20,13,11,20,16,11,15,18,27,24,28,36,41,35,41,55,63,66,72,80,90,104,109,
    115,127,135,147,162,163,186,194,200,216,223,241,249,258,275,277,299,302,300]
    tsi = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    if town == 1
        # Actual data up to day 55
        actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352,385,221,420,544,329,440,427,369,606,416,546,475,617,593,352,337,473,673,653,523,602,551,686,556,600]
        ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
        actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119,76,95,28,109,136,119,104,121,93,147,129,130,161,133,136,138,139,181,181,218,183,167,164,219,220]
        tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]
    end

    plot([start_day],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing)

    ratios = range(ratio_range[1], ratio_range[2], 15)
    betas = range(beta_range[1], beta_range[2], 40)
    coverages = range(coverage_range[1], coverage_range[2],15)

    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[]

    for c in coverages
        params2.p = c
        for r in ratios
            params.SIratio = r
            params2.SIratio = r
            previous_infected = Float64[]
            previous_seriously_infected = Float64[]


            for b in betas
                params.beta = b
                solution = solve_SIR(S0, I0, SI0, R0, days[1]-(start_day), params) # Solve the SIR model

                infected = Float64[]
                seriously_infected = Float64[]
                time = Float64[]

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i]+(start_day))
                end

                params2.beta = b
                solution = solve_SIR(solution.u[end][1], solution.u[end][2], solution.u[end][3], solution.u[end][4], days[2], params2) # Solve the SIR model

                for i = 1:length(solution.t)
                    push!(infected,solution.u[i][2])
                    push!(seriously_infected,solution.u[i][3])
                    push!(time, solution.t[i]+days[1])
                end

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

                # Determines whether to plot infected or seriously infected graph
                if infect == 1
                    plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels=nothing, fillrange=previous_infected, colour=:blue) # Plot the model
                else
                    plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels=nothing, fillrange=previous_seriously_infected, colour=:blue) # Plot the model
                end

                previous_seriously_infected = deepcopy(seriously_infected)
                previous_infected = deepcopy(infected)
            end
        end
    end

    if infect == 1
        if plot_actual == 1
            plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected", colour=:red) # Plot the model
        end
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Infected", labels="Infected", colour=:blue) # Plot the model
    else
        if plot_actual == 1
            plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected", colour=:red) # Plot the model
        end
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected", colour=:blue) # Plot the model
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