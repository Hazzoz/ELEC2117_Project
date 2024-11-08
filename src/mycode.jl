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

###############################################################
# This SIR! is a function that represents the differential
# equations that defines a more detailed force of infection
# version of the SIR model
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - SIratio = proportion of people being seriously infected
#   - delta = recovery rate of seriously infected
#   - alpha = number of recovered people losing immunity
#   - lambdas = recorded array of lambdas with timestamps
# - t = timespan
###############################################################
function SIR!(dP, P, params::SIRFoI, t)
    N = P[1] + P[2] + P[3] + P[4] # Total Population
    lambda = P[2]/N*params.beta*params.contacts # Force of Infection
    dP[1] = params.alpha*P[4] - lambda*P[1] # Change in Susceptible Population
    dP[2] = lambda*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*params.SIratio*P[2] - params.delta*P[3] # Change in Seriously Infected Population
    dP[4] = params.gamma*(1-params.SIratio)*P[2] + params.delta*P[3] - params.alpha*P[4] # Change in Recovered Population
end

###############################################################
# This SIR! is a function that represents the differential
# equations that defines a more detailed force of infection
# version of the SIR model
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - SIratio = proportion of people being seriously infected
#   - delta = recovery rate of seriously infected
#   - alpha = number of recovered people losing immunity
#   - lambdas = recorded array of lambdas with timestamps
# - t = timespan
###############################################################
function SIR!(dP, P, params::SIRFoIIntervention, t)
    N = P[1] + P[2] + P[3] + P[4] # Total Population
    lambda = P[2]/N*params.beta*params.contacts*(1-params.epsilon*params.p) # Force of Infection
    dP[1] = params.alpha*P[4] - lambda*P[1] # Change in Susceptible Population
    dP[2] = lambda*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*params.SIratio*P[2] - params.delta*P[3] # Change in Seriously Infected Population
    dP[4] = params.gamma*(1-params.SIratio)*P[2] + params.delta*P[3] - params.alpha*P[4] # Change in Recovered Population
end

"""
# solve_SIR is a driver function that chooses the required SIR model.
# Inputs:
# - S0 = Initial Susceptible Population
# - I0 = Initial Infected Population
# - SI0 = Initial Seriously Infected Population
# - R0 = Initial Recovered Population
# - days = No. of days modelled 
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - SIratio = proportion of people being seriously infected
#   - delta = recovery rate of seriously infected
#   - alpha = number of recovered people losing immunity
#   - lambdas = recorded array of lambdas with timestamps
"""
function solve_SIR(S0, I0, SI0, R0, days, params)
    P0 = [S0, I0, SI0, R0] # Initial populations vector

    tspan = (0, days) # Time span tuple

    solution = solve(ODEProblem(SIR!, P0, tspan, params)) # Solve the ODE with given parameters and timespan

    return solution # Return values
end

"""
# plot_SIR is a driver function that runs and then plots and SIR model.
# Inputs:
# - S0 = Initial Susceptible Population
# - I0 = Initial Infected Population
# - SI0 = Initial Seriously Infected Population
# - R0 = Initial Recovered Population
# - days = No. of days modelled 
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - SIratio = proportion of people being seriously infected
#   - delta = recovery rate of seriously infected
#   - alpha = number of recovered people losing immunity
"""
function plot_SIR(S0, I0, SI0, R0, days, params)
    solution = solve_SIR(S0, I0, SI0, R0, days, params) # Solve the SIR model

    plot(solution, xlabel="Time", ylabel="Population", title="Solution", labels=["Susceptible" "Infected" "Seriously Infected" "Recovered"]) # Plot the model
end

"""
# plot_infected is a function that plots either the infected or
# seriously infected data against the true data to allow for visual
# inspection of the models fit.
# Inputs:
# - S0 = Initial Susceptible Population
# - I0 = Initial Infected Population
# - SI0 = Initial Seriously Infected Population
# - R0 = Initial Recovered Population
# - days = No. of days modelled 
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - SIratio = proportion of people being seriously infected
#   - delta = recovery rate of seriously infected
#   - alpha = number of recovered people losing immunity
# - infect = determines whether to plot infected or seriously infected data
"""
function plot_infected(S0, I0, SI0, R0, days, params, infect, ratio_range)
    # Actual data up to day 30
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4]
    tsi = [21,22,23,24,25,26,27,28,29,30]

    plot([0],[0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing)

    ratios = range(ratio_range[1], ratio_range[2], 10)

    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[]
    previous_infected = Float64[]
    previous_seriously_infected = Float64[]

    for j in ratios
        params.SIratio = j
        solution = solve_SIR(S0, I0, SI0, R0, days, params) # Solve the SIR model

        infected = Float64[]
        seriously_infected = Float64[]
        time = Float64[]

        for i = 1:length(solution.t)
            push!(infected,solution.u[i][2])
            push!(seriously_infected,solution.u[i][3])
            push!(time, solution.t[i])
        end

        if isempty(previous_infected)
            previous_infected = deepcopy(infected)
        end
        if isempty(previous_seriously_infected)
            previous_seriously_infected = deepcopy(seriously_infected)
        end

        # Determines whether to plot infected or seriously infected graph
        if infect == 0
            plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels=nothing, fillrange=previous_infected, colour=:blue) # Plot the model
        else
            plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels=nothing, fillrange=previous_seriously_infected, colour=:blue) # Plot the model
        end

        previous_seriously_infected = deepcopy(seriously_infected)
        previous_infected = deepcopy(infected)
    end

    # Print associated R0 value
    R0 = params.contacts*params.beta/params.gamma
    println("R0: ", round(R0, digits=3))

    if infect == 0
        plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels="Infected", colour=:blue) # Plot the model
        plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected", colour=:red) # Plot the model
    else
        plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected", colour=:blue) # Plot the model
        plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected", colour=:red) # Plot the model
    end
end

"""
# error_beta is a function that plots the SME of a range of beta values
# using the predicted and actual data provided. It is placed around the estimate
# from level 2.
# Inputs:
# - S0 = Initial Susceptible Population
# - I0 = Initial Infected Population
# - SI0 = Initial Seriously Infected Population
# - R0 = Initial Recovered Population
# - days = No. of days modelled 
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - SIratio = proportion of people being seriously infected
#   - delta = recovery rate of seriously infected
#   - alpha = number of recovered people losing immunity
"""
function error_beta(S0, I0, SI0, R0, days, params, beta_range, ratio_range)
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
    ratios = range(ratio_range[1], ratio_range[2], 10)

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

function plot_intervention_with_error(S0, I0, SI0, R0, days, params, beta_range, params2, ratio_range)
    max_infected_peak = 0
    min_infected_peak = 0
    max_seriously_infected_peak = 0
    min_seriously_infected_peak = 0

    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title="Infected After Intervention", labels=nothing, colour=:black)

    ratios = range(ratio_range[1], ratio_range[2], 50)

    infected = Float64[]
    seriously_infected = Float64[]
    time = Float64[] 

    for r in ratios
        params.SIratio = r
        params2.SIratio = r
        previous_infected = Float64[]
        previous_seriously_infected = Float64[]
        betas = range(beta_range[1], beta_range[2], length=100)

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


function compare_intervention(S0, I0, SI0, R0, days, params, beta_range, params2, infect, ratio_range)
    max_infected_peak = 0
    min_infected_peak = 0
    max_seriously_infected_peak = 0
    min_seriously_infected_peak = 0
    max_infected_peak_intervention = 0
    min_infected_peak_intervention = 0
    max_seriously_infected_peak_intervention = 0
    min_seriously_infected_peak_intervention = 0

    ratios = range(ratio_range[1], ratio_range[2], length=20)
    betas = range(beta_range[1], beta_range[2], length=100)
    plot([31,31], [0,250], xlabel="Time (Days)", ylabel="Population", title=" ", labels=nothing, colour=:black)

    for r in ratios
        params.SIratio = r
        previous_infected = Float64[]
        previous_seriously_infected = Float64[]

        # Loops through range of betas
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

            if infect == 0
                plot!(solution.t, infected, xlabel="Time (Days)", ylabel="Population", title="Infected", labels=nothing, color=:red, fillrange=previous_infected)
            else
                plot!(solution.t, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected", labels=nothing, color=:red, fillrange=previous_seriously_infected)
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

            if infect == 0
                plot!(time, infected, xlabel="Time (Days)", ylabel="Population", title="Infected", labels=nothing, color=:blue, fillrange=previous_infected, alpha=0.1)
            else
                plot!(time, seriously_infected, xlabel="Time (Days)", ylabel="Population", title="Seriously Infected", labels=nothing, color=:blue, fillrange=previous_seriously_infected, alpha=0.1)
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

    if infect == 0
        println("Peak Infected No Intervention: ", round(min_infected_peak,digits=0), " - ", round(max_infected_peak,digits=0))
        println("Peak Infected Intervention: ", round(min_infected_peak_intervention,digits=0), " - ", round(max_infected_peak_intervention,digits=0))
        plot!([0], [0], xlabel="Time (Days)", ylabel="Population", title="Infected Intervention Effects", labels="No Intervention", color=:red)
        plot!([0], [0], xlabel="Time (Days)", ylabel="Population", title="Infected Intervention Effects", labels="Intervention", color=:blue)
    else
        println("Peak Seriously Infected No Intervention: ", round(min_seriously_infected_peak,digits=0), " - ", round(max_seriously_infected_peak,digits=0))
        println("Peak Seriously Infected Intervention: ", round(min_seriously_infected_peak_intervention,digits=0), " - ", round(max_seriously_infected_peak_intervention,digits=0))
        plot!([0], [0], xlabel="Time (Days)", ylabel="Population", title="Seriously Infected Intervention Effects", labels="No Intervention", color=:red)
        plot!([0], [0], xlabel="Time (Days)", ylabel="Population", title="Seriously Infected Intervention Effects", labels="Intervention", color=:blue)
    end
end

"""
# error_beta is a function that plots the SME of a range of beta values
# using the predicted and actual data provided. It is placed around the estimate
# from level 2.
# Inputs:
# - S0 = Initial Susceptible Population
# - I0 = Initial Infected Population
# - SI0 = Initial Seriously Infected Population
# - R0 = Initial Recovered Population
# - days = No. of days modelled 
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - SIratio = proportion of people being seriously infected
#   - delta = recovery rate of seriously infected
#   - alpha = number of recovered people losing immunity
"""
function error_coverage(S0, I0, SI0, R0, days, params, coverage_range, params2, beta_range, ratio_range)
    # Actual data up to day 55
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55,155,53,67,98,130,189,92,192,145,128,68,74,126,265,154,207,299,273,190,152,276,408,267,462,352]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4,22,0,15,48,38,57,9,18,20,0,41,15,35,36,27,38,24,40,34,57,18,29,63,66,119]
    tsi = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]

    # Create error array and range of coverages to test
    error_vals = Float64[]
    coverages = range(coverage_range[1], coverage_range[2], length=50)
    betas = range(beta_range[1], beta_range[2], length=50)
    ratios = range(ratio_range[1], ratio_range[2], length=50)

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
                push!(error_vals, sqrt(current_infected_error + current_seriously_infected_error)) # Calculate RSME at coverage value
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
    plot!(coverages, error_vals, xlabel="Coverage", ylabel="Error", title="Coverage vs Error", labels=nothing)
end

function plot_coverage(S0, I0, SI0, R0, days, params, infect, params2, ratio_range, beta_range)
    # Actual data up to day 30
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
            if infect == 0
                plot!(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels=nothing, fillrange=previous_infected, colour=:blue) # Plot the model
            else
                plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels=nothing, fillrange=previous_seriously_infected, colour=:blue) # Plot the model
            end

            previous_seriously_infected = deepcopy(seriously_infected)
            previous_infected = deepcopy(infected)
        end
    end

    if infect == 0
        plot!([0], [0], xlabel="Time", ylabel="Population", title="Infected", labels=nothing, colour=:blue) # Plot the model
        plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected", colour=:red) # Plot the model
    else
        plot!(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected", colour=:blue) # Plot the model
        plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected", colour=:red) # Plot the model
    end
end