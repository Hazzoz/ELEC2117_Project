# Structs for different versions of the SIR models
mutable struct SIRFoI
    beta::Float64  # Transmission chance of interaction
    gamma::Float64  # Recovery rate
    contacts::Float64  # Number of daily contacts
    SIratio::Float64  # Ratio of people becoming seriously infected
    delta::Float64  # Seriously infected recovery rate
    alpha::Float64  # Re-susceptible rate
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
function plot_infected(S0, I0, SI0, R0, days, params, infect)
    solution = solve_SIR(S0, I0, SI0, R0, days, params) # Solve the SIR model

    # Actual data up to day 30
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4]
    tsi = [21,22,23,24,25,26,27,28,29,30]

    infected = []
    seriously_infected = []
    time = []

    for i = 1:length(solution.t)
        push!(infected,solution.u[i][2])
        push!(seriously_infected,solution.u[i][3])
        push!(time, solution.t[i])
    end

    # Print associated R0 value
    R0 = params.contacts*params.beta/params.gamma
    println("R0: ", R0)

    # Determines whether to plot infected or seriously infected graph
    if infect == 0
        plot(time, infected, xlabel="Time", ylabel="Population", title="Infected", labels="Infected") # Plot the model
        plot!(ti, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected") # Plot the model
    else 
        plot(time, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected") # Plot the model
        plot!(tsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected") # Plot the model
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
function error_beta(S0, I0, SI0, R0, days, params)
    # Actual data up to day 30
    actual_infected = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55]
    ti = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    actual_seriously_infected = [0,0,1,2,5,5,5,2,9,4]
    tsi = [21,22,23,24,25,26,27,28,29,30]

    # Create error array and range of betas to test (around the initial guess of 0.037)
    error_vals = []
    betas = range(0.025, 0.042, length=50)

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

    println("Beta min: ", betas[argmin(error_vals)])
    plot(betas, error_vals, xlabel="Beta", ylabel="Error", title="Beta vs Error", labels=nothing)
end