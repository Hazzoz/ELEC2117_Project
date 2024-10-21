# Structs for different versions of the SIR models
struct SIRFoI
    beta::Float64  # Transmission chance of interaction
    gamma::Float64  # Recovery rate
    contacts::Float64  # Number of daily contacts
    SIratio::Float64  # Ratio of people becoming seriously infected
    delta::Float64  # Seriously infected recovery rate
    alpha::Float64  # Re-susceptible rate
    lambdas::Vector{Tuple{Float64,Float64}} # Recorded array of lambdas
end

struct SIRFoIIntervention
    beta::Float64  # Transmission chance of interaction
    gamma::Float64  # Recovery rate
    contacts::Float64  # Number of daily contacts
    SIratio::Float64  # Ratio of people becoming seriously infected
    delta::Float64  # Seriously infected recovery rate
    alpha::Float64  # Re-susceptible rate
    epsilon::Float64  # Intervention Efficacy Rate
    p::Float64 # Coverage rate
    lambdas::Vector{Tuple{Float64,Float64}} # Recorded array of lambdas
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
function SIR!(dP, P, params::SIRForceOfInfection, t)
    N = P[1] + P[2] + P[3] + P[4] # Total Population
    lambda = P[2]/N*params.beta*params.contacts # Force of Infection
    dP[1] = params.alpha*P[4] - lambda*P[1] # Change in Susceptible Population
    dP[2] = lambda*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*params.SIratio*P[2] - params.delta*P[3] # Change in Seriously Infected Population
    dP[4] = params.gamma*(1-params.SIratio)*P[2] + params.delta*P[3] - params.alpha*P[4] # Change in Recovered Population

    # Store lambda at this time step
    push!(params.lambdas, (lambda, t)) # Storing lambda values
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
function SIR!(dP, P, params::SIRForceOfInfectionIntervention, t)
    N = P[1] + P[2] + P[3] + P[4] # Total Population
    lambda = P[2]/N*params.beta*params.contacts # Force of Infection
    dP[1] = params.alpha*P[4] - lambda*P[1] # Change in Susceptible Population
    dP[2] = lambda*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*params.SIratio*P[2] - params.delta*P[3] # Change in Seriously Infected Population
    dP[4] = params.gamma*(1-params.SIratio)*P[2] + params.delta*P[3] - params.alpha*P[4] # Change in Recovered Population

    # Store lambda at this time step
    push!(params.lambdas, (lambda, t)) # Storing lambda values
end

"""
# solve_SIR is a driver function that chooses the required SIR
# model
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
#   - epsilon = intervention efficacy
#   - p = coverage rate
"""
function solve_SIR(S0, I0, SI0, R0, days, params)
    P0 = [S0, I0, SI0, R0] # Initial populations vector

    tspan = (0, days) # Time span tuple

    solution = solve(ODEProblem(SIR!, P0, tspan, params)) # Solve the ODE with given parameters and timespan

    return solution, params.lambdas # Return values
end

"""
# solve_SIR is a driver function that chooses the required SIR
# model
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
#   - epsilon = intervention efficacy
#   - p = coverage rate
"""
function plot_SIR(S0, I0, SI0, R0, days, params)
    solution, lambdas = solve_SIR(S0, I0, SI0, R0, days, params) # Solve the SIR model

    plot(solution, xlabel="Time", ylabel="Population", title="Solution", labels=["Susceptible" "Infected" "Seriously Infected" "Recovered"]) # Plot the model
end

"""
# solve_SIR is a driver function that chooses the required SIR
# model
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
#   - epsilon = intervention efficacy
#   - p = coverage rate
"""
function plot_infected(S0, I0, SI0, R0, days, params, infected)
    solution, lambdas = solve_SIR(S0, I0, SI0, R0, days, params) # Solve the SIR model

    actual_infected = [11,7,20,3,29,14,11,12,16,10,58]
    xi = [15,16,17,18,19,20,21,22,23,24,25]
    actual_seriously_infected = [0,0,1,2,5]
    xsi = [21,22,23,24,25]

    infected = []
    seriously_infected = []

    for i = 1:length(solution.t)
        push!(infected,solution.u[i][2])
        push!(seriously_infected,solution.u[i][3])
    end

    R0 = params.contacts*params.beta/params.gamma
    pc = 1 - 1/R0 # Herd immunity threshold
    #println("R0: ", R0)

    if infected == 0
        plot(solution.t, infected, xlabel="Time", ylabel="Population", title="Infected", labels="Infected") # Plot the model
        plot!(xi, actual_infected, xlabel="Time", ylabel="Population", title="Infected", labels="Actual Infected") # Plot the model
    else 
        plot(solution.t, seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Seriously Infected") # Plot the model
        plot!(xsi, actual_seriously_infected, xlabel="Time", ylabel="Population", title="Seriously Infected", labels="Actual Seriously Infected") # Plot the model
    end
end