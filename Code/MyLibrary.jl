using Plots, DifferentialEquations, LinearAlgebra, TernaryPlots, QuadGK, LinearSolve, DataFrames

####### Data Structures #######

@kwdef struct SimulationParameters
    m::Int64
    p::Int64
    L::Float64
    space::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    D::Float64
    S::Float64
    v::Float64
    δ::Float64
    EnzymeBudget::Float64
    MetabolicMatrix::Matrix{Float64}
    SupplyVector::Vector{Float64}
end

@kwdef mutable struct WorldData
    concentration_matrix::Matrix{Float64}
    A::Matrix{Float64}
    B::Matrix{Float64}
    n::Vector{Float64}
end

####### Core Functions #######

function BuildMetabolicMatrix(EnzymeBudget::Float64, m::Int, p::Int)
    MetabolicMatrix = rand(m,p)
    for i in 1:m
        MetabolicMatrix[i,:] = MetabolicMatrix[i,:]./sum(MetabolicMatrix[i,:]) * EnzymeBudget
    end
    return MetabolicMatrix
end

function BuildSupplyVector(S::Float64, p::Int, method::String="random")
    if method == "random"
        SupplyVector = rand(p)
        return SupplyVector / sum(SupplyVector) * S
    else
        print("Method not implemented")
    end
end

function PlotTernaryMetabolicMatrix(SimulationParameters,title::String = "Metabolic Rates", savepath::String="Latex/figures/"*title*".pdf")

    if SimulationParameters.p != 3
        raise("This function only works for 3 nutrients")
    end
    MetabolicMatrix = SimulationParameters.MetabolicMatrix
    m = SimulationParameters.m
    cartesian = zeros(m,2)
    for i in 1:m
        cartesian[i,:] = collect(tern2cart(MetabolicMatrix[i,:]))'
    end

    ternary_axes(
        title=title,
        xguide="α₁",
        yguide="α₂",
        zguide="α₃",
    )

    fig = scatter!(cartesian[:,1],cartesian[:,2], legend=false)
    # Plot Supply Vector
    cartesian_supply_vector = collect(tern2cart(SimulationParameters.SupplyVector))
    scatter!([cartesian_supply_vector[1]],[cartesian_supply_vector[2]],label="Supply Vector",color=:black,markersize=5)
    savepath = replace(savepath, " " => "_")
    savefig(savepath)
    return fig
end

function BuildInitialSpeciesTerritory(SimulationParameters, method::String="random")
    if method == "random"
        n = rand(SimulationParameters.m)
        return n / sum(n) * SimulationParameters.L
    elseif method == "uniform"
        n = ones(SimulationParameters.m)
        return n / sum(n) * SimulationParameters.L
    else
        print("Method not implemented")
    end
end

function SetupWorldData(SimulationParameters, InitialSpecies = "random")
    World = WorldData(
        concentration_matrix = zeros(length(SimulationParameters.space),SimulationParameters.p),
        A = zeros(SimulationParameters.m,SimulationParameters.p),
        B = zeros(SimulationParameters.m,SimulationParameters.p),
        n = BuildInitialSpeciesTerritory(SimulationParameters, InitialSpecies)
    )
    return World
end

function CalculateC(x, species_index, WorldData, SimulationParameters)
    # Calculate concentration for all nutrients in a given place for a given species
    SupplyVector = SimulationParameters.SupplyVector
    MetabolicMatrix = SimulationParameters.MetabolicMatrix
    D = SimulationParameters.D
    α = MetabolicMatrix[species_index,:]
    A = WorldData.A[species_index,:]
    B = WorldData.B[species_index,:]

    c = @. SupplyVector / α + A * exp(x * sqrt(α / D)) + B * exp(-x * sqrt(α / D))
    return c
end

function UpdateConcentration!(WorldData, SimulationParameters)
    # Extract relevant parameters
    space = SimulationParameters.space
    species_boundaries = cumsum(WorldData.n)
    species_boundaries[end] = SimulationParameters.L
    ComputeAllABCoefficients!(WorldData, SimulationParameters)
    species_index = 1
    for (i,x) in enumerate(space)
        if x > species_boundaries[species_index]
            species_index += 1
        end
        if species_index == 1
            c = CalculateC(x, species_index, WorldData, SimulationParameters)
        else
            c = CalculateC(x - species_boundaries[species_index - 1], species_index, WorldData, SimulationParameters)
        end
        WorldData.concentration_matrix[i,:] = c
    end
end

function ComputeABCoefficients!(nutrient_index, WorldData, SimulationParameters)
    m = SimulationParameters.m
    A = zeros(m*2,m*2)
    b = zeros(m*2)
    MetabolicMatrix = SimulationParameters.MetabolicMatrix
    SupplyVector = SimulationParameters.SupplyVector
    n = WorldData.n
    D = SimulationParameters.D


    for i in 1:m
        if i == m
            # From continuity of c
            A[i,i*2-1]   = exp(  n[i] * sqrt(MetabolicMatrix[i,nutrient_index] / D))
            A[i,i*2]     = exp(- n[i] * sqrt(MetabolicMatrix[i,nutrient_index] / D))
            A[i,1]       = -1
            A[i,2]       = -1
            b[i] = SupplyVector[nutrient_index] / MetabolicMatrix[1,nutrient_index] - SupplyVector[nutrient_index] / MetabolicMatrix[i,nutrient_index]

            # From continuity of dc
            A[m + i,i*2-1] =   sqrt(MetabolicMatrix[i,nutrient_index] / D) * exp(  n[i] * sqrt(MetabolicMatrix[i,nutrient_index] / D))
            A[m + i,i*2]   = - sqrt(MetabolicMatrix[i,nutrient_index] / D) * exp(- n[i] * sqrt(MetabolicMatrix[i,nutrient_index] / D))
            A[m + i,1]     = - sqrt(MetabolicMatrix[1,nutrient_index] / D)
            A[m + i,2]     = + sqrt(MetabolicMatrix[1,nutrient_index] / D)
            b[m + i] = 0
        else
            # From continuity of c
            A[i,i*2-1]   = exp(  n[i] * sqrt(MetabolicMatrix[i,nutrient_index] / D))
            A[i,i*2]     = exp(- n[i] * sqrt(MetabolicMatrix[i,nutrient_index] / D))
            A[i,i*2+1]   = -1
            A[i,i*2+2]   = -1
            b[i] = SupplyVector[nutrient_index] / MetabolicMatrix[i+1,nutrient_index] - SupplyVector[nutrient_index] / MetabolicMatrix[i,nutrient_index]

            # From continuity of dc
            A[m + i,i*2-1] =   sqrt(MetabolicMatrix[i,nutrient_index] / D) * exp(  n[i] * sqrt(MetabolicMatrix[i,nutrient_index] / D))
            A[m + i,i*2]   = - sqrt(MetabolicMatrix[i,nutrient_index] / D) * exp(- n[i] * sqrt(MetabolicMatrix[i,nutrient_index] / D))
            A[m + i,i*2+1] = - sqrt(MetabolicMatrix[i+1,nutrient_index] / D)
            A[m + i,i*2+2] = + sqrt(MetabolicMatrix[i+1,nutrient_index] / D)
            b[m + i] = 0
        end
    end
    prob = LinearProblem(A,b)
    AB = solve(prob)
    A = AB[1:2:end]
    B = AB[2:2:end]
    WorldData.A[:,nutrient_index] = A
    WorldData.B[:,nutrient_index] = B
    return nothing
end

function ComputeAllABCoefficients!(WorldData, SimulationParameters)
        ComputeABCoefficients!.(1:SimulationParameters.p, (WorldData,), (SimulationParameters,))
    return nothing
end

function CalculateGrowth(species_index, WorldData, SimulationParameters)
    MetabolicMatrix = SimulationParameters.MetabolicMatrix
    SupplyVector = SimulationParameters.SupplyVector
    α = MetabolicMatrix[species_index,:]
    A = WorldData.A[species_index,:]
    B = WorldData.B[species_index,:]
    n0 = WorldData.n[species_index]
    D = SimulationParameters.D

    return @. α*( SupplyVector / α*n0 + A * sqrt(D/α) * ( exp(n0*sqrt(α/D)) - 1 ) - B * sqrt(D/α) * ( exp(-n0*sqrt(α/D)) - 1 ) )
end

function SpatialEcologyDynamics!(du, u, p, t)
    WorldData = p[1]
    WorldData.n = u
    SimulationParameters = p[2]
    δ = SimulationParameters.δ
    UpdateConcentration!(WorldData, SimulationParameters)
    growths = CalculateGrowth.(1:SimulationParameters.m, (WorldData,), (SimulationParameters,))
    du .= sum.(growths) .- δ .* u
end

function AdvanceDynamics(t,WorldData, SimulationParameters, fast = false)
    u0 = WorldData.n
    tspan = (0.0, t)
    # Save only last state
    prob = ODEProblem(SpatialEcologyDynamics!, u0, tspan, (WorldData, SimulationParameters), save_everystep = !fast)
    sol = solve(prob)
    WorldData.n = sol.u[end]
    return sol
end

function ArrayToMatrix(u)
    matrix = zeros(length(u), length(u[1]))
    for i in eachindex(u)
        matrix[i,:] = u[i]
    end
    return matrix
end

function CompleteDynamics(t, m, p, SupplyVector, fast = false, D = 1.0)
    SimParam = SimulationParameters(
        m = m,                  # Number of Species
        p = p,                  # Number of Nutrients
        L = 1.0,                # Length of the domain (periodic boundary conditions)
        space = 0:0.01:1.0,     # Discretization of the space
        D = D,                # Diffusion Rate
        S = 1.0,                # Total Growth
        v = 1.0,                # Length-to-Population Conversion Factor
        δ = 1.0,                # Death Rate
        EnzymeBudget = 2.0,     # Total Enzyme Budget
        MetabolicMatrix = BuildMetabolicMatrix(2.0, m, p),# Metabolic Matrix
        SupplyVector = SupplyVector          # Supply Vector
    );

    World = SetupWorldData(SimParam);
    sol = AdvanceDynamics(t, World, SimParam, fast);
    return sol
end

function ExtractStatisticsOnDynamics(sol)

    # Calculate how many species remains
    SurvivingSpecies = sum(sol.u[end] .> 10^(-6))

    return Vector{Float32}([SurvivingSpecies])
end

function CollectStatisticsOnDynamics(realizations::Int64, m::Int64, p::Int64, SupplyVector::Vector{Float64}, D::Float64 )
    df = DataFrame(SurvivingSpecies = Float32[])
    for i in 1:realizations
        sol = CompleteDynamics(2000, m, p, SupplyVector, true, D)
        stats = ExtractStatisticsOnDynamics(sol)
        push!(df, stats)
    end
    return df
end

####### Plotting Functions #######

function PlotConcentrationInSpace(WorldData, SimulationParameters)
    UpdateConcentration!(WorldData, SimulationParameters)
    labels = ["Nutrient 1" "Nutrient 2" "Nutrient 3"]
    space = SimulationParameters.space
    concentration = WorldData.concentration_matrix
    n = WorldData.n
    plot(space, concentration, label=labels, xlabel="Space", ylabel="Concentration", title="Concentration of Nutrients")
    vline!(cumsum(n), label="Species", color=:black, linestyle=:dash, linewidth=1, alpha=0.5)
end

function PlotTerritoryInTime(sol, savepath::String="Latex/figures/TerritoryInTime.pdf")
    TerritoryInTime = ArrayToMatrix(sol.u)
    time = sol.t
    fig = portfoliocomposition(TerritoryInTime, sol.t, legend = false)
    plot!(framestyle = :box, grid = false, axis = false)
    plot!(xlabel = "Space (L) →", ylabel = "Time →", title = "Territory in Time: n₁, n₂, n₃...")
    savefig(savepath)
    return fig
end

function PlotPopulationInTime(sol, title="Population in Time", savepath="Latex/figures/"*title*".pdf")
    mat = ArrayToMatrix(sol.u)
    labels = Matrix{String}(undef, 1, size(mat,2))
    for i in 1:size(mat,2)
        labels[i] = "Species $i"
    end

    if minimum(mat) < 0.001
        ylims = (0.001,20)
    else
        ylims = (minimum(mat)/5,20)
    end

    fig = plot(sol.t,mat, yscale=:log10, ylims=ylims, legend=:bottomright, label=labels, title=title, xaxis="Time", yaxis="Population")
    savefig(fig, savepath)
    return fig
end
