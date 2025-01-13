using Plots, DifferentialEquations, LinearAlgebra, TernaryPlots, QuadGK, LinearSolve, NonlinearSolve, DataFrames

# Include a file
include("MyLibrary.jl");

# Setup Simulation Parameters
SimParam = SimulationParameters(
    m = 8,                  # Number of Species
    p = 3,                  # Number of Nutrients
    L = 1.0,                # Length of the domain (periodic boundary conditions)
    space = 0:0.01:1.0,     # Discretization of the space
    D = 0.01,                # Diffusion Rate
    S = 1.0,                # Total Growth
    v = 1.0,                # Length-to-Population Conversion Factor
    δ = 1.0,                # Death Rate
    EnzymeBudget = 2.0,     # Total Enzyme Budget
    MetabolicMatrix = BuildMetabolicMatrix(2.0, 8, 3),# Metabolic Matrix
    SupplyVector = [0.25,0.45,0.3]          # Supply Vector
);

# Setup World Data
World = SetupWorldData(SimParam, "uniform");

# Plot Metabolic Matrix in Ternary Space
PlotTernaryMetabolicMatrix(SimParam)
# Plot Concentration in Space
PlotConcentrationInSpace(World, SimParam)
# Run the dynamic
sol = AdvanceDynamics(1500, World, SimParam)
# Plot Concentration in Space
PlotConcentrationInSpace(World, SimParam)
# Plot territory vs Time
PlotTerritoryInTime(sol)
# Plot Population vs Time
PlotPopulationInTime(sol)


###### STATISTICS #######
df = CollectStatisticsOnDynamics(400, 10, 2, [0.4,0.6], 0.01)