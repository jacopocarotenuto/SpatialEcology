using Plots, DifferentialEquations, LinearAlgebra, TernaryPlots

# Global Parameters
const m::Int = 10 # Number of Species
const p::Int = 3  # Number of Nutrients
const L::Float32 = 1.0 # Length of the domain (periodic boundary conditions)

function BuildMetabolicMatrix(EnzymeBudget::Float64=1.0)
    MetabolicMatrix = rand(m,p)
    for i in 1:m
        MetabolicMatrix[i,:] = MetabolicMatrix[i,:]./sum(MetabolicMatrix[i,:]) * EnzymeBudget
    end
    return MetabolicMatrix
end

function PlotTernaryMetabolicMatrix(MetabolicMatrix::Matrix{Float64}, title::String = "Metabolic Rates", savepath::String="Latex/figures/"*title*".pdf")
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
    savepath = replace(savepath, " " => "_")
    savefig(savepath)
    return fig
end

# Parameters
SupplyVector = rand(p)
S = sum(SupplyVector)
MetabolicMatrix = BuildMetabolicMatrix()
PlotTernaryMetabolicMatrix(MetabolicMatrix, "Metabolic Rates")