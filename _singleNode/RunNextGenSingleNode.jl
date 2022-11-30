using DifferentialEquations,Plots,Parameters,LinearAlgebra,JLD,Statistics
include("functions.jl")
include("../_globalFunctions/_includes.jl")



NGp = NextGen2PopParams()

NGp.η_0I = -20.
NGp.κ=0.0203

tspan = (0.,30000.)
const AMP=0
const fr = 5.


u0 = rand(12)

u0[1] = 2.257
u0[2] = 0.036776
u0[3] = 1.6943
u0[4] = 0.98019
u0[5] = 3.3855
u0[6] = u0[5]
u0[7] = 2.257
u0[8] = u0[7]
u0[9] = 0.073552
u0[10] = u0[9]
u0[11] = 0.11033
u0[12] = u0[11]


#u0[39] = 100.

p = NGp
prob = ODEProblem(nextgen_stp_de,u0,tspan,p)
sol = solve(prob,maxiters=10e20,saveat=0:0.1:tspan[2])


plot(sol[1,end-1000:end])


