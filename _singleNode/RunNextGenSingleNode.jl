using DifferentialEquations,Plots,Parameters,LinearAlgebra,JLD,Statistics
include("functions.jl")
include("../_globalFunctions/_includes.jl")


NGp = NextGen2PopParams()

NGp.η_0I = -10.
NGp.η_0E = 20
NGp.κSII = 5.0
NGp.κSEI = 2.4
NGp.κSEE = 5.0
NGp.κSEI = 3.6
NGp.κ=0.0203

tspan = (0.,100000.)
const AMP=0. # or 10.
const fr = 2.

const U0 = 0.1



u0 =rand(14)

#u0[1] = 2.257
#u0[2] = 0.036776
#u0[3] = 1.6943
#u0[4] = -0.98019
#u0[5] = 3.3855
#u0[6] = u0[5]
#u0[7] = 2.257
#u0[8] = u0[7]
#u0[9] = 0.073552
#u0[10] = u0[9]
#u0[11] = 0.11033
#u0[12] = u0[11]


#u0[39] = 100.

u0[13] = NGp.κSEI

p = NGp
probSS = SteadyStateProblem(nextgen_stp_de,u0,p)
prob = ODEProblem(nextgen_stp_de,u0,tspan,p)

save_end = 30000
sol = solve(prob,maxiters=10e20,saveat=save_start = 0:0.01:save_end)



plot(sol.t[end-2000:end],sol[1,end-2000:end],label="rE")
plot!(sol.t[end-2000:end],sol[2,end-2000:end],label="rI")

